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
	  
	  
	  module ModDYNConvectivePhase
      !**********************************************************************
      !
      !    Function:  Contains the routines related to updating the mesh and
      !               particle data at the end of the Lagrangian Phase.
      !               This module is only used by the Dynamic MPM.
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               the updating of particle data.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr. 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMPMData
      use ModWriteTestData
      use ModLagrangianPhase
      use ModConvectivePhase
      use ModMPMDYNStresses
      use ModMeshAdjacencies
      use ModMPMMeshAdjustment
      use ModMeshInfo
      use ModRotBoundCond
      use ModEmptyElements
      use ModAdjustParticleDiscretisation
      use ModGlobalConstants
      use ModTwoLayerFormulation
      use ModString
      use ModRigidBody
      use ModLiquid
      use ModReadGeometryData
      
      implicit none

        real(REAL_TYPE), dimension(:),  &
          allocatable :: TemporaryMappingVector

      contains ! Routines of this module

      subroutine DYNConvectivePhase()
        !**********************************************************************
        !
        !    Function:  Calls the different subroutines required for updating the
        !               particle data. The basic steps are:
        !               
        !               - update particle velocity
        !               - update particle displacements and global position
        !               - map the new particles ve  locities to the nodes
        !               - update the nodal coordinates from the new nodal velocities
        !               - calculate particle strains 
        !               - calculate stresses for integration points
        !               - map stresses from Gauss points to particles for fully filled elements
        !               - reset the mesh
        !               - determine the elements that particles moved into
        !               - determine the new local particle coordinates
        !               - smoothen particle stresses within each element
        !               - update particle shape values
        !               - update the particle house-keeping data structure
        !               - ... further checks ...
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
         real(REAL_TYPE), dimension(Counters%N) :: NodalIncDisplacementMeshAdjust
         ! Local variables
         integer(INTEGER_TYPE) :: iOpt, I, EntityUsed

         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: Momentum
         !real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity) :: Momentum_VolLockSmooth
         
         
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: MomentumW
         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: MomentumG

         real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: NodalIncDisplacementG

         real(REAL_TYPE), dimension(Counters%N) :: DummyRotationVector1

         integer(INTEGER_TYPE) :: IEntity, StructureMaterialEntity, LoadedSides, IDOF
         
         ! Start: Bending elastic beam example ------------------------------------------------------------------------
         ! beam variables
          real(REAL_TYPE) :: LL
          real(REAL_TYPE) :: HH
          real(REAL_TYPE) :: Rho
          
          real(REAL_TYPE) :: kk
          
          real(REAL_TYPE) :: Y_max
          
          real(REAL_TYPE) :: C_1
          real(REAL_TYPE) :: C_2
          
          real(REAL_TYPE) :: ww
          
          real(REAL_TYPE) :: AA
          !integer(INTEGER_TYPE) :: ii
          real(REAL_TYPE) :: YoungsModulus !-> not the constitutive parameter
          real(REAL_TYPE) :: PoissonRatio !-> not the constitutive parameter
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          
          integer(INTEGER_TYPE) :: ii
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IControlPoint_Other
          integer(INTEGER_TYPE) :: IControlPoint_MovingMesh
          integer(INTEGER_TYPE) :: GlobalID
          
          
          ! Start: Bending elastic beam example ------------------------------------------------------------------------
          EntityUsed = 1   !!CC - each node has nEntity displacements, choose entity used for nodal update
          
          call GetRigidBodyAverageAcceleration()
          
         ! Start: Bending elastic beam example ------------------------------------------------------------------------
          ! ELASTIC BEAM PROBLEM - INITIAL VELOCITY
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !if ( (CalParams%TimeStep == 1) .and. (CalParams%IStep ==1) ) then 
          !      ! Bending of elastic beam 
          !      LL = 0.06 !m !length of beam 
          !      HH = 0.01 !m !height of beam 
          !      YoungsModulus = 318e6 !GPa ! need to get this from the GOM file 
          !      PoissonRatio = 0.054 !almost zero because 1D ! need to get this from the GOM file
          !      Rho = 1845 !kg/m3 !density of beam 
          !
          !      kk = 78.83 ! m^-1 
          !
          !      Y_max = 0.005 !m
          !
          !      C_1 = cosh(kk*LL) - cos(kk*LL)
          !      C_2 = sinh(kk*LL) - sin(kk*LL)
          !      ! 1e3 because YoungsModulus is in kPa
          !      ww = (kk**2) * ( (YoungsModulus*1e3*HH**2)/(12 * Rho * (1-PoissonRatio**2)) )**0.5 ! young's modulus need to be in Pa so put in the 1000 factor in 
          !
          !      AA = Y_max/(2*C_2)
          !
          !
          !      ! loop over material points 
          !      do ii = 1, Counters%NParticles
          !    
          !          !GlobPosArray 
          !          ! Y-coordinate 
          !          VelocityArray(ii,2) = AA * ww * &
          !    ( C_1 * (sinh(kk*(GlobPosArray(ii,1)-0.01)) + sin(kk*(GlobPosArray(ii,1)-0.01)) ) - C_2 * (cosh(kk*(GlobPosArray(ii,1)-0.01)) + cos(kk*(GlobPosArray(ii,1)-0.01))  )  )
          !    
          !          ! X-coordinate
          !          VelocityArray(ii,1) = 0.0
          !      end do 
          !
          !end if 
          ! End: Bending elastic beam example ------------------------------------------------------------------------
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Start: Shallow foundation example -------------------------------------------------------------------------------
          
          ! Update particles total velocities and Accelerations (accelerations added by CC)    
         
          ! 64 by 64 foundation 
            ! x-direction
            !VelocityArray(15361 , 1) = 0.0
            !VelocityArray(15362 , 1) = 0.0 
            !VelocityArray(15363 , 1) = 0.0
            !VelocityArray(15364 , 1) = 0.0 
            !VelocityArray(15365 , 1) = 0.0 
            !VelocityArray(15366 , 1) = 0.0
            !VelocityArray(15367 , 1) = 0.0
            !VelocityArray(15368 , 1) = 0.0
            !VelocityArray(15369 , 1) = 0.0
            !VelocityArray(15370 , 1) = 0.0
            !VelocityArray(15371 , 1) = 0.0
            !VelocityArray(15372 , 1) = 0.0
            !VelocityArray(15373 , 1) = 0.0
            !VelocityArray(15374 , 1) = 0.0
            !VelocityArray(15375 , 1) = 0.0
            !VelocityArray(15376 , 1) = 0.0
            !
            !
            !
            !! y-direction
            !VelocityArray(15361 , 2) = -1
            !VelocityArray(15362 , 2) = -1
            !VelocityArray(15363 , 2) = -1
            !VelocityArray(15364 , 2) = -1
            !VelocityArray(15365 , 2) = -1
            !VelocityArray(15366 , 2) = -1
            !VelocityArray(15367 , 2) = -1
            !VelocityArray(15368 , 2) = -1
            !VelocityArray(15369 , 2) = -1
            !VelocityArray(15370 , 2) = -1
            !VelocityArray(15371 , 2) = -1
            !VelocityArray(15372 , 2) = -1
            !VelocityArray(15373 , 2) = -1
            !VelocityArray(15374 , 2) = -1
            !VelocityArray(15375 , 2) = -1
            !VelocityArray(15376 , 2) = -1
          
        !VelocityArray(15873, 2) = -1     
        !VelocityArray(15874, 2) = -1
        !VelocityArray(15875, 2) = -1
        !VelocityArray(15876, 2) = -1
        !VelocityArray(15877, 2) = -1
        !VelocityArray(15878, 2) = -1
        !VelocityArray(15879, 2) = -1
        !VelocityArray(15880, 2) = -1
        !VelocityArray(15881, 2) = -1
        !VelocityArray(15882, 2) = -1
        !VelocityArray(15883, 2) = -1
        !VelocityArray(15884, 2) = -1
        !VelocityArray(15885, 2) = -1
        !VelocityArray(15886, 2) = -1
        !VelocityArray(15887, 2) = -1
        !VelocityArray(15888, 2) = -1
        !VelocityArray(15889, 2) = -1
        !VelocityArray(15890, 2) = -1
        !VelocityArray(15891, 2) = -1
        !VelocityArray(15892, 2) = -1
        !VelocityArray(15893, 2) = -1
        !VelocityArray(15894, 2) = -1
        !VelocityArray(15895, 2) = -1
        !VelocityArray(15896, 2) = -1
        !
        !VelocityArray(15873, 1) = 0  
        !VelocityArray(15874, 1) = 0
        !VelocityArray(15875, 1) = 0
        !VelocityArray(15876, 1) = 0
        !VelocityArray(15877, 1) = 0
        !VelocityArray(15878, 1) = 0
        !VelocityArray(15879, 1) = 0
        !VelocityArray(15880, 1) = 0
        !VelocityArray(15881, 1) = 0
        !VelocityArray(15882, 1) = 0
        !VelocityArray(15883, 1) = 0
        !VelocityArray(15884, 1) = 0
        !VelocityArray(15885, 1) = 0
        !VelocityArray(15886, 1) = 0
        !VelocityArray(15887, 1) = 0
        !VelocityArray(15888, 1) = 0
        !VelocityArray(15889, 1) = 0
        !VelocityArray(15890, 1) = 0
        !VelocityArray(15891, 1) = 0
        !VelocityArray(15892, 1) = 0
        !VelocityArray(15893, 1) = 0
        !VelocityArray(15894, 1) = 0
        !VelocityArray(15895, 1) = 0
        !VelocityArray(15896, 1) = 0

          
          ! x-direction
            !VelocityArray(769 , 1) = 0.0
            !VelocityArray(770 , 1) = 0.0 
            !VelocityArray(771 , 1) = 0.0
            !VelocityArray(772 , 1) = 0.0 
            !VelocityArray(773 , 1) = 0.0 
            !VelocityArray(774 , 1) = 0.0
            !VelocityArray(775 , 1) = 0.0
            !VelocityArray(776 , 1) = 0.0
            !VelocityArray(777 , 1) = 0.0
            !VelocityArray(778 , 1) = 0.0
            !VelocityArray(779 , 1) = 0.0
            !VelocityArray(780 , 1) = 0.0
            !VelocityArray(781 , 1) = 0.0
            !VelocityArray(782 , 1) = 0.0
            !VelocityArray(783 , 1) = 0.0
            !VelocityArray(784 , 1) = 0.0
            !VelocityArray(785 , 1) = 0.0
            !VelocityArray(786 , 1) = 0.0
            !VelocityArray(787 , 1) = 0.0
            !VelocityArray(788 , 1) = 0.0
            !VelocityArray(789 , 1) = 0.0
            !VelocityArray(790 , 1) = 0.0
            !VelocityArray(791 , 1) = 0.0
            !VelocityArray(792 , 1) = 0.0
            !VelocityArray(793 , 1) = 0.0
            !VelocityArray(794 , 1) = 0.0
            !VelocityArray(795 , 1) = 0.0
            !VelocityArray(796 , 1) = 0.0
            !VelocityArray(797 , 1) = 0.0
            !VelocityArray(798 , 1) = 0.0
            !VelocityArray(799 , 1) = 0.0
            !VelocityArray(800 , 1) = 0.0
            !
            !
            !! y-direction
            !VelocityArray(769 , 2) = -5
            !VelocityArray(770 , 2) = -5
            !VelocityArray(771 , 2) = -5
            !VelocityArray(772 , 2) = -5
            !VelocityArray(773 , 2) = -5
            !VelocityArray(774 , 2) = -5
            !VelocityArray(775 , 2) = -5
            !VelocityArray(776 , 2) = -5
            !VelocityArray(777 , 2) = -5
            !VelocityArray(778 , 2) = -5
            !VelocityArray(779 , 2) = -5
            !VelocityArray(780 , 2) = -5
            !VelocityArray(781 , 2) = -5
            !VelocityArray(782 , 2) = -5
            !VelocityArray(783 , 2) = -5
            !VelocityArray(784 , 2) = -5
            !VelocityArray(785 , 2) = -5
            !VelocityArray(786 , 2) = -5
            !VelocityArray(787 , 2) = -5
            !VelocityArray(788 , 2) = -5
            !VelocityArray(789 , 2) = -5
            !VelocityArray(790 , 2) = -5
            !VelocityArray(791 , 2) = -5
            !VelocityArray(792 , 2) = -5
            !VelocityArray(793 , 2) = -5
            !VelocityArray(794 , 2) = -5
            !VelocityArray(795 , 2) = -5
            !VelocityArray(796 , 2) = -5
            !VelocityArray(797 , 2) = -5
            !VelocityArray(798 , 2) = -5
            !VelocityArray(799 , 2) = -5
            !VelocityArray(800 , 2) = -5
         
         
          ! End: Shallow foundation example -------------------------------------------------------------------------------
         
          ! START: Two balls colliding example 
          !if ( (CalParams%TimeStep == 1) .and. (CalParams%IStep ==1) ) then
          !
          !    do ii=1, 4051
          !        
          !        VelocityArray(ii , 1) = 0.10!+4.0 !2.0
          !        VelocityArray(ii , 2) = 0.10!+4.0 !2.0
          !        
          !    end do 
          !    
          !    
          !    do ii=4052, 8102
          !        
          !        VelocityArray(ii , 1) = -0.10!-4.0 !-2.0
          !        VelocityArray(ii , 2) = -0.10!-4.0 !-2.0
          !        
          !    end do 
          !    
          !    
          !    
          !    
          !end if 
          ! END: Two balls colliding example 

          !--------------------------------------------------------------------------------------------------
          ! START: 
          
          
          
          
         ! Applying prescribed velocity
         call ApplyMPPrescribedVelocity_BoundLocus(TotalVelocitySys) 
         call ApplyMPPrescribedVelocity(TotalVelocitySys)
         call ApplyNodalPrescribedVelocity(DummyRotationVector1)

         if (IsMPMComputation()) then ! if MPM update particle velocity....
           call UpdateParticleVelocityAndMapMomentum(Momentum)!, Momentum_VolLockSmooth) ! From accelerations in global coordinate system
           
           
           if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3).or.(.not.(NFORMULATION==1))) then
             call UpdateParticleWaterVelocityAndMapMomentumW(MomentumW) ! From accelerations in global coordinate system
           end if
           if (CalParams%NumberOfPhases==3) then
             call UpdateParticleGasVelocityAndMapMomentumG(MomentumG)
           end if
         else !FEM  get velocity directly at the nodes
           ! Accelerations are in global coordinate system but TotalVelocitySoil, too?
           ! TotalVelocitySoil is global if contact formulation used and local without contact
           ! TotalVelocityWater is local with and without contact
           ! Rotate TotalVelocitySoil in case of no contact to global coordinate system
           ! Rotate TotalVelocityWater to global coordinate system
           if (IS3DCYLINDRIC) then
             do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityWater(:, IEntity), TotalVelocityWater(:, IEntity))
               if (.not.CalParams%ApplyContactAlgorithm.and.(.not.IsMPMSkipConvection())) then
                 call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocitySoil(:, IEntity), TotalVelocitySoil(:, IEntity))
               endif
             end do
           end if

           do IEntity = 1, Counters%nEntity
             do I = 1, Counters%N
                 
                 ! FEM and Generalized Alpha Scheme 
                 if ( (CalParams%ApplyGeneralizedAlphaTimeScheme%Switch == .TRUE.) .AND. (CalParams%TimeStep /= 1) ) then 
                     
                     TotalVelocitySoil(I, IEntity) = TotalVelocitySoil(I, IEntity) + &
                         ( CalParams%TimeIncrement * (CalParams%ApplyGeneralizedAlphaTimeScheme%Gamma*AccelerationSoil(I, IEntity) + &
                         (1-CalParams%ApplyGeneralizedAlphaTimeScheme%Gamma)*AccelerationSoilGeneralizedAlpha(I, IEntity)) ) * PBoundary(I)
                 
                 else 
                 
                     TotalVelocitySoil(I, IEntity) = TotalVelocitySoil(I, IEntity) + AccelerationSoil(I, IEntity) * CalParams%TimeIncrement  ! solid
                     
                 end if 
                 
               if (CalParams%NumberOfPhases==3) then
                 TotalVelocityGas (I, IEntity)= TotalVelocityGas(I, IEntity)+ AccelerationGas(I, IEntity) * CalParams%TimeIncrement
               end if
             end do
           end do

         end if
         
         
      
          ! Rotate Momentum from global to local coordinate system
          if (IS3DCYLINDRIC) then
            do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, Momentum(:, IEntity), Momentum(:, IEntity))
            end do
          end if

           if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
            ! Rotate water momentum from global to local coordinate system
            if (IS3DCYLINDRIC) then
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, MomentumW(:, IEntity), MomentumW(:, IEntity))
              end do
            end if ! rotation
          end if ! consolidation
      
          if (CalParams%NumberOfPhases==3) then
            if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity 
                call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, MomentumG(:, IEntity), MomentumG(:, IEntity))
              end do
            end if ! rotation
          end if ! unsat consolidation

         ! DoSystem = .false.  !not needed to map system velocity here
         ! Calculate new nodal velocities from the new momentums
         if (IsMPMComputation()) then !MPM
           call GetNodalVelocityFromNodalMomentumConv(Momentum)!, Momentum_VolLockSmooth) ! In local coordinate system
           if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3).or.(.not.(NFORMULATION==1))) then
             call GetNodalWaterVelocityFromNodalWaterMomentum(MomentumW) ! In local coordinate system
             
             !TotalVelocityWater(40,1) = TotalVelocityWater(2,1)
             !TotalVelocityWater(42,1) = TotalVelocityWater(2,1)
             !TotalVelocityWater(44,1) = TotalVelocityWater(2,1)
             !TotalVelocityWater(46,1) = TotalVelocityWater(2,1)
             !TotalVelocityWater(48,1) = TotalVelocityWater(2,1)
             
             if (CalParams%NumberOfPhases==3) then
               call GetNodalGasVelocityFromNodalGasMomentum(MomentumG)
             end if
           end if
         end if      
         
         
         ! I repeated this here to make sure the body behaves as a rigid body
         ! Applying prescribed velocity
         call ApplyMPPrescribedVelocity_BoundLocus(TotalVelocitySys) 
         call ApplyMPPrescribedVelocity(TotalVelocitySys)
         call ApplyNodalPrescribedVelocity(DummyRotationVector1)
         
         

          ! Rotate Momentum and TotalVelocitySoil from local to global coordinate system (but not for FEM which is already global)
          if (IS3DCYLINDRIC) then
            do IEntity = 1, Counters%nEntity
              call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, Momentum(:, IEntity), Momentum(:, IEntity))
              call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocitySoil(:, IEntity), TotalVelocitySoil(:, IEntity))
            end do
          end if

           ! Rotate momentum and velocity from local to global coordinate system
           if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
             if (IS3DCYLINDRIC) then
               do IEntity = 1, Counters%nEntity
                 call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, MomentumW(:, IEntity), MomentumW(:, IEntity))
                 if (IsMPMComputation()) then
                   call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityWater(:, IEntity), TotalVelocityWater(:, IEntity))
                 end if
               end do
             end if
           end if ! consolidation
          
          ! Rotate momentum and velocity from local to global coordinate system
          if (CalParams%NumberOfPhases==3) then ! rotation is needed
            if (IS3DCYLINDRIC) then ! rotation is needed
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, MomentumG(:, IEntity), MomentumG(:, IEntity))
                call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityGas(:, IEntity), TotalVelocityGas(:, IEntity))
              end do
            end if
          end if ! Unsat calculation 
          
          ! Overwrite nodal velocity and particle velocity if prescribed
           call ApplyMPPrescribedVelocity(TotalVelocitySys)
           call ApplyNodalPrescribedVelocity(DummyRotationVector1)    
        
           if (CalParams%BoundaryConditions%ApplyInfiltrationRate.or.CalParams%BoundaryConditions%ApplySeepageFace) then  
              call ApplyNodalInfiltrationRate(TotalVelocityWater,TotalVelocitySoil,LumpedNodalPorosityDegSat,NodalUnitMassGradient)
           end if
          
           ! incremental displacement soil for the higher order grid 
          call GetNodalIncrementalDisplacement(IncrementalDisplacementSoil, TotalVelocitySoil, AccelerationSoilGeneralizedAlpha, AccelerationSoil)
          !if (CalParams%ApplyBBar) then !--> volumetric locking mitigation
          !    call GetNodalIncrementalDisplacement_VolLockSmooth(IncrementalDisplacementSoil_VolLockSmooth, TotalVelocitySoil_VolLockSmooth)
          !end if   
          
          if (IsMPMSkipConvection()) then
            AccumulatedIncDisplacementSoil = AccumulatedIncDisplacementSoil + IncrementalDisplacementSoil
          end if

          if ((.not.(NFORMULATION==1))) then
            call GetNodalIncrementalDisplacement(IncrementalDisplacementWater, TotalVelocityWater, AccelerationWaterGeneralizedAlpha, AccelerationWater) ! Global coordinate system
          
            if (.not.IsMPMComputation()) then !FEM
             do I = 1, Counters%N
               TotalDisplacementWater(I)  = TotalDisplacementWater(I) + IncrementalDisplacementWater(I,1) ! Global coordinate system
             end do 
            end if
          end if

          if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
            call GetNodalIncrementalDisplacement(IncrementalDisplacementWater, TotalVelocityWater, AccelerationWaterGeneralizedAlpha, AccelerationWater) ! Global coordinate system
     
            if ( CalParams%ApplyAbsorbingBoundary.and.IsMPMComputation()) then
              !MPM ...need to get the accumulated displacement at the nodes
              call GetNodalAccumulatedDisplacementsWater(IncrementalDisplacementWater) ! Global coordinate system
            end if

            if ( CalParams%ApplyAbsorbingBoundary) then
              call CalcParticleDisplacementsWater()
              call UpdateParticleWaterAcceleration()
            end if

            if (.not.IsMPMComputation()) then !FEM
              do I = 1, Counters%N
                TotalDisplacementWater(I)  = TotalDisplacementWater(I) + IncrementalDisplacementWater(I,1) ! Global coordinate system
             end do 
            end if
          
          end if ! consolidation
      
          if (CalParams%NumberOfPhases==3) then
            call GetNodalIncrementalDisplacement(NodalIncDisplacementG, TotalVelocityGas, AccelerationGasGeneralizedAlpha, AccelerationGas)
     
            if ( CalParams%ApplyAbsorbingBoundary.and. IsMPMComputation()) then !MPM ...need to get the accumulated displacement at the nodes
              call GetNodalAccumulatedDisplacementsGas(NodalIncDisplacementG)
            end if 
                
            if ( CalParams%ApplyAbsorbingBoundary) then 
              call GetIncrementalDisplacementGas(NodalIncDisplacementG)
            end if 
 
            if (.not.IsMPMComputation()) then !FEM
              do I = 1, Counters%N
                TotalDisplacementGas(I)  = TotalDisplacementGas(I) + NodalIncDisplacementG(I,1)
              end do 
            end if
          
          end if ! Unsat Calculation
      
          ! Update the nodal coordinates in case of updated mesh analysis
          if (IsULFEMComputation()) then
            call UpdateNodes(EntityUsed)
            do IPatch = 1, Counters%NPatches ! loop over patches
            do I = 1, Counters%NEl(IPatch)!Counters%NEl
              call CoordLocalToGlobal(I, NodalCoordinatesUpd(:,:,IPatch), IPatch )
            end do ! elements
            end do ! patches
          end if
         
         !call ComputeStructureMoments()
         
         call UpdateNodalTotalDisplacement(EntityUsed)
         !call UpdateNodalTotalDisplacement_VolLockSmooth(EntityUsed)
         
         ! Determine nodal density field
          if(NFORMULATION==1) then ! 1 Constituent
              !call DetermineLiquidDensityField() !--> I commented this for now
          else ! 2 Constituents
              if(CalParams%NumberOfPhases==1) then ! 1 Phase
                call DetermineLiquidDensityField() !TEMPORARY
              end if
          end if
    
        if (.not.CalParams%ApplyStrainSmoothing) then 
           ! Update particle strain data
            call UpdateParticleStrains() ! update particle strains for solid and mixture
            
            
            if (CalParams%ApplyFBar .and. CalParams%ApplyFBarSmoothStrains) then ! only FBar smooth strains turned on
                
                ! Smooth the strain field 
                !--------------------------------------------------------------
                if (CalParams%ApplyFBarLowerOrder) then 
                    
                ! use lower order grid    
                call MapStrainIncrementTensor_VolLockSmooth_LowerOrder(WeightedStrainIncrementTensor_Numerator)
                call GetNodalStrainIncrement_VolLockSmooth_LowerOrder()
                
                call MapNodalStrainIncrementToMPs_VolLockSmooth_LowerOrder(NodalStrainIncrement_VolLockSmooth)
                
                else 
                    
                ! use same higher order grid
                ! we map strain increment * mass from the material point to the nodes 
                call MapStrainIncrementTensor_VolLockSmooth(WeightedStrainIncrementTensor_Numerator)
                ! we divide by the lumped mass to obtain the nodal strain increment
                call GetNodalStrainIncrement_VolLockSmooth()
                ! we map the nodal strain increment to the material points
                call MapNodalStrainIncrementToMPs_VolLockSmooth(NodalStrainIncrement_VolLockSmooth)
                
                end if
                !--------------------------------------------------------------
            end if     
        end if
        
        
        
        if ((.not.CalParams%ApplyStrainSmoothingLiquidTwoLayer)) then 
            call UpdateParticleStrainsLiquidTwoLayer() ! update particle strains for liquid MPs in TwoLayerForm
            if ((.not.(CalParams%NumberOfPhases==1)).and.(.not.(NFORMULATION==1))) then
              call ComputeTwoLayerVolumetricStrainLiquid(IncrementalDisplacementWater, IncrementalDisplacementSoil)
            end if
        end if
          
          
         if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1)) then
           call CalculateWaterVolumetricStrain(IncrementalDisplacementWater)
             if (IsMPMComputation()) then  ! MPM
               !update particle displacement (water)
               call UpdateParticleDisplacementsWater(IncrementalDisplacementWater)
            else !(.not.IsMPMComputation()) then !FEM
              call UpdateParticleVelocityWater(TotalVelocityWater)
             end if
         end if    
         
         if (CalParams%NumberOfPhases==3) then
            call CalculateGasVolumetricStrain(NodalIncDisplacementG)
            if (IsMPMComputation()) then  ! MPM
              !update particle displacement (gas)
              call UpdateParticleDisplacementsGas(NodalIncDisplacementG)
            else ! (.not.IsMPMComputation()) then !FEM
              call UpdateParticleVelocityGas(TotalVelocityGas)
           end if
         end if

         if ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation).and.(NFORMULATION==1)) then
            call CalculateParticleWaterAdvectiveFlux(IncrementalDisplacementWater)
         else if ((CalParams%NumberOfPhases==3).and.(NFORMULATION==1)) then
            call CalculateParticleWaterAdvectiveFlux(IncrementalDisplacementWater)
            call CalculateParticleAirAdvectiveFlux(NodalIncDisplacementG)
         end if

         iOpt = 0

         !Detect Liquid Free Surface
         call DetectLiquidFreeSurface();
         
         if(.not.(NFORMULATION==1)) then
            call TwoLayerData%DetermineTwoLayerStatus() 
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The balance equations are solved--> Incremental Water Pressure, Gas Pressure and Temperature are calculated, 
!and the new stresses are obtained applying the proper Constitutive equation
!Global coordinate system 
        ! The stresses, the water pressure and the gass pressure are updated (The mass balance equations are solved)
        call MPMDYNGetSig() ! TODO-VL: this would need to go last for USL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! Pore pressure smoothing 
        !if (CalParams%ApplyFBar) then  !.and. CalParams%ApplyFBarSmoothStresses) then !CalParams%ApplyFBarSmoothStrains
                
                ! Smooth the stress field 
                !--------------------------------------------------------------
                if (CalParams%ApplyFBarLowerOrder) then 
                    
                    if (CalParams%ApplyFBarSmoothStresses) then 
               
                        call MapStressIncrementTensor_VolLockSmooth_LowerOrder(WeightedStressIncrementTensor_Numerator)
                        call GetNodalStressIncrement_VolLockSmooth_LowerOrder()
                        call MapNodalStressIncrementToMPs_VolLockSmooth_LowerOrder(NodalStressIncrement_VolLockSmooth)
                
                    end if 
               
                    if ( (CalParams%NumberOfPhases==2) .and. CalParams%ApplyFBarSmoothWP) then 
                    
                        call MapWP_VolLockSmooth_LowerOrder(WeightedWP_Numerator)
                        call GetNodalWP_VolLockSmooth_LowerOrder()
                        call MapNodalWPToMPs_VolLockSmooth_LowerOrder(NodalWP_VolLockSmooth)
                    
               
                    end if 
                
                
                else 
                   
                    
                    if (CalParams%ApplyFBarSmoothStresses) then 
                
                        call MapStressIncrementTensor_VolLockSmooth(WeightedStressIncrementTensor_Numerator)
                        call GetNodalStressIncrement_VolLockSmooth()
                        call MapNodalStressIncrementToMPs_VolLockSmooth(NodalStressIncrement_VolLockSmooth)
               
                    end if 
               
                    if ( (CalParams%NumberOfPhases==2) .and. CalParams%ApplyFBarSmoothWP) then 
                    
                        call MapWP_VolLockSmooth(WeightedWP_Numerator)
                        call GetNodalWP_VolLockSmooth()
                        call MapNodalWPToMPs_VolLockSmooth(NodalWP_VolLockSmooth)
                
                    end if 
                
                end if
                !--------------------------------------------------------------
                
            !end if 
        
        
        
        
        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      ! Update particle stress data for fully filled elements
       if (IsMPMWithMixedIntegration() .and. .not.IsMPMSkipConvection()) then ! Not needed for pure MPM and FEM
         call AssignStressesToParticles() ! Need to check for parallelization
         call AssignStateParametersToParticles()  
       end if

       !if (CalParams%ApplyFBar) then ! if apply B-bar volumetric locking smoothing technique
           
           ! map stress from the material points to the nodes
           !call MapStressToNodes_VolLockSmooth(NodalStress_VolLockSmooth)
           !call UpdateNodalStressGradient()
           
           ! Smooth the stress field      
           !--------------------------------------------------------------
           !call MapStressTensor_VolLockSmooth(WeightedStressTensor_Numerator)
           !call GetNodalStress_VolLockSmooth()
           !call MapNodalStressToMPs_VolLockSmooth(NodalStress_VolLockSmooth)
           
           !--------------------------------------------------------------
           
           
           
           !call MapVolumetricStressGradient_VolLockSmooth(WeightedNodalVolumetricStressGradient_Numerator)
           !
           !call GetNodalVolumetricStressGradient_VolLockSmooth()     
           !call MapVolumetricStressGradientToMPs_VolLockSmooth(NodalVolumetricStressGradient_VolLockSmooth)
           ! 
           !call UpdateStress_VolLockSmooth()
           
           
           ! TODO: 
           
           
           
           
           !call MapStressToNodes() 
           !!call Get_StressGradient()    
           !call UpdateNodalStressGradient()
           !call FBar_UpdateMPStress()
       !end if 
       

!\\\\\\\ UPDATE PARTICLE PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      

       ! Update particle weights
      if (IsMPMComputation()) then
       call DynUpdateParticleWeights( )
       ! TODO-VL: this is where we would need to use the updated deformation gradient 
       ! to find the new particle volume. We need to store the initial volume because  
       ! the deformation gradient connects the initial 
       ! V_n+1_MP = V_0_MP * F_n_MP
      end if

      !Update porosity (Directly depends on: DEpsVol)
      if (CalParams%ApplyPorosityUpdate) then !Update porosity (Directly depends on: DEpsVol)
        call DynUpdateParticlePorosity( )
      end if
      
      ! Update permeability (depends on porosity)
      if (CalParams%ApplyDarcyPermeabilityUpdate) then 
          call DynUpdateParticleConductivity( )
      end if 

!\\\\\\\ END UPDATE PARTICLE PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      

       ! VAHID: This does not seem to be true! isMPMSkip..() must be swaped with %SkipConvection
       if (.not.IsMPMSkipConvection()) then
           
         if (CalParams%ApplyMeshSmoothing) then 
             
           NodalIncDisplacementMeshAdjust = 0.0 ! initialize
           
           if (.not.CalParams%SkipConvection) then
            
               if (CalParams%ApplyContactAlgorithm) then ! contact...use pile entity
                 
                   StructureMaterialEntity = getMaterialEntity(CalParams%MovingMesh%MovingMaterialID(IPatch_Temporary)) ! needs generalization 
               
                 do I = 1, Counters%N
                     NodalIncDisplacementMeshAdjust(I) = IncrementalDisplacementSoil(I,StructureMaterialEntity) !2 is the entity of the pile
                 end do
             
               else
               
                   ! no contact it is just one entity
                   do I = 1, Counters%N
                       NodalIncDisplacementMeshAdjust(I) = IncrementalDisplacementSoil(I,1)
                   end do
                   
               end if
               
           else
               
             do I = 1, Counters%N
               NodalIncDisplacementMeshAdjust(I) = IncrementalDisplacementSoil(I,1)
             end do
           
           end if

           ! Update the nodes for strain calculation 
           call UpdateNodesMeshAdjust(NodalIncDisplacementMeshAdjust)
           ! Update the nodes for strain calculation --> lower order for volumetric locking smoothing 
           !call UpdateNodesMeshAdjust_VolLockSmooth(IncrementalDisplacementSoil_VolLockSmooth(:,1)) !NodalIncDisplacementMeshAdjust)
            
           ! Perform mesh adjustment by modifying the initial nodal coordinates
           !call MeshAdjustment() ! if moving mesh ------------------------------------------
           ! This is currently only working for compressing mesh --> need to generalize for both compresisng and extending mesh
           call MeshAdjustment_MovingMesh_Multipatch(MovingMeshPatchID, CompressingMeshPatchID) ! Hardcoded patch number for now

           
           
           
           
           
           ! we need to loop the patches that do not belong to either of the moving nor the compressing mesh
           do IPatch = 1, Counters%NPatches
           
               ! 1- loop over every patch that is neither moving or compressing
               if ( (IPatch==MovingMeshPatchID) .or. (IPatch==CompressingMeshPatchID) ) cycle
              
               
               ! 2- loop over control points of "other" patch
               do IControlPoint_Other = 1, Counters%NodTot(IPatch) 
                   
                   GlobalID = Multipatch_Connecting_Local_To_Global_ControlPoints(IControlPoint_Other, IPatch)
                   
               
                   ! 3- loop over the control points of moving mesh 
                   do IControlPoint_MovingMesh = 1, Counters%NodTot(MovingMeshPatchID) 
                       
                   
                       if (Multipatch_Connecting_Local_To_Global_ControlPoints(IControlPoint_MovingMesh, MovingMeshPatchID) == GlobalID) then 
                           ! update the nodes
                           NodalCoordinates(IControlPoint_Other,:,IPatch) = NodalCoordinates(IControlPoint_MovingMesh,:,MovingMeshPatchID)
                       end if
                       
                   end do 
                   
                   
                   
                   ! 3- loop over the control points of compressing mesh 
                   do IControlPoint_MovingMesh = 1, Counters%NodTot(CompressingMeshPatchID) 
                       
                   
                       if (Multipatch_Connecting_Local_To_Global_ControlPoints(IControlPoint_MovingMesh, CompressingMeshPatchID) == GlobalID) then 
                           ! update the nodes
                           NodalCoordinates(IControlPoint_Other,:,IPatch) = NodalCoordinates(IControlPoint_MovingMesh,:,CompressingMeshPatchID)
                       end if
                       
                   end do 
                   
                   
               end do 
               
           end do 
           
           
           if (IsULFEMComputation()) then
             call UpdateMeshAdjacencyInformation(NodalCoordinatesUpd)!(:,:,IPatch_Temporary))
           else
             call UpdateMeshAdjacencyInformation(NodalCoordinates)!(:,:,IPatch_Temporary))
           end if
         end if 

         ! Update particle global positions, particle-element assignment, particle local positions and shape function values
         ! NOTE: EleParticles and Particle%ElementID get already updated while
         ! the remaining house-keeping data is updated in the succeeding routine!!
         if (IsMPMComputation()) then ! MPM
           call UpdateParticlePos()
         end if

         if ( CalParams%ApplyMeshSmoothing .and. ISAXISYMMETRIC ) then
            ! updating global position of GPs is only needed in case of moving mesh
            !call SetGPGlobalPositionElement(GPGlobalPositionElement)
         end if
         
         ! Update particle house-keeping lists and element switches
         if (IsMPMComputation()) then ! MPM
            !call UpdateParticleHouseKeeping() ! HARDCODED
             !if (MP_INIT_FILE) then 
             ! call UpdateParticleHouseKeeping_MPInitFile()
             !else ! MP_INIT_GP
                 call UpdateParticleHouseKeeping_MPInitGP()
             !end if
         end if
         ! Reset mesh (nodal coordinates)
         if (.not.IsULFEMComputation()) then
           NodalCoordinatesUpd = NodalCoordinates
         end if

         call SetActiveElement() 
         !if (CalParams%MPInitialization%MPInitFormat == MP_INIT_FILE) then 
         !    call SetParticleIndex_MPInitFile()
         !else if (CalParams%MPInitialization%MPInitFormat == MP_INIT_GAUSSPT) then 
             call SetParticleIndex()
         !end if 
         
         !if (CalParams%MPInitialization%MPInitFormat == MP_INIT_FILE) then 
         !    call SetParticleIndex_MPInitFile()
         !else if (CalParams%MPInitialization%MPInitFormat == MP_INIT_GAUSSPT) then 
         !    call SetParticleIndex()
         !end if 
         
         
         if (CalParams%ApplyEmptyElements) then
           ! Check whether elements are deactivated that should be considered fully filled
           call CheckEmptyElements()
           call AdjustParticleDiscretisation()
         end if

         call SetActiveElement() !hardcoded
         !if (CalParams%MPInitialization%MPInitFormat == MP_INIT_FILE) then 
         !    call SetParticleIndex_MPInitFile()
         !else if (CalParams%MPInitialization%MPInitFormat == MP_INIT_GAUSSPT) then 
             call SetParticleIndex()
         !end if 

         call CheckFillingOfElements()

         call SetUpEntityElements() !hardcoded
         call SetUpMaterialElements()
       endif
       
       
       !-----------------------------------------------------------------------------------------------------------
       ! Pore pressure smoothing algorithm 
       
       ! MG22: step 1 is to find how many and which material points are located in the local subzone
       if (IsGPSmoothFluidPressure()) then 
       call NumberOfMaterialPointsInSubElement()
       call StressAndPorePressureSmoothening_4GaussPoints()
       call AssignStressesToParticles_4GaussPoints()
       end if 
       

       ! Smoothing of particle stresses for fully filled elements
       if (IsMPMWithMixedIntegration() .and. .not.IsMPMSkipConvection()) then ! Not needed for pure MPM and FEM
         call StressAndPorePressureSmoothening()
         call StateParametersSmoothening ()
       else
         ! Store initial stresses of next step
         call SetInitialStressForNextLoadStep()
       end if
      
   !UPDATE particle properties that depends on stresses (NB: done after stress and pore pressure smoothening)    
       if ((CalParams%NumberOfPhases==3).or. &
          ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation))) then
      
        call DynUpdateParticleDegreeOfSaturation( ) !Update Degree of Saturation (Retention Curve) (Directly depends on: Pg, Pl) (Undirectly depends on: Smin, Smax,T,...)
        call DynUpdateParticleMixedWeight()!Update Mixed Weight according to degree of saturation  
        call DynUpdateParticleHydraulicConductivity( )!Update hydraulic condictivity according to degree of saturation
       end if
          


       if ((.not.(NFORMULATION==1)).and.(.not.(CalParams%NumberOfPhases==1))) then
         call TwoLayerData%DetermineContainedMaterialTypes()
         call DetermineDensityField(TwoLayerData%Nodes%DensityLiquidL, MaterialPointTypeLiquid, .false., .true.)
         call DetermineDensityField(TwoLayerData%Nodes%DensitySolidL, MaterialPointTypeSolid, .false., .true.)
         call DetermineDensityField(TwoLayerData%Nodes%DensityLiquidS, MaterialPointTypeLiquid, .false., .false.)
         call DetermineDensityField(TwoLayerData%Nodes%DensitySolidS, MaterialPointTypeSolid, .false., .false.)


         call DetermineConcentrationRatio(TwoLayerData%Nodes%DensitySolidL, TwoLayerData%Nodes%DensityLiquidL, MaterialPointTypeLiquid,.true.)
         call DetermineConcentrationRatio(TwoLayerData%Nodes%DensitySolidL, TwoLayerData%Nodes%DensityLiquidL, MaterialPointTypeSolid,.true.)
         call DetermineConcentrationRatio(TwoLayerData%Nodes%DensitySolidS, TwoLayerData%Nodes%DensityLiquidS, MaterialPointTypeSolid,.false.)

         call TwoLayerData%DetermineConcentrationRatioElementFND()

         call DetermineBoundaryElementSolidDomain2LayerForm()
         call DetermineBoundaryElementLiquidDomain2LayerForm()

         call DetermineFillingRatioField()
         call DetermineFillingRatioLiquid(TwoLayerData)

       end if

       ! calculate system initial kinetic energy
       if (CalParams%ConvergenceCheck%KineticEnergy0 == -1.0) then
         if (IsMPMComputation()) then ! MPM
           if (NFORMULATION==1) then
             call CalculateKineticEnergy(TotalVelocitySoilPrevious, TotalVelocityWaterPrevious)
           else
             call CalculateKineticEnergy2LayForm(TotalVelocitySoilPrevious, TotalVelocityWaterPrevious)
           end if
         else ! FEM
           call CalculateKineticEnergyFEM(TotalVelocitySoilPrevious, TotalVelocityWaterPrevious)
         end if

         CalParams%ConvergenceCheck%KineticEnergy0 = CalParams%ConvergenceCheck%KineticEnergy
         CalParams%ConvergenceCheck%KineticEnergySoil0 = CalParams%ConvergenceCheck%KineticEnergySoil
         CalParams%ConvergenceCheck%KineticEnergyWater0 = CalParams%ConvergenceCheck%KineticEnergyWater
       end if

       ! calculate system kinetic energy
       if (.not.IsMPMSkipConvection()) then
         if (IsMPMComputation()) then
           ! MPM
           if (NFORMULATION==1) then
             call CalculateKineticEnergy(TotalVelocitySoil, TotalVelocityWater)
           else
             call CalculateKineticEnergy2LayForm(TotalVelocitySoil, TotalVelocityWater)
           end if
         else
           ! FEM
           call CalculateKineticEnergyFEM(TotalVelocitySoil, TotalVelocityWater)
         end if
         ! calculate system internal and external works
       endif

       if (NFORMULATION==1) then
         call CalculateIntAndExtWorks(IncrementalDisplacementWater)
       else
         call CalculateIntAndExtWorks2LayForm(IncrementalDisplacementWater)
       end if
       
       
       
       
       
       
       
      
       if ((CalParams%ApplyFixedSolidSkeleton == .false.).and.(CalParams%PrescribedHead%HydraulicHead == .true.)) then
           HydraulicHeadLoadedElemID = .false.
           call GetAndUpdateHydraulicHeadLoad(HydraulicHeadLoadedElemID, LoadedSides, HydraulicHeadNodesConnectivities, HydraulicHeadLoad)
           Counters%HydraulicHeadSides = LoadedSides
           HydraulicHeadVector = 0.0
           if (NDIM == 3) then
               do I = 1, Counters%HydraulicHeadSides
                   ! TEMPORARILY COMMENTED
                   !call Load3D(HydraulicHeadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, N_BOUNDARY_NODES_HOE, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LOADTYPE_LIQUID)
               end do
           elseif (NDIM == 2) then
               do I = 1, Counters%HydraulicHeadSides
                   call Load2D(HydraulicHeadVector, ReducedDof, NodalCoordinates(:,:,IPatch_Temporary), 1, HydraulicHeadNodesConnectivities(1, I), HydraulicHeadLoad(1, 1, I), LOADTYPE_LIQUID)
               end do
           end if
           HydraulicHeadLoadTotal = 0.0
           do IEntity=1,Counters%NEntity
               do IDof=1,Counters%N
                 HydraulicHeadLoadTotal(IDof, IEntity) = HydraulicHeadVector(IDof)
               end do
            end do
       end if
       
       if (CalParams%SkipConvection) then
         if (.not.IsMPMSkipConvection()) then
           AccumulatedIncDisplacementSoil = 0.0
           call setElementDeterminant(DShapeValuesArray)
           call setIsDistorted(.false., DShapeValuesArray)
         else
           call setIsDistorted(.true., DShapeValuesArray)
         endif
       end if
       
       
       ! calculate material point energy
       call CalculateMPEnergy()

       
      end subroutine DYNConvectivePhase
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      
      subroutine NumberOfMaterialPointsInSubElement()
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! - counting the number of material points in each subzone for 
      ! allocation purposes. 
      ! - recording the material point indices that correspond to each subzone. 
      ! subzone 1 --> -1<xi<0, -1<eta<0
      ! subzone 2 --> 0<xi<1, -1<eta<0
      ! subzone 3 --> -1<xi<0, 0<eta<1
      ! subzone 4 --> 0<xi<1, 0<eta<1
      !                         ^ Eta
      !  (-1,+1)        4       |   (+1,+1)
      !                +---------------+ 3
      !                |        |      |
      !                |   SB3  |  SB4 |
      !                |        |      |
      !                |---------------|---> Xi
      !                |        |      |
      !                |    SB1 |  SB2 |
      !                |1       |      | 2
      !                +---------------+-
      !             (-1,-1)         (+1,-1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none 
      
      ! variables corresponding to the number of MPs in each subzone
      integer(INTEGER_TYPE) :: NumberOfMPsInSubZone1, NumberOfMPsInSubZone2
      integer(INTEGER_TYPE) :: NumberOfMPsInSubZone3, NumberOfMPsInSubZone4
      
      ! maximum number of MPs in a subzone in an element. This is for allocation purposes. 
      integer(INTEGER_TYPE) :: MaximumNumberOfMPsInSubZone
      
      ! loop counters
      integer(INTEGER_TYPE) :: iElement, iMP, NElemPart, IEl, ParticleIndex
      
      integer(INTEGER_TYPE) :: NumberOfElements  
      
      real(REAL_TYPE), dimension(NDIM) :: LocPos ! MP local position
      
      ! subzone counters
      integer(INTEGER_TYPE) :: MPCounter_subzone1, MPCounter_subzone2
      integer(INTEGER_TYPE) :: MPCounter_subzone3, MPCounter_subzone4
      
      ! multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary
      
      ! temporary variable 
      IPatch_Temporary = 1
      
      ! use this subroutine if you are in 'MPM-MIXED_MG22-NO-INTERPOLATION' or 'MPM-MIXED_MG22'
      ! otherwise not needed
      
      !if (CalParams%ComputationMethod == MPM_MIXED_MG22_INTEGRATION .or. &
      !    CalParams%ComputationMethod == MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER) then 
      
      ! total number of elements
      NumberOfElements = Counters%Sum_NEl
      
      
      
      MaximumNumberOfMPsInSubZone = 0 ! --> for allocation purposes
      
      !allocate(SubElementMPOrganization(NumberOfElements, Counters%NGaussPoints)) ! --> global matrix
      SubElementMPOrganization = 0
      
      ! PART A: finding the maximum counts of the materisl points in the subzones. 
      ! loop over active elements
      do iElement = 1, Counters%NAEl(IPatch_Temporary)!Counters%NAEl!NumberOfActiveElement
          
            ! initialize 
          NumberOfMPsInSubZone1 = 0 ! -1<xi<0, -1<eta<0
          NumberOfMPsInSubZone2 = 0 ! 0<xi<1, -1<eta<0
          NumberOfMPsInSubZone3 = 0 ! -1<xi<0, 0<eta<1
          NumberOfMPsInSubZone4 = 0 ! 0<xi<1, 0<eta<1
          
          ! obtain global element index
          IEl = ActiveElement(iElement, IPatch_Temporary)       
          
          ! obtain the number of MPs in element IEl
          NElemPart = NPartEle(IEl, IPatch_Temporary) !NumberOfIntegrationPoints(IEl)    ! Get the number of integration per element (Gaussian integration points or number of material points)
           
          ! loop over the material points in the elements
          do iMP = 1, NElemPart!TotalNumberOfMaterialPoints
          
              ! obtain MP global index
              ParticleIndex = GetParticleIndex(iMP, IEl, IPatch_Temporary)
              ! obtain the local position of the MP
              LocPos = Particles(ParticleIndex)%LocPos
              
              ! if in subzone 1
              if (LocPos(1)<0 .and. LocPos(2)<0) then 
                  
                  ! increase the count by 1
                  NumberOfMPsInSubZone1 = NumberOfMPsInSubZone1 + 1
          
              ! if in subzone 2    
              else if (LocPos(1)>0 .and. LocPos(2)<0) then 
                  
                  ! increase the count by 1
                  NumberOfMPsInSubZone2 = NumberOfMPsInSubZone2 + 1

              ! if in subzone 3
              else if (LocPos(1)<0 .and. LocPos(2)>0) then 
                  
                  ! increase the count by 1
                  NumberOfMPsInSubZone3 = NumberOfMPsInSubZone3 + 1
                  
              ! subzone 4    
              else if (LocPos(1)>0 .and. LocPos(2)>0) then 
                  
                  ! increae the count by 1
                  NumberOfMPsInSubZone4 = NumberOfMPsInSubZone4 + 1

              end if ! end checking which subzone we're in
      
          end do ! total number of material points
          
          ! find the overall maximum of all the zones and for all the elements
          MaximumNumberOfMPsInSubZone = max(MaximumNumberOfMPsInSubZone, &
                                            NumberOfMPsInSubZone1, NumberOfMPsInSubZone2, &
                                            NumberOfMPsInSubZone3, NumberOfMPsInSubZone4)
          
          SubElementMPOrganization(IEl, 1) = NumberOfMPsInSubZone1
          SubElementMPOrganization(IEl, 2) = NumberOfMPsInSubZone2
          SubElementMPOrganization(IEl, 3) = NumberOfMPsInSubZone3
          SubElementMPOrganization(IEl, 4) = NumberOfMPsInSubZone4
      end do ! loop over elements
      
      
      ! allocate the global matrix GetParticleIndexInSubElement using the total number of elements, number of gauss points = 4, 
      ! and the maximum number of material points to a subzone. 
      if (allocated(GetParticleIndexInSubElement)) then 
          deallocate(GetParticleIndexInSubElement)
      end if 
      
      allocate(GetParticleIndexInSubElement(NumberOfElements, Counters%NGaussPoints, MaximumNumberOfMPsInSubZone))
      ! initialize to zero initially
      GetParticleIndexInSubElement = 0
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PART B: find the indices of the material points in each element in each subzone
      ! loop over all the active elements
      do iElement = 1, Counters%NAEl(IPatch_Temporary)!Counters%NAEl!NumberOfActiveElement

          ! initialize the counter to 1
          MPCounter_subzone1 = 1
          MPCounter_subzone2 = 1
          MPCounter_subzone3 = 1
          MPCounter_subzone4 = 1

          ! obtain global element index corresponding to that active element count
          IEl = ActiveElement(iElement, IPatch_Temporary)            
          ! obtain the total number of particles within the active element
          NElemPart = NPartEle(IEl, IPatch_Temporary) !NumberOfIntegrationPoints(IEl)    ! Get the number of integration per element (Gaussian integration points or number of material points)
           
          ! we want to know whether the material points in the element
          do iMP = 1, NElemPart!TotalNumberOfMaterialPoints
          
              ! obtain global particle index corresponding the material point iMP
              ParticleIndex = GetParticleIndex(iMP, IEl, IPatch_Temporary)
              ! obtain the local position of ParticleIndex
              LocPos = Particles(ParticleIndex)%LocPos
              
              ! subzone 1
              if (LocPos(1)<0 .and. LocPos(2)<0) then 
                  
                  !NumberOfMPsInSubZone1 = NumberOfMPsInSubZone1 + 1
                  GetParticleIndexInSubElement(IEl, 1, MPCounter_subzone1) = ParticleIndex ! subzone 1
                  MPCounter_subzone1 = MPCounter_subzone1 + 1
          
              ! subzone 2    
              else if (LocPos(1)>0 .and. LocPos(2)<0) then 
                  
                  !NumberOfMPsInSubZone2 = NumberOfMPsInSubZone2 + 1
                  GetParticleIndexInSubElement(IEl, 2, MPCounter_subzone2) = ParticleIndex ! subzone 2
                  MPCounter_subzone2 = MPCounter_subzone2 + 1
              
              ! subzone 3
              else if (LocPos(1)<0 .and. LocPos(2)>0) then 
                  
                  !NumberOfMPsInSubZone3 = NumberOfMPsInSubZone3 + 1
                  GetParticleIndexInSubElement(IEl, 3, MPCounter_subzone3) = ParticleIndex ! subzone 3
                  MPCounter_subzone3 = MPCounter_subzone3 + 1
                  
              ! subzone 4    
              else if (LocPos(1)>0 .and. LocPos(2)>0) then 
                  
                  !NumberOfMPsInSubZone4 = NumberOfMPsInSubZone4 + 1
                  GetParticleIndexInSubElement(IEl, 4, MPCounter_subzone4) = ParticleIndex ! subzone 4
                  MPCounter_subzone4 = MPCounter_subzone4 + 1

              end if ! checking where the material point lies 
      
          end do ! total number of material points
          
          
          
      end do ! loop over elements
      
      
      
      !end if ! if mixed MG22 method
      
      end subroutine NumberOfMaterialPointsInSubElement
      
      
      
      
      
      
      
      
      
      
      
      
      
      
       
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine CalculateMPEnergy() 
      !-------------------------------------------------
      ! Calculate kinetic energy and strain energy 
      !-------------------------------------------------
       implicit none 
       
       integer(INTEGER_TYPE) :: IParticle
       integer(INTEGER_TYPE) :: IDim
       integer(INTEGER_TYPE) :: ITens
       
       real(REAL_TYPE) :: VelocityMagnitude
       
       ! need to calculate the material point kinetic 
       KineticEnergyTotal = 0.0
       StrainEnergyTotal = 0.0
       
       KineticEnergyMP = 0.0
       StrainEnergyMP = 0.0
       
       do IParticle = 1, Counters%NParticles 
           
           VelocityMagnitude = 0.0
           
           do IDim = 1, NDIM 
               VelocityMagnitude = VelocityMagnitude + ( VelocityArray(IParticle,IDim) )**2.0
           end do 
           
           VelocityMagnitude = sqrt(VelocityMagnitude)
           KineticEnergyMP(IParticle) = KineticEnergyMP(IParticle) + (0.50 * MassArray(IParticle) * ((VelocityMagnitude)**2.0))
       
           KineticEnergyTotal = KineticEnergyTotal + KineticEnergyMP(IParticle)
           
           do ITens = 1, NTENSOR
               StrainEnergyMP(IParticle) = StrainEnergyMP(IParticle) + ( 0.50 * SigmaEffArray(IParticle, ITens) * Particles(IParticle)%Eps(ITens) * Particles(IParticle)%IntegrationWeight )
           end do 
           

           StrainEnergyTotal = StrainEnergyTotal + StrainEnergyMP(IParticle)

       end do
       
       
      end subroutine CalculateMPEnergy
      
      
       
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine ComputeForceError()
        !**********************************************************************
        !
        !    Function:  Computes the out-of-balance forces error, for the
        !               mixture and for the water phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
       implicit none

         ! Local variables
         integer(INTEGER_TYPE) :: I, J, K, INode, IDof
         real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
         real(REAL_TYPE) :: OutOfBalanceNorm, ExternalLoadNorm
         real(REAL_TYPE) :: OutOfBalanceNormWater, ExternalLoadNormWater
         logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity
         
         ! Multipatch variables 
         integer(INTEGER_TYPE) :: IPatch
         !integer(INTEGER_TYPE) :: IPatch_Temporary

         integer(INTEGER_TYPE) :: GlobalNodeID

         ! Compute ForceError of the mixture

          ExternalLoadDueToVelocity = 0.0
          TemporaryMappingVector = 0.0

           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

          if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity
            
              do IPatch = 1, Counters%NPatches ! loop over patches
              do INode = 1, Counters%NodTot(IPatch)!Counters%Nodtot ! loop over control points
                  
                  GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
                  
              do J = 1, NVECTOR
                IDof = (INode-1) * NDOFL + J      ! global Dof
                do K = 1, Counters%NEntity
                  ExternalLoadDueToVelocity(IDof) = ExternalLoadDueToVelocity(IDof) + FReaction(IDoF, K)
                end do
              end do
              end do ! loop over control points 
              end do ! loop over patches
         end if
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
             if((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then    ! consolidation
               if(CalParams%ApplyAbsorbingBoundary) then  ! absorbing boundary
                 TemporaryMappingVector(I) = TemporaryMappingVector(I)+  &
                                         GravityLoadMixture(I, J) * PBoundary(I) + &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         DummyIntLoadWaterPrevious(I, J) * PBoundaryWater(I) &
                                         - VisDampForceSld(I, J) &
                                         - VisDampForceWat(I, J)
               else
                 TemporaryMappingVector(I) = TemporaryMappingVector(I)+   & ! no absorbing boundary
                                         GravityLoadMixture(I, J) * PBoundary(I) + &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         DummyIntLoadWaterPrevious(I, J) * PBoundaryWater(I)
               end if
             else     ! no consolidation
               if(CalParams%ApplyAbsorbingBoundary) then ! absorbing boundary
                 TemporaryMappingVector(I) = TemporaryMappingVector(I)+  &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) &
                                         - VisDampForceSld(I, J)
               else
                 TemporaryMappingVector(I) = TemporaryMappingVector(I) + &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I)
               end if
             end if
           end do
         end do
         
         
         !call Norm(Counters%NodTot, NDOFL, TemporaryMappingVector, ReducedDof, OutOfBalanceNorm)
         call Norm(NumberOfGlobalControlPointsUniqueMultipatch, NDOFL, TemporaryMappingVector, ReducedDof, OutOfBalanceNorm)

     
         TemporaryMappingVector = 0.0
          do J = 1, Counters%NEntity
             do I = 1, Counters%N
             if((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoadMixture(I, J) * PBoundary(I) + &
                                         ExtLoad(I, J) * PBoundary(I)
             else
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I)
             end if
           end do
         end do

           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

         if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity

         TemporaryMappingVector = 0.0
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) + ExternalLoadDueToVelocity (I)
           end do
         end do
        end if
    
         !call Norm(Counters%NodTot, NDOFL,  &
         !          TemporaryMappingVector, ReducedDof, &
         !          ExternalLoadNorm)
         call Norm(NumberOfGlobalControlPointsUniqueMultipatch, NDOFL,  &
                   TemporaryMappingVector, ReducedDof, &
                   ExternalLoadNorm)
     
        if (ExternalLoadNorm==0.0) then
          CalParams%ConvergenceCheck%ForceError = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceError = OutOfBalanceNorm / ExternalLoadNorm
        end if       
        
        if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then 


        CalParams%ConvergenceCheck%ForceErrorSoil = CONVERGENCE_ERROR_NOT_USED
       
            
         ! Compute ForceError of the water   
            
         TemporaryMappingVector = 0.0
          
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
             if(CalParams%ApplyAbsorbingBoundary) then  ! absorbing boundary
               TemporaryMappingVector(I) = TemporaryMappingVector(I)+  &
                                         GravityLoadWater(I, J) * PBoundaryWater(I) + &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         IntLoadWaterPrevious(I, J) * PBoundaryWater(I) + &
                                         QVW(I,J) - &
                                         VisDampForceWat(I, J)
             else
                TemporaryMappingVector(I) = TemporaryMappingVector(I)+   & ! no absorbing boundary 
                                         GravityLoadWater(I, J) * PBoundaryWater(I)+ &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         IntLoadWaterPrevious(I, J) * PBoundaryWater(I) + &
                                         QVW(I,J)
             end if
           end do
         end do          
         
         !call Norm(Counters%NodTot, NDOFL, TemporaryMappingVector, ReducedDof, OutOfBalanceNormWater)
         call Norm(NumberOfGlobalControlPointsUniqueMultipatch, NDOFL, TemporaryMappingVector, ReducedDof, OutOfBalanceNormWater)
     
         TemporaryMappingVector = 0.0
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoadWater(I, J) * PBoundaryWater(I) + &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) +   &
                                         QVW(I,J)
 
           end do
         end do

         !call Norm(Counters%NodTot, NDOFL, TemporaryMappingVector, ReducedDof, ExternalLoadNormWater)
         call Norm(NumberOfGlobalControlPointsUniqueMultipatch, NDOFL, TemporaryMappingVector, ReducedDof, ExternalLoadNormWater)
     
        if (ExternalLoadNormWater==0.0) then
          CalParams%ConvergenceCheck%ForceErrorWater = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceErrorWater = OutOfBalanceNormWater / ExternalLoadNormWater
        end if       
        
        end if
        
        end subroutine ComputeForceError
        
        
          subroutine ComputeForceError2LayForm()
        !**********************************************************************
        !
        !    Function:  For 2 Layer Formulation :: Computes the out-of-balance forces error, 
        !               for the solid and for the water phase.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
       implicit none

         ! Local variables
         integer(INTEGER_TYPE) :: I, J, K, INode, IDof
         real(REAL_TYPE) :: OutOfBalanceNormSolid, ExternalLoadNormSolid
         real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
         real(REAL_TYPE) :: OutOfBalanceNormWater, ExternalLoadNormWater
         logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity
         
         ! Multipatch variables 
         integer(INTEGER_TYPE) :: IPatch
         integer(INTEGER_TYPE) :: GlobalNodeID
 
        ! Compute ForceError of the SOLID
                  
          ExternalLoadDueToVelocity = 0.0
          TemporaryMappingVector = 0.0
          
           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do
           
           if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity

            do IPatch = 1, Counters%NPatches ! loop over patches
               do INode = 1, Counters%NodTot(IPatch)!Counters%Nodtot ! loop over control points
                   
                   GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
                   
              do J = 1, NVECTOR
                IDof = (INode-1) * NDOFL      ! global Dof
                do K = 1, Counters%NEntity
                  ExternalLoadDueToVelocity(IDof) = ExternalLoadDueToVelocity(IDof) + FReaction(IDoF, K)
                end do
              end do
               end do ! loop over control points
               end do ! loop over patches
          end if 
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
                 TemporaryMappingVector(I) = TemporaryMappingVector(I) + &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - &
                                         ! IntLoadPrevious gives the internal load at time level n, IntLoad at time level n+1
                                         IntLoadPrevious(I, J) * PBoundary(I) - QVW(I, J) + &
                                         TwoLayerData%InteractionForceSolid(I, J)
           end do
         end do        
         
         call Norm(Counters%Sum_NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof,  &
                   OutOfBalanceNormSolid)
     
         TemporaryMappingVector = 0.0
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoad(I, J) * PBoundary(I)+ &
                                         ExtLoad(I, J) * PBoundary(I) - QVW(I, J) +  &
                                         TwoLayerData%InteractionForceSolid(I, J)
           end do
         end do

            IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

          if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity
         TemporaryMappingVector = 0.0
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                    ExternalLoadDueToVelocity (I)
           end do
         end do
        end if
    
         call Norm(Counters%Sum_NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof, &
                   ExternalLoadNormSolid)
     
        if (ExternalLoadNormSolid==0.0) then
          CalParams%ConvergenceCheck%ForceErrorSoil = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceErrorSoil = OutOfBalanceNormSolid / ExternalLoadNormSolid
        end if       
       
        ! Compute ForceError of the LIQUID  

         TemporaryMappingVector = 0.0
          
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
                TemporaryMappingVector(I) = TemporaryMappingVector(I)+   & ! no absorbing boundary 
                                         GravityLoadWater(I, J) * PBoundaryWater(I)+ &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) - &
                                         ! IntLoadWaterPrevious gives the internal load at time level n,
                                         ! IntLoadWater at time level n+1
                                         IntLoadWaterPrevious(I, J) * PBoundaryWater(I) + &
                                         QVW(I,J) +  &
                                         TwoLayerData%InteractionForceLiquid(I, J)
           end do
         end do       
         
         call Norm(Counters%Sum_NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof,  &
                   OutOfBalanceNormWater)
     
         TemporaryMappingVector = 0.0
         
         do I = 1, Counters%N
           do J = 1, Counters%NEntity
               TemporaryMappingVector(I) = TemporaryMappingVector(I) +  &
                                         GravityLoadWater(I, J) * PBoundaryWater(I) + &
                                         ExtLoadWater(I, J) * PBoundaryWater(I) +   &
                                         QVW(I,J) +  &
                                         TwoLayerData%InteractionForceLiquid(I, J)
 
           end do
         end do        

         call Norm(Counters%Sum_NodTot, NDOFL,  &
                   TemporaryMappingVector, ReducedDof, &
                   ExternalLoadNormWater)
     
        if (ExternalLoadNormWater==0.0) then
          CalParams%ConvergenceCheck%ForceErrorWater = CONVERGENCE_ERROR_NOT_USED
        else
          CalParams%ConvergenceCheck%ForceErrorWater = OutOfBalanceNormWater / ExternalLoadNormWater
        end if       
        
        
          end subroutine ComputeForceError2LayForm
     
          
          
          
          
          
          
       subroutine UpdateParticleVelocityAndMapMomentum(Momentum)!, Momentum_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum       !CC added nEntity
          !real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(inout) :: Momentum_VolLockSmooth       !CC added nEntity
          
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
          
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof_VolLockSmooth
          
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleIncrementalVelocity
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: NodAcc, GeneralizedNodAcc
          
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          !real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape_VolLockSmooth
          
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          
          
          !!CC - changed to function call - needs to be set to zero here
          Momentum = 0.0
          !Momentum_VolLockSmooth = 0.0
          Time = CalParams%TimeIncrement
          NoEn = Counters%nEntity

          !do IAEl = 1, Counters%NAEl   ! Loop over all elements
              do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            IEl = ActiveElement(IAEl, IPatch)
            do I = 1, NVECTOR
                
              IDof(I, 1:ELEMENTNODES) = &
                  ReducedDof(Multipatch_Connecting_Local_To_Global_ControlPoints( &
                  ElementConnectivities(1:ELEMENTNODES, IEl, IPatch),IPatch)) + I
              
              !IDof_VolLockSmooth(I, 1:ELEMENTNODES_VolLockSmooth) = &
              !    ReducedDof_VolLockSmooth( & !Multipatch_Connecting_Local_To_Global_ControlPoints( &
              !    IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch)
              
              NodAcc(I, 1:ELEMENTNODES, 1:Counters%nEntity) = &
                  AccelerationSoil(IDof(I, 1:ELEMENTNODES), 1:NoEn) 
              
              GeneralizedNodAcc(I, 1:ELEMENTNODES, 1:Counters%nEntity) = AccelerationSoilGeneralizedAlpha(IDof(I, 1:ELEMENTNODES), 1:NoEn)
            
            end do
            
            do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 ParticleVelocity = VelocityArray(ParticleIndex,:)
                 PartilceMass = MassArray(ParticleIndex)
                 
                 ParticleShape = ShapeValuesArray(ParticleIndex,:)
                 !ParticleShape_VolLockSmooth = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
                 
                 ParticleIncrementalVelocity = 0.0
                 ParticleAcceleration = 0.0
                 if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) 
                 else
                  iEntity = 1
                 end if 
                  
                 do INode = 1, ELEMENTNODES!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES  ! loop over element nodes
                      
                  do I = 1, NVECTOR   ! loop over each direction 
                    ! Particle accelerations
                    ParticleAcceleration(I) = ParticleAcceleration(I) + ( ParticleShape(INode) * NodAcc(I, INode, iEntity) )
                    
                    ! Particle velocity
                    if ( (CalParams%ApplyGeneralizedAlphaTimeScheme%Switch == .TRUE.) .and. (CalParams%TimeStep /=1) ) then 
                    
                        !Mapping the particle velocity using the generalized alpha time scheme
                        ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + &
                            ( CalParams%TimeIncrement * ParticleShape(INode) * (CalParams%ApplyGeneralizedAlphaTimeScheme%Gamma*NodAcc(I, INode, iEntity) + &
                            (1-CalParams%ApplyGeneralizedAlphaTimeScheme%Gamma)*GeneralizedNodAcc(I, INode, iEntity)) )
                        
                    else 
                        
                        ! Particle x-velocity
                        ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + ( Time * ParticleShape(INode) * NodAcc(I, INode, iEntity) )
                        
                    end if 
                    
                  end do 
                  
                end do !Loop over nodes
               
                ParticleVelocity = ParticleVelocity + ParticleIncrementalVelocity

                do I = 1, NVECTOR ! nodal i-momentum
                  !Momentum(IDof(I,1:(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)), iEntity) = &
                  !    Momentum(IDof(I,1:ELEMENTNODES), iEntity) + ( PartilceMass * ParticleShape * ParticleVelocity(I) )
                  
                    Momentum(IDof(I,1:ELEMENTNODES), iEntity) = &
                      Momentum(IDof(I,1:ELEMENTNODES), iEntity) + ( PartilceMass * ParticleShape * ParticleVelocity(I) )
                  
                  
                    !Momentum_VolLockSmooth(IDof_VolLockSmooth(I,1:ELEMENTNODES_VolLockSmooth), iEntity) = &
                    !  Momentum_VolLockSmooth(IDof_VolLockSmooth(I,1:ELEMENTNODES_VolLockSmooth), iEntity) + ( PartilceMass * ParticleShape_VolLockSmooth * ParticleVelocity(I) )
                  
                end do     

                VelocityArray(ParticleIndex,:) = ParticleVelocity
                AccelerationArray(ParticleIndex,:) =  ParticleAcceleration
                 
              end if ! NumbOfLayers = 1 or SOLID MatPoint
            end do !Loop over particles
                end do !elements    
                
                end do ! patches
          
       end subroutine UpdateParticleVelocityAndMapMomentum  
       
       
       
       
       
       
       
       
       
       
       
       !subroutine MapFJacDet_VolLockSmooth(WeightedFJacDet_Numerator) !(Momentum)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       !
       ! implicit none
       ! 
       ! real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(out) :: WeightedFJacDet_Numerator       !CC added nEntity
       ! 
       ! ! Local variables
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       !   integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   real(REAL_TYPE) :: ParticleVolStrain
       !   real(REAL_TYPE) :: Time, PartilceMass
       !   
       !   ! Multipatch variables 
       !   integer(INTEGER_TYPE) :: IPatch
       !   integer(INTEGER_TYPE) :: aa 
       !   
       !   !!CC - changed to function call - needs to be set to zero here
       !   WeightedFJacDet_Numerator = 0.0
       !   NoEn = Counters%nEntity
       !
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch) ! Loop over all elements 
       !           
       !           ! the background grid should have the same exact active element discretization 
       !           IEl = ActiveElement(IAEl, IPatch)
       !     
       !           do I = 1, NVECTOR
       !             
       !               ! degree of freedom for active element 
       !               IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !                   ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
       !       
       !           end do
       !     
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !               if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !
       !                   PartilceMass = MassArray(ParticleIndex) ! particle mass
       !                   ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! particle shape values array on the lower order grid
       !          
       !                   ! reset volumetric strain 
       !                   !Particles(ParticleIndex)%VolumetricStrain = 0.0 
       !          
       !                   if (CalParams%ApplyContactAlgorithm) then
       !                      iEntity = EntityIDArray(ParticleIndex) 
       !                   else
       !                      iEntity = 1
       !                   end if 
       !          
       !                   
       !                   
       !                   do I = 1, NVECTOR ! nodal i-momentum
       !                       WeightedFJacDet_Numerator( IDof(I,1:ELEMENTNODES_VolLockSmooth), iEntity) = &
       !                           WeightedFJacDet_Numerator( IDof(I,1:ELEMENTNODES_VolLockSmooth), iEntity ) &
       !                           + ( PartilceMass * ParticleShape * Particles(ParticleIndex)%FJacDet )
       !                   end do     
       !
       !               end if ! NumbOfLayers = 1 or SOLID MatPoint
       !     
       !           end do !Loop over particles
       !       
       !       end do !elements    
       !   
       !   end do ! patches
       !
       !end subroutine MapFJacDet_VolLockSmooth  
       
       
       
       

       !call MapStrainTensor_VolLockSmooth(WeightedStrainTensor_Numerator)
       !---------------------------------------------------------------------------
       subroutine MapStrainIncrementTensor_VolLockSmooth(WeightedStrainIncrementTensor_Numerator) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        ! Hardcoded NControlPoints for the higher order mesh of patch 1
        real(REAL_TYPE), dimension(NTENSOR, Counters%NodTot(1),Counters%nEntity), intent(out) :: WeightedStrainIncrementTensor_Numerator       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedStrainIncrementTensor_Numerator = 0.0
          NoEn = Counters%nEntity
       
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, Counters%NAEl(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !do I = 1, NVECTOR
                    
                      ! degree of freedom for active element 
                      !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
                      !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
              
                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES) = &
                          ElementConnectivities(1:ELEMENTNODES, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          
                          
                          do I = 1, NTENSOR ! nodal i-momentum
                              WeightedStrainIncrementTensor_Numerator( I, INode(1:ELEMENTNODES), iEntity) = &
                                  WeightedStrainIncrementTensor_Numerator( I, INode(1:ELEMENTNODES), iEntity ) &
                                  + ( PartilceMass * ParticleShape * Particles(ParticleIndex)%EpsStep(I) )
                          end do     
       
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
       end subroutine MapStrainIncrementTensor_VolLockSmooth  
       
       !---------------------------------------------------------------------------
       subroutine MapStressIncrementTensor_VolLockSmooth(WeightedStressIncrementTensor_Numerator) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        ! Hardcoded NControlPoints for the higher order mesh of patch 1
        real(REAL_TYPE), dimension(NTENSOR, Counters%NodTot(1),Counters%nEntity), intent(out) :: WeightedStressIncrementTensor_Numerator       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedStressIncrementTensor_Numerator = 0.0
          NoEn = Counters%nEntity
       
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, Counters%NAEl(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !do I = 1, NVECTOR
                    
                      ! degree of freedom for active element 
                      !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
                      !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
              
                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES) = &
                          ElementConnectivities(1:ELEMENTNODES, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          
                          
                          do I = 1, NTENSOR ! nodal i-momentum
                              WeightedStressIncrementTensor_Numerator( I, INode(1:ELEMENTNODES), iEntity) = &
                                  WeightedStressIncrementTensor_Numerator( I, INode(1:ELEMENTNODES), iEntity ) &
                                  + ( PartilceMass * ParticleShape * SigmaEffArray(ParticleIndex,I) )
                          end do     
       
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
       end subroutine MapStressIncrementTensor_VolLockSmooth  
       
       !--------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       subroutine MapStrainIncrementTensor_VolLockSmooth_LowerOrder(WeightedStrainIncrementTensor_Numerator) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        real(REAL_TYPE), dimension(NTENSOR, NURBS%NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(out) :: WeightedStrainIncrementTensor_Numerator       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
          real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedStrainIncrementTensor_Numerator = 0.0
          NoEn = Counters%nEntity
            
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, Counters%NAEl(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !do I = 1, NVECTOR
                    
                      ! degree of freedom for active element 
                      !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
                      !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
              
                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES_VolLockSmooth) = &
                          NURBS%ElementConnectivities_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          
                          
                          do I = 1, NTENSOR ! nodal i-momentum
                              WeightedStrainIncrementTensor_Numerator( I, INode(1:ELEMENTNODES_VolLockSmooth), iEntity) = &
                                  WeightedStrainIncrementTensor_Numerator( I, INode(1:ELEMENTNODES_VolLockSmooth), iEntity ) &
                                  + ( PartilceMass * ParticleShape * Particles(ParticleIndex)%EpsStep(I) )
                          end do     
       
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
       end subroutine MapStrainIncrementTensor_VolLockSmooth_LowerOrder  
       
       !---------------------------------------------------------------------------
       subroutine MapStressIncrementTensor_VolLockSmooth_LowerOrder(WeightedStressIncrementTensor_Numerator) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        real(REAL_TYPE), dimension(NTENSOR, NURBS%NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(out) :: WeightedStressIncrementTensor_Numerator       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
          !real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedStressIncrementTensor_Numerator = 0.0
          NoEn = Counters%nEntity
            
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, Counters%NAEl(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !do I = 1, NVECTOR
                    
                      ! degree of freedom for active element 
                      !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
                      !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
              
                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES_VolLockSmooth) = &
                          NURBS%ElementConnectivities_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          
                          
                          do I = 1, NTENSOR ! nodal i-momentum
                              WeightedStressIncrementTensor_Numerator( I, INode(1:ELEMENTNODES_VolLockSmooth), iEntity) = &
                                  WeightedStressIncrementTensor_Numerator( I, INode(1:ELEMENTNODES_VolLockSmooth), iEntity ) &
                                  + ( PartilceMass * ParticleShape * SigmaEffArray(ParticleIndex,I) ) !Particles(ParticleIndex)%EpsStep(I)
                          end do     
       
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
       end subroutine MapStressIncrementTensor_VolLockSmooth_LowerOrder  
       
       !---------------------------------------------------------------------------
       subroutine MapWP_VolLockSmooth(WeightedWP_Numerator) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        ! Hardcoded NControlPoints for the higher order mesh of patch 1
        real(REAL_TYPE), dimension(Counters%NodTot(1),Counters%nEntity), intent(out) :: WeightedWP_Numerator       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedWP_Numerator = 0.0
          NoEn = Counters%nEntity
       
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, Counters%NAEl(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !do I = 1, NVECTOR
                    
                      ! degree of freedom for active element 
                      !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
                      !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
              
                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES) = &
                          ElementConnectivities(1:ELEMENTNODES, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          !if ( abs(Particles(ParticleIndex)%WaterPressure) > 1e-5) then 
                          
                          !do I = 1, NTENSOR ! nodal i-momentum
                              WeightedWP_Numerator( INode(1:ELEMENTNODES), iEntity) = & !I, 
                                  WeightedWP_Numerator( INode(1:ELEMENTNODES), iEntity ) & !I, 
                                  + ( PartilceMass * ParticleShape * Particles(ParticleIndex)%WaterPressure )
                          !end do     
       
                         ! end if 
                          
                              
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
       end subroutine MapWP_VolLockSmooth  
       
       
       
       !---------------------------------------------------------------------------
       subroutine MapWP_VolLockSmooth_LowerOrder(WeightedWP_Numerator) !(Momentum)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
       
        implicit none
        
        ! Hardcoded NControlPoints for the higher order mesh of patch 1
        real(REAL_TYPE), dimension(NURBS%NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(out) :: WeightedWP_Numerator       !CC added nEntity
        
        ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
          !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
          integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: INode !NVECTOR, 
          
          real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
          real(REAL_TYPE) :: ParticleVolStrain
          real(REAL_TYPE) :: Time, PartilceMass
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: aa 
          
          !!CC - changed to function call - needs to be set to zero here
          WeightedWP_Numerator = 0.0
          NoEn = Counters%nEntity
       
          do IPatch = 1, Counters%NPatches ! Loop over patches
          
              do IAEl = 1, Counters%NAEl(IPatch) ! Loop over all elements 
                  
                  ! the background grid should have the same exact active element discretization 
                  IEl = ActiveElement(IAEl, IPatch)
            
                  !do I = 1, NVECTOR
                    
                      ! degree of freedom for active element 
                      !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
                      !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
              
                      ! degree of freedom for active element 
                      INode(1:ELEMENTNODES_VolLockSmooth) = &
                          NURBS%ElementConnectivities_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
              
                      
                      
                  !end do
            
                  do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                      ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                      if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       
                          PartilceMass = MassArray(ParticleIndex) ! particle mass
                          ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! particle shape values array on the lower order grid
                 
                          ! reset volumetric strain 
                          !Particles(ParticleIndex)%VolumetricStrain = 0.0 
                 
                          if (CalParams%ApplyContactAlgorithm) then
                             iEntity = EntityIDArray(ParticleIndex) 
                          else
                             iEntity = 1
                          end if 
                 
                          
                          
                          !do I = 1, NTENSOR ! nodal i-momentum
                              WeightedWP_Numerator( INode(1:ELEMENTNODES_VolLockSmooth), iEntity) = & !I, 
                                  WeightedWP_Numerator( INode(1:ELEMENTNODES_VolLockSmooth), iEntity ) & !I, 
                                  + ( PartilceMass * ParticleShape * Particles(ParticleIndex)%WaterPressure )
                          !end do     
       
                      end if ! NumbOfLayers = 1 or SOLID MatPoint
            
                  end do !Loop over particles
              
              end do !elements    
          
          end do ! patches
       
       end subroutine MapWP_VolLockSmooth_LowerOrder
       
       
       
       !------------------------------------------------------------------------
       
       
       
       
       !subroutine MapStressTensor_VolLockSmooth(WeightedStressTensor_Numerator) !(Momentum)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       !
       ! implicit none
       ! 
       ! real(REAL_TYPE), dimension(NTENSOR, NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(out) :: WeightedStressTensor_Numerator       !CC added nEntity
       ! 
       ! ! Local variables
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, iEntity, ParticleIndex, NoEn ! INode,
       !   !integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: INode !NVECTOR, 
       !   
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   real(REAL_TYPE) :: ParticleVolStrain
       !   real(REAL_TYPE) :: Time, PartilceMass
       !   
       !   ! Multipatch variables 
       !   integer(INTEGER_TYPE) :: IPatch
       !   integer(INTEGER_TYPE) :: aa 
       !   
       !   !!CC - changed to function call - needs to be set to zero here
       !   WeightedStressTensor_Numerator = 0.0
       !   NoEn = Counters%nEntity
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch) ! Loop over all elements 
       !           
       !           ! the background grid should have the same exact active element discretization 
       !           IEl = ActiveElement(IAEl, IPatch)
       !     
       !           !do I = 1, NVECTOR
       !             
       !               ! degree of freedom for active element 
       !               !IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !               !    ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
       !       
       !               ! degree of freedom for active element 
       !               INode(1:ELEMENTNODES_VolLockSmooth) = &
       !                   IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) !ReducedDof_VolLockSmooth( ) + I
       !       
       !               
       !               
       !           !end do
       !     
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !               if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !
       !                   PartilceMass = MassArray(ParticleIndex) ! particle mass
       !                   ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! particle shape values array on the lower order grid
       !          
       !                   ! reset volumetric strain 
       !                   !Particles(ParticleIndex)%VolumetricStrain = 0.0 
       !          
       !                   if (CalParams%ApplyContactAlgorithm) then
       !                      iEntity = EntityIDArray(ParticleIndex) 
       !                   else
       !                      iEntity = 1
       !                   end if 
       !          
       !                   
       !                   
       !                   do I = 1, NTENSOR ! nodal i-momentum
       !                       WeightedStressTensor_Numerator( I, INode(1:ELEMENTNODES_VolLockSmooth), iEntity) = &
       !                           WeightedStressTensor_Numerator( I, INode(1:ELEMENTNODES_VolLockSmooth), iEntity ) &
       !                           + ( PartilceMass * ParticleShape * SigmaEffArray(ParticleIndex, I) ) !Particles(ParticleIndex)%EpsStep(I)
       !                   end do     
       !
       !               end if ! NumbOfLayers = 1 or SOLID MatPoint
       !     
       !           end do !Loop over particles
       !       
       !       end do !elements    
       !   
       !   end do ! patches
       !
       !end subroutine MapStressTensor_VolLockSmooth  
       
       
       
       
       
       
       
       
       
       
       
       !------------------------------------------------------------------------
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       !subroutine MapVolumetricStrain_VolLockSmooth(WeightedNodalVolumetricStrainIncrement_Numerator) !(Momentum)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       !
       ! implicit none
       ! 
       ! real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(out) :: WeightedNodalVolumetricStrainIncrement_Numerator       !CC added nEntity
       ! 
       ! ! Local variables
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       !   integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   real(REAL_TYPE) :: ParticleVolStrain
       !   real(REAL_TYPE) :: Time, PartilceMass
       !   
       !   ! Multipatch variables 
       !   integer(INTEGER_TYPE) :: IPatch
       !   integer(INTEGER_TYPE) :: aa 
       !   
       !   !!CC - changed to function call - needs to be set to zero here
       !   WeightedNodalVolumetricStrainIncrement_Numerator = 0.0
       !   NoEn = Counters%nEntity
       !
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch) ! Loop over all elements 
       !           
       !           ! the background grid should have the same exact active element discretization 
       !           IEl = ActiveElement(IAEl, IPatch)
       !     
       !           do I = 1, NVECTOR
       !             
       !               ! degree of freedom for active element 
       !               IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !                   ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch) ) + I
       !       
       !           end do
       !     
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !               if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !
       !                   PartilceMass = MassArray(ParticleIndex) ! particle mass
       !                   ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! particle shape values array on the lower order grid
       !          
       !                   ! reset volumetric strain 
       !                   Particles(ParticleIndex)%VolumetricStrain = 0.0 
       !          
       !                   if (CalParams%ApplyContactAlgorithm) then
       !                      iEntity = EntityIDArray(ParticleIndex) 
       !                   else
       !                      iEntity = 1
       !                   end if 
       !          
       !                   do I = 1, NVECTOR  
       !                      ! Calculate volumetric strain
       !                      Particles(ParticleIndex)%VolumetricStrain = Particles(ParticleIndex)%VolumetricStrain + &
       !                                                             GetEpsStepI(Particles(ParticleIndex), I)!(Particles(ParticleIndex)%EpsStep(I))         
       !                       !Particles(ParticleIndex)%VolumetricStrain = 1
       !                   end do    
       !        
       !                   if (Particles(ParticleIndex)%VolumetricStrain>0) then 
       !                       aa = 1
       !                   end if 
       !                   
       !                   
       !                   do I = 1, NVECTOR ! nodal i-momentum
       !                       WeightedNodalVolumetricStrainIncrement_Numerator( IDof(I,1:ELEMENTNODES_VolLockSmooth), iEntity) = &
       !                           WeightedNodalVolumetricStrainIncrement_Numerator( IDof(I,1:ELEMENTNODES_VolLockSmooth), iEntity ) &
       !                           + ( PartilceMass * ParticleShape * Particles(ParticleIndex)%VolumetricStrain )
       !                   end do     
       !
       !               end if ! NumbOfLayers = 1 or SOLID MatPoint
       !     
       !           end do !Loop over particles
       !       
       !       end do !elements    
       !   
       !   end do ! patches
       !
       !end subroutine MapVolumetricStrain_VolLockSmooth  
       
       
       
       
       
       
       
       
       
       
       !subroutine MapVolumetricStressGradient_VolLockSmooth(WeightedNodalVolumetricStressGradient_Numerator) !(Momentum)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       !
       ! implicit none
       ! 
       !   !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum       !CC added nEntity
       !
       ! 
       ! real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(out) :: WeightedNodalVolumetricStressGradient_Numerator       !CC added nEntity
       ! 
       ! ! Local variables
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       !   integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   real(REAL_TYPE) :: ParticleVolStrain
       !   real(REAL_TYPE) :: Time, PartilceMass
       !   
       !   ! Multipatch variables 
       !   integer(INTEGER_TYPE) :: IPatch
       !   integer(INTEGER_TYPE) :: aa
       !   
       !   !!CC - changed to function call - needs to be set to zero here
       !   WeightedNodalVolumetricStressGradient_Numerator = 0.0
       !   NoEn = Counters%nEntity
       !
       !
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !     
       !           
       !           ! the background grid should have the same exact active element discretization 
       !           IEl = ActiveElement(IAEl, IPatch)
       !     
       !           do I = 1, NVECTOR
       !       
       !               IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !                   ReducedDof_VolLockSmooth( IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch 
       !       
       !     end do
       !     
       !     do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !
       !          PartilceMass = MassArray(ParticleIndex)
       !          ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !          
       !          Particles(ParticleIndex)%VolumetricStressGradient = 0.0
       !          
       !          if (CalParams%ApplyContactAlgorithm) then
       !           iEntity = EntityIDArray(ParticleIndex) 
       !          else
       !           iEntity = 1
       !          end if 
       !           
       !          !do INode = 1, ELEMENTNODES_VolLockSmooth!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES  ! loop over element nodes
       !               
       !           do I = 1, NVECTOR  
       !               
       !             Particles(ParticleIndex)%VolumetricStressGradient = Particles(ParticleIndex)%VolumetricStressGradient + &
       !                                                             SigmaEffArray(ParticleIndex, I)
       !                                                             !(Particles(ParticleIndex)%StressGradient(I))  
       !             
       !             !Particles(ParticleIndex)%VolumetricStressGradient = 1.0
       !             
       !           end do    
       !           
       !           if (Particles(ParticleIndex)%VolumetricStressGradient>0) then  
       !               aa = 1    
       !           end if 
       !           
       !         !end do !Loop over nodes
       !        
       !         !ParticleVelocity = ParticleVelocity + ParticleIncrementalVelocity
       !
       !         do I = 1, NVECTOR ! nodal i-momentum
       !
       !           WeightedNodalVolumetricStressGradient_Numerator(IDof(I,1:ELEMENTNODES_VolLockSmooth), iEntity) = &
       !               WeightedNodalVolumetricStressGradient_Numerator(IDof(I,1:ELEMENTNODES_VolLockSmooth), iEntity) + &
       !               ( PartilceMass * ParticleShape * Particles(ParticleIndex)%VolumetricStressGradient )
       !           
       !         end do     
       !
       !         !VelocityArray(ParticleIndex,:) = ParticleVelocity
       !         !AccelerationArray(ParticleIndex,:) =  ParticleAcceleration
       !          
       !       end if ! NumbOfLayers = 1 or SOLID MatPoint
       !     end do !Loop over particles
       !         end do !elements    
       !         
       !         end do ! patches
       !   
       !end subroutine MapVolumetricStressGradient_VolLockSmooth  
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       !subroutine MapStressToNodes_VolLockSmooth(NodalStress_VolLockSmooth)!(WeightedNodalVolumetricStrainIncrement_Numerator) !(Momentum)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       !
       ! implicit none
       ! 
       !   !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum       !CC added nEntity
       !
       ! real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: NodalStress_VolLockSmooth ! to the actual mesh
       ! 
       ! ! Local variables
       ! 
       ! integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       ! 
       ! integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof !_VolLockSmooth
       !   
       ! real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape !_VolLockSmooth
       ! 
       ! real(REAL_TYPE), dimension(NVECTOR) :: ParticleStress !NTENSOR
       !     
       !   !real(REAL_TYPE) :: ParticleVolStrain
       ! 
       ! real(REAL_TYPE) :: Time!, PartilceMass
       !   
       !   ! Multipatch variables 
       ! 
       ! integer(INTEGER_TYPE) :: IPatch
       !   
       !   !!CC - changed to function call - needs to be set to zero here
       ! 
       ! NodalStress_VolLockSmooth = 0.0
       ! 
       ! NoEn = Counters%nEntity
       !
       ! do IPatch = 1, Counters%NPatches ! Loop over patches
       !    
       !     do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !     
       !         ! the background grid should have the same exact active element discretization 
       !         IEl = ActiveElement(IAEl, IPatch)
       !     
       !         
       !         do I = 1, NVECTOR    
       !             IDof(I, 1:ELEMENTNODES) = & !_VolLockSmooth
       !                   ReducedDof( IEN(1:ELEMENTNODES, IEl, IPatch)) + I !,IPatch 
       !         end do
       !     
       !         do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !   
       !             ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !       
       !             if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !          
       !                 ParticleStress = SigmaEffArray(ParticleIndex,1:NVECTOR) ! only volumetric stress needed
       !                 ParticleShape = ShapeValuesArray(ParticleIndex,:) !_VolLockSmooth
       !         
       !                 if (CalParams%ApplyContactAlgorithm) then
       !                     iEntity = EntityIDArray(ParticleIndex) 
       !                 else
       !                     iEntity = 1
       !                 end if 
       !
       !                   
       !                 do I = 1, NVECTOR ! nodal i-momentum
       !                     NodalStress_VolLockSmooth(IDof(I,1:ELEMENTNODES), iEntity) = & !_VolLockSmooth
       !                     NodalStress_VolLockSmooth(IDof(I,1:ELEMENTNODES), iEntity) + ( ParticleShape * ParticleStress(I) )
       !                 end do
       !
       !             end if ! NumbOfLayers = 1 or SOLID MatPoint
       !     
       !         end do !Loop over particles
       !      
       !     end do !elements    
       !   
       ! end do ! patches
       !   
       !end subroutine MapStressToNodes_VolLockSmooth  
       
       
       
       
       
       
       
       
       
       
       
       
       !subroutine MapVolumetricStrainToMPs_VolLockSmooth(NodalVolumetricStrain_VolLockSmooth)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       ! 
       ! implicit none
       ! 
       ! real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(in) :: NodalVolumetricStrain_VolLockSmooth       !CC added nEntity
       !   
       !   
       ! ! Local variables
       ! integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       ! integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   
       ! real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodStrain
       ! real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   
       ! 
       ! ! Multipatch variables
       ! integer(INTEGER_TYPE) :: IPatch
       !
       ! !!CC - changed to function call - needs to be set to zero here
       ! NoEn = Counters%nEntity
       !
       !   
       ! do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !     do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !     
       !         IEl = ActiveElement(IAEl, IPatch)
       !       
       !         do I = 1, NVECTOR
       !       
       !             IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !                 ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
       !                 IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       !
       !               NodStrain(I, 1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = &
       !                   NodalVolumetricStrain_VolLockSmooth(IDof(I, 1:ELEMENTNODES_VolLockSmooth), 1:NoEn) 
       !     
       !         end do
       !    
       !           
       !         do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !             ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !             if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !          
       !                 Particles(ParticleIndex)%BBarVolumetricStrain = 0.0                
       !                 ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !          
       !                 if (CalParams%ApplyContactAlgorithm) then
       !                     iEntity = EntityIDArray(ParticleIndex)
       !                 else
       !                     iEntity = 1  
       !                 end if 
       !          
       !                 do INode = 1, ELEMENTNODES_VolLockSmooth
       !
       !                     do I = 1, 1!NVECTOR  
       !
       !                         ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
       !                         ! Particle strain 
       !                         Particles(ParticleIndex)%BBarVolumetricStrain = Particles(ParticleIndex)%BBarVolumetricStrain + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
       !                     end do    
       !                 
       !                 end do !Loop over nodes
       !               
       !             end if ! NumbOfLayers = 1 or SOLID MatPoint
       !           
       !         end do !Loop over particles
       !       
       !     end do !elements    
       !   
       ! end do ! patches
       ! 
       !
       !end subroutine MapVolumetricStrainToMPs_VolLockSmooth  
       
       !subroutine MapFJacDetToMPs_VolLockSmooth(NodalFJacDet_VolLockSmooth)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       ! 
       ! implicit none
       ! 
       ! real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(in) :: NodalFJacDet_VolLockSmooth       !CC added nEntity
       !   
       !   
       ! ! Local variables
       ! integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       ! integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   
       ! real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodFJacDet
       ! real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   
       ! 
       ! ! Multipatch variables
       ! integer(INTEGER_TYPE) :: IPatch
       !
       ! !!CC - changed to function call - needs to be set to zero here
       ! NoEn = Counters%nEntity
       !
       !   
       ! do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !     do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !     
       !         IEl = ActiveElement(IAEl, IPatch)
       !       
       !         do I = 1, NVECTOR
       !       
       !             IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !                 ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
       !                 IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       !
       !               NodFJacDet(I, 1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = &
       !                   NodalFJacDet_VolLockSmooth(IDof(I, 1:ELEMENTNODES_VolLockSmooth), 1:NoEn) 
       !     
       !         end do
       !    
       !           
       !         do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !             ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !             if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !          
       !                 !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
       !                 Particles(ParticleIndex)%FBarJacDet = 0.0
       !                 ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !          
       !                 if (CalParams%ApplyContactAlgorithm) then
       !                     iEntity = EntityIDArray(ParticleIndex)
       !                 else
       !                     iEntity = 1  
       !                 end if 
       !          
       !                 do INode = 1, ELEMENTNODES_VolLockSmooth
       !
       !                     do I = 1, 1!NVECTOR  
       !
       !                         ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
       !                         ! Particle strain 
       !                         Particles(ParticleIndex)%FBarJacDet = Particles(ParticleIndex)%FBarJacDet + ( ParticleShape(INode) * NodFJacDet(I, INode, iEntity) )
       !                     end do    
       !                 
       !                 end do !Loop over nodes
       !               
       !             end if ! NumbOfLayers = 1 or SOLID MatPoint
       !           
       !         end do !Loop over particles
       !       
       !     end do !elements    
       !   
       ! end do ! patches
       ! 
       !
       !end subroutine MapFJacDetToMPs_VolLockSmooth !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       !--------------------------------------------------------------------------------------------------------------------
       
       subroutine MapNodalStrainIncrementToMPs_VolLockSmooth(NodalStrainIncrement_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(NTENSOR,Counters%NodTot(1),Counters%nEntity), intent(in) :: NodalStrainIncrement_VolLockSmooth       !CC added nEntity
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: Nodes_Element
          
        real(REAL_TYPE), dimension(NTENSOR, ELEMENTNODES, Counters%nEntity) :: NodStrain
        real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES) = &
                        ElementConnectivities(1:ELEMENTNODES, IEl, IPatch)
                
                do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodStrain(I, 1:ELEMENTNODES, 1:Counters%nEntity) = &
                          NodalStrainIncrement_VolLockSmooth(I, Nodes_Element(1:ELEMENTNODES), 1:NoEn) 
            
                end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        Particles(ParticleIndex)%EpsStep = 0.0
                        ParticleShape = ShapeValuesArray(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES
       
                            do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                            end do    
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapNodalStrainIncrementToMPs_VolLockSmooth !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       !--------------------------------------------------------------------------------------------------------------------
       
       subroutine MapNodalStressIncrementToMPs_VolLockSmooth(NodalStressIncrement_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(NTENSOR,Counters%NodTot(1),Counters%nEntity), intent(in) :: NodalStressIncrement_VolLockSmooth       !CC added nEntity
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: Nodes_Element
          
        real(REAL_TYPE), dimension(NTENSOR, ELEMENTNODES, Counters%nEntity) :: NodStress
        real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES) = &
                        ElementConnectivities(1:ELEMENTNODES, IEl, IPatch)
                
                do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodStress(I, 1:ELEMENTNODES, 1:Counters%nEntity) = &
                          NodalStressIncrement_VolLockSmooth(I, Nodes_Element(1:ELEMENTNODES), 1:NoEn) 
            
                end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        !Particles(ParticleIndex)%EpsStep = 0.0
                        SigmaEffArray(ParticleIndex,:) = 0.0
                        ParticleShape = ShapeValuesArray(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES
       
                            do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                !Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                                SigmaEffArray(ParticleIndex,I) = SigmaEffArray(ParticleIndex,I) + ( ParticleShape(INode) * NodStress(I, INode, iEntity) )
                            end do    
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    !call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapNodalStressIncrementToMPs_VolLockSmooth !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       !--------------------------------------------------------------------------------------------------------------------

       
       
       subroutine MapNodalStrainIncrementToMPs_VolLockSmooth_LowerOrder(NodalStrainIncrement_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(NTENSOR,NURBS%NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(in) :: NodalStrainIncrement_VolLockSmooth       !CC added nEntity
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: Nodes_Element
          
        real(REAL_TYPE), dimension(NTENSOR, ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodStrain
        real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES_VolLockSmooth) = &
                        NURBS%ElementConnectivities_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)
                
                do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodStrain(I, 1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = &
                          NodalStrainIncrement_VolLockSmooth(I, Nodes_Element(1:ELEMENTNODES_VolLockSmooth), 1:NoEn) 
            
                end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        Particles(ParticleIndex)%EpsStep = 0.0
                        ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES_VolLockSmooth
       
                            do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                            end do    
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapNodalStrainIncrementToMPs_VolLockSmooth_LowerOrder !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       !--------------------------------------------------------------------------------------------------------------------

       
       
       subroutine MapNodalStressIncrementToMPs_VolLockSmooth_LowerOrder(NodalStressIncrement_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(NTENSOR,NURBS%NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(in) :: NodalStressIncrement_VolLockSmooth       !CC added nEntity
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: Nodes_Element
          
        real(REAL_TYPE), dimension(NTENSOR, ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodStress
        real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES_VolLockSmooth) = &
                        NURBS%ElementConnectivities_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)
                
                do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodStress(I, 1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = &
                          NodalStressIncrement_VolLockSmooth(I, Nodes_Element(1:ELEMENTNODES_VolLockSmooth), 1:NoEn) 
            
                end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        !Particles(ParticleIndex)%EpsStep = 0.0
                        SigmaEffArray(ParticleIndex,:)=0.0
                        ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES_VolLockSmooth
       
                            do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                SigmaEffArray(ParticleIndex,I) = SigmaEffArray(ParticleIndex,I) + (ParticleShape(INode) * NodStress(I, INode, iEntity))
                                !Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                            end do    
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    !call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapNodalStressIncrementToMPs_VolLockSmooth_LowerOrder !MapVolumetricStrainToMPs_VolLockSmooth  
       
       !--------------------------------------------------------------------------------------------------------------------
       
       subroutine MapNodalWPToMPs_VolLockSmooth(NodalWP_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(Counters%NodTot(1),Counters%nEntity), intent(in) :: NodalWP_VolLockSmooth       !CC added nEntity !NTENSOR,
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: Nodes_Element
          
        real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: NodWP !NTENSOR, 
        real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES) = &
                        ElementConnectivities(1:ELEMENTNODES, IEl, IPatch)
                
                !do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodWP( 1:ELEMENTNODES, 1:Counters%nEntity) = & !I,
                          NodalWP_VolLockSmooth( Nodes_Element(1:ELEMENTNODES), 1:NoEn) !I,
            
                !end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        Particles(ParticleIndex)%WaterPressure = 0.0
                        ParticleShape = ShapeValuesArray(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES
       
                            !do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                !Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                                Particles(ParticleIndex)%WaterPressure = Particles(ParticleIndex)%WaterPressure + &
                                    ( ParticleShape(INode) * NodWP(INode, iEntity) ) !I, 
                            !end do    
                                
                                
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    !call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    !if (Particles(ParticleIndex)%WaterPressure>0.0) then ! hardcoded 
                    !    Particles(ParticleIndex)%WaterPressure = 0.0
                    !end if
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapNodalWPToMPs_VolLockSmooth !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       !--------------------------------------------------------------------------------------------------------------------

       
       
       subroutine MapNodalWPToMPs_VolLockSmooth_LowerOrder(NodalWP_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        ! hardcoded on the first patch
        real(REAL_TYPE), dimension(NURBS%NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(in) :: NodalWP_VolLockSmooth       !CC added nEntity
          
          
        ! Local variables
        integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
        integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: Nodes_Element
          
        real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodWP !NTENSOR, 
        real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
          
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch
       
        !!CC - changed to function call - needs to be set to zero here
        NoEn = Counters%nEntity
       
          
        do IPatch = 1, Counters%NPatches ! Loop over patches
          
            do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            
                IEl = ActiveElement(IAEl, IPatch)
              
                Nodes_Element(1:ELEMENTNODES_VolLockSmooth) = &
                        NURBS%ElementConnectivities_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)
                
                !do I = 1, NTENSOR
              
                    
                        !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
                        !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       
                    ! start here 
                      NodWP(1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = & !I, 
                          NodalWP_VolLockSmooth(Nodes_Element(1:ELEMENTNODES_VolLockSmooth), 1:NoEn) !I, 
            
                !end do
           
                  
                do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
              
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
              
                    if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
                 
                        !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
                        Particles(ParticleIndex)%WaterPressure = 0.0
                        ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
                 
                        if (CalParams%ApplyContactAlgorithm) then
                            iEntity = EntityIDArray(ParticleIndex)
                        else
                            iEntity = 1  
                        end if 
                 
                        do INode = 1, ELEMENTNODES_VolLockSmooth
       
                            !do I = 1, NTENSOR!NVECTOR  
       
                                ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
                                ! Particle strain 
                                !Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
                                Particles(ParticleIndex)%WaterPressure = Particles(ParticleIndex)%WaterPressure + ( ParticleShape(INode) * NodWP(INode, iEntity) ) !I, 
                            !end do    
                        
                        end do !Loop over nodes
                      
                    end if ! NumbOfLayers = 1 or SOLID MatPoint
                  
                    
                    !call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
                    !if (Particles(ParticleIndex)%WaterPressure>0.0) then ! hardcoded 
                    !    Particles(ParticleIndex)%WaterPressure = 0.0
                    !end if
                end do !Loop over particles
              
            end do !elements    
          
        end do ! patches
        
       
       end subroutine MapNodalWPToMPs_VolLockSmooth_LowerOrder !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       
       
       
       
       !-----------------------------------------------------------------------
       !subroutine MapNodalStressToMPs_VolLockSmooth(NodalStress_VolLockSmooth)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       ! 
       ! implicit none
       ! 
       ! ! hardcoded on the first patch
       ! real(REAL_TYPE), dimension(NTENSOR,NControlPoints_VolLockSmooth(1),Counters%nEntity), intent(in) :: NodalStress_VolLockSmooth       !CC added nEntity
       !   
       !   
       ! ! Local variables
       ! integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       ! integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: Nodes_Element
       !   
       ! real(REAL_TYPE), dimension(NTENSOR, ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodStress
       ! real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   
       ! 
       ! ! Multipatch variables
       ! integer(INTEGER_TYPE) :: IPatch
       !
       ! !!CC - changed to function call - needs to be set to zero here
       ! NoEn = Counters%nEntity
       !
       !   
       ! do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !     do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !     
       !         IEl = ActiveElement(IAEl, IPatch)
       !       
       !         Nodes_Element(1:ELEMENTNODES_VolLockSmooth) = &
       !                 IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)
       !         
       !         do I = 1, NTENSOR
       !       
       !             
       !                 !ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
       !                 !IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       !
       !             ! start here 
       !               NodStress(I, 1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = &
       !                   NodalStress_VolLockSmooth(I, Nodes_Element(1:ELEMENTNODES_VolLockSmooth), 1:NoEn) 
       !     
       !         end do
       !    
       !           
       !         do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !             ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !             if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !          
       !                 !Particles(ParticleIndex)%BBarVolumetricStrain = 0.0    
       !                 !Particles(ParticleIndex)%EpsStep = 0.0
       !                 SigmaEffArray(ParticleIndex, :) = 0.0
       !                 ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !          
       !                 if (CalParams%ApplyContactAlgorithm) then
       !                     iEntity = EntityIDArray(ParticleIndex)
       !                 else
       !                     iEntity = 1  
       !                 end if 
       !          
       !                 do INode = 1, ELEMENTNODES_VolLockSmooth
       !
       !                     do I = 1, NTENSOR!NVECTOR  
       !
       !                         ! --> scalar that is the same for both dofs in this case. Not intersted in mapping it NVECTOR times
       !                         ! Particle strain 
       !                         !Particles(ParticleIndex)%EpsStep(I) = Particles(ParticleIndex)%EpsStep(I) + ( ParticleShape(INode) * NodStrain(I, INode, iEntity) )
       !                     
       !                         SigmaEffArray(ParticleIndex, I) = SigmaEffArray(ParticleIndex, I) + ( ParticleShape(INode) * NodStress(I, INode, iEntity) )
       !                         
       !                     end do    
       !                 
       !                 end do !Loop over nodes
       !               
       !             end if ! NumbOfLayers = 1 or SOLID MatPoint
       !           
       !             
       !             !call IncreaseEps(Particles(ParticleIndex), Particles(ParticleIndex)%EpsStep)
       !             
       !         end do !Loop over particles
       !       
       !     end do !elements    
       !   
       ! end do ! patches
       ! 
       !
       !end subroutine MapNodalStressToMPs_VolLockSmooth !MapVolumetricStrainToMPs_VolLockSmooth  
       
       
       
       
       
       
       
       
       
       
       
       
       
       !------------------------------------------------------------------------
       !subroutine CalculateFBar()
       !
       !implicit none 
       !
       !integer(INTEGER_TYPE), dimension(NDIM,NDIM) :: IdentityMatrix
       !real(REAL_TYPE), dimension(NDIM,NDIM) :: FMatrix, FMatrixInv
       !real(REAL_TYPE), dimension(NDIM,NDIM) :: FBarFMatrixInv
       !real(REAL_TYPE), dimension(NDIM,NDIM) :: F_vol
       !real(REAL_TYPE), dimension(NDIM,NDIM) :: F_dev
       !real(REAL_TYPE), dimension(NDIM,NDIM) :: FBar
       !
       !real(REAL_TYPE), dimension(NDIM,NDIM) :: StrainFBar
       !
       !integer(INTEGER_TYPE) :: IParticle, IDim
       !real(REAL_TYPE) :: DetFMatrix
       !real(REAL_TYPE) :: DetFMatrixInv
       !
       !! initialise
       !IdentityMatrix = 0 
       !
       !do IDim = 1, NDIM 
       !    IdentityMatrix(IDim,IDim) = 1     
       !end do 
       !
       !
       !do IParticle = 1, Counters%NParticles
       ! 
       !    if (NDIM == 2) then 
       !        FMatrix(1,1) = Particles(IParticle)%FDeformationGradient(1)
       !        FMatrix(1,2) = Particles(IParticle)%FDeformationGradient(3)
       !        FMatrix(2,1) = Particles(IParticle)%FDeformationGradient(4)
       !        FMatrix(2,2) = Particles(IParticle)%FDeformationGradient(2)
       !
       !    else if (NDIM == 3) then
       !        
       !        FMatrix(1,1) = Particles(IParticle)%FDeformationGradient(1)
       !        FMatrix(2,2) = Particles(IParticle)%FDeformationGradient(2)
       !        FMatrix(3,3) = Particles(IParticle)%FDeformationGradient(3)
       !        
       !        FMatrix(1,2) = Particles(IParticle)%FDeformationGradient(4)
       !        FMatrix(2,1) = Particles(IParticle)%FDeformationGradient(4)
       !        
       !        FMatrix(1,3) = Particles(IParticle)%FDeformationGradient(6)
       !        FMatrix(3,1) = Particles(IParticle)%FDeformationGradient(6)
       !
       !        FMatrix(2,3) = Particles(IParticle)%FDeformationGradient(5)
       !        FMatrix(3,2) = Particles(IParticle)%FDeformationGradient(5)
       !        
       !    end if 
       !    
       !        
       !    F_vol = ( Particles(IParticle)%FBarJacDet **(1.0/2.0) ) * IdentityMatrix
       !
       !    F_dev = ( Particles(IParticle)%FJacDet **(-1.0/2.0) ) * FMatrix
       !
       !    FBar = matmul(F_vol, F_dev)
       !
       !    
       !
       !    if (NDIM == 2) then 
       !
       !    Particles(IParticle)%FDeformationGradient(1) = FBar(1,1) 
       !    Particles(IParticle)%FDeformationGradient(3) = 0.0!FBar(1,2)
       !    Particles(IParticle)%FDeformationGradient(4) = FBar(2,1)
       !    Particles(IParticle)%FDeformationGradient(2) = FBar(2,2)
       !    
       !    elseif (NDIM == 3) then 
       !        
       !    Particles(IParticle)%FDeformationGradient(1) = FBar(1,1) 
       !    Particles(IParticle)%FDeformationGradient(2) = FBar(2,2)  
       !    Particles(IParticle)%FDeformationGradient(3) = FBar(3,3)  
       !    
       !    !Particles(IParticle)%FDeformationGradient(3) = FBar(1,2)
       !    Particles(IParticle)%FDeformationGradient(4) = FBar(2,1)
       !    Particles(IParticle)%FDeformationGradient(5) = FBar(2,3)
       !    Particles(IParticle)%FDeformationGradient(6) = FBar(3,1)
       !      
       !    
       !    end if 
       !    
       !    
       !
       !    !call RJacInv(NDIM, FMatrix, FMatrixInv, DetFMatrix, DetFMatrixInv)
       !    !
       !    !FBarFMatrixInv = matmul(FBar, FMatrixInv)
       !    
       !    StrainFBar = 0.5*(FBar - IdentityMatrix) !0.5*
       !    
       !    if (NDIM == 2) then 
       !    Particles(IParticle)%EpsStep(1) = StrainFBar(1,1)
       !    Particles(IParticle)%EpsStep(3) = StrainFBar(1,2)*2.0 ! shear component
       !    Particles(IParticle)%EpsStep(4) = 0.0!StrainFBar(2,1)*2.0 ! shear component
       !    Particles(IParticle)%EpsStep(2) = StrainFBar(2,2)
       !    
       !    else if (NDIM == 3) then 
       !    
       !    Particles(IParticle)%EpsStep(1) = StrainFBar(1,1)
       !    Particles(IParticle)%EpsStep(3) = StrainFBar(3,3) !StrainFBar(1,2)*2.0 ! shear component
       !    Particles(IParticle)%EpsStep(4) = StrainFBar(1,2)*2.0  !0.0!StrainFBar(2,1)*2.0 ! shear component
       !    Particles(IParticle)%EpsStep(2) = StrainFBar(2,2)
       !    Particles(IParticle)%EpsStep(5) = StrainFBar(2,3)*2.0
       !    Particles(IParticle)%EpsStep(6) = StrainFBar(1,3)*2.0
       !        
       !    end if 
       !    
       !    
       !
       !    !call IncreaseEps(Particles(IParticle), Particles(IParticle)%EpsStep)
       !    
       !end do
       !
       !end subroutine CalculateFBar
       
       !subroutine MapVolumetricStressGradientToMPs_VolLockSmooth(NodalVolumetricStressGradient_VolLockSmooth)
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       ! 
       ! implicit none
       ! 
       !   real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(in) :: NodalVolumetricStressGradient_VolLockSmooth       !CC added nEntity
       !   ! Local variables
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       !   integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth) :: IDof
       !   !real(REAL_TYPE), dimension(NVECTOR) :: ParticleIncrementalVelocity
       !   !real(REAL_TYPE), dimension(NVECTOR) :: ParticleAcceleration
       !   !real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
       !   real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth, Counters%nEntity) :: NodStressGradient
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape
       !   !real(REAL_TYPE) :: Time, PartilceMass
       !   
       !   ! Multipatch variables 
       !   !integer(INTEGER_TYPE) :: IPatch_Temporary
       !   integer(INTEGER_TYPE) :: IPatch
       !   
       !   
       !   !!CC - changed to function call - needs to be set to zero here
       !   NoEn = Counters%nEntity
       !
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !     
       !           IEl = ActiveElement(IAEl, IPatch)
       !     
       !           do I = 1, NVECTOR
       !       
       !               IDof(I, 1:ELEMENTNODES_VolLockSmooth) = &
       !                   ReducedDof_VolLockSmooth(&!Multipatch_Connecting_Local_To_Global_ControlPoints( &
       !                   IEN_VolLockSmooth(1:ELEMENTNODES_VolLockSmooth, IEl, IPatch)) + I !,IPatch) !ElementConnectivities
       !
       !               NodStressGradient(I, 1:ELEMENTNODES_VolLockSmooth, 1:Counters%nEntity) = &
       !                   NodalVolumetricStressGradient_VolLockSmooth(IDof(I, 1:ELEMENTNODES_VolLockSmooth), 1:NoEn) 
       !     
       !           end do
       !    
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !       
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !       
       !               if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !          
       !                   Particles(ParticleIndex)%FBarVolumetricStressGradient = 0.0                
       !                   ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !          
       !          
       !          
       !
       !                   if (CalParams%ApplyContactAlgorithm) then
       !                       iEntity = EntityIDArray(ParticleIndex) 
       !                   else
       !                       iEntity = 1
       !                   end if 
       !           
       !
       !                   do INode = 1, ELEMENTNODES_VolLockSmooth
       !
       !                       do I = 1, 1!NVECTOR  !1
       !                           ! Particle strain 
       !                           Particles(ParticleIndex)%FBarVolumetricStressGradient = Particles(ParticleIndex)%FBarVolumetricStressGradient + ( ParticleShape(INode) * NodStressGradient(I, INode, iEntity) )
       !                       end do    
       !                   end do !Loop over nodes
       !       
       !               end if ! NumbOfLayers = 1 or SOLID MatPoint
       !     
       !           end do !Loop over particles
       !     
       !       end do !elements    
       !       
       !   end do ! patches
       !   
       !end subroutine MapVolumetricStressGradientToMPs_VolLockSmooth  
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       ! VolLockSmooth_BOOKMARK
       !subroutine BBar_StrainInterpolation()
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       ! implicit none
       ! 
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       !
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape ! shape of this is dependent on lower order control net 
       !   real(REAL_TYPE), dimension(Counters%NParticles) :: VolumetricStrainIncrement ! volumetric (i.e., dilational) strain !, dimension(NTENSOR)
       !   
       !   ! Multipatch variables 
       !   integer(INTEGER_TYPE) :: IPatch_Temporary
       !   integer(INTEGER_TYPE) :: IPatch
       !   integer(INTEGER_TYPE) :: NodeID
       !   
       !   ! Volumetric locking variables 
       !   real(REAL_TYPE) :: ParticleWeight
       !   real(REAL_TYPE) :: MeanStrainIncrement
       !   real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR) :: DeviatoricStrainIncrement ! storing those temporarily here
       !   real(REAL_TYPE), dimension(Counters%NParticles) :: NewVolumetricStrainIncrement
       !   real(REAL_TYPE), dimension(NTENSOR) :: ParticleStrainIncrement
       !   real(REAL_TYPE), allocatable, dimension(:) :: WeightedNodalVolumetricStrainIncrement_Numerator !,:
       !   real(REAL_TYPE), allocatable, dimension(:) :: NodalVolume_Denominator
       !   real(REAL_TYPE), allocatable, dimension(:) :: NodalVolumetricStrainIncrement !, :
       !   
       !   integer(INTEGER_TYPE) :: IErrorTemp
       !   integer(INTEGER_TYPE) :: IDim
       !   
       !   integer(INTEGER_TYPE), dimension(NTENSOR) :: IdentityMatrix
       ! 
       !   if (NDIM == 2) then 
       !       IdentityMatrix = 0.0
       !       IdentityMatrix(1) = 1.0
       !       IdentityMatrix(2) = 1.0
       !       
       !   else if (NDIM ==3) then 
       !       IdentityMatrix = 0.0
       !       IdentityMatrix(1) = 1.0
       !       IdentityMatrix(2) = 1.0
       !       IdentityMatrix(3) = 1.0
       !       
       !   end if
       !   
       !   
       !   IPatch_Temporary = 1
       !   
       !   ! we need to allocate these local variables with the number of control points of the lower order control net
       !   allocate(WeightedNodalVolumetricStrainIncrement_Numerator(NControlPoints_VolLockSmooth(IPatch_Temporary) ), stat=IErrorTemp ) !, NTENSOR
       !   WeightedNodalVolumetricStrainIncrement_Numerator = 0.0
       !   allocate(NodalVolume_Denominator(NControlPoints_VolLockSmooth(IPatch_Temporary) ), stat=IErrorTemp )
       !   NodalVolume_Denominator = 0.0
       !   allocate(NodalVolumetricStrainIncrement(NControlPoints_VolLockSmooth(IPatch_Temporary) ), stat=IErrorTemp ) !, NTENSOR 
       !   NodalVolumetricStrainIncrement = 0.0
       !   
       !   ! initialize 
       !   DeviatoricStrainIncrement = 0.0 ! deviatoric strain increment for each particle (NParticles, NTENSOR)
       !   NewVolumetricStrainIncrement = 0.0 ! new volumetric strain increment for each particle (NParticles, NTENSOR)
       !   
       !   NoEn = Counters%nEntity
       !   
       !   
       !   
       !   
       !   
       !   
       !   
       !   
       !   
       !   
       !   
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch) ! Loop over all active elements 
       !       
       !           IEl = ActiveElement(IAEl, IPatch) ! Obtain active element index
       !           
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !               
       !               if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !                   
       !                   ParticleWeight = Particles(ParticleIndex)%IntegrationWeight   ! obtain particle weight
       !                   ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! obtain lower order shape
       !                   ParticleStrainIncrement = Particles(ParticleIndex)%EpsStep ! NTENSOR sized strain voigt vector
       !                   
       !                   if (CalParams%ApplyContactAlgorithm) then         
       !                       iEntity = EntityIDArray(ParticleIndex) 
       !                   else
       !                       iEntity = 1
       !                   end if
       !                     
       !                   VolumetricStrainIncrement = 0.0
       !                   
       !                   ! Step 1 --> calculate volumetric strains (equivalent to div velocity * TimeIncrement)
       !                   do IDim = 1, NDim
       !                      VolumetricStrainIncrement(ParticleIndex) = VolumetricStrainIncrement(ParticleIndex) + ParticleStrainIncrement(IDim)
       !                   end do
       !                   
       !                   ! calculate deviatoric strain increment
       !                   !DeviatoricStrainIncrement(ParticleIndex,:) = ParticleStrainIncrement - ( (VolumetricStrainIncrement/3.0) * IdentityMatrix)
       !                   
       !                   do INode = 1, ELEMENTNODES_VolLockSmooth ! loop over element nodes
       !                       
       !                       NodeID = IEN_VolLockSmooth(INode, IEl, IPatch) 
       !                       NodalVolume_Denominator(NodeID) = NodalVolume_Denominator(NodeID) + ( ParticleShape(INode) * ParticleWeight )
       !                       
       !                       WeightedNodalVolumetricStrainIncrement_Numerator(NodeID) = WeightedNodalVolumetricStrainIncrement_Numerator(NodeID) + ( ParticleShape(INode) * VolumetricStrainIncrement(ParticleIndex) * ParticleWeight )
       !                   
       !                   end do !Loop over nodes
       !                         
       !               end if ! NumbOfLayers = 1 or SOLID MatPoint
       !             
       !           end do !Loop over particles
       !             
       !       end do !elements                        
       !       
       !   end do ! patches
       !       
       !       
       !       
       !       
       !       
       !       
       !   
       !   ! Calculating the nodal strain ---------------------------------------------------------
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !       
       !           IEl = ActiveElement(IAEl, IPatch)
       !             
       !           do INode = 1, ELEMENTNODES_VolLockSmooth ! loop over element nodes
       !               
       !               NodeID = IEN_VolLockSmooth(INode, IEl, IPatch) 
       !             
       !               !do I = 1, NTENSOR       
       !               NodalVolumetricStrainIncrement(NodeID) = WeightedNodalVolumetricStrainIncrement_Numerator(NodeID) / NodalVolume_Denominator(NodeID)
       !               !end do
       !             
       !           end do
       !         
       !       end do
       !       
       !   end do
       !       
       !   
       !   !-------------------------------------------------------------------------------------------
       !       
       !       
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !       
       !           IEl = ActiveElement(IAEl, IPatch)
       !           
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !                 
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !               ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !                 
       !               do INode = 1, ELEMENTNODES_VolLockSmooth
       !               
       !                   NodeID = IEN_VolLockSmooth(INode, IEl, IPatch) 
       !               
       !                   !do I = 1, NTENSOR
       !
       !                    
       !                   NewVolumetricStrainIncrement(ParticleIndex) = NewVolumetricStrainIncrement(ParticleIndex) + ( ParticleShape(INode) * NodalVolumetricStrainIncrement(NodeID) )
       !                     !, I                                            !, I                                                            !, I
       !                   !end do 
       !                 
       !               end do 
       !                       
       !           end do 
       !           
       !       end do 
       !                   
       !   end do 
       !   
       !       
       !   
       !   do IPart = 1, Counters%NParticles
       !       !Particles(IPart)%EpsStep = ((NewVolumetricStrainIncrement(IPart)/3)*IdentityMatrix) + DeviatoricStrainIncrement(IPart,:)
       !       Particles(IPart)%EpsStep = Particles(IPart)%EpsStep + ((1.0/3.0) * (NewVolumetricStrainIncrement(IPart)-VolumetricStrainIncrement(IPart)))
       !       call IncreaseEps(Particles(IPart), Particles(IPart)%EpsStep)
       !   end do
       !   
       !   
       !end subroutine BBar_StrainInterpolation
       
       
       
       !subroutine UpdateStrainIncrement_VolLockSmooth()
       !
       !implicit none
       !
       !integer(INTEGER_TYPE) :: IPart, IDim
       !integer(INTEGER_TYPE), dimension(NTENSOR) :: IdentityMatrix
       !
       !
       !if (NDIM == 2) then 
       !
       !    IdentityMatrix(1) = 1
       !    IdentityMatrix(2) = 1
       !    IdentityMatrix(3) = 0
       !    IdentityMatrix(4) = 0
       !    
       !elseif (NDIM ==3) then 
       !    
       !    IdentityMatrix(1) = 1
       !    IdentityMatrix(2) = 1
       !    IdentityMatrix(3) = 1
       !    IdentityMatrix(4) = 0
       !    IdentityMatrix(5) = 0
       !    IdentityMatrix(6) = 0
       !end if 
       !
       !do IPart = 1, Counters%NParticles
       !
       !    do IDim = 1, NVECTOR 
       !        !Particles(IPart)%EpsStep(IDim) = Particles(IPart)%EpsStep_VolLockSmooth(IDim)
       !        Particles(IPart)%EpsStep(IDim) = Particles(IPart)%EpsStep(IDim) + &
       !        ( (1.0/3.0) * (Particles(IPart)%VolumetricStrainIncrement_VolLockSmooth - Particles(IPart)%VolumetricStrainIncrement)  ) !* IdentityMatrix
       !         
       !    end do 
       !    
       !    !Particles(IPart)%EpsStep + &
       !        !( (1.0/3.0) * (Particles(IPart)%VolumetricStrainIncrement_VolLockSmooth - Particles(IPart)%VolumetricStrainIncrement) * IdentityMatrix )
       !         
       !    call IncreaseEps(Particles(IPart), Particles(IPart)%EpsStep)
       !
       !end do
       !
       !end subroutine UpdateStrainIncrement_VolLockSmooth
       
       
       
       
       
       
       
       
       
       
       !subroutine UpdateStress_VolLockSmooth
       !
       !implicit none 
       !
       !integer(INTEGER_TYPE) :: IPart, IDim
       !integer(INTEGER_TYPE) :: aa
       !integer(INTEGER_TYPE), dimension(NTENSOR) :: IdentityMatrix
       !
       !if (NDIM == 2) then 
       !IdentityMatrix(1) = 1
       !IdentityMatrix(2) = 1
       !IdentityMatrix(3) = 0
       !IdentityMatrix(4) = 0
       !elseif (NDIM ==3) then 
       !IdentityMatrix(1) = 1
       !IdentityMatrix(2) = 1
       !IdentityMatrix(3) = 1
       !IdentityMatrix(4) = 0
       !IdentityMatrix(5) = 0
       !IdentityMatrix(6) = 0
       !end if 
       !
       !do IPart = 1, Counters%NParticles
       !
       !    
       !    !if (abs(Particles(IPart)%FBarVolumetricStressGradient) > 0.0) then 
       !    !
       !    !    aa = 1
       !    !    
       !    !end if 
       !    !
       !    !
       !    !if (abs(Particles(IPart)%VolumetricStressGradient) > 0.0) then 
       !    !
       !    !    aa = 1
       !    !    
       !    !end if 
       !    
       !    
       !    !SigmaEffArray(IPart,:) = SigmaEffArray(IPart,:) + &
       !    !    ( (1.0/3.0) * ( Particles(IPart)%FBarVolumetricStressGradient - Particles(IPart)%VolumetricStressGradient ) * IdentityMatrix )
       !    do IDim = 1, NDIM 
       !    
       !        SigmaEffArray(IPart,IDim) = SigmaEffArray(IPart,IDim) + &
       !        ( (1.0/3.0) * ( Particles(IPart)%FBarVolumetricStressGradient - Particles(IPart)%VolumetricStressGradient )  ) !* IdentityMatrix
       !    
       !        !SigmaEffArray(IPart,IDim) = SigmaEffArray(IPart,IDim) + &
       !        !( (1.0/3.0) * ( Particles(IPart)%FBarVolumetricStressGradient - Particles(IPart)%VolumetricStressGradient )  ) !* IdentityMatrix
       !        
       !    end do 
       !    
       !end do
       !
       !
       !end subroutine UpdateStress_VolLockSmooth
       
       
       
       
       
       
       
       ! VolLockSmooth_BOOKMARK
       !subroutine FBar_UpdateMPStress()
       ! !**********************************************************************
       ! !
       ! !    Function:  To update particles total velocities and Accelerations
       ! !
       ! !
       ! ! Implemented in the frame of the MPM project.
       ! !
       ! !**********************************************************************
       ! implicit none
       ! 
       !   integer(INTEGER_TYPE) :: I, IEl, IAEl, IPart, INode, iEntity, ParticleIndex, NoEn
       !
       !   real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ParticleShape ! shape of this is dependent on lower order control net 
       !   real(REAL_TYPE), dimension(Counters%NParticles)  :: VolumetricStressGradient ! volumetric (i.e., dilational) strain !, dimension(NTENSOR)
       !   
       !   ! Multipatch variables 
       !   integer(INTEGER_TYPE) :: IPatch_Temporary
       !   integer(INTEGER_TYPE) :: IPatch
       !   integer(INTEGER_TYPE) :: NodeID
       !   
       !   ! Volumetric locking variables 
       !   real(REAL_TYPE) :: ParticleWeight
       !   real(REAL_TYPE) :: MeanStrainIncrement
       !   real(REAL_TYPE), dimension(Counters%NParticles, NTENSOR) :: DeviatoricStrainIncrement ! storing those temporarily here
       !   real(REAL_TYPE), dimension(Counters%NParticles) :: NewVolumetricStressGradient
       !   real(REAL_TYPE), dimension(NDIM) :: ParticleStressGradient
       !   real(REAL_TYPE), allocatable, dimension(:) :: WeightedNodalVolumetricStressGradient_Numerator !,:
       !   real(REAL_TYPE), allocatable, dimension(:) :: NodalVolume_Denominator
       !   real(REAL_TYPE), allocatable, dimension(:) :: NodalVolumetricStressGradient !, :
       !   
       !   integer(INTEGER_TYPE) :: IErrorTemp
       !   integer(INTEGER_TYPE) :: IDim
       !   
       !   integer(INTEGER_TYPE), dimension(NTENSOR) :: IdentityMatrix
       ! 
       !   if (NDIM == 2) then 
       !       IdentityMatrix = 0.0
       !       IdentityMatrix(1) = 1.0
       !       IdentityMatrix(2) = 1.0
       !       
       !   else if (NDIM ==3) then 
       !       IdentityMatrix = 0.0
       !       IdentityMatrix(1) = 1.0
       !       IdentityMatrix(2) = 1.0
       !       IdentityMatrix(3) = 1.0
       !       
       !   end if
       !   
       !   
       !   IPatch_Temporary = 1
       !   
       !   ! we need to allocate these local variables with the number of control points of the lower order control net
       !   allocate(WeightedNodalVolumetricStressGradient_Numerator(NControlPoints_VolLockSmooth(IPatch_Temporary) ), stat=IErrorTemp ) !, NTENSOR
       !   WeightedNodalVolumetricStressGradient_Numerator = 0.0
       !   allocate(NodalVolume_Denominator(NControlPoints_VolLockSmooth(IPatch_Temporary) ), stat=IErrorTemp )
       !   NodalVolume_Denominator = 0.0
       !   allocate(NodalVolumetricStressGradient(NControlPoints_VolLockSmooth(IPatch_Temporary) ), stat=IErrorTemp ) !, NTENSOR 
       !   NodalVolumetricStressGradient = 0.0
       !   
       !   ! initialize 
       !   !DeviatoricStrainIncrement = 0.0 ! deviatoric strain increment for each particle (NParticles, NTENSOR)
       !   NewVolumetricStressGradient = 0.0 ! new volumetric strain increment for each particle (NParticles, NTENSOR)
       !   
       !   NoEn = Counters%nEntity
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch) ! Loop over all active elements 
       !       
       !           IEl = ActiveElement(IAEl, IPatch) ! Obtain active element index
       !           
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !               
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !               
       !               if ( (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) .or. (NFORMULATION==1) ) then ! NumbOfLayers = 1 or SOLID MatPoint
       !                   
       !                   ParticleWeight = Particles(ParticleIndex)%IntegrationWeight   ! obtain particle weight
       !                   ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:) ! obtain lower order shape
       !                   ParticleStressGradient = Particles(ParticleIndex)%SigmaGradient ! NTENSOR sized strain voigt vector
       !                   
       !                   if (CalParams%ApplyContactAlgorithm) then         
       !                       iEntity = EntityIDArray(ParticleIndex) 
       !                   else
       !                       iEntity = 1
       !                   end if
       !                     
       !                   VolumetricStressGradient(ParticleIndex) = 0.0
       !                   
       !                   ! Step 1 --> calculate volumetric strains (equivalent to div velocity * TimeIncrement)
       !                   do IDim = 1, NDim
       !                      VolumetricStressGradient(ParticleIndex) = VolumetricStressGradient(ParticleIndex) + ParticleStressGradient(IDim)
       !                   end do
       !                   
       !                   ! calculate deviatoric strain increment
       !                   !DeviatoricStrainIncrement(ParticleIndex,:) = ParticleStrainIncrement - ( (VolumetricStrainIncrement/3.0) * IdentityMatrix)
       !                   
       !                   do INode = 1, ELEMENTNODES_VolLockSmooth ! loop over element nodes
       !                       NodeID = IEN_VolLockSmooth(INode, IEl, IPatch) 
       !                       NodalVolume_Denominator(NodeID) = NodalVolume_Denominator(NodeID) + ( ParticleShape(INode) * ParticleWeight )
       !                       
       !                       WeightedNodalVolumetricStressGradient_Numerator(NodeID) = WeightedNodalVolumetricStressGradient_Numerator(NodeID) + ( ParticleShape(INode) * VolumetricStressGradient(ParticleIndex) * ParticleWeight )
       !                   end do !Loop over nodes
       !                         
       !               end if ! NumbOfLayers = 1 or SOLID MatPoint
       !             
       !           end do !Loop over particles
       !             
       !       end do !elements                        
       !       
       !   end do ! patches
       !       
       !       
       !       
       !       
       !       
       !       
       !   
       !   ! Calculating the nodal strain ---------------------------------------------------------
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !       
       !           IEl = ActiveElement(IAEl, IPatch)
       !             
       !           do INode = 1, ELEMENTNODES_VolLockSmooth ! loop over element nodes
       !               
       !               NodeID = IEN_VolLockSmooth(INode, IEl, IPatch) 
       !             
       !               !do I = 1, NTENSOR       
       !               NodalVolumetricStressGradient(NodeID) = WeightedNodalVolumetricStressGradient_Numerator(NodeID) / NodalVolume_Denominator(NodeID)
       !               !end do
       !             
       !           end do
       !         
       !       end do
       !       
       !   end do
       !       
       !   
       !   !-------------------------------------------------------------------------------------------
       !       
       !       
       !   
       !   do IPatch = 1, Counters%NPatches ! Loop over patches
       !   
       !       do IAEl = 1, nael_NURBS(IPatch)!Counters%NEl ! Loop over all elements 
       !       
       !           IEl = ActiveElement(IAEl, IPatch)
       !           
       !           do IPart = 1, NPartEle(IEl,IPatch)   ! Loop over all particles in element
       !                 
       !               ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
       !               ParticleShape = ShapeValuesArray_VolLockSmooth(ParticleIndex,:)
       !                 
       !               do INode = 1, ELEMENTNODES_VolLockSmooth
       !               
       !                   NodeID = IEN_VolLockSmooth(INode, IEl, IPatch) 
       !               
       !                   !do I = 1, NTENSOR
       !
       !                    
       !                   NewVolumetricStressGradient(IPart) = NewVolumetricStressGradient(IPart) + ( ParticleShape(INode) * NodalVolumetricStressGradient(NodeID) )
       !                     !, I                                            !, I                                                            !, I
       !                   !end do 
       !                 
       !               end do 
       !                       
       !           end do 
       !           
       !       end do 
       !                   
       !   end do 
       !   
       !       
       !   
       !   do IPart = 1, Counters%NParticles
       !       !Particles(IPart)%EpsStep = ((NewVolumetricStrainIncrement(IPart)/3)*IdentityMatrix) + DeviatoricStrainIncrement(IPart,:)
       !       !Particles(IPart)%EpsStep = Particles(IPart)%EpsStep + ((1.0/3.0) * (NewVolumetricStrainIncrement-VolumetricStrainIncrement))
       !       !call IncreaseEps(Particles(IPart), Particles(IPart)%EpsStep)
       !       SigmaEffArray(IPart,:) = SigmaEffArray(IPart,:) + ((1.0/3.0) * (NewVolumetricStressGradient(IPart) - VolumetricStressGradient(IPart)) * IdentityMatrix)
       !   end do
       !   
       !   
       !   end subroutine FBar_UpdateMPStress
       
       
       
       
       
       
       
       
       
       
       
       
       
       
         
         subroutine GetNodalIncrementalDisplacement(NodalIncDisplacement, NodalVelocities, AccelerationSoilGeneralizedAlpha_Input, AccelerationSoil_Input)
                                                   
        !**********************************************************************
        !
        !    Function:  To calculate the nodal incremental displacements from the nodal velocities
        !
        !    NodalIncDisplacement: The nodal incremental displacements vector
        !    NodalVelocities : Nodal velocities vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: NodalIncDisplacement
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: NodalVelocities
          
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoilGeneralizedAlpha_Input
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: AccelerationSoil_Input
     
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J
          
          do IDOF = 1, Counters%N       !Loop over all degrees of freedom
            do J = 1,Counters%nEntity   !loop over all entities
                
                if ( (CalParams%ApplyGeneralizedAlphaTimeScheme%Switch == .TRUE.) .AND. (.not.IsMPMComputation()) .AND. (CalParams%TimeStep /= 1) ) then ! FEM with Generalized Alpha scheme 
                    NodalIncDisplacement(IDOF,J) = & !NodalIncDisplacement(IDOF,J) + &
                                             (NodalVelocities(IDOF,J)*CalParams%TimeIncrement) + &
                                             ( (CalParams%TimeIncrement**2) * &
                                             ( (0.5-CalParams%ApplyGeneralizedAlphaTimeScheme%Beta)*AccelerationSoilGeneralizedAlpha_Input(IDOF,J) + CalParams%ApplyGeneralizedAlphaTimeScheme%Beta*AccelerationSoil_Input(IDOF,J) ) ) 
                else 

                    NodalIncDisplacement(IDOF,J) = NodalVelocities(IDOF,J) * CalParams%TimeIncrement
              
              end if 
              
              
              
              
              
            end do
          end do
                   
         end subroutine GetNodalIncrementalDisplacement    
         
         
         
        ! subroutine GetNodalIncrementalDisplacement_VolLockSmooth(NodalIncDisplacement, NodalVelocities)
        !                                           
        !!**********************************************************************
        !!
        !!    Function:  To calculate the nodal incremental displacements from the nodal velocities
        !!
        !!    NodalIncDisplacement: The nodal incremental displacements vector
        !!    NodalVelocities : Nodal velocities vector
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(inout) :: NodalIncDisplacement
        !  real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(in) :: NodalVelocities
        !
        !  ! Local variables
        !  integer(INTEGER_TYPE) :: IDOF, J
        !  
        !  do IDOF = 1, Counters%N_VolLockSmooth       !Loop over all degrees of freedom
        !    do J = 1,Counters%nEntity   !loop over all entities
        !      NodalIncDisplacement(IDOF,J) = NodalVelocities(IDOF,J) * CalParams%TimeIncrement
        !    end do
        !  end do
        !           
        ! end subroutine GetNodalIncrementalDisplacement_VolLockSmooth    
         
         
         subroutine UpdateNodalTotalDisplacement(EntityUsed)
                                                   
        !**********************************************************************
        !
        !    Function:  To add the nodal incremental displacements to the total displacemnts
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, EntityUsed
          
          do IDOF = 1, Counters%N    !Loop over all degrees of freedom
            TotalDisplacementSoil(IDOF)= TotalDisplacementSoil(IDOF)+ IncrementalDisplacementSoil(IDOF,EntityUsed)
          end do
          
          if(.not.(NFORMULATION==1)) then
          do IDOF = 1, Counters%N    !Loop over all degrees of freedom
            TotalDisplacementWater(IDOF)= TotalDisplacementWater(IDOF)+ IncrementalDisplacementWater(IDOF,EntityUsed)
          end do
          end if
          
         end subroutine UpdateNodalTotalDisplacement    
         
         
         
         
        ! subroutine UpdateNodalTotalDisplacement_VolLockSmooth(EntityUsed)
        !                                           
        !!**********************************************************************
        !!
        !!    Function:  To add the nodal incremental displacements to the total displacemnts
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !  ! Local variables
        !  integer(INTEGER_TYPE) :: IDOF, EntityUsed
        !  
        !  do IDOF = 1, Counters%N_VolLockSmooth    !Loop over all degrees of freedom
        !    TotalDisplacementSoil_VolLockSmooth(IDOF)= TotalDisplacementSoil_VolLockSmooth(IDOF)+ IncrementalDisplacementSoil_VolLockSmooth(IDOF,EntityUsed)
        !  end do
        !  
        !  if(.not.(NFORMULATION==1)) then
        !  do IDOF = 1, Counters%N    !Loop over all degrees of freedom
        !    TotalDisplacementWater(IDOF)= TotalDisplacementWater(IDOF)+ IncrementalDisplacementWater(IDOF,EntityUsed)
        !  end do
        !  end if
        !  
        !end subroutine UpdateNodalTotalDisplacement_VolLockSmooth    
        
         
        subroutine CalculateKineticEnergy(TotalVelocitySoilLoc, TotalVelocityWaterLoc)
        !**********************************************************************
        !
        !   Function:  Calculate system kinetic energy
        !
        !   Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: TotalVelocitySoilLoc, TotalVelocityWaterLoc
          ! local variables
          integer(INTEGER_TYPE) :: I, INode, IDof
          real(REAL_TYPE) :: SquaredVelocityS, KES, KESn, SquaredVelocityW, KEW, KEWn
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          ! Calculation of kinetic energy of the system using the velocity fields given by nodal values
          KES = 0.0
          KEW = 0.0  
          KESn = 0.0
          KEWn = 0.0
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do INode = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over all nodes
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
          
            IDof = ReducedDof(GlobalNodeID)!INode)

            SquaredVelocityS = 0.0
            do I = 1, NVECTOR
              SquaredVelocityS = SquaredVelocityS + (TotalVelocitySoilLoc(IDof+I, 1)) * (TotalVelocitySoilLoc(IDof+I, 1))
            end do
            KESn = 0.5 * LumpedMassDry(IDof+1, 1) * SquaredVelocityS
            
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation
              SquaredVelocityW = 0.0
              do I = 1, NVECTOR
                SquaredVelocityW = SquaredVelocityW + (TotalVelocityWaterLoc(IDof+I, 1)) * (TotalVelocityWaterLoc(IDof+I, 1))
              end do
              KEWn = 0.5 * LumpedMassNWater(IDof+1, 1) * SquaredVelocityW
            end if ! end if consolidation
            KES = KES + KESn
            KEW = KEW + KEWn

          end do ! end loop over nodes (i.e., control points)
          end do ! loop over patches
          
          CalParams%ConvergenceCheck%KineticEnergy = KES + KEW
          CalParams%ConvergenceCheck%KineticEnergySoil = KES
          CalParams%ConvergenceCheck%KineticEnergyWater = KEW
          
        end subroutine CalculateKineticEnergy
        
        
        subroutine CalculateKineticEnergy2LayForm(TotalVelocitySoilLoc, TotalVelocityWaterLoc)
        !**********************************************************************
        !
        !   Function:  Calculate system kinetic energy in 2 Layer formulation
        !
        !   Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: TotalVelocitySoilLoc, TotalVelocityWaterLoc    
          ! local variables
          integer(INTEGER_TYPE) :: I, INode, IDof
          real(REAL_TYPE) :: SquaredVelocityS, KES, KESn, SquaredVelocityW, KEW, KEWn
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          
          ! Calculation of kinetic energy of the system
          KES = 0.0
          KEW = 0.0 
          KESn = 0.0
          KEWn = 0.0
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do INode = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over all nodes
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
          
            IDof = ReducedDof(GlobalNodeID)!INode)
            
            ! Kinetic energy of solid phase
            SquaredVelocityS = 0.0
            do I = 1, NVECTOR
              SquaredVelocityS = SquaredVelocityS + (TotalVelocitySoilLoc(IDof+I, 1)) * (TotalVelocitySoilLoc(IDof+I, 1))
            end do
            KESn = 0.5 * LumpedMassDry(IDof+1, 1) * SquaredVelocityS
            
            ! Kinetic energy of water phase
            SquaredVelocityW = 0.0
            do I = 1, NVECTOR
                SquaredVelocityW = SquaredVelocityW + (TotalVelocityWaterLoc(IDof+I, 1)) * (TotalVelocityWaterLoc(IDof+I, 1))
            end do
            KEWn = 0.5 * LumpedMassWater(IDof+1, 1) * SquaredVelocityW
            KES = KES + KESn
            KEW = KEW + KEWn

          end do ! end loop over nodes
          
          end do ! loop over control points

          
          CalParams%ConvergenceCheck%KineticEnergy = KES + KEW       
          CalParams%ConvergenceCheck%KineticEnergySoil = KES          
          CalParams%ConvergenceCheck%KineticEnergyWater = KEW
            
        end subroutine CalculateKineticEnergy2LayForm
        
        
        subroutine CalculateKineticEnergyFEM(TotalVelocitySoilLoc, TotalVelocityWaterLoc)
        !**********************************************************************
        !
        !   Function:  Calculate system kinetic energy in FEM
        !
        !**********************************************************************
        implicit none 
        
          ! arguments
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: TotalVelocitySoilLoc, TotalVelocityWaterLoc
          
          !local variable
          integer(INTEGER_TYPE) :: I, INode, IDof
          real(REAL_TYPE) :: SquaredVelocityS, KES, SquaredVelocityW, KEW
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          KES = 0.0
          KEW = 0.0          
          
          ! Calculation of kinetic energy of the system

          do IPatch = 1, Counters%NPatches ! loop over patches
          do INode = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over all nodes
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
              
            IDof = ReducedDof(GlobalNodeID)!INode)
            
            SquaredVelocityS = 0.0
            do I = 1, NVECTOR
              SquaredVelocityS = SquaredVelocityS + (TotalVelocitySoilLoc(IDof+I, 1)) * (TotalVelocitySoilLoc(IDof+I, 1))
            end do
            KES  = KES + 0.5 * LumpedMassDry(IDof+1, 1) * SquaredVelocityS
          
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation
              SquaredVelocityW = 0.0
              do I = 1, NVECTOR
                SquaredVelocityW = SquaredVelocityW + (TotalVelocityWaterLoc(IDof+I, 1)) * (TotalVelocityWaterLoc(IDof+I, 1))
              end do
              KEW  = KEW + 0.5 * LumpedMassNWater(IDof+1, 1) * SquaredVelocityW
            end if ! end if consolidation
          
          end do ! end loop over nodes
          
          end do ! end loop for patches
         
          CalParams%ConvergenceCheck%KineticEnergy = KES + KEW       
          CalParams%ConvergenceCheck%KineticEnergySoil = KES          
          CalParams%ConvergenceCheck%KineticEnergyWater = KEW
        
        end subroutine CalculateKineticEnergyFEM

        subroutine CalculateIntAndExtWorks(IncrementalDisplacementWater)
        !**********************************************************************
        !
        !    Function:  Calculate internal and external works
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IncrementalDisplacementWater
          
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, I, J, K, INode
          real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
          real(REAL_TYPE) :: IncrementalExtWorkSoil
          real(REAL_TYPE) :: IncrementalIntWorkSoil
          real(REAL_TYPE) :: IncrementalExtWorkWater
          real(REAL_TYPE) :: IncrementalIntWorkWater
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: DummyIntLoadWater!, DummyIntLoadWaterPrevious
          logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID

          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then ! if consolidation

            IntLoadPrevious = IntLoad
            IntLoadWaterPrevious = IntLoadWater
            IntLoadWaterPorosityPrevious = IntLoadWaterPorosity
            QVWPorosityPrevious = QVWPorosity
            
            if(CalParams%ApplyPartialSaturation) then
             DummyIntLoadWaterPrevious =  BishopIntLoad  !account for partial saturation to be used in mixture equation
            else
             DummyIntLoadWaterPrevious = IntLoadWater
            end if
            
            call GetNodalIntForces()
            call ConsolidationIntForces(IntLoadWater) ! this is where we calculate the IntLoadWater
            call ConsolidationIntForcesPorosity(IntLoadWaterPorosity)
            call GetQVWArrayPorosity(QVWPorosity)
            
            if((CalParams%NumberOfPhases==3).or.(CalParams%ApplyPartialSaturation)) then
              call ConsolidationForcesBishop(BishopIntLoad)  
              DummyIntLoadWater =  BishopIntLoad  !account for partial saturation
            else
              DummyIntLoadWater = IntLoadWater
            end if
              
            IncrementalIntWorkSoil = 0.0
            IncrementalExtWorkSoil = 0.0
            IncrementalIntWorkWater = 0.0
            IncrementalExtWorkWater = 0.0
            
            do IDOF = 1, Counters%N ! Loop over all degrees of freedom
                do J = 1, Counters%nEntity
                    
                    ! calculate incremental internal and external work of soil
                    IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                             0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J) + DummyIntLoadWater(IDOF,J) + &
                                             DummyIntLoadWaterPrevious(IDOF,J) - IntLoadWaterPorosity(IDOF,J) - &
                                             IntLoadWaterPorosityPrevious(IDOF,J)) * &
                                             IncrementalDisplacementSoil(IDOF,J)                        ! mid point rule
                    IncrementalExtWorkSoil = IncrementalExtWorkSoil + &
                                             (ExtLoad(IDOF,J) - ExtLoadWaterPorosity(IDOF,J) + &
                                             GravityLoadMixture(IDOF,J) - GravityLoadWaterPorosity(IDOF,J) - &
                                             0.5 * (QVWPorosity(IDOF,J) + QVWPorosityPrevious(IDOF,J))) *  &
                                             IncrementalDisplacementSoil(IDOF,J)                        ! left/mid point rule
                    
                    ! calculate incremental internal and external work of water
                    IncrementalIntWorkWater = IncrementalIntWorkWater +  &
                                              0.5 * (IntLoadWaterPorosity(IDOF,J) + IntLoadWaterPorosityPrevious(IDOF,J)) *   &
                                              IncrementalDisplacementWater(IDOF,J)                      ! mid point rule
                    IncrementalExtWorkWater = IncrementalExtWorkWater +  &
                                              (ExtLoadWaterPorosity(IDOF,J) + &
                                              GravityLoadWaterPorosity(IDOF,J) + &
                                              0.5 * (QVWPorosity(IDOF,J) + QVWPorosityPrevious(IDOF,J))) *  &
                                              IncrementalDisplacementWater(IDOF,J)                      ! left/mid point rule
                end do
            end do
          
            CalParams%ConvergenceCheck%InternalWorkSoil = CalParams%ConvergenceCheck%InternalWorkSoil + IncrementalIntWorkSoil
            CalParams%ConvergenceCheck%ExternalWorkSoil = CalParams%ConvergenceCheck%ExternalWorkSoil + IncrementalExtWorkSoil
            
            CalParams%ConvergenceCheck%InternalWorkWater = CalParams%ConvergenceCheck%InternalWorkWater + IncrementalIntWorkWater 
            CalParams%ConvergenceCheck%ExternalWorkWater = CalParams%ConvergenceCheck%ExternalWorkWater + IncrementalExtWorkWater
            
            CalParams%ConvergenceCheck%InternalWork = CalParams%ConvergenceCheck%InternalWorkSoil + CalParams%ConvergenceCheck%InternalWorkWater
            CalParams%ConvergenceCheck%ExternalWork = CalParams%ConvergenceCheck%ExternalWorkSoil + CalParams%ConvergenceCheck%ExternalWorkWater
           
          else ! 1-phase calculation
          
            IntLoadPrevious = IntLoad
            call GetNodalIntForces()  
              
            IncrementalIntWorkSoil = 0.0
            IncrementalExtWorkSoil = 0.0
            ExternalLoadDueToVelocity = 0.0
          
            do IPatch = 1, Counters%NPatches ! loop over patches
            do INode = 1, Counters%NodTot(IPatch)!Counters%Nodtot ! loop over control points
                
                GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
                
              do J = 1, NVECTOR
                IDof = (INode-1) * NDOFL + J      ! global Dof
                do K = 1, Counters%NEntity
                    ExternalLoadDueToVelocity(IDof) = ExternalLoadDueToVelocity(IDof) + FReaction(IDoF, K)
                end do
              end do
            end do ! loop over control points
            end do ! loop over patches

           IsPrescribedVelocity = .false.
           IsPrescribedSurfaceVelocity = .false.
           do I = 1, NVECTOR
             IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
             IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
           end do

           if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity 
              do IDOF = 1, Counters%N ! Loop over all degrees of freedom
                do J = 1, Counters%nEntity   
                  IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                           0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                           IncrementalDisplacementSoil(IDOF,J)                  ! mid point rule
                  IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                           ExternalLoadDueToVelocity(IDOF) *  &
                                           TotalVelocitySoil(IDOF,J) * CalParams%TimeIncrement  ! left point rule
                end do
              end do
            else
            ! load
              do IDOF = 1, Counters%N ! Loop over all degrees of freedom
                do J = 1, Counters%nEntity
                  IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                           0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                           IncrementalDisplacementSoil(IDOF,J)                  ! mid point rule
                  IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                           (ExtLoad(IDOF,J) + GravityLoad(IDOF,J)) *  &
                                           IncrementalDisplacementSoil(IDOF,J)                  ! left point rule
                end do
              end do
            endif

            CalParams%ConvergenceCheck%InternalWork = CalParams%ConvergenceCheck%InternalWork + IncrementalIntWorkSoil 
            CalParams%ConvergenceCheck%ExternalWork = CalParams%ConvergenceCheck%ExternalWork + IncrementalExtWorkSoil

          end if
          
        end subroutine CalculateIntAndExtWorks
        
        subroutine CalculateIntAndExtWorks2LayForm(IncrementalDisplacementWater)
        !**********************************************************************
        !
        !    Function:  Calculate internal and external works
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none       
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: IncrementalDisplacementWater
          
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, I, J, K, INode
          real(REAL_TYPE), dimension(Counters%N) :: ExternalLoadDueToVelocity
          real(REAL_TYPE) :: IncrementalExtWorkSoil
          real(REAL_TYPE) :: IncrementalIntWorkSoil
          real(REAL_TYPE) :: IncrementalExtWorkWater
          real(REAL_TYPE) :: IncrementalIntWorkWater
          logical :: IsPrescribedVelocity, IsPrescribedSurfaceVelocity
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
   
          ! Calculate Internal and External Work of the SOLID and WATER

          IntLoadPrevious = IntLoad
          IntLoadWaterPrevious = IntLoadWater
            
          call GetNodalIntForces()
          call ConsolidationIntForces(IntLoadWater)
            
          IncrementalIntWorkSoil = 0.0
          IncrementalExtWorkSoil = 0.0
          IncrementalIntWorkWater = 0.0
          IncrementalExtWorkWater = 0.0
          ExternalLoadDueToVelocity = 0.0
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do INode = 1, Counters%NodTot(IPatch)!Counters%Nodtot ! loop over control points
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(INode,IPatch)
              
            do J = 1, NVECTOR
              IDof = ReducedDof(GlobalNodeID)+J!(INode-1) * NDOFL      ! global Dof
              do K = 1, Counters%NEntity
                ExternalLoadDueToVelocity(IDof) =  &
                  ExternalLoadDueToVelocity(IDof) +  &
                  FReaction(IDoF, K)
              end do
            end do
          end do !loop over control points
          end do !loop over patches
          
           do IDOF = 1, Counters%N ! Loop over all degrees of freedom
            do J = 1, Counters%nEntity   
 
              IsPrescribedVelocity = .false.
              IsPrescribedSurfaceVelocity = .false.
              do I = 1, NVECTOR
                IsPrescribedVelocity = IsPrescribedVelocity .or. CalParams%ApplyPrescribedVelocity(I)
                IsPrescribedSurfaceVelocity = IsPrescribedSurfaceVelocity .or. CalParams%ApplySurfacePrescribedVelocity(I)
              end do

              if ( IsPrescribedVelocity .or. ( IsPrescribedSurfaceVelocity .and. (.not.IsMPMComputation()) ) ) then ! prescribedvelocity 
            
                IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                         0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                         IncrementalDisplacementSoil(IDOF,J)                  ! mid point rule
                IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                         ExternalLoadDueToVelocity(IDOF) *  &
                                         TotalVelocitySoil(IDOF,J) * CalParams%TimeIncrement

                IncrementalIntWorkWater = IncrementalIntWorkWater +  &
                                          0.5 * (IntLoadWater(IDOF,J) + IntLoadWaterPrevious(IDOF,J)) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! mid point rule
                IncrementalExtWorkWater = IncrementalExtWorkWater +  &
                                          (ExtLoadWater(IDOF,J) + GravityLoadWater(IDOF,J) +  &
                                          QVW(IDOF,J) + TwoLayerData%InteractionForceLiquid(IDOF, J) ) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! left point rule
     
             else ! load
                
                IncrementalIntWorkSoil = IncrementalIntWorkSoil +  &
                                         0.5 * (IntLoad(IDOF,J) + IntLoadPrevious(IDOF,J)) *  &
                                         IncrementalDisplacementSoil(IDOF,J)                        ! mid point rule
                IncrementalExtWorkSoil = IncrementalExtWorkSoil +  &
                                         (ExtLoad(IDOF,J) + GravityLoad(IDOF,J) - &
                                         QVW(IDOF, J) + TwoLayerData%InteractionForceSolid(IDOF, J)) *  &
                                         IncrementalDisplacementSoil(IDOF,J)                        ! left point rule
                    
                IncrementalIntWorkWater = IncrementalIntWorkWater +  &
                                          0.5 * (IntLoadWater(IDOF,J) + IntLoadWaterPrevious(IDOF,J)) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! mid point rule 
                IncrementalExtWorkWater = IncrementalExtWorkWater +  &
                                          (ExtLoadWater(IDOF,J) + GravityLoadWater(IDOF,J) +  &
                                          QVW(IDOF, J) + TwoLayerData%InteractionForceLiquid(IDOF, J)) *  &
                                          IncrementalDisplacementWater(IDOF,J)                      ! left point rule
     
             end if
            end do   
          end do
          
          CalParams%ConvergenceCheck%InternalWorkSoil = CalParams%ConvergenceCheck%InternalWorkSoil + IncrementalIntWorkSoil 
          CalParams%ConvergenceCheck%ExternalWorkSoil = CalParams%ConvergenceCheck%ExternalWorkSoil + IncrementalExtWorkSoil

          CalParams%ConvergenceCheck%InternalWorkWater = CalParams%ConvergenceCheck%InternalWorkWater + IncrementalIntWorkWater 
          CalParams%ConvergenceCheck%ExternalWorkWater = CalParams%ConvergenceCheck%ExternalWorkWater + IncrementalExtWorkWater
            
        end subroutine CalculateIntAndExtWorks2LayForm

        subroutine  GetAndUpdateHydraulicHeadLoad(HHElemID, NumberOfSides, HHNodesCon, HHDummyLoad)
        !HH= HydraulicHead  
 
        use ModMeshAdjacencies
        use ModReadCalculationData
        use ModMeshInfo
                
        logical, dimension(:), allocatable, intent(out) :: HHElemID
        integer(INTEGER_TYPE), dimension(:,:), allocatable, intent(out) :: HHNodesCon
        integer(INTEGER_TYPE),intent(out) ::NumberOfSides
        real(REAL_TYPE), dimension(:, :, :), allocatable,intent(out)  :: HHDummyLoad
        !local
        integer(INTEGER_TYPE) :: IElement, IError,ISide, IsBoundarySide
        integer(INTEGER_TYPE) :: NodeID, J, I, K
            
        
        integer(INTEGER_TYPE), dimension(ELEMENTBOUNDARYNODES) :: LocalNodes
        logical, dimension (ELEMENTBOUNDARYNODES) :: NodePair
        logical, dimension (Counters%Sum_NEl, ELEMENTSIDES) :: HHElemSide
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary
        integer(INTEGER_TYPE) :: IPatch
           
       ! initialization
        NumberOfSides = 0
        allocate (HHElemID(Counters%Sum_NEl), stat = IError)
        HHElemID = .false.
        
        HHElemSide = .false.
       
        do IPatch = 1, Counters%NPatches ! loop over patches
        do IElement = 1, Counters%NEl(IPatch)!Counters%NEl ! loop all elements
            if (IsActiveElement(IElement,IPatch)) then
                do ISide =1,ELEMENTSIDES !loop over side nodes
                    IsBoundarySide = BoundaryElementSurface(IElement,ISide,IsActiveElement(:,IPatch), Counters%Sum_NEl) !Give 1 if ISide is adiacent to inactive element
                    if (IsBoundarySide==1) then !side is on a free surface
                        call DetermineSideNodes(ISide,LocalNodes) !Local ID (1, 2, 3...) of nodes at the boundary of ISide, TRI3,TRI6 = 2 nodes, TETRA = 3 nodes
                        do J=1,ELEMENTBOUNDARYNODES

                            NodeID = ElementConnectivities(LocalNodes(J),IElement, IPatch_Temporary) !Global name of boundary node
                            call DetectHydraulicHeadNode(IsHydraulicHeadNode)
                            NodePair(J)= IsHydraulicHeadNode(NodeID)
                        end do
                        
                        if (count(NodePair) == ELEMENTBOUNDARYNODES) then !If both nodes of the ISide are inside infiltration area then store name of element and side

                            NumberOfSides = NumberOfSides + 1      !count sides
                            HHElemID(IElement) = .true.
                            HHElemSide(IElement,ISide) = .true. !matrix which store loaded side for the elements
                        end if
                    end if
                end do
            end if
        end do ! elements
        end do ! patches
        
        NumberOfSides = count(HHElemSide)
              
        if (allocated(HHNodesCon)) then
            deallocate(HHNodesCon, stat = IError)
        end if
                
        if (NDIM == 3) then
            allocate(HHNodesCon(N_BOUNDARY_NODES_HOE,NumberOfSides), stat = IError)
        elseif (NDIM == 2) then
            allocate(HHNodesCon(ELEMENTBOUNDARYNODES,NumberOfSides), stat = IError)
        end if
      
        K = 1
        do IPatch = 1, Counters%NPatches ! loop over patches
        do IElement = 1, Counters%NEl(IPatch)!Counters%NEl
            do ISide = 1,ELEMENTSIDES
                if (HHElemSide(IElement,ISide)) then
                    call DetermineSideNodes(ISide,LocalNodes)
                    do I=1,ELEMENTBOUNDARYNODES
                        NodeID = ElementConnectivities(LocalNodes(I), IElement, IPatch_Temporary) !Global name of boundary node
                        HHNodesCon(I, K) =  NodeID
                    end do
                    K = K+1
                end if
            end do
        end do ! elements
        end do ! patches
                 
        if (allocated(HHDummyLoad)) then
            deallocate(HHDummyLoad, stat = IError)
        end if
        
         if (NDIM == 3) then
            allocate(HHDummyLoad(N_BOUNDARY_NODES_HOE, NDOFL,NumberOfSides), stat = IError)
        elseif (NDIM == 2) then
            allocate(HHDummyLoad(ELEMENTBOUNDARYNODES, NDOFL,NumberOfSides), stat = IError)
        end if
        HHDummyLoad = 1.0
                 
        end subroutine GetAndUpdateHydraulicHeadLoad       
        
        logical function ConvergenceCheck() 
        !**********************************************************************
        !
        !    Function:  Check Convergence
        !
        !**********************************************************************
        
        implicit none        
        
          logical :: ConvergenceSoil, ConvergenceWater, Convergence
        
          if(NFORMULATION==1) then
              call ComputeForceError() ! hardcoded
          else
              call ComputeForceError2LayForm()
          end if
 
          if (CalParams%ConvergenceCheck%ExternalWork==0.0) then
            CalParams%ConvergenceCheck%KineticError = CONVERGENCE_ERROR_NOT_USED
          else
            CalParams%ConvergenceCheck%KineticError =  &
              abs(CalParams%ConvergenceCheck%KineticEnergy / CalParams%ConvergenceCheck%ExternalWork)
          end if
          
          if (CalParams%ConvergenceCheck%ExternalWorkSoil==0.0) then
            CalParams%ConvergenceCheck%KineticErrorSoil = CONVERGENCE_ERROR_NOT_USED
          else
            CalParams%ConvergenceCheck%KineticErrorSoil =  &
              abs(CalParams%ConvergenceCheck%KineticEnergySoil / CalParams%ConvergenceCheck%ExternalWorkSoil)
          end if

          if (CalParams%ConvergenceCheck%ExternalWorkWater==0.0) then
            CalParams%ConvergenceCheck%KineticErrorWater = CONVERGENCE_ERROR_NOT_USED
          else
            CalParams%ConvergenceCheck%KineticErrorWater =  &
              abs(CalParams%ConvergenceCheck%KineticEnergyWater / CalParams%ConvergenceCheck%ExternalWorkWater)            
          end if
          
          if (CalParams%ApplyQuasiStatic) then
            Convergence = &
              (CalParams%ConvergenceCheck%ForceError<=CalParams%ToleratedErrorForce.or. &
               CalParams%ConvergenceCheck%ForceError==CONVERGENCE_ERROR_NOT_USED).and. &
              (CalParams%ConvergenceCheck%KineticError<=CalParams%ToleratedErrorEnergy.or. &
               CalParams%ConvergenceCheck%KineticError==CONVERGENCE_ERROR_NOT_USED)
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
              ConvergenceSoil =  &
                (CalParams%ConvergenceCheck%ForceErrorSoil<=CalParams%ToleratedErrorForce.or. &
                 CalParams%ConvergenceCheck%ForceErrorSoil==CONVERGENCE_ERROR_NOT_USED).and. &
                (CalParams%ConvergenceCheck%KineticErrorSoil<=CalParams%ToleratedErrorEnergy.or. &
                 CalParams%ConvergenceCheck%KineticErrorSoil==CONVERGENCE_ERROR_NOT_USED)
              ConvergenceWater =  &
                (CalParams%ConvergenceCheck%ForceErrorWater<=CalParams%ToleratedErrorForceWater.or. &
                 CalParams%ConvergenceCheck%ForceErrorWater==CONVERGENCE_ERROR_NOT_USED).and. &
                (CalParams%ConvergenceCheck%KineticErrorWater<=CalParams%ToleratedErrorEnergyWater.or. &
                 CalParams%ConvergenceCheck%KineticErrorWater==CONVERGENCE_ERROR_NOT_USED)
            else
              ConvergenceSoil = .true.
              ConvergenceWater = .true.
            end if
            ConvergenceCheck = ((Convergence.and.ConvergenceSoil.and.ConvergenceWater).or. &
                                 (CalParams%TimeStep>CalParams%MaxTimeSteps))
          else 
            ! Use the real time instead of number of increments
            ConvergenceCheck = (CalParams%TotalRealTime>=CalParams%TotalTime)
          end if
                    
        end function ConvergenceCheck
        
        logical function DivergenceCheck() 
        !**********************************************************************
        !
        !    Function:  Check Divergence
        !
        !    Divergence in case of negative energy dissipation
        !    Divergence tolerance is -ToleratedDivergence (in J)
        !
        !**********************************************************************
        
        implicit none        
        
          logical :: DivergenceSoil, DivergenceWater, Divergence
 
          CalParams%ConvergenceCheck%EnergyDissipation = CalParams%ConvergenceCheck%ExternalWork -  &
                                                         CalParams%ConvergenceCheck%InternalWork - &
                                                         CalParams%ConvergenceCheck%KineticEnergy + &
                                                         CalParams%ConvergenceCheck%KineticEnergy0
          
          Divergence = (CalParams%ConvergenceCheck%EnergyDissipation<=(-1 * CalParams%ToleratedDivergence))
          
          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
            CalParams%ConvergenceCheck%EnergyDissipationSoil = CalParams%ConvergenceCheck%ExternalWorkSoil -  &
                                                               CalParams%ConvergenceCheck%InternalWorkSoil - &
                                                               CalParams%ConvergenceCheck%KineticEnergySoil + &
                                                               CalParams%ConvergenceCheck%KineticEnergySoil0
            
            CalParams%ConvergenceCheck%EnergyDissipationWater = CalParams%ConvergenceCheck%ExternalWorkWater -  &
                                                                CalParams%ConvergenceCheck%InternalWorkWater - &
                                                                CalParams%ConvergenceCheck%KineticEnergyWater + &
                                                                CalParams%ConvergenceCheck%KineticEnergyWater0
            
            DivergenceSoil = (CalParams%ConvergenceCheck%EnergyDissipationSoil<=(-1 * CalParams%ToleratedDivergence))
            
            DivergenceWater = (CalParams%ConvergenceCheck%EnergyDissipationWater<=(-1 * CalParams%ToleratedDivergence))
          else
            DivergenceSoil = .false.
            DivergenceWater = .false.
          end if
          
          DivergenceCheck = (Divergence.or.DivergenceSoil.or.DivergenceWater)
          
        end function DivergenceCheck
                
         subroutine DynUpdateParticleWeights( )
        !**********************************************************************
        !
        !    Function:  Update the weights of the material points.
        ! 
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: I, IAEl, IEl, Int, IntGlo, MaterialID, J
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
          real(REAL_TYPE) :: VolumetricStrain, ConstDensitySolid, ConstDensityLiquid
          real(REAL_TYPE) :: LiquidPressure, Density, ConcRatioSolid
          real(REAL_TYPE) :: RatioDensity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ElementNodalDensity, ElementNodalConcRatio
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          logical :: PressureGreaterThreshold, PressureLowerThreshold, FullyFilled, FreeSurface 
          logical :: UpdDensityLiquid
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
        
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 1-layer formulation or 2-layer form with 1 Phase
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==1).or. &
               ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1))) then
                                
          
              !do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
                  IEl = ActiveElement(IAEl, IPatch)   
            
            do Int = 1, NPartEle(IEl,IPatch) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex(Int, IEl, IPatch)
              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI(Particles(IntGlo),I)
              end do
          
             if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
                 
                 ! Consider stress state 
                 LiquidPressure = (SigmaEffArray(IntGlo,1) + SigmaEffArray(IntGlo,2) + SigmaEffArray(IntGlo,3))/3.0 ! this is valid for 2D and 3D
     
                 MaterialID = MaterialIDArray(IntGlo)
     
                 ! Determine interpolated density
                 call CalculateRatioDensity(IEl,IntGlo,RatioDensity)
      
                 PressureGreaterThreshold = (Particles(IntGlo)%WaterPressure>=CalParams%LiquidPressureCavitationThreshold)
                 PressureLowerThreshold = (Particles(IntGlo)%WaterPressure<CalParams%LiquidPressureCavitationThreshold)
                 FreeSurface = (Particles(IntGlo)%LiquidFreeSurface==1.0)
                 FullyFilled = (RatioDensity>1.0)
               

                 UpdDensityLiquid = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
                    ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))
                 
                 if (.not.UpdDensityLiquid) then
                    Particles(IntGlo)%Density = MatParams(MaterialID)%FluidThresholdDensity/1000         
                 else 
                    Particles(IntGlo)%LiquidFreeSurface = 0.0 ! Particle doesn't belong to the free surface
                    ! hardcoded to remove below and keep it the same as updated in the stress subroutine 
                    Particles(IntGlo)%Density = Particles(IntGlo)%Density /  (VolumetricStrain + 1.0)
                    
                    ! Duan Zhang et al. (2011)
                    !Particles(IntGlo)%Density = 2*Particles(IntGlo)%ShearModulus[1-]
                 end if
                 
                 Particles(IntGlo)%IntegrationWeight = MassArray(IntGlo) / Particles(IntGlo)%Density

             else
                 ! Duan Zhang et al. (2011)
                 ! hardcoded to remove below and keep it the same as updated in the stress subroutine 
                 ! for SOLID or MIXTURE Material Point               
                 ! In Duang Zhang et al. (2011), the density should be a function of GetEpsStepI(Particles(IntGlo),2)
                 !Particles(IntGlo)%IntegrationWeight
                 Particles(IntGlo)%IntegrationWeight = Particles(IntGlo)%IntegrationWeight * (1.0 + VolumetricStrain)
             end if
             end do !loop over material points
                end do ! loop over elements
                end do ! loop over patches
            
          end if ! 1-layer formulation or 2-layer form with 1 Phase
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! 2-layer formulation with 2 Phases
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if((NFORMULATION==2).and.(CalParams%NumberOfPhases==2)) then
              
          !--- Calculate  density  ----
              ConstDensityLiquid = 0.0
              ConstDensitySolid = 0.0
                
              do J = 1, Counters%NLayers ! loop over all material points in element
                if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)    
                end if
              end do 
          
          
              !do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
            do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
                  IEl = ActiveElement(IAEl, IPatch)   
            
            do Int = 1, NPartEle(IEl,IPatch) ! Loop over all material points of the element
                
              ! Determine global ID of Material Point 
              IntGlo = GetParticleIndex(Int, IEl, IPatch)
              
              ! Determine volumetric strain of Material Point 
              VolumetricStrain = 0.0
              do I = 1, NVECTOR ! for 2D and 3D all three components have to be taken into account as eps_tt in axisymmetric \=0
                 VolumetricStrain = VolumetricStrain + GetEpsStepI(Particles(IntGlo),I)
              end do
          
             if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
     
              PressureGreaterThreshold = (Particles(IntGlo)%WaterPressure>=CalParams%LiquidPressureCavitationThreshold)
              PressureLowerThreshold = (Particles(IntGlo)%WaterPressure<CalParams%LiquidPressureCavitationThreshold)
              FreeSurface = (Particles(IntGlo)%LiquidFreeSurface==1.0)
              if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
               FullyFilled = (Particles(IntGlo)%FillingRatio>0.98)
              else
               FullyFilled = (Particles(IntGlo)%FillingRatio>1.0)
              end if

              UpdDensityLiquid = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
              ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))
                 
                 if (.not.UpdDensityLiquid) then
                  if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
                    ! Compute density from nodal density field
                    NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl, IPatch_Temporary)    
                    ElementNodalConcRatio = TwoLayerData%Nodes(NodeIDs)%DensitySolidS / ConstDensitySolid
                 
                    ConcRatioSolid = 0.0
                    ParticleShape = ShapeValuesArray(IntGlo,:)
                    Particles(IntGlo)%Density = 0.0
                    do I = 1, ELEMENTNODES
                      ConcRatioSolid = ConcRatioSolid + ParticleShape(I) * ElementNodalConcRatio(I) 
                    end do
                  else
                    ConcRatioSolid = 0.0
                  end if
                    Particles(IntGlo)%Density = Particles(IntGlo)%ConstDensity * (1.0 - ConcRatioSolid)        
                 else ! For output purposes, possible future use, update density
                    Particles(IntGlo)%LiquidFreeSurface = 0.0d0 ! Particle doesn't belong to the free surface
                    Particles(IntGlo)%Density = Particles(IntGlo)%Density / (VolumetricStrain + 1.0)
                 end if
                 
                 Particles(IntGlo)%IntegrationWeight = MassArray(IntGlo) / Particles(IntGlo)%Density
                 Particles(IntGlo)%EffConcentrationRatioLiquid = Particles(IntGlo)%Density / Particles(IntGlo)%ConstDensity
                 
              end if ! LIQUID Material Point
                 
             if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) then
                 
                 if(Particles(IntGlo)%PhaseStatus==PhaseStatusLIQUID) then 
                    ! Compute density from nodal density field
                    NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl, IPatch_Temporary)
                    ElementNodalDensity = TwoLayerData%Nodes(NodeIDs)%DensitySolidS
                                        
                    Density = 0.0
                    ParticleShape = ShapeValuesArray(IntGlo,:)
                    Particles(IntGlo)%Density = 0.0
                    do I = 1, ELEMENTNODES
                      Density = Density + ParticleShape(I) * ElementNodalDensity(I) 
                    end do
                    Particles(IntGlo)%Density = Density
                    
                    ! Calculate new effective porosity
                    Particles(IntGlo)%EffConcentrationRatioSolid = Particles(IntGlo)%Density / Particles(IntGlo)%ConstDensity
                    Particles(IntGlo)%EffPorosity = 1.0 - Particles(IntGlo)%EffConcentrationRatioSolid
                    
                    if(Particles(IntGlo)%EffPorosity<CalParams%LimitPorosity) then
                      Particles(IntGlo)%EffPorosity = CalParams%LimitPorosity        
                      Particles(IntGlo)%EffConcentrationRatioSolid = 1.0 - Particles(IntGlo)%EffPorosity
                      Particles(IntGlo)%Density = Particles(IntGlo)%EffConcentrationRatioSolid * Particles(IntGlo)%ConstDensity
                    end if
                    
                 else  ! Not Fluidized soil
                    Particles(IntGlo)%Density = Particles(IntGlo)%Density / (VolumetricStrain + 1.0)
                 end if
                 
                 Particles(IntGlo)%IntegrationWeight = MassArray(IntGlo) / Particles(IntGlo)%Density
                 Particles(IntGlo)%EffConcentrationRatioSolid = Particles(IntGlo)%Density / Particles(IntGlo)%ConstDensity
                 Particles(IntGlo)%EffPorosity = 1.0 - Particles(IntGlo)%EffConcentrationRatioSolid
                 
              end if   ! SOLID Material Point              
                 
             end do !loop over material points
                end do ! loop over elements
                
                end do ! loop over patches
            
          end if !2-layer formulation with 2 Phases
         
        end subroutine DynUpdateParticleWeights
        

        
        subroutine DynUpdateParticlePorosity( )
        !**********************************************************************
        !
        !    Function:  Updates particle porosity
        !
        !    - The "mass balance of the solid" has been considered. 
        !    - The variation of the solid density is neglected.          
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: I, IParticle, ParticleIndex
          real(REAL_TYPE) :: DEpsVol

          do IParticle = 1, Counters%NParticles  ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            DEpsVol = 0.0
             do I = 1, NVECTOR
              DEpsVol = DEpsVol + GetEpsStepI(Particles(ParticleIndex),I)
             end do

            Particles(ParticleIndex)%Porosity = Particles(ParticleIndex)%Porosity + (1.0d0 - Particles(ParticleIndex)%Porosity) * DEpsVol
            
            if (Particles(ParticleIndex)%Porosity<0.0) then
                Particles(ParticleIndex)%Porosity = 0.0
            end if
            
          end do ! loop over particles
         
        end subroutine DynUpdateParticlePorosity
        
        
        
        
        subroutine DynUpdateParticleConductivity( )
      !**********************************************************************
      !
      !    Function:  Updates particle permeability based on Xie and Leo (2004) 
      !               Large strain oedometer is the benchmark. 
      !
      !                - need to update porosity 
      !
      ! Implemented in the frame of the MPM project.
      !
      !**********************************************************************

      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: I, IParticle, ParticleIndex
      real(REAL_TYPE) :: DEpsVol
      
      real(REAL_TYPE) :: two
      
      ! defining
      two = 2.0
      
      

      do IParticle = 1, Counters%NParticles  ! Loop over particles
          ParticleIndex = GetParticleIndexFromList(IParticle)
         
          !! hardcoded
          !Particles(IParticle)%InitialConductivity = MatParams(1)%HydraulicConductivityLiquid ! --> stored for updating ..%Conductivity later on
          
         Particles(ParticleIndex)%Conductivity = MatParams(1)%HydraulicConductivityLiquid * &!Particles(ParticleIndex)%InitialConductivity * &
             ((1-MatParams(1)%InitialPorosity)/(1-Particles(ParticleIndex)%Porosity))**two
         
         
         if (Particles(ParticleIndex)%Conductivity<0.0) then
            Particles(ParticleIndex)%Conductivity = 0.0
         end if

      end do ! loop over particles

        end subroutine DynUpdateParticleConductivity
        
        
        

        
        subroutine DynUpdateParticleDegreeOfSaturation( )
        !**********************************************************************
        !
        !    Function:  Updates Degree of Saturation (Sr) of the particle 
        !               The expression of the Degree of Saturation is the RETENTION CURVE
        !
        !**********************************************************************
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: ParticleIndex, IParticle, ISet
        real(REAL_TYPE) :: Sr
        
          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            ISet = MaterialIDArray(ParticleIndex)
        
              if (MatParams(ISet)%RetentionCurve== SWRC_VANGENUCHTEN) then
                  call DynUpdateParticleDegreeOfSaturationVanGenuchten(ParticleIndex,ISet,Sr)
              else if (MatParams(ISet)%RetentionCurve==SWRC_LINEAR) then
                  call DynUpdateParticleDegreeOfSaturationLinear(ParticleIndex,ISet,Sr)
              end if
          
          end do
        
        end subroutine DynUpdateParticleDegreeOfSaturation
        
        
        subroutine DynUpdateParticleHydraulicConductivity( )
         !**********************************************************************
        !
        !    Function:  Updates hydraulic conductivity (K) of the particle 
        !    as a function of Degree of Saturation (Sr)
        !               
        !
        !**********************************************************************
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: ParticleIndex, IParticle, ISet
        
          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            ISet = MaterialIDArray(ParticleIndex)
            
            if (MatParams(ISet)%HydraulicConductivityCurve==HCC_CONSTANT) then
            ! do nothing, Kw doesn't vary
            else if (MatParams(ISet)%HydraulicConductivityCurve==HCC_HILLEL) then
                  call DynUpdateParticleHydraulicConductivityHillel(ParticleIndex,ISet)
            else if (MatParams(ISet)%HydraulicConductivityCurve==HCC_MUALEM) then
                  call DynUpdateParticleHydraulicConductivityMualem(ParticleIndex,ISet)  
            end if
          
          end do
        
        end subroutine DynUpdateParticleHydraulicConductivity
        
        subroutine DynUpdateParticleHydraulicConductivityHillel(ParticleIndex,ISet)
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
        if (Sr >=1.0)  RETURN

        r = MatParams(ISet)%rexponentHillel
        ksat = MatParams(ISet)%HydraulicConductivityLiquid !% this is the saturated permeability
        
        Particles(ParticleIndex)%Conductivity = ksat * (Sr**(r))
        
        end subroutine DynUpdateParticleHydraulicConductivityHillel
        
        subroutine DynUpdateParticleHydraulicConductivityMualem(ParticleIndex,ISet)
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
        
        end subroutine DynUpdateParticleHydraulicConductivityMualem
       
        
        subroutine DynUpdateParticleLiquidWeight()
        !**********************************************************************
        !
        !    Function: updates the liquid density (WD)
        !               
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        real(REAL_TYPE) ::  Pw, T, Kw, g
        real(REAL_TYPE) ::  Alpha, Beta, Pw0, WD0

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)

            Pw = -Particles(ParticleIndex)%WaterPressure
            T = Particles(ParticleIndex)%Temperature   !Temperature (�C)
            Kw = Particles(ParticleIndex)%BulkWater
            g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
            Alpha = -0.00034d0  !Default=-0.00034, volumetric thermal expansion coefficient for water (1/�C)
            Beta = 1/Kw         !Default=0.00045, water compressibility (1/MPa) ---> Inverse of the bulk modulus
            Pw0 =  0.0d0        !Default=0.1, Reference water pressure (MPa) 
            WD0 = 1.0026d0      !Default=1002.6, Reference water density (kg/m3)
        
            Particles(ParticleIndex)%WaterWeight = WD0*exp(Beta*(Pw-Pw0) + Alpha*T)*g
     
          end do ! loop over particles

        end subroutine DynUpdateParticleLiquidWeight


        subroutine DynUpdateParticleGasWeight()
        !**********************************************************************
        !
        !    Function: Updates the gas density (GD)
        !               with respect de Pw, Pg, and T
        !               The expression of the Gas Density is based
        !               on the IDEAL GASES LAW (PV=nRT)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        real(REAL_TYPE) ::  Pw, Pg, T, WD
        real(REAL_TYPE) ::  g, Suc, Ma, Mw, R, Pv0, F, Pv
        real(REAL_TYPE) ::  n1, n2

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)

            g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
            Pw = Particles(ParticleIndex)%WaterPressure
            Pg = Particles(ParticleIndex)%GasPressure
            T = Particles(ParticleIndex)%Temperature           !Temperature (�C)
            WD = Particles(ParticleIndex)%WaterWeight/g        !Water density of the particle
      
            Suc = Pw-Pg     !Suction

            if (Suc<=0.0) RETURN

            T = T+273.15d0              !Temperature (�K)
            Ma = 0.028d0                !Default=0.02895, Molecular gas of dry air (kg/mol)
            Mw = 0.018d0                !Default=0.01801528, Molecular mass of water (=vapour) (kg/mol)
            R = 0.008314d0              !Default=8.3144621, Gas constant value (J/mol*K)

            !Calculation of the vapour Pressure (Pg = Pv + Pa) using the psychrometric law 
            n1 = 136075.0d0
            n2 = -5239.7d0
            Pv0 = n1*exp(n2/T)
            F = exp(-(Mw*Suc)/(R*T*WD))     !Psychrometric Law
            Pv = Pv0*F                      !Vapour Pressure

            Particles(ParticleIndex)%GasWeight = (Pv*(Mw-Ma) + Pg*Ma)/(R*T)*g
        
          end do ! loop over particles

        end subroutine DynUpdateParticleGasWeight


        subroutine DynUpdateParticleMixedWeight
        !**********************************************************************
        !
        !    Function: Updates the density of the mixture
        !               with respect de new densities, porosity and degree of saturation
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        ! Local variables
        integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        real(REAL_TYPE) ::  MW, WW, GW, Sr, Sg, N, MixedWeight
        real(REAL_TYPE), dimension(NVECTOR) :: GravityVector
        
        GravityVector = CalParams%GravityData%GravityVector

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
 
          MW = Particles(ParticleIndex)%MaterialWeight      ! Material weight (dry) of the particle       
          WW = Particles(ParticleIndex)%WaterWeight         ! Water weight of the particle
          GW = Particles(ParticleIndex)%GasWeight           ! Gas weight of the particle
          Sr = Particles(ParticleIndex)%DegreeSaturation    ! Degree of Saturation (liquid) of the particle
          Sg = 1.0d0 - Sr                                   ! Degree of Satutation (gas) of the particle
          N = Particles(ParticleIndex)%Porosity             ! Porosity of the particle
          MixedWeight = MW + N*Sr*WW + N*Sg*GW
          

          Particles(ParticleIndex)%MixedWeight = MixedWeight
          Particles(ParticleIndex)%FbodyMixed(:) = MixedWeight*GravityVector(:)*Particles(ParticleIndex)%IntegrationWeight
    
          end do ! loop over particles
        
        end subroutine DynUpdateParticleMixedWeight
        
        

        subroutine DynUpdateParticleAirInWaterMassFraction( )
        !**********************************************************************
        !
        !    Function:  Update Mass Fraction (Henry's Law)
        !
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
          real(REAL_TYPE) :: H, Mw, Ma, HM     ! Henry's constant, molecular water mass, molecular air mass   
     
          H = 100.0              ! Taking Henry's coefficient as a constant. Default=10000000 (KPa)
          Ma = 0.028d0                ! Default=0.02895, Molecular gas of dry air (kg/mol)
          Mw = 0.018d0                ! Default=0.01801528, Molecular mass of water (kg/mol)
          HM = Ma/(H*Mw)
            
          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            
            Particles(ParticleIndex)%AirInWaterMassFraction = Particles(ParticleIndex)%GasPressure*HM
     
           end do ! loop over particles
         
        end subroutine DynUpdateParticleAirInWaterMassFraction
        
        
        
        subroutine DynUpdateParticleVapourInGasMassFraction( )
        !**********************************************************************
        !
        !    Function:  Update Mass Fraction (Psychrometric's Law)
        !
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
 
        implicit none
       
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
          real(REAL_TYPE) ::  Pw, Pg, T, wVG, WD
          real(REAL_TYPE) ::  g, Mw, R, Suc

          g = CalParams%GravityData%GAccel  !Gravity (m/s2)
          Mw = 0.018d0                !Default=0.01801528, Molecular mass of water (kg/mol)
          R = 0.008314d0              !Default= 8.3144621, Gas constant value (J/mol*K)

          do IParticle = 1, Counters%NParticles ! Loop over particles
            ParticleIndex = GetParticleIndexFromList(IParticle)
            
            Pw = Particles(ParticleIndex)%WaterPressure
            Pg = Particles(ParticleIndex)%GasPressure
            T = Particles(ParticleIndex)%Temperature           !Temperature (�C)
            wVG = Particles(ParticleIndex)%VapourInGasMassFraction
            WD = Particles(ParticleIndex)%WaterWeight/g        !Water density of the particle
      
            T = T+273.15d0  !Temperature (�K)
            Suc = Pw-Pg     !Suction

           Particles(ParticleIndex)%VapourInGasMassFraction = exp(-Mw*Suc/(R*T*WD))
     
           end do ! loop over particles
         
        end subroutine DynUpdateParticleVapourInGasMassFraction



        subroutine UpdateNodes(EntityUsed)
        !**********************************************************************
        !
        !    Function:  Updating of nodal coordinates from displacements for
        !               updated mesh analysis.
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE) :: EntityUsed
          
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, J0
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot ! loop over control points
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            J0 = ReducedDof(GlobalNodeID)
            do J = 1, NVECTOR
              
              NodalCoordinatesUpd(I, J,IPatch) = NodalCoordinatesUpd(I, J,IPatch) +  IncrementalDisplacementSoil(J0 + J,EntityUsed)
                
            end do
          end do ! loop over control points 
          end do ! loop over patches
          
        end subroutine UpdateNodes

        
        subroutine UpdateNodesMeshAdjust(IncUTot)
        !**********************************************************************
        !
        !    Function:  Updating of nodal coordinates from displacements for
        !               updated mesh analysis.
        ! 
        !     IncUTot : Incremental nodal displacements of load step
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N),intent(in) :: IncUTot  
     
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, J0
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: GlobalNodeID
          
          do IPatch = 1, Counters%NPatches ! loop over patches
          do I = 1, Counters%NodTot(IPatch)!Counters%NodTot
              
              GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
              
            J0 = ReducedDof(GlobalNodeID)!I)
            do J = 1, NVECTOR
              
              NodalCoordinatesUpd(I, J, IPatch) = NodalCoordinates(I, J, IPatch) + IncUTot(J0 + J)  
             end do
          end do ! loop over control points
          end do ! loop over patches
                   
        end subroutine UpdateNodesMeshAdjust
        
        
        
        !subroutine UpdateNodesMeshAdjust_VolLockSmooth(IncUTot)
        !!**********************************************************************
        !!
        !!    Function:  Updating of nodal coordinates from displacements for
        !!               updated mesh analysis.
        !! 
        !!     IncUTot : Incremental nodal displacements of load step
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  real(REAL_TYPE), dimension(Counters%N_VolLockSmooth),intent(in) :: IncUTot  
        !
        !  ! Local variables
        !  integer(INTEGER_TYPE) :: I, J, J0
        !  
        !  ! Multipatch variables 
        !  !integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        !  integer(INTEGER_TYPE) :: IPatch
        !  integer(INTEGER_TYPE) :: GlobalNodeID
        !  
        !  do IPatch = 1, Counters%NPatches ! loop over patches
        !  do I = 1, NControlPoints_VolLockSmooth(IPatch)!Counters%NodTot
        !      
        !      !GlobalNodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(I,IPatch)
        !      
        !    J0 = ReducedDof(I)!I)GlobalNodeID)
        !    do J = 1, NVECTOR
        !      
        !      NodalCoordinatesUpd_VolLockSmooth(I, J, IPatch) = NodalCoordinates_VolLockSmooth(I, J, IPatch) + IncUTot(J0 + J)  
        !       
        !    end do
        !    
        !  end do ! loop over control points
        !  end do ! loop over patches
        !           
        !end subroutine UpdateNodesMeshAdjust_VolLockSmooth
    

        subroutine StressAndPorePressureSmoothening()
        !**********************************************************************
        !
        !    Function: Smoothing of stress and pore pressure in case of GP integration.
        !
        ! Implemented in the frame of the MPM project.
        !
        ! Note : This subroutine works only if NumberOfLayers = 1 !!! Temprary Solution
        !
        !**********************************************************************
        
        implicit none
      
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, IEl, NElemPart, IPart, ParticleIndex
          integer(INTEGER_TYPE) :: CheckMaterialID, IAEl
          real(REAL_TYPE) :: ParticleStress
          real(REAL_TYPE) :: WPAverage, DPBVDAverage
          real(REAL_TYPE) :: GPAverage
          real(REAL_TYPE), dimension(NTENSOR) :: StressAverage
          logical :: IsMixedElement
          real(REAL_TYPE) :: ParticleVol, ParticlesVols, AverageWeight
          real(REAL_TYPE), dimension(NTENSOR) :: SigmaEffStep, SigmaEff
          
          integer(INTEGER_TYPE) :: IPatch

          if(.not.(NFORMULATION==1)) RETURN

          !do IAEl = 1, Counters%NAEl
              do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            IEl = ActiveElement(IAEl, IPatch)
              NElemPart = NPartEle(IEl,IPatch) 
              StressAverage = 0.0 ! Initialise average stress of element IEl
              WPAverage = 0.0 ! Initialise average water pressure of element IEl
              GPAverage = 0.0 ! Initialise average gas pressure of element IEl
              DPBVDAverage = 0.0
              
              ParticleIndex = GetParticleIndex(1, IEl, IPatch)
              CheckMaterialID = MaterialIDArray(ParticleIndex)
              IsMixedElement = .false.
              
              ParticlesVols = 0.0
              ParticleVol = 0.0

              if (CalParams%ApplyEmptyElements) then
                AverageWeight = 0.0
                if (NMaterialParticlesEle(IEl)>0) then
                  do IPart = 1, NElemPart
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
                    if (IsMaterialParticle(ParticleIndex)) then
                      AverageWeight = AverageWeight + Particles(ParticleIndex)%IntegrationWeight
                    end if
                  end do
                  AverageWeight = AverageWeight / NMaterialParticlesEle(IEl)
                else ! Only virtual particle(s)
                  AverageWeight = 1.0
                end if
              end if

              do IPart = 1, NElemPart ! Loop over particles
                ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
                
                if (MaterialIDArray(ParticleIndex)/= CheckMaterialID) then
                  IsMixedElement = .true.
                end if 
                
                if (CalParams%ApplyEmptyElements .and. (.not.IsMaterialParticle(ParticleIndex))) then
                  ParticleVol = AverageWeight ! Use average weight
                else  
                  ParticleVol = Particles(ParticleIndex)%IntegrationWeight
                end if
                ! Sum up stress components, water pressure and gas pressure
                do I = 1, NTENSOR
                  ParticleStress = SigmaEffArray(ParticleIndex,I) 
                  StressAverage(I) = StressAverage(I) + ParticleStress * ParticleVol
                end do
                
                WPAverage = WPAverage + Particles(ParticleIndex)%WaterPressure * ParticleVol
                DPBVDAverage = DPBVDAverage + Particles(ParticleIndex)%DBulkViscousPressure * ParticleVol
                GPAverage = GPAverage + Particles(ParticleIndex)%GasPressure * ParticleVol
                ParticlesVols = ParticlesVols + ParticleVol
              end do

              ! Take the average of each stress component, the water pressure and the gas pressure
              do J = 1, NTENSOR
                StressAverage(J) = StressAverage(J) / ParticlesVols
              end do
              WPAverage = WPAverage / ParticlesVols
              GPAverage = GPAverage / ParticlesVols
              DPBVDAverage = DPBVDAverage / ParticlesVols

              ! Assign the average stress, water pressure and gas pressure to particles
              do IPart = 1, NElemPart
                ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
                if (.not.IsMixedElement) then
                    SigmaEffArray(ParticleIndex, :) = StressAverage(:)
                    Particles(ParticleIndex)%WaterPressure = WPAverage
                    Particles(ParticleIndex)%DBulkViscousPressure = DPBVDAverage
                    Particles(ParticleIndex)%GasPressure = GPAverage
                end if
              end do

              do IPart = 1, NPartEle(IEl,IPatch) ! Loop over particles in IElement
                ParticleIndex = GetParticleIndex(IPart, IEl, IPatch)
 
                SigmaEffStep = SigmaEffArray(ParticleIndex,:) - SigmaEff0Array(ParticleIndex,:)
                call SetSigmaEffStep(Particles(ParticleIndex), SigmaEffStep)

                SigmaEff =SigmaEffArray(ParticleIndex,:)
                SigmaEff0Array(ParticleIndex,:)=SigmaEff

                Particles(ParticleIndex)%WaterPressure0 = Particles(ParticleIndex)%WaterPressure

                Particles(ParticleIndex)%GasPressure0 = Particles(ParticleIndex)%GasPressure

              end do
                end do ! elements
                
                end do ! patches

        end subroutine StressAndPorePressureSmoothening
        
        !---------------------------------------------------------------------------------------------------
        
        
        
         subroutine StressAndPorePressureSmoothening_4GaussPoints()
      !**********************************************************************
      !
      ! Function: Smoothing of stress and pore pressure in case of GP integration.
      !
      ! Implemented in the frame of the MPM project.
      !
      ! Note : This subroutine works only if NumberOfLayers = 1 !!! Temprary Solution
      !
      !**********************************************************************

      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: I, J, IEl, NElemPart, IPart, ParticleIndex
      integer(INTEGER_TYPE) :: CheckMaterialID, IAEl
      real(REAL_TYPE) :: ParticleStress
      real(REAL_TYPE) :: WPAverage, DPBVDAverage
      real(REAL_TYPE) :: GPAverage
      real(REAL_TYPE), dimension(NTENSOR) :: StressAverage
      logical :: IsMixedElement
      real(REAL_TYPE) :: ParticleVol, ParticlesVols, AverageWeight
      real(REAL_TYPE), dimension(NTENSOR) :: SigmaEffStep, SigmaEff
      
      integer(INTEGER_TYPE) :: IDim
      integer(INTEGER_TYPE) :: IGaussPoint
      
      integer(INTEGER_TYPE), dimension(4) :: IdentityMatrix
      
      real(REAL_TYPE) :: MeanEffStress
      real(REAL_TYPE) :: MeanEffStressAverage
      real(REAL_TYPE) :: ParticlesVolsTotalElement
      real(REAL_TYPE) :: MeanEffStressAverageElement
      real(REAL_TYPE) :: MeanEffArrayGP
      
      ! Multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary
      
      ! initialize temporary variables 
      IPatch_Temporary = 1
      
      ! initialize identity matrix
      IdentityMatrix = 0
      
      ! hardcode indentity matrix (3x3 in voight notation)
      ! 3D
      !IdentityMatrix(1) = 1
      !IdentityMatrix(2) = 1
      !IdentityMatrix(3) = 1
      !IdentityMatrix(4) = 0
      !IdentityMatrix(5) = 0
      !IdentityMatrix(6) = 0
      ! 2D --> only 2D for now
      !IdentityMatrix(1) = 1
      !IdentityMatrix(2) = 1
      !IdentityMatrix(3) = 0
      !IdentityMatrix(4) = 0

      if(.not.(NFORMULATION==1)) RETURN

      ! loop over active elements
      do IAEl = 1, Counters%NAEl(IPatch_Temporary)!Counters%NAEl
          
          ParticlesVolsTotalElement = 0.0
          MeanEffStressAverage = 0.0
          !MeanEffStressAverage 
          !StressAverage = 0.0
          
          ! obtain IEl global index for IAEl
          IEl = ActiveElement(IAEl,IPatch_Temporary)
         
          ! loop over each of the 4 gauss points
          do IGaussPoint = 1, Counters%NGaussPoints
         
              ! obtain the number of material points corresponding to the gauss point
              NElemPart = SubElementMPOrganization(IEl, IGaussPoint) !NPartEle(IEl)
            
             !StressAverage = 0.0 ! Initialise average stress of element IEl
             WPAverage = 0.0 ! Initialise average water pressure of element IEl
             GPAverage = 0.0 ! Initialise average gas pressure of element IEl
             DPBVDAverage = 0.0

             !ParticleIndex = GetParticleIndexInSubElement(IEl, ISubZone, 1) !GetParticleIndex(1, IEl)
             IsMixedElement = .false.

             ParticlesVols = 0.0
             ParticleVol = 0.0

             
             ! Apply Empty Elements---------------------------
         if (CalParams%ApplyEmptyElements) then
            AverageWeight = 0.0
            if (NMaterialParticlesEle(IEl)>0) then
               do IPart = 1, NElemPart
                  ParticleIndex = GetParticleIndex(IPart, IEl, IPatch_Temporary) ! Get the particle ID
                  if (IsMaterialParticle(ParticleIndex)) then
                     AverageWeight = AverageWeight + Particles(ParticleIndex)%IntegrationWeight
                  end if
               end do
               AverageWeight = AverageWeight / NMaterialParticlesEle(IEl)
            else ! Only virtual particle(s)
               AverageWeight = 1.0
            end if
         end if !---------------------------------------------------

             ! loop over particles corresponding to the gauss point
             do IPart = 1, NElemPart ! Loop over particles
                 
                 ! obtain the particle index corresponding to IPart
                 ParticleIndex = GetParticleIndexInSubElement(IEl, IGaussPoint, IPart)!GetParticleIndex(IPart, IEl) ! Get the particle ID

                 CheckMaterialID = MaterialIDArray(ParticleIndex)

                 if (MaterialIDArray(ParticleIndex)/= CheckMaterialID) then
                     IsMixedElement = .true.
                 end if

                 if (CalParams%ApplyEmptyElements .and. (.not.IsMaterialParticle(ParticleIndex))) then
                     ParticleVol = AverageWeight ! Use average weight
                 else
                     ! obtain the particle volume (i.e., integration weight in the physical domain)
                     ParticleVol = Particles(ParticleIndex)%IntegrationWeight
                 end if
            
            
            
            
                 ! Sum up stress solid components
                 do I = 1, NTENSOR
                     ! obtain the particle stress for ParticleIndex
                     !ParticleStress = SigmaEffArray(ParticleIndex,I)
                     !
                     !! sum up all the (stress*integration weight) for all the MPs in the 
                     !StressAverage(I) = StressAverage(I) + ParticleStress * ParticleVol
                 end do

                ! Sum up water pressure and gas pressure
                 WPAverage = WPAverage + Particles(ParticleIndex)%WaterPressure * ParticleVol ! water pressure 
                 DPBVDAverage = DPBVDAverage + Particles(ParticleIndex)%DBulkViscousPressure * ParticleVol
                 GPAverage = GPAverage + Particles(ParticleIndex)%GasPressure * ParticleVol ! gas pressure
            
            
                 ! sum up all the particle volumes
                 ParticlesVols = ParticlesVols + ParticleVol
         
             end do! Loop over particles

         ! Take the average of each stress component, the water pressure and the gas pressure
             if (ParticlesVols>0) then ! avoid dividing by zero
             
             !do J = 1, NTENSOR
             !StressAverage(J) = StressAverage(J) / ParticlesVols !solid
             !
             !end do
             
             WPAverage = WPAverage / ParticlesVols ! water
             GPAverage = GPAverage / ParticlesVols ! gas
             DPBVDAverage = DPBVDAverage / ParticlesVols ! bulk viscous pressure
             
             else 
            !StressAverage = 0.0
             WPAverage = 0.0
             GPAverage = 0.0
             DPBVDAverage = 0.0
             end if 
         
         

         ! we do not want to overwrite at this stage so I decided not include these here
         ! Assign the average stress, water pressure and gas pressure to particles
         !do IPart = 1, NElemPart
         !   ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
         !   if (.not.IsMixedElement) then
         !      SigmaEffArray(ParticleIndex, :) = StressAverage(:)
         !      Particles(ParticleIndex)%WaterPressure = WPAverage
         !      Particles(ParticleIndex)%DBulkViscousPressure = DPBVDAverage
         !      Particles(ParticleIndex)%GasPressure = GPAverage
         !   end if
         !end do

         ! if we update the stresses we need to do these again
         !do IPart = 1, NElemPart!NPartEle(IEl) ! Loop over particles in in subzone in IElement
         !   ParticleIndex = GetParticleIndexInSubElement(IEl, IGaussPoint, IPart)!GetParticleIndex(IPart, IEl)
         !
         !   SigmaEffStep = SigmaEffArray(ParticleIndex,:) - SigmaEff0Array(ParticleIndex,:)
         !   call SetSigmaEffStep(Particles(ParticleIndex), SigmaEffStep)
         !
         !   SigmaEff =SigmaEffArray(ParticleIndex,:)
         !   SigmaEff0Array(ParticleIndex,:)=SigmaEff
         !
         !   Particles(ParticleIndex)%WaterPressure0 = Particles(ParticleIndex)%WaterPressure
         !
         !   Particles(ParticleIndex)%GasPressure0 = Particles(ParticleIndex)%GasPressure
         !end do ! end loop over particles in subzones
         
         
         ! assign the average stress to the gaussian point
         !SigmaEffArrayGaussPoints(IEl, IGaussPoint, :) = StressAverage ! tensor 
         SigmaEffArrayGaussPointsWaterPressure(IEl, IGaussPoint) = WPAverage ! scalar
         SigmaEffArrayGaussPointsGasPressure(IEl, IGaussPoint) = GPAverage ! scalar
         SigmaEffArrayGaussPointsBulkViscousPressure(IEl, IGaussPoint) = DPBVDAverage ! scalar
         
         
         
         ! Mean effective stress correction for the entire element ---------------------------------------------------------------
         ! calculate the mean effective stress for the gaussian point p_GP
         !                      (Sigma_X                                  +       Sigma_Y                                 )/3
         !MeanEffStress = 0.0
         !do IDim = 1, NDIM 
         !!MeanEffStress = ( SigmaEffArrayGaussPoints(IEl, IGaussPoint, 1) + SigmaEffArrayGaussPoints(IEl, IGaussPoint, 2)  ) / 3.0
         !MeanEffStress = MeanEffStress + (SigmaEffArrayGaussPoints(IEl, IGaussPoint, IDim) )
         !end do 
         !MeanEffStress = MeanEffStress / 3.0
         !
         !
         !! Numerator = p_GP * Vol_GP
         !MeanEffStressAverage = MeanEffStressAverage + (MeanEffStress * ParticlesVols)
         !
         !
         !! find total particle volume in the element (summing accross all the subzones))
         !ParticlesVolsTotalElement = ParticlesVols + ParticlesVolsTotalElement
         
          end do ! loop over all gauss points
          
          
         ! obtain mean effective stress in the element
         ! if (ParticlesVolsTotalElement>0) then ! avoid dividing by zero
         !MeanEffStressAverageElement = MeanEffStressAverage / ParticlesVolsTotalElement
         ! else 
         !MeanEffStressAverageElement = 0.0
         !end if 
         
         !do IGaussPoint = 1, Counters%NGaussPoints ! loop over all gauss points
         !    !                               (Sigma_X                                  +       Sigma_Y                       )/3
         !    MeanEffArrayGP =  (SigmaEffArrayGaussPoints(IEl, IGaussPoint, 1)+SigmaEffArrayGaussPoints(IEl, IGaussPoint, 2) )/3.0!+SigmaEffArrayGaussPoints(IEl, IGaussPoint, 3) )
         !    
         !    ! subtract the effective stress from the gauss point and add in a stabilized mean effective stress for the ENTIRE element
         !    SigmaEffArrayGaussPoints(IEl, IGaussPoint, :) = SigmaEffArrayGaussPoints(IEl, IGaussPoint, :) - ((MeanEffArrayGP - MeanEffStressAverageElement)*IdentityMatrix)
         !end do 
         
         
      end do  ! loop over all active elements
      

   end subroutine StressAndPorePressureSmoothening_4GaussPoints
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !---------------------------------------------------------------------------------------------------
        
        subroutine StateParametersSmoothening()
        !**********************************************************************
        !
        !    Function: Smoothing of state parameters in case of GP integration.
        !
        ! Implemented in the frame of the MPM project.
        !
        ! Note : This subroutine works only if NumberOfLayers = 1 !! Temprary Solution
        !
        !**********************************************************************
        implicit none
      
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, NElemPart, IPart, ParticleIndex
          integer(INTEGER_TYPE) :: CheckMaterialID, IAEl
          real(REAL_TYPE) :: ParticleStateVariable
          real(REAL_TYPE) :: AverageHPStateVariables (2), &
                              AverageModifiedHPStateVariables(2), &
                              AverageHPIGStateVariables (7), &
                              AverageEpsP(NTENSOR), &
                              AverageSigmaPrin(NTENSOR)
          real(REAL_TYPE) :: AverageCohesionStSoft, &
                              AveragePhiStSoft, &
                              AveragePsiStSoft, &
                              AveragePP
          real(REAL_TYPE) :: AverageESMstatev(NSTATEVAR) !user defined model
          logical :: IsMixedElement
          real(REAL_TYPE) :: ParticleVol, ParticlesVols, AverageWeight
          
          integer(INTEGER_TYPE) :: IPatch
     
         if(.not.(NFORMULATION==1)) RETURN
          
         !if (.not.IsMPMComputation()) RETURN ! if FEM no need for mapping...already one particle
         
          !do IAEl = 1, Counters%NAEl
              do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            IEl = ActiveElement(IAEl, IPatch)
              NElemPart = NPartEle(IEl,IPatch)
              AverageHPStateVariables         = 0.0
              AverageModifiedHPStateVariables = 0.0
              AverageHPIGStateVariables       = 0.0
              AverageEpsP = 0.0
              AverageSigmaPrin = 0.0
              AverageCohesionStSoft = 0.0
              AveragePhiStSoft = 0.0
              AveragePsiStSoft = 0.0
              AveragePP = 0.0
              
              AverageESMstatev = 0.0
              
              ParticleIndex = GetParticleIndex(1, IEl, IPatch)
              CheckMaterialID = MaterialIDArray(ParticleIndex) ! the same as material ID
              IsMixedElement = .false.
              
              ParticlesVols = 0.0
              ParticleVol = 0.0

              if (CalParams%ApplyEmptyElements) then
                AverageWeight = 0.0
                if (NMaterialParticlesEle(IEl)>0) then
                  do IPart = 1, NElemPart
                    ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
                    if (IsMaterialParticle(ParticleIndex)) then
                      AverageWeight = AverageWeight + Particles(ParticleIndex)%IntegrationWeight
                    end if
                  end do
                  AverageWeight = AverageWeight / NMaterialParticlesEle(IEl)
                else ! Only virtual particle(s)
                  AverageWeight = 1.0
                end if
              end if

              do IPart = 1, NElemPart ! Loop over particles
                ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
                
                if (MaterialIDArray(ParticleIndex) /= CheckMaterialID) then
                  IsMixedElement = .true.
                  goto 100 ! skip this element 
                end if 
              
                if (CalParams%ApplyEmptyElements .and. (.not.IsMaterialParticle(ParticleIndex))) then
                  ParticleVol = AverageWeight ! Use average weight
                else  
                  ParticleVol = Particles(ParticleIndex)%IntegrationWeight
                end if

               ! Sum up state variables for MCStrainSoftening model
                do I = 1, NTENSOR
                  ParticleStateVariable = GetEpsPI(Particles(ParticleIndex), I)  
                  AverageEpsP(I) = AverageEpsP(I) + ParticleStateVariable * ParticleVol
                  ParticleStateVariable = GetSigmaPrinI(Particles(ParticleIndex), I)  
                  AverageSigmaPrin(I) = AverageSigmaPrin(I) + ParticleStateVariable * ParticleVol
                end do
                AverageCohesionStSoft = AverageCohesionStSoft + Particles(ParticleIndex)%CohesionStSoft * ParticleVol
                AveragePhiStSoft = AveragePhiStSoft + Particles(ParticleIndex)%PhiStSoft * ParticleVol
                AveragePsiStSoft = AveragePsiStSoft + Particles(ParticleIndex)%PsiStSoft * ParticleVol
                
              ! Sum up state variables for MCC model
                AveragePP = AveragePP + Particles(ParticleIndex)%PP * ParticleVol
                

                do I = 1, 2  ! HP model 
                  ParticleStateVariable = GetHPStateVariablesI(Particles(ParticleIndex), I)
                  AverageHPStateVariables(I) = AverageHPStateVariables(I) + ParticleStateVariable * ParticleVol
                end do

                do I = 1, 7  ! HPIG model 
                  ParticleStateVariable = GetHPIGStateVariablesI(Particles(ParticleIndex), I)
                  AverageHPIGStateVariables(I) = AverageHPIGStateVariables(I) + ParticleStateVariable * ParticleVol
                end do

                ! Modified HP model
                do I = 1, 2
                  ParticleStateVariable = GetModifiedHPStateVariablesI(Particles(ParticleIndex), I)
                  AverageModifiedHPStateVariables(I) = AverageModifiedHPStateVariables(I) + ParticleStateVariable * ParticleVol
                end do

                AverageESMstatev=  AverageESMstatev +  ESMstatevArray(ParticleIndex,:) * ParticleVol
   
                ParticlesVols = ParticlesVols + ParticleVol
              end do ! loop over particles


              ! Take the average of each state parameter

              ! HP model
              AverageHPStateVariables = AverageHPStateVariables / ParticlesVols

              ! HPIG model
              AverageHPIGStateVariables = AverageHPIGStateVariables / ParticlesVols
          

              ! modified HP model
              AverageModifiedHPStateVariables = AverageModifiedHPStateVariables / ParticlesVols

              AverageESMstatev =  AverageESMstatev / ParticlesVols
 
              ! Strain softening model
              AverageEpsP = AverageEpsP / ParticlesVols
              AverageSigmaPrin = AverageSigmaPrin / ParticlesVols

              AverageCohesionStSoft = AverageCohesionStSoft / ParticlesVols
              AveragePhiStSoft = AveragePhiStSoft / ParticlesVols
              AveragePsiStSoft = AveragePsiStSoft / ParticlesVols
              
              ! MCC model
              AveragePP = AveragePP / ParticlesVols
              
              ! Assign the average state parameters to particles
              if (.not.IsMixedElement) then ! double check for mixed element
                do IPart = 1, NElemPart
                  ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
            
                  call SetHPStateVariables(Particles(ParticleIndex), AverageHPStateVariables)
                  call SetHPIGStateVariables(Particles(ParticleIndex), AverageHPIGStateVariables)
                  call SetModifiedHPStateVariables(Particles(ParticleIndex), AverageModifiedHPStateVariables)
     
                   ! Strain softening model
                   call SetEpsP(Particles(ParticleIndex), AverageEpsP)
                   call SetSigmaPrin(Particles(ParticleIndex), AverageSigmaPrin)
                   Particles(ParticleIndex)%CohesionStSoft = AverageCohesionStSoft
                   Particles(ParticleIndex)%PhiStSoft = AveragePhiStSoft
                   Particles(ParticleIndex)%PsiStSoft = AveragePsiStSoft
                   
                   ! MCC model
                   Particles(ParticleIndex)%PP = AveragePP
                   !user-defined model
                   ESMstatevArray(ParticleIndex,:) = AverageESMstatev
                end do
              end if
                         
100             end do ! elements
                
                end do ! patches

          end subroutine StateParametersSmoothening
                
        subroutine MapVeloToParticle(NodalVelocities)
                              
        !**********************************************************************
        !
        !    Function:  To update particles total velocities and Accelerations
        !
        !
        !    NodalVelocities : Nodal accelerations vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(in) :: NodalVelocities
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IPart, INode, IDof, ParticleIndex, NodeID
          real(REAL_TYPE), dimension (NVECTOR) :: ParticleIncrementalVelocity
          real(REAL_TYPE), dimension (NVECTOR) :: ParIncVel
          
          ! Multipatch variables 
          !integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE) :: IPatch

          do IPatch = 1, Counters%NPatches ! loop over patches
          do IEl = 1, Counters%NEl(IPatch)!Counters%NEl   ! Loop over all elements
              do IPart = 1, NPartEle(IEl,IPatch) ! Loop over all particles in element
                
                ParticleIndex = GetParticleIndex(IPart, IEl, IPatch) ! Get the particle ID
                ParticleIncrementalVelocity = 0.0
                ParIncVel = 0.0
                
                  do INode = 1, ELEMENTNODES

                    NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(ElementConnectivities(INode, IEl, IPatch),IPatch)    ! Global node ID
                    IDof = ReducedDof(NodeID)

                    do I = 1, NVECTOR
                      ParIncVel(I) = ShapeValuesArray(ParticleIndex,INode) * NodalVelocities(IDof+1, 1)
                      ParticleIncrementalVelocity(I) = ParticleIncrementalVelocity(I) + ParIncVel(I)
                    end do  
               
               end do !Loop over nodes
                   
                 VelocityArray(ParticleIndex,:) = ParticleIncrementalVelocity
     
            end do !Loop over particles
          end do !elements
          end do ! patches
          
         end subroutine MapVeloToParticle

        subroutine GetNodalVelocityFromNodalMomentumConv(Momentum)!, Momentum_VolLockSmooth)
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Momentum
          
          !real(REAL_TYPE), dimension(Counters%N_VolLockSmooth,Counters%nEntity), intent(in) :: Momentum_VolLockSmooth
          
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J

         do IDOF = 1, Counters%N                                       !Loop over all degrees of freedom
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
                if(LumpedMassDry(IDOF,J)/=0) then
                    
                    ! According to Tran and Solowski (2020), this time intergration to obtain velocity was done at the nodal level...
                    ! However, in the FLIP approach, we do not do this! <-- why is it unstable 
                    if ( (CalParams%ApplyGeneralizedAlphaTimeScheme%Switch== 1) .and. (CalParams%TimeStep /=1) ) then  
                        
                        TotalVelocitySoil(IDOF, J) = TotalVelocitySoil(IDOF, J) + &
                         ( CalParams%TimeIncrement * (CalParams%ApplyGeneralizedAlphaTimeScheme%Gamma*AccelerationSoil(IDOF, J) + &
                         (1-CalParams%ApplyGeneralizedAlphaTimeScheme%Gamma)*AccelerationSoilGeneralizedAlpha(IDOF, J)) ) * PBoundary(IDOF)
                    else 
                    
                        TotalVelocitySoil(IDOF,J) = ( Momentum(IDOF,J) / LumpedMassDry(IDOF,J) ) *  &
                                                PBoundary(IDOF)
                    
                    end if 
                    
                else
                    TotalVelocitySoil(IDOF,J) = 0.0
                end if
            end do
         end do
         
         
         
         ! solve for the lower order velocity field --> volumetric locking mitigation
         !if (CalParams%ApplyBBar) then 
         !
         !    do IDOF = 1, Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
         !        do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
         !            if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
         !                TotalVelocitySoil_VolLockSmooth(IDOF,J) = ( Momentum_VolLockSmooth(IDOF,J) / LumpedMassDry_VolLockSmooth(IDOF,J) ) *  &
         !                    PBoundary_VolLockSmooth(IDOF)
         !            else
         !                TotalVelocitySoil_VolLockSmooth(IDOF,J) = 0.0
         !            end if
         !        end do
         !    end do
         !    
         !end if 
         
          
        end subroutine GetNodalVelocityFromNodalMomentumConv
        
        
        
        
        !subroutine GetNodalVolumetricStrain_VolLockSmooth()
        !!**********************************************************************
        !!
        !!    Function:  To calculate the nodal velocities from nodal mass and momentum
        !!
        !!    Momentum : Nodal momentum vector
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !
        !integer(INTEGER_TYPE) :: IDOF, J
        !
        !
        !do IDOF = 1, Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        !
        !    do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
        !    
        !        if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
        !        
        !            NodalVolumetricStrain_VolLockSmooth(IDOF,J) = ( WeightedNodalVolumetricStrainIncrement_Numerator(IDOF,J) / LumpedMassDry_VolLockSmooth(IDOF,J) ) !*  &
        !                                        !PBoundary(IDOF)
        !        
        !        else
        !        
        !            NodalVolumetricStrain_VolLockSmooth(IDOF,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
        !        
        !        end if
        !    
        !    end do
        !  
        !end do
        !  
        !
        !end subroutine GetNodalVolumetricStrain_VolLockSmooth
        
        
        !subroutine GetNodalFJacDet_VolLockSmooth()
        !!**********************************************************************
        !!
        !!    Function:  To calculate the nodal velocities from nodal mass and momentum
        !!
        !!    Momentum : Nodal momentum vector
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !
        !integer(INTEGER_TYPE) :: IDOF, J
        !
        !
        !do IDOF = 1, Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        !
        !    do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
        !    
        !        if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
        !        
        !            NodalFJacDet_VolLockSmooth(IDOF,J) = ( WeightedFJacDet_Numerator(IDOF,J) / LumpedMassDry_VolLockSmooth(IDOF,J) ) !*  &
        !                                        !PBoundary(IDOF)
        !        
        !        else
        !        
        !            NodalFJacDet_VolLockSmooth(IDOF,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
        !        
        !        end if
        !    
        !    end do
        !  
        !end do
        !  
        !
        !end subroutine GetNodalFJacDet_VolLockSmooth
        
        !----------------------------------------------------------------------------------
        
        subroutine GetNodalStrainIncrement_VolLockSmooth()
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        IPatch_Temporary = 1
        
        do INode = 1, Counters%NodTot(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                do ITens = 1, NTENSOR
                
                    if(LumpedMassDry(IDOF,J)/=0) then
                
                    NodalStrainIncrement_VolLockSmooth(ITens,INode,J) = ( WeightedStrainIncrementTensor_Numerator(ITens,INode,J) &
                                                                            / LumpedMassDry(IDOF,J) ) !*  &
                                                !PBoundary(IDOF)
                
                    else
                
                    NodalStrainIncrement_VolLockSmooth(ITens,INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
                
                
                    end if
                
                end do 
            
            end do
          
        end do
          
        
        end subroutine GetNodalStrainIncrement_VolLockSmooth
        
        !----------------------------------------------------------------------------------
        
        subroutine GetNodalStressIncrement_VolLockSmooth()
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        IPatch_Temporary = 1
        
        do INode = 1, Counters%NodTot(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                do ITens = 1, NTENSOR
                
                    if(LumpedMassDry(IDOF,J)/=0) then
                
                    NodalStressIncrement_VolLockSmooth(ITens,INode,J) = ( WeightedStressIncrementTensor_Numerator(ITens,INode,J) &
                                                                            / LumpedMassDry(IDOF,J) ) !*  &
                                                !PBoundary(IDOF)
                
                    else
                
                    NodalStressIncrement_VolLockSmooth(ITens,INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
                
                
                    end if
                
                end do 
            
            end do
          
        end do
          
        
        end subroutine GetNodalStressIncrement_VolLockSmooth
        
        
        
        !----------------------------------------------------------------------------------
        
        subroutine GetNodalWP_VolLockSmooth()
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        IPatch_Temporary = 1
        
        do INode = 1, Counters%NodTot(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                !do ITens = 1, NTENSOR
                
                    if(LumpedMassDry(IDOF,J)/=0) then
                
                    NodalWP_VolLockSmooth(INode,J) = ( WeightedWP_Numerator(INode,J) & !ITens, ITens,
                                                                            / LumpedMassDry(IDOF,J) ) !*  &
                                                !PBoundary(IDOF)
                
                    else
                
                    NodalWP_VolLockSmooth(INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0 !ITens,
                
                
                    end if
                
                !end do 
            
            end do
          
        end do
          
        
        end subroutine GetNodalWP_VolLockSmooth
        
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------
        
        
        
        
        
        
        
        subroutine GetNodalWP_VolLockSmooth_LowerOrder()
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        IPatch_Temporary = 1
        
        do INode = 1, NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof_VolLockSmooth(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                !do ITens = 1, NTENSOR
                
                    if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
                
                    NodalWP_VolLockSmooth(INode,J) = ( WeightedWP_Numerator(INode,J) & !ITens, ITens,
                                                                            / LumpedMassDry_VolLockSmooth(IDOF,J) ) !*  &
                                                !PBoundary(IDOF)
                
                    else
                
                    NodalWP_VolLockSmooth(INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0 !ITens,
                    
                    
                    end if
                
                !end do 
            
            end do
          
        end do
          
        
        end subroutine GetNodalWP_VolLockSmooth_LowerOrder
        
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------
        
        
        
        subroutine GetNodalStrainIncrement_VolLockSmooth_LowerOrder()
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        IPatch_Temporary = 1
        
        do INode = 1, NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof_VolLockSmooth(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                do ITens = 1, NTENSOR
                
                    if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
                
                    NodalStrainIncrement_VolLockSmooth(ITens,INode,J) = ( WeightedStrainIncrementTensor_Numerator(ITens,INode,J) &
                                                                            / LumpedMassDry_VolLockSmooth(IDOF,J) ) !*  &
                                                !PBoundary(IDOF)
                
                    else
                
                    NodalStrainIncrement_VolLockSmooth(ITens,INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
                
                
                    end if
                
                end do 
            
            end do
          
        end do
          
        
        end subroutine GetNodalStrainIncrement_VolLockSmooth_LowerOrder
        
        
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------
        
        
        
        subroutine GetNodalStressIncrement_VolLockSmooth_LowerOrder()
        !**********************************************************************
        !
        !    Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
        
        integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        integer(INTEGER_TYPE) :: INode, ITens
        
        IPatch_Temporary = 1
        
        do INode = 1, NURBS%NControlPoints_VolLockSmooth(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        
            do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
            
                IDOF = ReducedDof_VolLockSmooth(INode)+1 ! obtain the x dof so that we can find the mass. 
                                           ! It does not matter what dof you choose (here x-direction), 
                                           ! because mass is scalar. 
                
                do ITens = 1, NTENSOR
                
                    if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
                
                    NodalStressIncrement_VolLockSmooth(ITens,INode,J) = ( WeightedStrainIncrementTensor_Numerator(ITens,INode,J) &
                                                                            / LumpedMassDry_VolLockSmooth(IDOF,J) ) !*  &
                                                !PBoundary(IDOF)
                
                    else
                
                    NodalStressIncrement_VolLockSmooth(ITens,INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
                
                
                    end if
                
                end do 
            
            end do
          
        end do
          
        
        end subroutine GetNodalStressIncrement_VolLockSmooth_LowerOrder
        
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------------------
        
        !-----------------------------------------------------------------
        !subroutine GetNodalStress_VolLockSmooth()
        !!**********************************************************************
        !!
        !!    Function:  To calculate the nodal velocities from nodal mass and momentum
        !!
        !!    Momentum : Nodal momentum vector
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !
        !integer(INTEGER_TYPE) :: IDOF, J, IPatch_Temporary
        !integer(INTEGER_TYPE) :: INode, ITens
        !
        !IPatch_Temporary = 1
        !
        !do INode = 1, NControlPoints_VolLockSmooth(IPatch_Temporary) !Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        !
        !    do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
        !    
        !        IDOF = ReducedDof_VolLockSmooth(INode)+1 ! obtain the x dof so that we can find the mass. 
        !                                   ! It does not matter what dof you choose (here x-direction), 
        !                                   ! because mass is scalar. 
        !        
        !        do ITens = 1, NTENSOR
        !        
        !            if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
        !        
        !            NodalStress_VolLockSmooth(ITens,INode,J) =  WeightedStressTensor_Numerator(ITens,INode,J) &
        !                                                         / LumpedMassDry_VolLockSmooth(IDOF,J)  !*  &
        !                                        !PBoundary(IDOF)
        !        
        !            else
        !        
        !            NodalStress_VolLockSmooth(ITens,INode,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
        !        
        !            end if
        !        
        !        end do 
        !    
        !    end do
        !  
        !end do
        !  
        !
        !end subroutine GetNodalStress_VolLockSmooth
        
        
        
        !-----------------------------------------------------------------
        
        
        
        !subroutine GetNodalVolumetricStressGradient_VolLockSmooth()
        !!**********************************************************************
        !!
        !!    Function:  To calculate the nodal velocities from nodal mass and momentum
        !!
        !!    Momentum : Nodal momentum vector
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !  !real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Momentum
        !  ! Local variables
        !  integer(INTEGER_TYPE) :: IDOF, J
        !
        ! do IDOF = 1, Counters%N_VolLockSmooth                                       !Loop over all degrees of freedom
        !    do J = 1, Counters%nEntity !loop through all entities       CC - added dimension nEntity (J)
        !        if(LumpedMassDry_VolLockSmooth(IDOF,J)/=0) then
        !            NodalVolumetricStressGradient_VolLockSmooth(IDOF,J) = WeightedNodalVolumetricStressGradient_Numerator(IDOF,J) / LumpedMassDry_VolLockSmooth(IDOF,J) !*  &
        !                                        !PBoundary(IDOF)
        !        else
        !            NodalVolumetricStressGradient_VolLockSmooth(IDOF,J) = 0.0!TotalVelocitySoil(IDOF,J) = 0.0
        !        end if
        !    end do
        !  end do
        !  
        !end subroutine GetNodalVolumetricStressGradient_VolLockSmooth
        
        
        
        
        
        

         subroutine UpdateParticleDisplacementsWater(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (water)
        !
        !     DUTot : Obtained incremental displacements (water)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
          ! Multipatch variables
          !integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: PatchID
          
          do IParticle = 1, Counters%NParticles ! Loop over particles
            
            
            PatchID = PatchIDArray(IParticle)  
              
            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, PatchID)), PatchID) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle displacement
            Particles(IParticle)%UW = Particles(IParticle)%UW + ParticleDisplacement
        
         end do ! particles
          
         end subroutine UpdateParticleDisplacementsWater
         
         subroutine UpdateParticleDisplacementsGas(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (gas)
        !
        !     DUTot : Obtained incremental displacements (gas)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: PatchID
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

               
            PatchID = PatchIDArray(IParticle)    
               
            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle)
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, PatchID)), IPatch_Temporary) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle displacement
            Particles(IParticle)%UG = Particles(IParticle)%UG + ParticleDisplacement
        
         end do ! particles
          
         end subroutine UpdateParticleDisplacementsGas
         
          subroutine UpdateParticleVelocityWater(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (water)
        !
        !     DUTot : Obtained incremental displacements (water)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
     
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: IPatch
          integer(INTEGER_TYPE) :: PatchID
          
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

            PatchID = PatchIDArray(IParticle)    
               
            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, IPatch_Temporary)), IPatch_Temporary) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle velocity
            
            VelocityWaterArray(IParticle,:) = ParticleDisplacement 
            
                   
         end do ! particles
          
         end subroutine UpdateParticleVelocityWater
         
         subroutine UpdateParticleVelocityGas(DUTot)
        !**********************************************************************
        !
        !    Function:  Calculate particle displacements from nodal displacements (gas)
        !
        !     DUTot : Obtained incremental displacements (gas)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none        
     
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: DUTot
          ! Local variables
          integer(INTEGER_TYPE) :: IDof, INode, iEntity
          integer(INTEGER_TYPE) :: NodeID, DofID, IParticle
          integer(INTEGER_TYPE) :: ElementID
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleDisplacement
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE) :: PatchID
          
           do IParticle = 1, Counters%NParticles ! Loop over particles

            PatchID = PatchIDArray(IParticle)    
               
            ElementID = ElementIDArray(IParticle)
          
            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IParticle) 
            else
              iEntity = 1
            end if 
              
            ParticleDisplacement = 0.0
                    
            do IDof = 1, NVECTOR ! Loop over all degrees of freedom (x, y, z)
              do INode = 1, ELEMENTNODES  ! Loop over all nodes in element               
                NodeID = Multipatch_Connecting_Local_To_Global_ControlPoints(iabs(ElementConnectivities(INode, ElementID, IPatch_Temporary)), IPatch_Temporary) 
                DofID = ReducedDof(NodeID) + IDof
                ParticleDisplacement(IDof) = ParticleDisplacement(IDof) + ShapeValuesArray(IParticle,INode) * DUTot(DofID,iEntity)                             
              end do                                             
            end do 
                          
            ! Total particle velocity
            
            VelocityGasArray(IParticle,:) = ParticleDisplacement
                   
         end do ! particles
          
         end subroutine UpdateParticleVelocityGas
         
         
 
         subroutine InitialiseConvectivePhaseData()   
        !**********************************************************************
        !    SUBROUTINE: InitialiseConvectivePhaseData
        ! 
        !    DESCRIPTION: 
        !>   Initialises the arrays globally defined in this module
        !
        !**********************************************************************     
           implicit none
         
           ! Local variable
           integer(INTEGER_TYPE) :: IError
           
           call DestroyConvectivePhaseData()
           
           allocate(TemporaryMappingVector(Counters%N), stat = IError)
           TemporaryMappingVector = 0.0
         
         end subroutine InitialiseConvectivePhaseData
         
   
         subroutine DestroyConvectivePhaseData()    
         !**********************************************************************
         !    SUBROUTINE: DestroyConvectivePhaseData
         ! 
         !    DESCRIPTION: 
         !>   Destroys the arrays globally defined in this module
         !
         !**********************************************************************   
           implicit none
         
           ! Local variable
           integer(INTEGER_TYPE) :: IError
           
           if (allocated(TemporaryMappingVector) ) then
             deallocate(TemporaryMappingVector, stat = IError)
           end if
         
    end subroutine DestroyConvectivePhaseData

 
      subroutine CalculateCriticalTimeStep()
      !**********************************************************************
      !
      !  Function: calculate critical time step size
      !            1-phase: based on CFL condition
      !            2-phase: based on MSc Thesis Mieremet (TU Delft, 2015)
      !            mass scaling is taken into account
      !
      !**********************************************************************
      
      implicit none
      
        !local variable
        integer(INTEGER_TYPE) :: IAEl, IEl, NParticles, ParticleIndex, IParticle, J, MaterialIndex
        real(REAL_TYPE) :: WaveSpeed, MaxWaveSpeed, TimeStep, MinTimeStep, CourantFactor
        real(REAL_TYPE) :: ConstDensityLiquid, ConstDensitySolid, BulkW, Nu, G
        real(REAL_TYPE) :: WaveSpeed_c1, WaveSpeed_c2, WaveSpeed_c3
        real(REAL_TYPE) :: rho_tilde, Timestep_3
        real(REAL_TYPE) :: TimeIncrementNew
        real(REAL_TYPE) :: Chi = 0.0
        
        ! Multipatch variables 
        integer(INTEGER_TYPE) :: IPatch_Temporary = 1
        integer(INTEGER_TYPE) :: IPatch
        
        !--------------------------------------------------------------------
        !-------------------- 1 layer formulation ---------------------------
        !--------------------------------------------------------------------
        
        if((NFORMULATION==1).or. &
          ((NFORMULATION==2).and.(CalParams%NumberOfPhases==1)))then
          TimeStep = 0.0
          MinTimeStep = huge(MinTimeStep)
          CourantFactor = CalParams%CourantNumber
         
          !do IAEl = 1, Counters%NAEl ! loop over all active elements
              
              do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
              
            IEl = ActiveElement(IAEl, IPatch)
            NParticles = NPartEle(IEl,IPatch) 
            
            Chi = 0.0
            WaveSpeed = 0.0
            MaxWaveSpeed = 0.0
            
            do IParticle = 1, NParticles ! loop over all particles in element
              ParticleIndex = GetParticleIndex(IParticle, IEl, IPatch)
			  MaterialIndex = MaterialIDArray(ParticleIndex)
			  if (.not.MatParams(MaterialIndex)%MaterialModel==ESM_RIGID_BODY) then ! do not compute wavespeed for a rigid body
              call GetWaveSpeedForTimeStepSize(ParticleIndex, IEl, WaveSpeed)
              
              if (IsNan(WaveSpeed).or.(abs(WaveSpeed) > huge(WaveSpeed))) then
                call GiveError('WaveSpeed is not a number or infinity')
              end if
              
              if (WaveSpeed>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed
			  end if
            end do ! end loop over all particles in element
            
            if (MaxWaveSpeed>(0.0)) then
              TimeStep = ElementLMin(IEl,IPatch) / MaxWaveSpeed !hardcoded
              
              if (CalParams%ApplyBulkViscosityDamping) then
                Chi = CalParams%BulkViscosityDamping1 - CalParams%BulkViscosityDamping2**2 * TimeStep * RateVolStrain(IEl)
                TimeStep = (sqrt(1 + Chi**2) - Chi) * TimeStep
              end if
              
			end if
			
			if (TimeStep>0.0) then
				if (TimeStep<MinTimeStep) MinTimeStep = TimeStep
			end if
               
                end do ! end loop over all active elements
                
                end do ! patches
          
        end if ! end 1 layer formulation
        
        !------------------------------------------------------------------------- 
        !--------------- 2 layer formulation and 2 phases ------------------------
        !-------------------------------------------------------------------------
        
        if ((NFORMULATION==2).and.(CalParams%NumberOfPhases==2)) then
            
          !--- Calculate  Soil and Water parameters  ----
          ConstDensityLiquid = 0.0
          ConstDensitySolid = 0.0
          BulkW  = 0.0
          Nu = 0.0
          G = 0.0
                
          do J = 1, Counters%NLayers ! loop over all materials
            if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
              ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
              ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
              BulkW  = MatParams(J)%BulkModulusLiquid
              Nu  = MatParams(J)%PoissonRatio
              G = MatParams(J)%ShearModulus
            else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
              ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
              BulkW  = MatParams(J)%BulkModulusLiquid 
            else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
              ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
              Nu  = MatParams(J)%PoissonRatio
              G = MatParams(J)%ShearModulus
            end if
          end do 
             
          TimeStep = 0.0
          MinTimeStep = huge(MinTimeStep)
          CourantFactor = CalParams%CourantNumber

          !do IAEl = 1, Counters%NAEl ! loop over all active elements
              do IPatch = 1, Counters%NPatches ! Loop over patches
                do IAEl = 1, Counters%NAEl(IPatch)!Counters%NEl ! Loop over all elements 
            IEl = ActiveElement(IAEl, IPatch)
            NParticles = NPartEle(IEl,IPatch)
            
            WaveSpeed_c1 = 0.0
            WaveSpeed_c2 = 0.0
            WaveSpeed_c3 = 0.0
            MaxWaveSpeed = 0.0
            rho_tilde = 0.0
            Chi = 0.0

            do IParticle = 1, NParticles ! loop over all particles in element
              ParticleIndex = GetParticleIndex(IParticle, IEl, IPatch)
              call GetWaveSpeed2LayerForm (ParticleIndex, IEl, ConstDensityLiquid, ConstDensitySolid, &
                                           BulkW, Nu, G, WaveSpeed_c1, WaveSpeed_c2)
              if (IsNan(WaveSpeed_c1).or.(abs(WaveSpeed_c1) > huge(WaveSpeed_c1))) then
                call GiveError('WaveSpeed_c1 is not a number or infinity'       //NEW_LINE('A')// &
                               'ParticleIndex: ' // trim(String(ParticleIndex)) //NEW_LINE('A')// &
                               'Element      : ' // trim(String(IEl)))
              end if

              if (IsNan(WaveSpeed_c2).or.(abs(WaveSpeed_c2) > huge(WaveSpeed_c2))) then
                call GiveError('WaveSpeed_c2 is not a number or infinity'       //NEW_LINE('A')// &
                               'ParticleIndex: ' // trim(String(ParticleIndex)) //NEW_LINE('A')// &
                               'Element      : ' // trim(String(IEl)))
              end if

              if(WaveSpeed_c1>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed_c1
              if(WaveSpeed_c2>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed_c2

              ! ---------- Permeability criterion ------------!
              if((TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID).and. &
                 (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid))then

                rho_tilde = (Particles(ParticleIndex)%EffConcentrationRatioSolid * ConstDensitySolid +  &
                            Particles(ParticleIndex)%EffPorosity * ConstDensityLiquid + &
                            (1.0 / Particles(ParticleIndex)%EffPorosity - 2.0) * ConstDensityLiquid)*CalParams%ScalingMassFactor

                Timestep_3 = 2.0 * rho_tilde * Particles(ParticleIndex)%Conductivity /  &
                             (ConstDensityLiquid * CalParams%GravityData%GAccel)

                WaveSpeed_c3 = ElementLMin(IEl,IPatch_Temporary) / Timestep_3

                if (IsNan(WaveSpeed_c3).or.(abs(WaveSpeed_c3) > huge(WaveSpeed_c3))) then
                  call GiveError('WaveSpeed_c3 is not a number or infinity'       //NEW_LINE('A')// &
                                 'ParticleIndex: ' // trim(String(ParticleIndex)) //NEW_LINE('A')// &
                                 'Element      : ' // trim(String(IEl)))

                end if

                if(WaveSpeed_c3>MaxWaveSpeed) MaxWaveSpeed = WaveSpeed_c3

              end if
            end do ! end loop over all particles in element
            
            if (MaxWaveSpeed>(0.0)) then
               TimeStep = ElementLMin(IEl,IPatch) / MaxWaveSpeed
               if (CalParams%ApplyBulkViscosityDamping) then
                 Chi = CalParams%BulkViscosityDamping1 - CalParams%BulkViscosityDamping2**2 * TimeStep * RateVolStrain(IEl)
                 TimeStep = (sqrt(1 + Chi**2) - Chi) * TimeStep
               end if
            end if

            if (TimeStep<MinTimeStep) MinTimeStep = TimeStep

                end do ! end loop over all active elements
                
                end do ! patches

        end if ! 2 layer formulation and 2 phases

        TimeIncrementNew = MinTimeStep * CourantFactor

        if (.not.CalParams%ApplyQuasiStatic) then
          if ((CalParams%TotalRealTime + TimeIncrementNew)>CalParams%TotalTime) then
            TimeIncrementNew = CalParams%TotalTime - CalParams%TotalRealTime
          end if
        end if

        ! hardcoded minimum time step size, to prevent time steps smaller than 1d-10 s
        if (TimeIncrementNew <= 1d-10) then
          TimeIncrementNew = 1d-10
        end if

        if (CalParams%TimeStep == 1 .or. (CalParams%IStep == 1.and.CalParams%TimeStep == 2)) then
          ! only for the first time step
          if (abs(TimeIncrementNew - CalParams%TimeIncrement) > SMALL) then
            call GiveMessage(' MaxWaveSpeed : ' // trim(String(MaxWaveSpeed)))
            call GiveMessage(' MinTimeStep  : ' // trim(String(MinTimeStep)))
            call GiveMessage(' TimeIncrement: ' // trim(String(TimeIncrementNew)))
          endif
        endif

        CalParams%TimeIncrement = TimeIncrementNew ! hardcoded ! Abdelrahman Alsardi !0.0001 

      end subroutine CalculateCriticalTimeStep




      subroutine GetWaveSpeed2LayerForm(IParticle, IEl, ConstDensityLiquid, ConstDensitySolid, &
                                        BulkW, Nu, G, WaveSpeed_c1, WaveSpeed_c2)
      !**********************************************************************
      !
      !  Function: calculate wave speed for particle IParticle in element IEl
      !
      !********************************************************************** 
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IParticle, IEl
        real(REAL_TYPE), intent(in) :: ConstDensityLiquid, ConstDensitySolid, BulkW, Nu, G
        real(REAL_TYPE), intent(inout) :: WaveSpeed_c1, WaveSpeed_c2

        ! local variables
        real(REAL_TYPE) :: Eoed_Skel, Eoed, rho_sat, Beta_s

        WaveSpeed_c1 = 0.0
        WaveSpeed_c2 = 0.0
        rho_sat = 0.0
        Beta_s = 0.0

        !!=========== ELEMENT CONTAINS = SOLID+LIQUID =====================!!
        if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID) then

          if (MaterialPointTypeArray(IParticle)==MaterialPointTypeSolid) then !this takes a lot of time!

            if (Particles(IParticle)%PhaseStatus==PhaseStatusSOLID) then

              Eoed_Skel = (2 * G * (1 - Nu)) / (1 - 2 * Nu) ! kPa  ! Eoed

              ! Assumed to be fully saturated :: %EffPorosity = %ConcRatioLiquid
              Eoed = Eoed_Skel + BulkW / Particles(IParticle)%EffPorosity
              rho_sat = Particles(IParticle)%EffConcentrationRatioSolid * ConstDensitySolid +  &
                        Particles(IParticle)%EffPorosity * ConstDensityLiquid

              Beta_s = sqrt((Particles(IParticle)%EffPorosity * Eoed_Skel / BulkW) /  &
                            (Particles(IParticle)%EffConcentrationRatioSolid +  &
                             Particles(IParticle)%EffPorosity * Eoed_Skel / BulkW))
            end if

            if (Particles(IParticle)%PhaseStatus==PhaseStatusLIQUID) then
              Beta_s = 1.0
            end if

          end if ! Solid Material point

          if (MaterialPointTypeArray(IParticle)==MaterialPointTypeLiquid) then 
            Beta_s = 1.0
          end if ! Liquid Material point
        end if
           
        !!=========== ELEMENT CONTAINS = LIQUID =====================!!
        if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
          Beta_s = 1.0
        end if

        !!=========== ELEMENT CONTAINS = SOLID =====================!!
        if (TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) then
          Eoed_Skel = ( 2 * G * (1 - Nu))/ (1 - 2 * Nu) ! kPa  ! Eoed
          Eoed = Eoed_Skel 
          rho_sat = Particles(IParticle)%EffConcentrationRatioSolid * ConstDensitySolid
          Beta_s = 0.0
        end if

        !--------- Calculate wave speed --------------!
        if(rho_sat>0.0) then
          WaveSpeed_c1 = sqrt((Eoed / rho_sat)/CalParams%ScalingMassFactor)   
        else
          WaveSpeed_c1 = 0.0
        end if
   
        if(ConstDensityLiquid>0.0) then
          WaveSpeed_c2 = Beta_s * sqrt((BulkW / ConstDensityLiquid)/CalParams%ScalingMassFactor)
        else
          WaveSpeed_c2 = 0.0
        end if
         
      end subroutine GetWaveSpeed2LayerForm

        
      end module ModDYNConvectivePhase
