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
	  
	  
	  module ModLagrangianPhase
      !**********************************************************************
      !
      ! Function:  To initialise the Dynamic MPM time step i.e. to map load and velocities from 
      !            particles to nodes, fill the mass matrix and solve the equations
      !            is necessary to be applied.
      !
      !     $Revision: 9908 $
      !     $Date: 2023-03-06 09:54:29 +0100 (ma, 06 mrt 2023) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters 
      use ModReadCalculationData
      use ModMPMData
      use ModParticle
      use ModMPMStresses
      use ModMPMDYNBTSigma
      use ModMPMDYN2PhaseSP
      use ModMPMDYN3PhaseSP
      use ModMeshInfo
      use ModRotBoundCond
      use ModDynViscousBoudary
      use ModTwoLayerFormulation
      use ModRigidBody
      use ModMPMInit
      use ModMPMDynContact
      
      
      implicit none

      contains ! Routines of this module

        subroutine LagrangianPhase()
        !**********************************************************************
        !
        ! Function:  Initialise the Lagrangian Phase for the dynamic MPM calculation:
        !            - Loads and velocities are transferred from particles to nodes
        !            - Fill the mass matrix by mapping of mass from particles to nodes
        !            - Integrate the momentum equation to find the accelerations,...velocities..., and displacements
        ! 
        !**********************************************************************

        implicit none

          
        ! Local variables
        integer(INTEGER_TYPE) :: IEntity, I
        integer(INTEGER_TYPE) :: iDof, J, ILoadSystem
        real(REAL_TYPE), dimension(Counters%N,Counters%nEntity) :: &
            Momentum, MomentumW, RateofMomentumW, WaterInertia, &
            MomentumG, RateofMomentumG, GasInertia

        real(REAL_TYPE), dimension(Counters%NodTot,NVECTOR) :: PrescribedVelocityNodes

        logical DoSystem                              !CC - flag used to do system velocity too...


          Momentum        = 0.0
          MomentumW       = 0.0
          RateofMomentumW = 0.0
          WaterInertia    = 0.0
          MomentumG       = 0.0
          RateofMomentumG = 0.0
          GasInertia      = 0.0
          
          ! Form the mass matrix and map momentum to nodes
          if (.not.IsMPMSkipConvection()) then
            call MapMomentumAndMassP2N(Momentum, LumpedMassDry) !  Map mass and momentum from particle to node (global coord system)
          endif

          ! Form the total lumped mass by adding the lumped masses of all entities
          if (CalParams%ApplyContactAlgorithm) then
              if (.not.IsMPMSkipConvection()) then
                  call FormSystemLumpedMass()
              endif
          end if ! contact algorithm

          ! Rotate Momentum from global to local coordinate system
          if (IS3DCYLINDRIC) then ! only 3D
            do IEntity = 1, Counters%nEntity
              call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, Momentum(:, IEntity), Momentum(:, IEntity))
            end do
          end if ! rotate BC

          ! Calculate nodal velocities
          if (IsMPMComputation()) then 
            ! If MPM then map velocity from partcles to nodes otherwise it is already at the nodes (FEM)
            if (.not.IsMPMSkipConvection()) then
              DoSystem = .true.     !if contact model is used, get TotalVelocitySys too
              call GetNodalVelocityFromNodalMomentum(Momentum, DoSystem) ! Local coordinate system
            endif
          endif

          TotalVelocitySoilPrevious = TotalVelocitySoil
 
  ! **  ******************************************************************************************************************************
  ! **  ******************************************************************************************************************************
  ! **  ********************************* SINGLE LAYER FORMULATION *******************************************************************
  ! **  ******************************************************************************************************************************
  ! **  ******************************************************************************************************************************
      if (NFORMULATION==1) then
  ! **  *********************************************************************************************************
  ! **  *** Dynamic Equilibrium of the liquid phase ----> acceleration of the liquid ****************************
  ! **  *********************************************************************************************************

          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then

              call FormConsolidationMatrices(TotalVelocitySoil, TotalVelocityWater) !To extrapolate loads and masses from particles to nodes

              if (IsMPMComputation()) then
                  if (.not.IsMPMSkipConvection()) then
                      ! If MPM then map velocity from material points to nodes otherwise it is already at the nodes (FEM)
                      call MapWaterMomentumFromMaterialPointsToNodes(MomentumW) ! Global coordinate system
                  endif

                  ! Rotate momentum from global to local coordinate system
                  if (IS3DCYLINDRIC) then ! only 3D
                      do IEntity = 1, Counters%nEntity
                          call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, MomentumW(:, IEntity), MomentumW(:, IEntity))
                      end do
                  end if ! only for RotBoundCond
                  if (.not.IsMPMSkipConvection()) then
                      call GetNodalWaterVelocityFromNodalWaterMomentum(MomentumW) ! in local coordinate system
                  endif
              end if ! only for MPMComputation

              call GetQVWArray(QVW)  ! (Drag force) does not need rotation because V and W are already rotated

     
                  ! v and w are in global c.s. so rotate QVW to local c.s.
                  if (IS3DCYLINDRIC) then ! only 3D
                      do IEntity = 1, Counters%nEntity
                          call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, QVW(:, IEntity), QVW(:, IEntity))
                      end do
                  end if ! only for RotBoundCond
          

              if ( CalParams%ApplyAbsorbingBoundary) then
                  call GetVisDampingForceWat(AccelerationWater) ! rotated inside the subroutine to local coordinate system
              end if ! only for AbsorbingBoundary

              if (CalParams%ApplySubmergedCalculation) then
                  if (CalParams%IStep <= CalParams%NumberSubmergedCalculation) then  ! for submerged calculation,
                      QVW = 0.0  ! the velocity and acceleration of wave are zero at gravity phase
                      QVWPorosity = 0.0
                  end if
              end if ! only for SubmergedCalculation

              RateofMomentumW = ExtLoadWater + GravityLoadWater - IntLoadWater ! momentum of liquid (local coordinate system)

              if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
                  if (CalParams%DampingFactor>0.0) then ! apply damping (liquid)
                      call ApplyDampingForceWater (RateofMomentumW) ! apply the local damping on liquid phase everywhere, local coordinate system
                  end if
              else ! use the damping from GOM file

                  call MapWaterDampingFromParticles(RateofMomentumW) ! damping applied on layers

              end if ! damping on liquid

              RateofMomentumW = (RateofMomentumW + QVW - VisDampForceWat)/CalParams%ScalingMassFactor ! momentum of liquid

              call CalculateWaterIncrementalNodalAcceleration(RateofMomentumW) ! local coordinate system

              ! Rotate liquid velocity and acceleration from local to global coordinate system
              if (IS3DCYLINDRIC) then ! only 3D
                  do IEntity = 1, Counters%nEntity
                      call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, AccelerationWater(:,IEntity), AccelerationWater(:,IEntity))
                  end do
                  if (.not.IsMPMSkipConvection()) then
                      do IEntity = 1, Counters%nEntity
                          call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocityWater(:,IEntity), &
                              TotalVelocityWater(:,IEntity))
                      end do
                  endif
              endif ! only for RotBoundCond

              TotalVelocityWaterPrevious = TotalVelocityWater ! In global coordinate system

              ! Update nodal TotalVelocityWater needed for contact formulation (in global coordinate system)
              do IEntity = 1, Counters%nEntity
                  do I = 1, Counters%N
                      TotalVelocityWater(I, IEntity) = TotalVelocityWater(I, IEntity) &
                          + AccelerationWater(I, IEntity) * CalParams%TimeIncrement
                  end do
              end do


              if(CalParams%ApplyContactAlgorithm) then
                  call GetNodePrescribedVelocity(PrescribedVelocityNodes)
                  call CorrectWaterForPrescribedVelocity(InterfaceNodes,ContactNodeNormals,PrescribedVelocityNodes)
              end if


              ! rotate liquid velocity and acceleration back from global to local coordinate system
              if (IS3DCYLINDRIC) then ! 3D function
                  do IEntity = 1, Counters%nEntity
                      call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, AccelerationWater(:,IEntity), AccelerationWater(:,IEntity))
                      call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, TotalVelocityWater(:,IEntity), TotalVelocityWater(:,IEntity))
                  end do
              end if ! only for RotBoundCond

              call GetWaterInertiaArray(WaterInertia) ! (Mass*Acceleration) In local coordinate system

              ! rotate the liquid acceleration from local to global coordinate system
              if (IS3DCYLINDRIC) then ! 3D function
                  do IEntity = 1, Counters%nEntity
                      call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, AccelerationWater(:,IEntity), AccelerationWater(:,IEntity))
                  end do
              end if ! only for RotBoundCond

          end if


! **  *********************************************************************************************************
! **  *** Dynamic Equilibrium of the gas phase ----> acceleration of the gas ********************************** 
! **  *********************************************************************************************************

          if (CalParams%NumberOfPhases==3) then

            call FormConsolidationMatricesGas()

            if (IsMPMComputation()) then
              ! If MPM then map velocity from partcles to nodes otherwise it is already at the nodes (FEM)
              call MapGasMomentumFromParticlesToNodes(MomentumG) 
              if (IS3DCYLINDRIC) then ! 3D function
                do IEntity = 1, Counters%nEntity
                  call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, MomentumG(:, IEntity), MomentumG(:, IEntity))
                end do
              end if ! only for RotBoundCond
              call GetNodalGasVelocityFromNodalGasMomentum(MomentumG)
            end if ! only for MPMComputation

            call GetQVWgArray()  ! (Drag force) does not need rotation because V and W are already rotated
          
            if ( CalParams%ApplyAbsorbingBoundary) then
              call GetVisDampingForceGas(AccelerationGas) ! rotated inside the subroutine 
            end if ! only for AbsorbingBoundary

            if (CalParams%ApplySubmergedCalculation) then  
              if (CalParams%IStep <= CalParams%NumberSubmergedCalculation) then
                ! for submerged calculation,
                ! the velocity and acceleration of wave are zero at gravity phase
                QVWg = 0.0
                ! the velocity and acceleration of wave are zero at gravity phase * porosity * degree of saturation
                QVWgPorosityDegreeSat = 0.0
              end if
            end if ! only for SubmergedCalculation
                
            RateofMomentumG = ExtLoadGas + GravityLoadGas - IntLoadGas ! momentum of gas

            if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
              if (CalParams%DampingFactor>0.0) then ! apply damping (gas)
                call ApplyDampingForceGas (RateofMomentumG) ! apply the local damping on gas phase everywhere
              end if
            else ! use the damping from GOM file

                call MapGasDampingFromParticles(RateofMomentumG) ! damping applied on layers 

            end if ! damping on gas

            RateofMomentumG = (RateofMomentumG + QVWg - VisDampForceGas)/CalParams%ScalingMassFactor ! momentum of gas

            call CalculateGasIncrementalNodalAcceleration(RateofMomentumG)

            call GetGasInertiaArray(GasInertia) ! (Mass*Acceleration) no rotation because the acceleration is already rotated

            if (IS3DCYLINDRIC) then ! 3D function
              ! rotation is needed (rotate the gas acceleration back to xyz because it will not be needed here again)
              do IEntity = 1, Counters%nEntity
                call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, AccelerationGas(:,IEntity), AccelerationGas(:,IEntity))
              end do
            end if ! only for RotBoundCond

          end if

! **  *********************************************************************************************************
! **  *** Dynamic Equilibrium of the mixture ----> acceleration of the solid skeleton**************************  
! **  *********************************************************************************************************
          
          ! Increase load by a load step is moved into forces calculation code
          ! Transfer external and gravity loads from particles to nodes and calculate the internal force vector (BT*Sig)
          if(CalParams%TimeStep==1) then ! first time step of load step
            call GetNodalExtAndIntForces() ! rotated to local coordinate system
          else
            call GetNodalExtForces() ! rotated to local coordinate system
          end if
          
          if ( CalParams%ApplyAbsorbingBoundary) then
            call GetVisDampingForceSld() ! rotated inside the subroutine
          end if ! only for AbsorbingBoundary
             
          
         !!!!!!!!!!!!!! START : Calculate the nodal rate of momentum !!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!! 2-phase calculation (solid+liquid) !!!!!!!!!!!!!!!!!!!!!!!!!!! 
         if (CalParams%NumberOfPhases==2) then
           
              
           if (CalParams%TimeStep==1) then ! first time step of load step
             ! get the real liquid info, taking into account the porosity and (eventually) degree of saturation
             call ConsolidationForcesPorosity(ExtLoadWaterPorosity, GravityLoadWaterPorosity, IntLoadWaterPorosity)
           else
             ! get the real liquid info, taking into account the porosity and (eventually) degree of saturation
             call ConsolidationExtForcesPorosity(ExtLoadWaterPorosity, GravityLoadWaterPorosity)
           end if
           call MapMPPorosityDegSatToNodes(LumpedNodalPorosityDegSat)
         
           do ILoadSystem = 1, Counters%NWaterLoadSystems    
               ExtLoadWaterPorosity = ExtLoadWaterPorosity + (ExtLoadWaterTotal(:,:,ILoadSystem) * &
                   CalParams%Multipliers%WaterACurrent(ILoadSystem) ) *  LumpedNodalPorosityDegSat
           end do
           ExtLoadWaterPorosity = ExtLoadWaterPorosity +  SpaceTimeExTLoadWater *  LumpedNodalPorosityDegSat

           ! liquid, local coordinate system
           RateofMomentumW = (ExtLoadWaterPorosity + GravityLoadWaterPorosity - IntLoadWaterPorosity)/CalParams%ScalingMassFactor

           if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
             if (CalParams%DampingFactor>0.0) then ! apply damping (liquid)
               call ApplyDampingForceWater (RateofMomentumW) ! apply the local damping on liquid phase everywhere, local coordinate system
             end if
           else ! use the damping from GOM file

               call MapWaterDampingFromParticles(RateofMomentumW) ! damping applied on layers

           end if ! damping on liquid

           if (CalParams%ApplyFixedSolidSkeleton) then ! Solid is fixed, compresibility of solid skeleton is neglected.
               RateofMomentum = 0.0
           else !solids are free to move, update the momentum equations
               if(CalParams%ApplyPartialSaturation) then
                   !calculate internal load water accounting for partial saturation
                   if (CalParams%TimeStep==1) call ConsolidationForcesBishop(BishopIntLoad) !if IsFollowUpPhase it is already computed in CalculateIntAndExtWorks

                   RateofMomentum = ((ExtLoad - ExtLoadWaterPorosity) + & ! (this term corresponds to the effective loading, partial saturation is included in ExtLoadWaterPorosity)
                       (GravityLoadMixture - GravityLoadWaterPorosity) - & ! (this term corresponds to the gravity of dry soil, partial saturation is included in GravityLoadWaterPorosity)
                       (IntLoad + BishopIntLoad - IntLoadWaterPorosity))/CalParams%ScalingMassFactor! solid, local coordinate system
               else
                   RateofMomentum = ((ExtLoad - ExtLoadWaterPorosity) + & ! (this term corresponds to the effective loading)
                       (GravityLoadMixture - GravityLoadWaterPorosity) - & ! (this term corresponds to the gravity of dry soil)
                       (IntLoad + (IntLoadWater - IntLoadWaterPorosity)))/CalParams%ScalingMassFactor! solid, local coordinate system
               end if

               if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
                   if (CalParams%DampingFactor>0.0) then ! apply damping (solid)
                       call ApplyDampingForce () ! use the damping from the con file (applied for all domain)
                   end if
               else ! use the damping from GOM file

                   call MapDampingFromParticles() ! solid

               end if ! damping on solid

               if (CalParams%ApplyBulkViscosityDamping) then
                   RateofMomentum = RateofMomentum - BulkViscLoad / CalParams%ScalingMassFactor
               end if

               RateofMomentum = RateofMomentum + RateofMomentumW &
                   - WaterInertia &
                   - VisDampForceSld - VisDampForceWat ! momentum of mixture, local coordinate system
           end if
           

         !!!!!!!!!!!!!!!!!!!! 3-phase calculation (solid+liquid+gas) !!!!!!!!!!!!!!!!!!!!!!!!!!!
         else if (CalParams%NumberOfPhases==3) then

           call ConsolidationForcesBishop(BishopIntLoad)
           call ConsolidationForcesPorosityDegreeSat(ExtLoadWaterPorosityDegreeSat, ExtLoadGasPorosityDegreeSat,  &
                             GravityLoadWaterPorosityDegreeSat,GravityLoadGasPorosityDegreeSat,  &
                             IntLoadWaterPorosityDegreeSat, IntLoadGasPorosityDegreeSat) ! get the real liquid and gas info, taking into account the porosity and the degree of saturation
           
           if(CalParams%NumberOfMaterials==1) then
               do ILoadSystem=1, Counters%NWaterLoadSystems
                   ExtLoadWaterPorosityDegreeSat = ExtLoadWaterPorosityDegreeSat +  (ExtLoadWaterTotal(:,:,ILoadSystem) * &
                       CalParams%Multipliers%WaterACurrent(ILoadSystem) ) *  &
                       (MatParams(1)%InitialPorosity * MatParams(1)%InitialDegreeOfSaturation)
               end do
               ExtLoadWaterPorosityDegreeSat = ExtLoadWaterPorosityDegreeSat +  &
                   SpaceTimeExTLoadWater * (MatParams(1)%InitialPorosity * MatParams(1)%InitialDegreeOfSaturation)
           else
           end if


           if(CalParams%NumberOfMaterials==1) then
               do ILoadSystem=1, Counters%NWaterLoadSystems
                   ExtLoadGasPorosityDegreeSat = ExtLoadGasPorosityDegreeSat + (ExtLoadWaterTotal(:,:,ILoadSystem) * &
                       CalParams%Multipliers%WaterACurrent(ILoadSystem)) *  &
                       (MatParams(1)%InitialPorosity *(1.0 - MatParams(1)%InitialDegreeOfSaturation))
               end do
               ExtLoadGasPorosityDegreeSat = ExtLoadGasPorosityDegreeSat + &
                   SpaceTimeExTLoadWater * (MatParams(1)%InitialPorosity *(1.0 - MatParams(1)%InitialDegreeOfSaturation))
           else
           end if
          
           
           RateofMomentumW = (ExtLoadWaterPorosityDegreeSat + GravityLoadWaterPorosityDegreeSat - IntLoadWaterPorosityDegreeSat)/ &
                              CalParams%ScalingMassFactor ! liquid, local coordinate system
           
           if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
             if (CalParams%DampingFactor>0.0) then ! apply damping (liquid)
               call ApplyDampingForceWater (RateofMomentumW) ! apply the local damping on water phase everywhere
             end if
           else 

               call MapWaterDampingFromParticles(RateofMomentumW) ! damping applied on layers 

           end if ! damping on liquid

!         
           RateofMomentumG = (ExtLoadGasPorosityDegreeSat + GravityLoadGasPorosityDegreeSat - IntLoadGasPorosityDegreeSat)/ &
                              CalParams%ScalingMassFactor  ! gas

           if (CalParams%ApplyCPSDamping) then  ! use the damping from the CPS file (applied for complete domain)
             if (CalParams%DampingFactor>0.0) then ! apply damping (gas)
               call ApplyDampingForceGas (RateofMomentumG)  ! apply the local damping on gas phase everywhere
             end if
           else 

               call MapGasDampingFromParticles(RateofMomentumG) ! damping applied on layers 

           end if ! damping on gas
           if (CalParams%ApplyFixedSolidSkeleton) then ! Solid is fixed, compresibility of solid skeleton is neglected.
               RateofMomentum = 0.0
           else !solids are free to move, update the momentum equations
               RateofMomentum = ((ExtLoad - ExtLoadWaterPorosityDegreeSat - ExtLoadGasPorosityDegreeSat) + & ! (this term corresponds to the effective loading)
                   (GravityLoadMixture - GravityLoadWaterPorosityDegreeSat - GravityLoadGasPorosityDegreeSat) - & ! (this term corresponds to the gravity of dry soil)
                   ((IntLoad + BishopIntLoad)- IntLoadWaterPorosityDegreeSat - IntLoadGasPorosityDegreeSat))/ &
                   CalParams%ScalingMassFactor ! solid, local coordinate system

               if (CalParams%ApplyCPSDamping) then  ! use the damping from the CPS file (applied for complete domain)
                   if (CalParams%DampingFactor>0.0) then  ! apply damping (solid)
                       call ApplyDampingForce ()  ! use the damping from the con file (applied for all domain)
                   end if
               else ! use the damping from GOM file

                   call MapDampingFromParticles() ! solid

               end if ! damping on solid

               RateofMomentum = RateofMomentum + RateofMomentumW + RateofMomentumG &
                   - WaterInertia - GasInertia &
                   - VisDampForceSld - VisDampForceWat - VisDampForceGas ! momentum of mixture
           end if
     
     
         !!!!!!!!!!!!!!!!!!!! 1-phase calculation (solid) !!!!!!!!!!!!!!!!!!!!!!!!!!               
         else if (CalParams%NumberOfPhases==1) then
          if (CalParams%ApplyFixedSolidSkeleton) then ! Solid is fixed, compresibility of solid skeleton is neglected.
               RateofMomentum = 0.0
            else !solids are free to move, update the momentum equations 
           RateofMomentum = (ExtLoad + GravityLoad - IntLoad)/CalParams%ScalingMassFactor
 
           if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (used for complete domain)
             if (CalParams%DampingFactor>0.0) then 
                 call ApplyDampingForce() ! normal case
               end if
                   else ! apply the damping of GID input

               call MapDampingFromParticles() ! normal case

           end if ! damping on solid

           if (CalParams%ApplyBulkViscosityDamping) then
             RateofMomentum = RateofMomentum - BulkViscLoad / CalParams%ScalingMassFactor
           end if

           if ( CalParams%ApplyAbsorbingBoundary) then
             RateofMomentum = RateofMomentum - VisDampForceSld
           end if ! only for AbsorbingBoundary
      end if !free/fixed solid skeleton                  
         end if ! difference between number of phases
       call GetInternalLoadRigidBody()
       call OverWriteRateofMomentumRigidBody()
         !!!!!!!!!!!!!! END : Calculate the nodal rate of momentum !!!!!!!!!!!!!!!!

          !!!!!!!!!!!!!! START : Calculate nodal acceleration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         call CalculateIncrementalNodalAcceleration() ! local coordinate system      
           
         if (IS3DCYLINDRIC) then ! 3D function
               ! Rotate acceleration and velocity of soil from local to global coordinate system
             do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, AccelerationSoil(:, IEntity), AccelerationSoil(:, IEntity))
               call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocitySoil(:,IEntity), TotalVelocitySoil(:,IEntity))
             end do
         end if ! only for RotBoundCond
         
         !!!!!!!!!!!!!! END : Calculate nodal acceleration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (CalParams%ApplyContactAlgorithm) then
           
           if (IS3DCYLINDRIC) then ! 3D function
             if (.not.IsMPMSkipConvection()) then
               ! Rotate total velocity from local to global coordinate system
               call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, TotalVelocitySys, TotalVelocitySys) ! back to xy
               ! Rotate acceleration from local to global coordinate system
               call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, AccelerationSys, AccelerationSys) ! back to xy
             endif
           end if ! only for RotBoundCond
           
           do IDOF = 1, Counters%N ! loop over all degrees of freedom

             TotalVelocitySys(IDOF) = TotalVelocitySys(IDOF) + AccelerationSys(IDOF) * CalParams%TimeIncrement ! Global coordinate system
           end do
         end if
         
 
         ! store velocity before correction and update nodal velocity
         ! Needed for predictor-correction scheme (e.g. contact algorithm, infiltration, seepage face...)
           do IDOF = 1, Counters%N ! loop over all degrees of freedom
             do J = 1, Counters%nEntity ! loop over all entities
               TotalVelocitySoilPrevious(IDOF,J) = TotalVelocitySoil(IDOF,J) ! global coordinate system
               TotalVelocitySoil(IDOF,J) = TotalVelocitySoil(IDOF,J) + AccelerationSoil(IDOF,J) *  CalParams%TimeIncrement ! Global coordinate system
             end do
           end do
           
         call ApplyMPPrescribedVelocity(TotalVelocitySys)
         call ApplyNodalPrescribedVelocity(TotalVelocitySys)           

         ! Calculate Non Advective and Conductive Fluxes (this is needed in the nodes)
         if((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation)) then
             call GetAdvectiveFluxes(AdvectiveFluxDarcyWater, AdvectiveFluxDarcyAir) !
         else if (CalParams%NumberOfPhases==3) then ! Unsat Calculation (gas)
             call GetNonAdvectiveFluxes(NonAdvectiveFluxAirInWater, NonAdvectiveFluxVapourInGas) !
             call GetAdvectiveFluxes(AdvectiveFluxDarcyWater, AdvectiveFluxDarcyAir) !
         end if
 
      
     !Correct Liquid and Solid acceleration if there is an infiltration/seepage boundary condition
         if((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation)) then
             if (CalParams%BoundaryConditions%UseInfiltration) then
                 call DetectInfiltrationNodes(IsInfiltrationNode)
             end if
             if (CalParams%PrescribedHead%HydraulicHead.or.CalParams%BoundaryConditions%UseSeepageFace) then
             call DetectSeepageFaceNodes(IsSeepageFaceNode)
         end if
             if (CalParams%BoundaryConditions%ApplyInfiltrationRate.or.CalParams%BoundaryConditions%ApplySeepageFace) then
                 call GetNodalUnitMassGradient(NodalUnitMassGradient,PBoundaryWater)
                 call ApplyNodalInfiltrationRate(TotalVelocityWater,TotalVelocitySoil,LumpedNodalPorosityDegSat,NodalUnitMassGradient)
                 call RecalculateAccelerationFromCorrectedVelocity(TotalVelocityWater,TotalVelocityWaterPrevious,AccelerationWater)
                 call RecalculateAccelerationFromCorrectedVelocity(TotalVelocitySoil,TotalVelocitySoilPrevious,AccelerationSoil)
             end if
         end if
         end if
         if (NFORMULATION==2) then
             ! **  ******************************************************************************************************************************
             ! **  ******************************************************************************************************************************
             ! **  *********************************   TWO LAYER FORMULATION   ******************************************************************
             ! **  ******************************************************************************************************************************
             ! **  ******************************************************************************************************************************



             if(CalParams%NumberOfPhases==1) then

                 !!!!!! -- LIQUID PHASE -- !!!!!!

                 !!!!!!!!!!!!!! START : Calculate the nodal rate of momentum !!!!!!!!!!!!!!!!

                 call FormConsolidationMatrices(TotalVelocitySoil, TotalVelocityWater)

                 call MapWaterMomentumFromMaterialPointsToNodes(MomentumW) ! Global coordinate system

                 call GetNodalWaterVelocityFromNodalWaterMomentum(MomentumW)

                 RateofMomentumW = ExtLoadWater + GravityLoadWater - IntLoadWater ! momentum of liquid (local coordinate system)

                 if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
                     if (CalParams%DampingFactor>0.0) then ! apply damping (liquid)
                         call ApplyDampingForceWater (RateofMomentumW) ! apply the local damping on liquid phase everywhere, local coordinate system
                     end if
                 else ! use the damping from GOM file

                     call MapWaterDampingFromParticles(RateofMomentumW) ! damping applied on layers

                 end if ! damping on liquid


                 call CalculateWaterIncrementalNodalAcceleration(RateofMomentumW) ! local coordinate system


                 TotalVelocityWaterPrevious = TotalVelocityWater ! In global coordinate system

                 ! Update nodal TotalVelocityWater needed for contact formulation (in global coordinate system)
                 do IEntity = 1, Counters%nEntity
                     do I = 1, Counters%N
                         TotalVelocityWater(I, IEntity) = TotalVelocityWater(I, IEntity) &
                             + AccelerationWater(I, IEntity) * CalParams%TimeIncrement
                     end do
                 end do


                 !!!!!! -- SOLID PHASE -- !!!!!!
                 if (CalParams%ApplyFixedSolidSkeleton) then ! Solid is fixed, compresibility of solid skeleton is neglected.
                     RateofMomentum = 0.0
                 else !update the momentum equations
                     ! Increase load by a load step is moved into forces calculation code
                     ! Transfer external and gravity loads from particles to nodes and calculate the internal force vector (BT*Sig)
                     call GetNodalExtAndIntForces() ! rotated to local coordinate system

                     RateofMomentum = ExtLoad  + & ! (this term corresponds to the effective loading)
                         GravityLoad  - & ! (this term corresponds to the gravity of dry soil)
                         IntLoad ! solid, local coordinate system

                     if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
                         if (CalParams%DampingFactor>0.0) then ! apply damping (solid)
                             call ApplyDampingForce () ! use the damping from the con file (applied for all domain)
                         end if
                     else ! use the damping from GOM file

                         call MapDampingFromParticles() ! solid

                     end if ! damping on solid

                 end if

                 !!!!!!!!!!!!!! END : Calculate the nodal rate of momentum !!!!!!!!!!!!!!!!

                 !!!!!!!!!!!!!! START : Calculate nodal acceleration !!!!!!!!!!!!!!!!!!!!!!


                 call CalculateIncrementalNodalAcceleration() ! local coordinate system


                 !!!!!!!!!!!!!! END : Calculate nodal acceleration !!!!!!!!!!!!!!!!!!!!!!!!

             end if ! if NumbOfPhases=1

             !!!!!!!!!!!!!!!!!!!!!!!!! Number Of Phases = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if(CalParams%NumberOfPhases==2) then

                 !!!!!! -- LIQUID PHASE -- !!!!!!
                 call FormConsolidationMatrices(TotalVelocitySoil, TotalVelocityWater)

                 call MapWaterMomentumFromMaterialPointsToNodes(MomentumW) ! Global coordinate system

                 call GetNodalWaterVelocityFromNodalWaterMomentum(MomentumW)

                 call GetQVWArray(QVW)  ! does not need rotation because V and W are already rotated

                 RateofMomentumW =  ExtLoadWater + GravityLoadWater - IntLoadWater

                 call TwoLayerData%DetermineInteractionForces()
                 RateofMomentumW =  RateofMomentumW - TwoLayerData%InteractionForceSolid ! momentum of liquid (local coordinate system)

                 if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
                     if (CalParams%DampingFactor>0.0) then ! apply damping (liquid)
                         call ApplyDampingForceWater (RateofMomentumW) ! apply the local damping on liquid phase everywhere, local coordinate system
                     end if
                 else ! use the damping from GOM file

                     call MapWaterDampingFromParticles(RateofMomentumW) ! damping applied on layers

                 end if ! damping on liquid

                 RateofMomentumW = (RateofMomentumW + QVW)/CalParams%ScalingMassFactor

                 call CalculateWaterIncrementalNodalAcceleration(RateofMomentumW) ! local coordinate system


                 !!!!!! -- SOLID PHASE or MIXTURE -- !!!!!!
                 if (CalParams%ApplyFixedSolidSkeleton) then ! Solid is fixed, compresibility of solid skeleton is neglected.
                     RateofMomentum = 0.0
                 else
                     ! Increase load by a load step is moved into forces calculation code
                     ! Transfer external and gravity loads from particles to nodes and calculate the internal force vector (BT*Sig)
                     call GetNodalExtAndIntForces() ! rotated to local coordinate system


                     RateofMomentum =  ExtLoad + GravityLoad  - IntLoad

                     RateofMomentum = RateofMomentum +  TwoLayerData%InteractionForceSolid ! solid, local coordinate system

                     if (CalParams%ApplyBulkViscosityDamping) then
                         RateofMomentum = RateofMomentum- BulkViscLoad !apply viscous damping
                     end if


                     if (CalParams%ApplyCPSDamping) then ! use the damping from the CPS file (applied for complete domain)
                         if (CalParams%DampingFactor>0.0) then ! apply damping (solid)
                             call ApplyDampingForce () ! use the damping from the con file (applied for all domain)
                         end if
                     else ! use the damping from GOM file

                         call MapDampingFromParticles() ! solid

                     end if ! damping on solid

                     RateofMomentum = (RateofMomentum - QVW)/CalParams%ScalingMassFactor

                 end if
                 !!!!!!!!!!!!!! END : Calculate the nodal rate of momentum !!!!!!!!!!!!!!!!

                 !!!!!!!!!!!!!! START : Calculate nodal acceleration !!!!!!!!!!!!!!!!!!!!!!


                 call CalculateIncrementalNodalAcceleration() ! local coordinate system


             end if ! if NumbOfPhases=2

         end if ! if NumbOfLayers=2

          
 
        end subroutine LagrangianPhase 

   
        subroutine GetNonAdvectiveFluxes(NonAdvectiveFluxAirInWater, NonAdvectiveFluxVapourInGas)
        !**********************************************************************
        !
        !  Function:  Calculation of the non advective flux of the air in the liquid because 
        !             of the diffusion
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) ::  &
               NonAdvectiveFluxAirInWater,  &
               NonAdvectiveFluxVapourInGas
          
          ! Local variables
          real(REAL_TYPE), dimension(NDIM, ELEMENTNODES) :: B
          
          integer(INTEGER_TYPE) :: I, IntGlo, IEl, Int, NN, IAEl, NElemPart, iEntityID, iNode
          
          real(REAL_TYPE) :: WtN, Det, S1, S2
          real(REAL_TYPE) :: g,N,Pg,T,WD,GD,Sr,Sg,R,D1,D2,Q,num,Tort
          real(REAL_TYPE) :: DiffusionAirInWater,DiffusionVapourInGas
          
          NonAdvectiveFluxVapourInGas = 0.0
          NonAdvectiveFluxAirInWater = 0.0

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
              
                g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
                N = Particles(IntGlo)%Porosity              !Porosity of the particle
                Pg = Particles(IntGlo)%GasPressure          !Gas Pressure (MPa)
                T = Particles(IntGlo)%Temperature           !Temperature (ºC)
                WD = Particles(IntGlo)%WaterWeight/g        !Water density of the particle
                GD = Particles(IntGlo)%GasWeight/g          !Gas density of the particle
                Sr = Particles(IntGlo)%DegreeSaturation     !Degree of Saturation (liquid) of the particle
                Sg = 1.0d0 - Sr                             !Degree of Satutation (gas) of the particle
                    
                T = T+273.15d0  !Temperature (ºK)
                R = 0.008314d0  !Default=8.3144621, Gas constant value (J/mol*K)
                D1 = 0.00011d0  !Default=0.00011, Coefficient (m2/s)
                Q = 24.53d0     !Default=24530, Coefficient (J/mol)
                D2 = 0.0000000059d0   !999999 Default=0.0000059, Coefficient (m2*Pa*K**(-num)/s)
                num = 2.3d0     !Default=2.3
                Tort = 1.0d0        !999999 Tortuosity reduction coefficient

                DiffusionAirInWater = D1*exp(-Q/(R*T))
                DiffusionVapourInGas = D2*(T**num)/Pg
                    
                if (Pg<0.000001.or.Pg>-0.000001) then
                  DiffusionVapourInGas = 0.0d0
                end if
            
                S1 = - N * WD * Sr * Tort * DiffusionAirInWater * Particles(IntGlo)%AirInWaterMassFraction   !///////// Fick's Law
                S2 = - N * GD * Sg * Tort * DiffusionVapourInGas * Particles(IntGlo)%VapourInGasMassFraction !///////// Fick's Law
        
                !get particle entity
                if (.not.CalParams%ApplyContactAlgorithm) then      
                  iEntityID = 1
                else                                           
                  iEntityID = EntityIDArray(IntGlo)
                end if        
                
                do INode=1,ELEMENTNODES  !loop through element nodes

                  nn=ElementConnectivities(iNode,iel) ! get global node number
                
                  !///////// Fick's Law: calculation of the nodal Non Advective flux of air in the water ( iAW = - D * Grad(wAW) ) ///////////
                  do I = 1, NVECTOR
                    NonAdvectiveFluxAirInWater(ReducedDof(nn)+I,iEntityID) = NonAdvectiveFluxAirInWater(ReducedDof(nn)+I,iEntityID) + B(I,INode)*S1 
                  end do

                  !///////// Fick's Law: calculation of the nodal Non Advective flux of vapour in the gas ( iVG = - D * Grad(wVG) ) ///////////
                  do I = 1, NVECTOR
                    NonAdvectiveFluxVapourInGas(ReducedDof(nn)+I,iEntityID) = NonAdvectiveFluxVapourInGas(ReducedDof(nn)+I,iEntityID) + B(I,INode)*S2 
                  end do  
                                             
                end do ! node loop             
              
              end do ! Loop over 1 .. NElemPart
              !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
          end do ! Loop over elements
          
        end subroutine GetNonAdvectiveFluxes

       
        subroutine GetAdvectiveFluxes(AdvectiveFluxDarcyWater,AdvectiveFluxDarcyAir)
        !**********************************************************************
        !
        ! Function:  To map momentum from particles to the grid points (nodes)
        !
        ! O  Momentum : Nodal momentum vector, the output of this subroutine
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: AdvectiveFluxDarcyWater,AdvectiveFluxDarcyAir
     
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IPart, INode, Ni, ParticleIndex, NodeID, iEntity
          real(REAL_TYPE) :: g, N, WD, GD, Sr, mass, shape
          
          AdvectiveFluxDarcyWater = 0.0
          AdvectiveFluxDarcyAir = 0.0
          
          do IEl = 1, Counters%NEl ! loop over all elements
            do IPart = 1, NPartEle(IEl) ! loop over all material points in element
              
              ParticleIndex = GetParticleIndex(IPart,IEl) ! get the particle ID
              g = CalParams%GravityData%GAccel  !Gravity (m/s2)
              N = Particles(ParticleIndex)%Porosity              !Porosity of the particle
              WD = Particles(ParticleIndex)%WaterWeight/g        !Water density of the particle
              GD = Particles(ParticleIndex)%GasWeight/g        !Water density of the particle
              Sr = Particles(ParticleIndex)%DegreeSaturation     !Degree of Saturation (liquid) of the particle
              mass = MassArray(ParticleIndex) !Mass of particle
             
                
              if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex)! entity to which particle belongs
              else
                iEntity = 1
              end if ! only for ContactAlgorithm
                  
              do INode = 1, ELEMENTNODES ! loop over nodes
                
                NodeID = ElementConnectivities(INode, IEl) ! global node ID
                shape = ShapeValuesArray(ParticleIndex,INode)
                !!!!!!!!!!!!!!!!!!!Water!!!!!!!!!!!!!!!!!!!!!
                do I = 1, NVECTOR
                  Ni = ReducedDof(NodeID) + I  ! global storage coordinate
                  AdvectiveFluxDarcyWater(Ni,iEntity) = AdvectiveFluxDarcyWater(Ni,iEntity) + (Mass/ LumpedMassDry(Ni,iEntity)) * Shape * N*WD*Sr
                end do
                !!!!!!!!!!!!!!!!!!!Air!!!!!!!!!!!!!!!!!!!
                if (CalParams%NumberOfPhases==3) then
                  do I = 1, NVECTOR
                    Ni = ReducedDof(NodeID) + I  ! global storage coordinate
                    AdvectiveFluxDarcyAir(Ni,iEntity) = AdvectiveFluxDarcyAir(Ni,iEntity) + (Mass/ LumpedMassDry(Ni,iEntity)) * Shape* N*GD*(1-Sr)
                  end do
                end if

               
              end do ! loop over nodes
            end do ! loop over material points
          end do ! loop over elements
          
        end subroutine GetAdvectiveFluxes
        
        

        real(REAL_TYPE) function GetParticlePrescribedVelocityI(IDim, ParticleIndex)
        !**********************************************************************
        !
        ! Function:  Modify particles velocity considering the current multiplier
        ! in the context of moving mesh feature
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IDim, ParticleIndex
         
            GetParticlePrescribedVelocityI = Particles(ParticleIndex)%PrescrVelo(IDim) * CalParams%Multipliers%VelocitySolidCurrent
             

        end function GetParticlePrescribedVelocityI
        
        real(REAL_TYPE) function GetParticlePrescribedAccelerationI(IDim, ParticleIndex)
        !**********************************************************************
        !
        ! Function:  in case of prescribed velocity, if it changes in time, the acceleration is computed.
        !
         !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IDim, ParticleIndex
          

            GetParticlePrescribedAccelerationI = Particles(ParticleIndex)%PrescrVelo(IDim) * CalParams%Multipliers%AccelerationSolid  
              

        end function GetParticlePrescribedAccelerationI
        
        real(REAL_TYPE) function GetNodalPrescribedVelocityI(IDim, INode)
        !**********************************************************************
        !
        !This function gives the value of prescribed velocity along coordinate IDim of node INode
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IDim, INode
          
          GetNodalPrescribedVelocityI = &
              CalParams%PrescribedVelo%NodalPrescribedVelocityValue(INode,IDim)* CalParams%Multipliers%VelocitySolidCurrent
          
        end function GetNodalPrescribedVelocityI
        
        
        real(REAL_TYPE) function GetNodalPrescribedAccelerationI(IDim, INode)
        !**********************************************************************
        !
        !This function gives the value of prescribed acceleration along coordinate IDim of node INode
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IDim, INode
          
          GetNodalPrescribedAccelerationI = &
              CalParams%PrescribedVelo%NodalPrescribedVelocityValue(INode,IDim)* CalParams%Multipliers%AccelerationSolid
          
        end function GetNodalPrescribedAccelerationI

        subroutine GetNodePrescribedVelocity(PrescribedVelocityNodes)
        !**********************************************************************
        !
        !    Function:  get nodal prescribed velocity, if Contact Algorithm
        !
        !*********************************************************************     

        implicit none
        
          real(REAL_TYPE), dimension(Counters%NodTot,NVECTOR), intent(out) :: PrescribedVelocityNodes
          ! Local variables
          real(REAL_TYPE), dimension(Counters%N) :: DofPrescribedVelocity
          real(REAL_TYPE), dimension(Counters%N) :: ISPrescribedDof
          integer(INTEGER_TYPE) :: IDoF, IDim, NodeID
          logical :: IsPrescribedVelocity
          
                     
          IsPrescribedVelocity = CalParams%PrescribedVelo%ApplyPrescribedVelo
          if (IsPrescribedVelocity) then
              if (CalParams%PrescribedVelo%NNodePrescribedVelo>0) then ! prescribed velocity at the nodes

                  do NodeID = 1, Counters%NodTot !loop nodes
                      do IDim=1,NVECTOR
                          PrescribedVelocityNodes(NodeID,IDim) = GetNodalPrescribedVelocityI(IDim,NodeID)
                      end do
                  end do

              else !prescribed velocity at MP

                  call MapMaterialPointPrescribedVelocityToNodes(DofPrescribedVelocity,ISPrescribedDof)
                  do NodeID = 1, Counters%NodTot !loop nodes

                      IDoF = ReducedDoF(NodeID)
                      do IDim=1,NVECTOR
                          PrescribedVelocityNodes(NodeID,IDim) = DofPrescribedVelocity(IDof+IDim)
                      end do

                  end do
              end if

          end if
        end subroutine GetNodePrescribedVelocity
      
          
        subroutine GetNodalExtAndIntForces()
        !**********************************************************************
        !
        !    Function:  To extrapolate loads from material points to nodes using 
        !               the shape function values evaluated at the material points local position.
        !
        !**********************************************************************

        implicit none


          integer(INTEGER_TYPE) :: IEntity, I, IDoF, J
          integer(INTEGER_TYPE) :: IDim, ILoadSystem
          
          logical :: DoConsiderReactionForces, IsPrescribedVelocity
          IDim = NVECTOR
          
           ! Internal and External forces are calculated inside "BtSig"
           ! Calculate internal nodal forces
           call MPMDYNBTSig(ExtLoad, IntLoad, GravityLoad, FReaction, FReactionWater, BulkViscLoad)
           
           if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
            do ILoadSystem = 1, Counters%NSoilLoadSystems    
              ExtLoad = ExtLoad + ExtLoadTotal(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
            end do

           end if
           
            if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
           ExtLoad = ExtLoad + SpaceTimeExTLoadWater
           end if
          
           IsPrescribedVelocity = CalParams%PrescribedVelo%ApplyPrescribedVelo
           
           DoConsiderReactionForces =  IsPrescribedVelocity  
           
           if (DoConsiderReactionForces) then ! Compute in global coordinate system for output
             do I = 1, Counters%NodTot
               IDoF = ReducedDof(I)
               do J = 1, IDim
                 FReaction(IDoF + J, 1:Counters%NEntity) =  &
                    FReaction(IDoF + J, 1:Counters%NEntity) - GravityLoad(IDoF + J, 1:Counters%NEntity)
                  if (CalParams%NumberOfPhases==2) then
                 FReactionWater(IDoF + J, 1:Counters%NEntity) = &
                    FReactionWater(IDoF + J, 1:Counters%NEntity) - GravityLoadWater(IDoF + J, 1:Counters%NEntity)
                  else if (CalParams%NumberOfPhases==3) then
                 FReactionGas(IDoF + J, 1:Counters%NEntity) = &
                    FReactionGas(IDoF + J, 1:Counters%NEntity) - GravityLoadGas(IDoF + J, 1:Counters%NEntity)
                  end if
               end do
             end do 
           end if ! only for ConsiderReactionForeces

           if (IS3DCYLINDRIC) then ! Rotate vectors from global to local coordinate system ! 3D function
             do IEntity = 1, Counters%nEntity 
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoad(:, IEntity),     ExtLoad(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoad(:, IEntity), GravityLoad(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IntLoad(:, IEntity),     IntLoad(:, IEntity))
             end do
           end if ! only for RotBoundCond
          
        end subroutine GetNodalExtAndIntForces
        
        subroutine GetNodalExtForces()
        !**********************************************************************
        !
        !    Function:  To extrapolate loads from material points to nodes using 
        !               the shape function values evaluated at the material points local position.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE) :: IEntity
          integer(INTEGER_TYPE) :: IDim, ILoadSystem

          IDim = NVECTOR
           
          call MPMDYNLoad(ExtLoad, GravityLoad)

           if ((Counters%NLoadedElementSidesSolidNodes+Counters%NLoadedElementSidesSolidNodesB) > 0) then
             do ILoadSystem =1, Counters%NSoilLoadSystems
                 ExtLoad = ExtLoad + ExtLoadTotal(:,:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
             end do
           end if
           
            if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
           ExtLoad = ExtLoad + SpaceTimeExTLoadWater
          end if
           
           if (IS3DCYLINDRIC) then ! 3D function
             do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, ExtLoad(:, IEntity),     ExtLoad(:, IEntity))
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, GravityLoad(:, IEntity), GravityLoad(:, IEntity))
             end do
           end if

        end subroutine GetNodalExtForces

        subroutine GetNodalIntForces()
        !**********************************************************************
        !
        !    Function:  get nodal internal forces
        !
        !*********************************************************************     
        implicit none

          integer(INTEGER_TYPE) :: IEntity, I, IDoF, J
          logical :: DoConsiderReactionForces, IsPrescribedVelocity

           call MPMDYNBTSigOnly(IntLoad, FReaction, BulkViscLoad)

           IsPrescribedVelocity = CalParams%PrescribedVelo%ApplyPrescribedVelo
           
           DoConsiderReactionForces =  IsPrescribedVelocity 

           if (DoConsiderReactionForces) then
             do I = 1, Counters%NodTot
               IDoF = ReducedDof(I)
               do J = 1, NDOFL
                 FReaction(IDoF + J, 1:Counters%NEntity) =  &
                    FReaction(IDoF + J, 1:Counters%NEntity) - GravityLoad(IDoF + J, 1:Counters%NEntity)
                  if (CalParams%NumberOfPhases==2) then
                  end if
               end do
             end do
           end if

           if (IS3DCYLINDRIC) then ! 3D function
             do IEntity = 1, Counters%nEntity
               call RotVec(IRotation, NRotNodes, RotMat, ReducedDof, IntLoad(:, IEntity), IntLoad(:, IEntity))
             end do
           end if

        end subroutine GetNodalIntForces

        subroutine ApplyDampingForce()
        !**********************************************************************
        !
        !   Function:  To calculate the damping force 
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: I, J
          real(REAL_TYPE) :: DampedOutOffBalance ! change the vector below to this scalar
                             
          DampedOutOffBalance = 0.0
           
          do I = 1, Counters%N ! loop over all degrees-of-freedom
            do J = 1, Counters%nEntity ! loop over all entities
         
              if (TotalVelocitySoil(I,J)/=0.0) then ! needed, beacuse sign() will return (+1) if NodalVelocities(I) = 0, and not zero as expected
                DampedOutOffBalance = - CalParams%DampingFactor * abs(RateofMomentum(I,J)) * sign(1.d0, TotalVelocitySoil(I,J))
                RateofMomentum(I,J) = RateofMomentum(I,J) + DampedOutOffBalance ! apply damping force
              end if ! only for TotalVelocitySoil not zero 

            end do ! loop over all entities
          end do ! loop over all degrees-of-freedom
        
        end subroutine ApplyDampingForce
         

        subroutine ApplyDampingForceWater(RateofMomentumW)
        !**********************************************************************
        !
        !   Function:  To calculate the damping force for liquid phase
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: RateofMomentumW 
        
          ! Local variables
          integer(INTEGER_TYPE) :: I, J
          real(REAL_TYPE) :: DampedOutOffBalance                      
                             
          DampedOutOffBalance = 0.0
           
         do I = 1, Counters%N ! loop over all degrees-of-freedom
            do J = 1, Counters%nEntity ! loop over all entities
              
              if (TotalVelocityWater(I,J)/=0.0) then ! needed, beacuse sign() will return (+1) if NodalVelocities(I) = 0, and not zero as expected
                DampedOutOffBalance = - CalParams%DampingFactor * abs(RateofMomentumW(I,J)) * sign(1.d0, TotalVelocityWater(I,J))
                    RateofMomentumW(I,J) = RateofMomentumW(I,J) + DampedOutOffBalance ! apply damping force
               end if ! only for TotalVelocityWater not zero  

            end do ! loop over all entities
          end do ! loop over all degrees-of-freedom
        
        end subroutine ApplyDampingForceWater
        
         
        subroutine ApplyDampingForceGas(RateofMomentumG)
        !**********************************************************************
        !
        !   Function:  To calculate the damping force for gas phase
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: RateofMomentumG 
        
          ! Local variables
          integer(INTEGER_TYPE) :: I,J
          real(REAL_TYPE) :: DampedOutOffBalance                      
                             
          DampedOutOffBalance = 0.0
           
          do I = 1, Counters%N ! loop over all degrees-of-freedom
            do J = 1, Counters%nEntity ! loop over all entities
            
              if (TotalVelocityGas(I,J)/=0.0) then ! needed, beacuse sign() will return (+1) if NodalVelocities(I) = 0, and not zero as expected
                DampedOutOffBalance = - CalParams%DampingFactor * abs(RateofMomentumG(I,J)) * sign(1.d0, TotalVelocityGas(I,J))
                RateofMomentumG(I,J) = RateofMomentumG(I,J) + DampedOutOffBalance ! apply damping force
              end if ! only for TotalVelocityGas not zero 
            
            end do ! loop over all entities
          end do ! loop over all degrees-of-freedom
        
        end subroutine ApplyDampingForceGas
         

        subroutine MapMomentumAndMassP2N(Momentum, LumpedMass)
        !**********************************************************************
        !
        ! Function:  to form the lumped mass vector and map momentum from material points to nodes
        !
        ! O  Momentum : nodal momentum vector, the output of this subroutine
        ! O  LumpedMass : lumped soil mass
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity), intent(inout) ::   &
              Momentum,  &
              LumpedMass
          
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(ELEMENTNODES) :: PartShape
          real(REAL_TYPE), dimension(ELEMENTNODES, Counters%nEntity) :: LumMassX
          real(REAL_TYPE) :: ParticleMass
          
          integer(INTEGER_TYPE) :: I, IEl, IPart, IAEl, ParticleIndex, iEntity, NoEn

          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: Ni ! number of columns is equal to the number of nodes of the element
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES, Counters%nEntity) :: Mom
         
          LumpedMass = 0.0
          Momentum = 0.0
          NoEn = Counters%nEntity
          
          do IAEl = 1, Counters%NAEl ! loop over all elements

              IEl = ActiveElement(IAEl)

              do I = 1, NVECTOR
                  Ni(I, 1:ELEMENTNODES) = ReducedDof( ElementConnectivities(1:ELEMENTNODES, IEl) ) + I
                  Mom(I, 1:ELEMENTNODES, 1:Counters%nEntity) = Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn)
              end do
              LumMassX = LumpedMass( Ni(1, 1:ELEMENTNODES), 1:NoEn)

              do IPart = 1, NPartEle(IEl) ! loop over all material points in element
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID

                  if((NFORMULATION==1).or.(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid)) then ! Only if NumberOfLayers = 1 or SOLID Material Point

                      ParticleVelocity = VelocityArray(ParticleIndex,:) !get particle velocity vector
                      ParticleMass = MassArray(ParticleIndex) !get particle's mass
                      PartShape = ShapeValuesArray(ParticleIndex,:) !get particle shape fucntions

                      if (CalParams%ApplyContactAlgorithm) then
                          iEntity = EntityIDArray(ParticleIndex)
                      else
                          iEntity = 1
                      end if

                      ! nodal x-mass
                      LumMassX(:, iEntity) = LumMassX(:, iEntity) + ParticleMass * PartShape(:)

                      ! nodal momentum in directions from 1 to NDIM
                      do I=1,NVECTOR
                          Mom(I, :, iEntity) = Mom(I, :, iEntity) + ParticleMass * PartShape(:) * ParticleVelocity(I)
                      end do

                  end if ! Only if NumberOfLayers = 1 or SOLID Material Point

              end do ! loop over material points in element

              do I = 1, NVECTOR
                  LumpedMass( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = LumMassX( 1:ELEMENTNODES, 1:NoEn )
                  Momentum( Ni(I, 1:ELEMENTNODES), 1:NoEn ) = Mom( I, 1:ELEMENTNODES, 1:NoEn )
              end do

          end do ! loop over Eelements
          
        end subroutine MapMomentumAndMassP2N
         

        subroutine FormSystemLumpedMass()
        !**********************************************************************
        !
        ! Function:  To Form the system Lumped mass vector
        !
        !**********************************************************************

        implicit none
           
          ! Local variables
          integer(INTEGER_TYPE) :: IEntity, IDof
          
          LumpedMassSys = 0.0
         
          do IDof = 1, Counters%N ! loop over degrees-of-freedom
            do IEntity = 1, Counters%NEntity ! loop over entities       
              LumpedMassSys(IDof) = LumpedMassSys(IDof) + LumpedMassDry(IDof,IEntity) 
            end do ! loop over entities
          end do ! loop over degrees-of-freedom
          
        end subroutine FormSystemLumpedMass

         
        subroutine MapMomentumFromParticlesToNodes(Momentum)
        !**********************************************************************
        !
        ! Function:  To map momentum from particles to the grid points (nodes)
        !
        ! O  Momentum : Nodal momentum vector, the output of this subroutine
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: Momentum
     
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IPart, INode, ParticleIndex, NodeID, iEntity
          real(REAL_TYPE), dimension (NVECTOR) :: ParticleVelocity
          integer(INTEGER_TYPE), dimension(NVECTOR) :: Ni 
          
          Momentum = 0.0
          
          do IEl = 1, Counters%NEl ! loop over all elements
            do IPart = 1, NPartEle(IEl) ! loop over all material points in element
              
              ParticleIndex = GetParticleIndex(IPart,IEl) ! get the particle ID
              ParticleVelocity = VelocityArray(ParticleIndex,:) ! get particle velocity vector
                
              if (CalParams%ApplyContactAlgorithm) then
                iEntity = EntityIDArray(ParticleIndex)
              else
                iEntity = 1
              end if ! only for ContactAlgorithm
                  
              do INode = 1, ELEMENTNODES ! loop over nodes
                
                NodeID = ElementConnectivities(INode, IEl) ! global node ID
                do I = 1, NVECTOR
                  Ni = ReducedDof(NodeID) + I ! global storage coordinate of x-val
                  Momentum(Ni,iEntity) = Momentum(Ni,iEntity) + MassArray(ParticleIndex) * ShapeValuesArray(ParticleIndex,INode) * ParticleVelocity(I)
                end do
               
              end do ! loop over nodes
            end do ! loop over material points
          end do ! loop over elements
          
        end subroutine MapMomentumFromParticlesToNodes


        subroutine GetNodalVelocityFromNodalMomentum(Momentum, DoSystem)
        !**********************************************************************
        !
        ! Function:  To calculate the nodal velocities from nodal mass and momentum
        !
        !    Momentum : Nodal momentum vector
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(in) :: Momentum
          
          logical DoSystem
                   
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF,J
          real(REAL_TYPE), dimension(Counters%N) :: MomentumSystem
          
          if (CalParams%ApplyContactAlgorithm .and. DoSystem) then
            MomentumSystem = 0.0
          end if ! only for ContactAlgorithm

          do IDOF = 1, Counters%N ! loop over all degrees of freedom
            do J = 1, Counters%nEntity ! loop over all entities
                
                if (LumpedMassDry(IDOF,J)/=0) then
                  TotalVelocitySoil(IDOF,J) = ( Momentum(IDOF,J) / LumpedMassDry(IDOF,J) ) * PBoundary(IDOF)
                else
                  TotalVelocitySoil(IDOF,J) = 0.0
                end if
                
                if (CalParams%ApplyContactAlgorithm .and. DoSystem) then ! for contact model, add entity momentum to get system momentum
                  MomentumSystem(idof) = MomentumSystem(idof) + Momentum(idof,J)
                end if ! only for ContactAlgorithm
                              
            end do ! loop over all entities
            
            if (CalParams%ApplyContactAlgorithm .and. DoSystem) then ! for contact model, calculate system nodal velocities
              if (LumpedMassSys(IDOF)/=0) then
                TotalVelocitySys(IDOF) = (MomentumSystem(IDOF)/ LumpedMassSys(IDOF)) * PBoundary(IDOF)
              else
                TotalVelocitySys(IDOF) = 0.0
              end if
            end if ! only for ContactAlgorithm
            
          end do ! loop over all degrees of freedom
          
        end subroutine GetNodalVelocityFromNodalMomentum
         
        
        subroutine CalculateIncrementalNodalAcceleration()
        !**********************************************************************
        !
        !  Function:  To calculate the incremental nodal accelerations from nodal mass and rateofmomentum
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IDOF, J, INode, C
          real(REAL_TYPE), dimension(Counters%N) :: MomentumSystem
          
          if (CalParams%ApplyContactAlgorithm) then
            MomentumSystem = 0.0
          end if ! only for ContactAlgorithm
                   

          do IDOF = 1, Counters%N ! loop over all degrees of freedom
            do J = 1, Counters%nEntity ! loop over all entities
          
                if(LumpedMassDry(IDOF,J)/=0) then
                  AccelerationSoil(IDOF,J)=( RateofMomentum(IDOF,J) / LumpedMassDry(IDOF,J) ) * PBoundary(IDOF)
                else
                  AccelerationSoil(IDOF,J) = 0.0
                end if
                
                if (CalParams%ApplyContactAlgorithm) then 
                  MomentumSystem(idof) = MomentumSystem(idof) + RateofMomentum(idof,J) ! sum nodal rates of momentum
                end if ! only for ContactAlgorithm               
                
            end do ! loop over all entities
            
            if (CalParams%ApplyContactAlgorithm) then 
              if(LumpedMassSys(IDOF)/=0) then
                AccelerationSys(IDOF)=(MomentumSystem(IDOF) / LumpedMassSys(IDOF) ) * PBoundary(IDOF)
                    if (CalParams%RigidBody%IsRigidBody) then !If rigid body is active the system acceleration must be constrained
                        if (MOD(idof,NDIM)>0) then
                            INode=(idof/NDIM)+1! Node number                
					    else
                           INode=(idof/NDIM)
                        end if
                        if (RigdBodyInterface(INode)) then ! Node belongs to interface
                            C=IDOF-NDIM*(INode-1)
                            if (CalParams%RigidBody%Constrains(C)==1) then
                                AccelerationSys(IDOF)=0
                            else
                                AccelerationSys(IDOF)=	AccelerationSoil(IDOF,HARD_ENTITY)
                            end if
                        end if
                    end if
                else
                AccelerationSys(IDOF) = 0.0
              end if
            end if ! only for ContactAlgorithm
            
          end do ! loop over all degrees of freedom
         
        end subroutine CalculateIncrementalNodalAcceleration

        subroutine InitialiseMaterialPointsForK0Stresses()                                         
        !**********************************************************************
        !
        !  Function:  Initialises material point stresses by means of K0.
        !             Also updates material point density in case of liquid material
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: ParticleIndex, MaterialIndex
          real(REAL_TYPE) :: SigY, SigYeff, SigXeff, SigWP, K0Value, UpdatedFluidDensity
          real(REAL_TYPE) :: TopLayerDepth, BottomLayerDepth
          real(REAL_TYPE), dimension(NTENSOR) :: SigmaEff0
          real(REAL_TYPE), dimension(NVECTOR) :: GaussPointGlobalCoord
          integer(INTEGER_TYPE) :: IGaussPoint
          integer(INTEGER_TYPE) :: IAEl, IEl, IElTyp, IPart
          integer(INTEGER_TYPE) :: NGaussPoints
          integer(INTEGER_TYPE) :: ISoilLayer
          real(REAL_TYPE) :: WeiGP, Depth, xMax, yMin, zMin
          real(REAL_TYPE) :: PosGP(NVECTOR)
          real(REAL_TYPE) :: PointCoord(2), Sr
          real(REAL_TYPE), dimension(ELEMENTNODES) :: GaussPointShapeValues
          real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: DShapeValues
          integer(INTEGER_TYPE), dimension(CalParams%NumberSoilLayers, 1) :: MaterialSetLayer
          logical :: FindMaterialID, IsUndrTotalStress
            
          if (IsFollowUpPhase().or.(.not.CalParams%ApplyK0Procedure)) RETURN
          
          MaterialSetLayer = 0.0
          IGaussPoint  =   1
          IElTyp       =   ELEMENTNODES
          NGaussPoints =   ELEMENTGAUSSPOINTS
          TopLayerDepth = 0.0
          BottomLayerDepth = 0.0
          Sr=1.0d0
          if (CalParams%NumberSoilLayers>1) then !using multilayer K0 for CPT
              call GenerateBoreholeInfo()
          endif
          
          if (Counters%SoilSurfaceNumberofSides == 0) then ! if a soil surface has been specified in GID
              allocate(SoilSurfaceNodeCoordMatrix(NDIM,NDIM))
              SoilSurfaceNodeCoordMatrix = 0.0
          elseif (Counters%SoilSurfaceNumberofSides > 0) then
              call SurfaceNodesCoordinates(SoilSurfaceNodeCoordMatrix, SoilSurfaceNodesConnectivities, Counters%SoilSurfaceNumberofSides)
          end if
          
        xMax= maxval(NodalCoordinates(:,1))   
        yMin= minval(NodalCoordinates(:,2))
        if (NDIM == 3) then
            zMin= minval(NodalCoordinates(:,2))
        end if
        
        if (Counters%PhreaticSurfaceNumberofSides == 0)then
            allocate(PhreaticSurfaceNodeCoordMatrix(NDIM,NDIM))
        end if
        
        if ((Counters%PhreaticSurfaceNumberofSides == 0).and.(Counters%SoilSurfaceNumberofSides == 0))then
            PhreaticSurfaceNodeCoordMatrix = 0.0
            !
        else if ((Counters%PhreaticSurfaceNumberofSides == 0).and.(Counters%SoilSurfaceNumberofSides>0))then
            !
            If (CalParams%ApplyPartialSaturation) then
                if (NDIM == 3) then
                    PhreaticSurfaceNodeCoordMatrix(:,1)=xMax
                    PhreaticSurfaceNodeCoordMatrix(:,2)=yMin
                    PhreaticSurfaceNodeCoordMatrix(:,3)=zMin
                else if (NDIM == 2) then
                    PhreaticSurfaceNodeCoordMatrix(:,1)=xMax
                    PhreaticSurfaceNodeCoordMatrix(:,2)=yMin
                end if
                !
            else
                call SurfaceNodesCoordinates(PhreaticSurfaceNodeCoordMatrix, SoilSurfaceNodesConnectivities, Counters%SoilSurfaceNumberofSides)
            end if
            
        else if (Counters%PhreaticSurfaceNumberofSides > 0) then
            call SurfaceNodesCoordinates(PhreaticSurfaceNodeCoordMatrix, PhreaticSurfaceNodesConnectivities, Counters%PhreaticSurfaceNumberofSides)
        end if
          
          do IAEl = 1, Counters%NAEl
              IEl = ActiveElement(IAEl)
              
              do IPart = 1, NPartEle(IEl) ! Loop over all particles in element
                  ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
                  
                  ! initialize GaussPointPosition in the element
                  GaussPointGlobalCoord = 0
                  MaterialIndex = MaterialIDArray(ParticleIndex)
                  K0Value = MatParams(MaterialIndex)%K0Value
                
                         
                  if (K0Value<SKIP_K0_THRESHOLD) CYCLE ! skip material point
                     
                  ! Determine the location of IGaussPoint
                  call GaussPointLocalCoordinates(IGaussPoint, WeiGP, PosGP)
                     
                  ! Determine the shape function of the GaussPoint
                  call ShapeFunctionData(PosGP, IElTyp, GaussPointShapeValues, DShapeValues)
     
                  ! Determine the global location of IGaussPoint
                  call CoordGaussPointLocalToGlobal(IEl, NodalCoordinates, GaussPointShapeValues, GaussPointGlobalCoord)
     
                  !calculate total vertical stress, water pressure and for liquids also the density at Gauss Point location
                  if(NFORMULATION==1) then
                     call ComputeVerticalStressForK0(MaterialIndex, GaussPointGlobalCoord(1:2), UpdatedFluidDensity, SigY, SigWP, SoilSurfaceNodeCoordMatrix, PhreaticSurfaceNodeCoordMatrix)
                  else
                  
                    if(MaterialPointTypeArray(ParticleIndex)/=MaterialPointTypeLiquid) then
                      MaterialSetLayer = MaterialIDArray(ParticleIndex) !Assign material ID
                      FindMaterialID = .false. 
                    end if

                    call ComputeVerticalStressForK02Layers(GaussPointGlobalCoord(1:2), UpdatedFluidDensity, SigY, SigWP, ParticleIndex, MaterialSetLayer) !Double Point
                  end if
                    
                  if(NFORMULATION==1) then 
                  IsUndrTotalStress = trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_TOTAL
                  if (IsUndrTotalStress) then ! effectivestress=totalstress
                     SigYeff = SigY 
                     SigXeff = K0Value * SigYeff
                     SigWP = 0.0
                  elseif ((CalParams%NumberOfPhases==3).or. &
                     ((CalParams%NumberOfPhases==2).and.(CalParams%ApplyPartialSaturation))) then
                     !Calculate degree of saturation (needed for initialization considering Bishop stress principle)
                          if (MatParams(MaterialIndex)%RetentionCurve== SWRC_VANGENUCHTEN) then
                              call DynUpdateParticleDegreeOfSaturationVanGenuchten(ParticleIndex,MaterialIndex,Sr)
                          else if (MatParams(MaterialIndex)%RetentionCurve==SWRC_LINEAR) then
                              call DynUpdateParticleDegreeOfSaturationLinear(ParticleIndex,MaterialIndex,Sr)
                          end if

                      !calculate effective stresses at the Gauss Point location
                      SigYeff = SigY - SigWP*Sr !Bishop effective stress principle
                      SigXeff = K0Value * SigYeff
                   else
                      SigYeff = SigY - SigWP
                      SigXeff = K0Value * SigYeff
                   end if
              
                    ! assign initial and current stresses to the particle
                    if (NDIM == 3) then ! 3D case
                      SigmaEff0Array(ParticleIndex, 1) = SigXeff
                      SigmaEff0Array(ParticleIndex, 2) = SigYeff
                      SigmaEff0Array(ParticleIndex, 3) = SigXeff
                      SigmaEff0Array(ParticleIndex, 4) = 0.d0
                      SigmaEff0Array(ParticleIndex, 5) = 0.d0
                      SigmaEff0Array(ParticleIndex, 6) = 0.d0
                    else if (NDIM == 2) then ! 2D case
                      SigmaEff0Array(ParticleIndex, 1) = SigXeff
                      SigmaEff0Array(ParticleIndex, 2) = SigYeff
                      SigmaEff0Array(ParticleIndex, 3) = SigXeff
                      SigmaEff0Array(ParticleIndex, 4) = 0.d0
                    else
                      call GiveError('Dimension not defined. [subroutine InitialiseMaterialPointsForK0Stresses()].')
                    end if
                    
                    SigmaEff0 = SigmaEff0Array(ParticleIndex,:)
                    SigmaEffArray(ParticleIndex,:) =  SigmaEff0
                    Particles(ParticleIndex)%WaterPressure = SigWP
                    Particles(ParticleIndex)%WaterPressure0 = SigWP
                    

                    
                    ! assign current density to the liquid material point  
                    if (MatParams(MaterialIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') then
                      Particles(ParticleIndex)%Density = UpdatedFluidDensity
                      MassArray(ParticleIndex) = UpdatedFluidDensity * Particles(ParticleIndex)%IntegrationWeight
                      Particles(ParticleIndex)%FBody = MassArray(ParticleIndex) *  CalParams%GravityData%GAccel * CalParams%GravityData%GravityVector
                    end if ! only for liquid material points
                    
                  end if
                      
                  if(NFORMULATION==2) then       
                    if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then
                      !calculate effective stresses at the Gauss Point location
                      SigYeff = SigY - SigWP
                      SigXeff = K0Value * SigYeff
              
                      ! assign initial and current stresses to the particle
                      if (NDIM == 3) then ! 3D case
                        SigmaEff0Array(ParticleIndex, 1) = SigXeff
                        SigmaEff0Array(ParticleIndex, 2) = SigYeff
                        SigmaEff0Array(ParticleIndex, 3) = SigXeff
                        SigmaEff0Array(ParticleIndex, 4) = 0.d0
                        SigmaEff0Array(ParticleIndex, 5) = 0.d0
                        SigmaEff0Array(ParticleIndex, 6) = 0.d0
                      else if (NDIM == 2) then ! 2D case
                        SigmaEff0Array(ParticleIndex, 1) = SigXeff
                        SigmaEff0Array(ParticleIndex, 2) = SigYeff
                        SigmaEff0Array(ParticleIndex, 3) = SigXeff
                        SigmaEff0Array(ParticleIndex, 4) = 0.d0
                      else
                        call GiveError('Dimension not defined. [subroutine InitialiseMaterialPointsForK0Stresses()].')
                      end if    
                      SigmaEff0 = SigmaEff0Array(ParticleIndex,:)
                      SigmaEffArray(ParticleIndex,:) =  SigmaEff0
                    end if
                      
                    if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then 
                      ! assign initial and current stresses to the particle
                      SigmaEff0 = 0.0
                      if (NDIM == 3) then ! 3D case
                        SigmaEff0Array(ParticleIndex, 1:3) = SigWP
                        SigmaEff0Array(ParticleIndex, 4:6) = 0.d0                    
                        SigmaEffArray(ParticleIndex, 1:3) = SigWP
                        SigmaEffArray(ParticleIndex, 4:6) = 0.d0
                      else if (NDIM == 2) then ! 2D case
                        SigmaEff0Array(ParticleIndex, 1:3) = SigWP
                        SigmaEff0Array(ParticleIndex, 4) = 0.d0                    
                        SigmaEffArray(ParticleIndex, 1:3) = SigWP
                        SigmaEffArray(ParticleIndex, 4) = 0.d0
                      else
                        call GiveError('Dimension not defined. [subroutine InitialiseMaterialPointsForK0Stresses()].')
                      end if    
                      Particles(ParticleIndex)%WaterPressure = SigWP
                      Particles(ParticleIndex)%WaterPressure0 = SigWP
                    end if
                        
                  end if
                     
              end do ! Loop over the particles within the element
            
          end do ! Loop over all elements

        end subroutine InitialiseMaterialPointsForK0Stresses 


              
        subroutine AssignTractionToEntity() 
        !**********************************************************************
        !
        ! Function:  To assign the traction to the corresponding entity
        !
        ! NOTE : This function works ONLY if NumberOfLayers = 1
        !
        !**********************************************************************
        
        implicit none

          ! local variables
          integer(INTEGER_TYPE) :: IDof, IParticle, LoadedEntityID, I
          logical :: Terminate
          logical, dimension(Counters%NParticles) :: LoadedParticle
          integer(INTEGER_TYPE), dimension(:), allocatable :: ConnectivitiesSide
          real(REAL_TYPE), dimension(:, :), allocatable :: LoadValuesSide
          
  
          
          LoadedParticle = .false. 
          ExtLoad = 0.0
          ExtLoadTotal = 0.0
          Terminate = .false.
          LoadedEntityID = 1
          
          if (CalParams%ApplyContactAlgorithm) then !if contact is used then there are 2 entities!
            !SOLID TRACTION
            if ( Counters%NLoadedElementSidesSolid > 0 ) then ! solid   
             !Load system A   
              allocate( ConnectivitiesSide( size(LoadOnNodesConnectivitiesSolid,1) ) )
              allocate( LoadValuesSide( size(LoadValuesOnNodesSolid,1), size(LoadValuesOnNodesSolid,2) ) )  
              do I = 1, Counters%NLoadedElementSidesSolidNodes 
                ConnectivitiesSide = LoadOnNodesConnectivitiesSolid(:, I) ! extract values for each loaded element side I
                LoadValuesSide = LoadValuesOnNodesSolid(:, :, I)
                call DetermineParticlesInLoadedElement(ConnectivitiesSide, LoadValuesSide, LoadedParticle)
              end do
              
              do IParticle =1, Counters%Nparticles    ! Entity of loaded elements is determined by the first particle in the array LoadedParticle   
                if (LoadedParticle(IParticle)) then   ! in future loads applied to different entities should be considered
                  LoadedEntityID = EntityIDArray(IParticle) !Particles(IParticle)%EntityID
                  Terminate = .true.
                  if (Terminate) goto 1
                end if  
              end do
              
1           do IDof = 1, Counters%N
              ExtLoad(IDof, LoadedEntityID) = ExtLoadVector(IDof)
              ExtLoadTotal(IDof, LoadedEntityID,1) = ExtLoadVector(IDof)
            end do 
            deallocate( ConnectivitiesSide )
            deallocate( LoadValuesSide)
            
            !Load system B
            if (Counters%NLoadedElementSidesSolidNodesB>0 ) then
            allocate( ConnectivitiesSide( size(LoadOnNodesConnectivitiesSolidB,1) ) )
             allocate( LoadValuesSide( size(LoadValuesOnNodesSolidB,1), size(LoadValuesOnNodesSolidB,2) ) )  
              do I = 1, Counters%NLoadedElementSidesSolidNodesB
                ConnectivitiesSide = LoadOnNodesConnectivitiesSolidB(:, I) ! extract values for each loaded element side I
                LoadValuesSide = LoadValuesOnNodesSolidB(:, :, I)
                call DetermineParticlesInLoadedElement(ConnectivitiesSide, LoadValuesSide, LoadedParticle)
              end do
              
              do IParticle =1, Counters%Nparticles    ! Entity of loaded elements is determined by the first particle in the array LoadedParticle   
                if (LoadedParticle(IParticle)) then   ! in future loads applied to different entities should be considered
                  LoadedEntityID = EntityIDArray(IParticle) !Particles(IParticle)%EntityID
                  Terminate = .true.
                  if (Terminate) goto 11
                end if  
              end do
              
11          do IDof = 1, Counters%N
              ExtLoadTotal(IDof, LoadedEntityID,2) = ExtLoadVectorB(IDof)
            end do 
            deallocate( ConnectivitiesSide )
            deallocate( LoadValuesSide)
            end if
            end if !solid traction
            
            !LIQUID TRACTION
            if ( Counters%NLoadedElementSidesWater > 0 ) then ! liquid
            !Load system water A    
               allocate( ConnectivitiesSide( size(LoadOnNodesConnectivitiesWater,1) ) )
               allocate( LoadValuesSide( size(LoadValuesOnNodesWater,1), size(LoadValuesOnNodesWater,2) ) )  
              do I = 1, Counters%NLoadedElementSidesWaterNodes                
                ConnectivitiesSide = LoadOnNodesConnectivitiesWater(:, I) ! extract values for each loaded element side I
                LoadValuesSide = LoadValuesOnNodesWater(:, :, I)
                call DetermineParticlesInLoadedElement(ConnectivitiesSide, LoadValuesSide, LoadedParticle)
              end do
              
              do IParticle =1, Counters%Nparticles    ! Entity of loaded elements is determined by the first particle in the array LoadedParticle   
                if (LoadedParticle(IParticle)) then   ! in future loads applied to different entities should be considered
                  LoadedEntityID = EntityIDArray(IParticle) !Particles(IParticle)%EntityID
                  Terminate = .true.
                  if (Terminate) goto 2
                end if  
              end do
              
2             do IDof = 1, Counters%N
                  ExtLoadWaterTotal (IDof, LoadedEntityID,1) = ExtLoadVectorWater(IDof)
              end do
    
              deallocate( ConnectivitiesSide )
              deallocate( LoadValuesSide )
             
              !Load system water B
              if (Counters%NLoadedElementSidesWaterNodesB>0 ) then
              allocate( ConnectivitiesSide( size(LoadOnNodesConnectivitiesWaterB,1) ) )
              allocate( LoadValuesSide( size(LoadValuesOnNodesWaterB,1), size(LoadValuesOnNodesWaterB,2) ) )  
              do I = 1, Counters%NLoadedElementSidesWaterNodesB                
                ConnectivitiesSide = LoadOnNodesConnectivitiesWaterB(:, I) ! extract values for each loaded element side I
                LoadValuesSide = LoadValuesOnNodesWaterB(:, :, I)
                call DetermineParticlesInLoadedElement(ConnectivitiesSide, LoadValuesSide, LoadedParticle)
              end do
              
              do IParticle =1, Counters%Nparticles    ! Entity of loaded elements is determined by the first particle in the array LoadedParticle   
                if (LoadedParticle(IParticle)) then   ! in future loads applied to different entities should be considered
                  LoadedEntityID = EntityIDArray(IParticle) !Particles(IParticle)%EntityID
                  Terminate = .true.
                  if (Terminate) goto 21
                end if  
              end do
              
21            do IDof = 1, Counters%N
                  ExtLoadWaterTotal (IDof, LoadedEntityID,2) = ExtLoadVectorWaterB(IDof)
               end do
    
              deallocate( ConnectivitiesSide )
              deallocate( LoadValuesSide )
              end if
            end if !Liquid traction
            
            !Total head
            if ( Counters%HydraulicHeadSides > 0 ) then ! liquid total head
                
             do LoadedEntityID=1,Counters%NEntity
               do IDof=1,Counters%N
                 HydraulicHeadLoadTotal(IDof, LoadedEntityID) = HydraulicHeadVector(IDof)
               end do
             end do
              
            end if !total head
            
            !GAS TRACTION
            if ( Counters%NLoadedElementSidesGas > 0 ) then ! gas
                !Load System A
              allocate( ConnectivitiesSide( size(LoadOnNodesConnectivitiesGas,1) ) )
              allocate( LoadValuesSide( size(LoadValuesOnNodesGas,1), size(LoadValuesOnNodesGas,2) ) )
              do I = 1, Counters%NLoadedElementSidesGasNodes                
                ConnectivitiesSide = LoadOnNodesConnectivitiesGas(:, I) ! extract values for each loaded element side I 
                LoadValuesSide = LoadValuesOnNodesGas(:, :, I)
                call DetermineParticlesInLoadedElement(ConnectivitiesSide, LoadValuesSide, LoadedParticle)
              end do
              do IParticle =1, Counters%Nparticles    ! Entity of loaded elements is determined by the first particle in the array LoadedParticle   
                if (LoadedParticle(IParticle)) then   ! in future loads applied to different entities should be considered
                  LoadedEntityID = EntityIDArray(IParticle) !Particles(IParticle)%EntityID
                  Terminate = .true.
                  if (Terminate) goto 4
                end if  
              end do
              
4             do IDof = 1, Counters%N
                  ExtLoadGasTotal (IDof, LoadedEntityID,1) = ExtLoadVectorGas(IDof)
              end do
              
              deallocate( ConnectivitiesSide )
              deallocate( LoadValuesSide )
              
              !Load system B
              if (Counters%NLoadedElementSidesGasNodesB>0 ) then
              allocate( ConnectivitiesSide( size(LoadOnNodesConnectivitiesGasB,1) ) )
              allocate( LoadValuesSide( size(LoadValuesOnNodesGasB,1), size(LoadValuesOnNodesGasB,2) ) )
              do I = 1, Counters%NLoadedElementSidesGasNodesB                
                ConnectivitiesSide = LoadOnNodesConnectivitiesGasB(:, I) ! extract values for each loaded element side I 
                LoadValuesSide = LoadValuesOnNodesGasB(:, :, I)
                call DetermineParticlesInLoadedElement(ConnectivitiesSide, LoadValuesSide, LoadedParticle)
              end do
              do IParticle =1, Counters%Nparticles    ! Entity of loaded elements is determined by the first particle in the array LoadedParticle   
                if (LoadedParticle(IParticle)) then   ! in future loads applied to different entities should be considered
                  LoadedEntityID = EntityIDArray(IParticle) !Particles(IParticle)%EntityID
                  Terminate = .true.
                  if (Terminate) goto 41
                end if  
              end do
              
41            do IDof = 1, Counters%N
                  ExtLoadGasTotal (IDof, LoadedEntityID,2) = ExtLoadVectorGasB(IDof)
              end do
              
              deallocate( ConnectivitiesSide )
              deallocate( LoadValuesSide )
                         
              end if
            end if !gas traction

            else !No contact algorithm  

                do IDof = 1, Counters%N
                    if (CalParams%NumberOfPhases>=2) then
                        ExtLoadWaterTotal (IDof, LoadedEntityID,1) = ExtLoadVectorWater(IDof)
                        if (Counters%NWaterLoadSystems>1)  then
                            ExtLoadWaterTotal (IDof, LoadedEntityID,2) = ExtLoadVectorWaterB(IDof) !Load system B
                        end if
                        HydraulicHeadLoadTotal(IDof, LoadedEntityID) = HydraulicHeadVector(IDof)
                    end if

                    if (CalParams%NumberOfPhases==3) then  ! consolidation
                        ExtLoadGasTotal (IDof, LoadedEntityID,1) = ExtLoadVectorGas(IDof)
                        if (Counters%NGasLoadSystems>1)  then
                            ExtLoadGasTotal (IDof, LoadedEntityID,2) = ExtLoadVectorGasB(IDof) !Load system B
                        end if
                    end if

                    ExtLoad(IDof, LoadedEntityID) = ExtLoadVector(IDof)
                    ExtLoadTotal(IDof, LoadedEntityID,1) = ExtLoadVector(IDof)
                    if (Counters%NSoilLoadSystems>1)  then
                        ExtLoadTotal(IDof, LoadedEntityID,2) = ExtLoadVectorB(IDof) !Load system B
                    end if
                end do
           end if !contact or not

        
        end subroutine AssignTractionToEntity
        

        subroutine MapDampingFromParticles()
        !**********************************************************************
        !
        ! Function:  To calculate the damping force solid
        !
        !**********************************************************************

        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: I,J, NElemPart, IPart, ParticleIndex, IEl, INode, NN, EntityID, IEnt
          integer(INTEGER_TYPE), dimension(Counters%nEntity) :: NPart
          integer(INTEGER_TYPE), dimension(Counters%NodTot, Counters%nEntity) :: NodeTimes
          integer(INTEGER_TYPE), dimension(NVECTOR) :: Ni
          
          real(REAL_TYPE) :: DampedOutOffBalance, Damping
          real(REAL_TYPE), dimension(Counters%NEl, Counters%nEntity) :: AverageDamping
          real(REAL_TYPE), dimension(Counters%NodTot, Counters%nEntity) :: NodalDamping
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: Dampingfactors  
          
          NodeTimes = 0
          AverageDamping = 0.0
          NodalDamping = 0.0

          do IEl = 1, Counters%NEl ! loop over elements
            
            if (IsActiveElement(IEl)) then ! only for active elements
              NElemPart = NPartEle(IEl) ! number of material points in element
              NPart = 0
              
              do IPart = 1, NElemPart ! loop over material points
                ParticleIndex = GetParticleIndex(IPart, IEl)

                if (.not.CalParams%ApplyContactAlgorithm) then ! NO CONTACT
                  EntityID = 1
                else ! CONTACT IS APPLIED
                  EntityID = EntityIDArray(ParticleIndex) 
                end if ! wheter or not ContactAlgorithm

                NPart(EntityID) = NPart(EntityID) + 1
                Damping = Particles(ParticleIndex)%Damping
                AverageDamping(IEl,EntityID) = AverageDamping(IEl,EntityID) + Damping
              end do ! loop over material points
              
              do IEnt = 1, Counters%nEntity ! loop over entities
                if (NPart(IEnt)/=0) then
                  AverageDamping(IEl,IEnt) = AverageDamping(IEl,IEnt) / NPart(IEnt)
                end if
              end do ! loop over entities
               
              do INode=1,ELEMENTNODES ! loop over element nodes
                NN=ElementConnectivities(INode,IEl) ! get global node number
                do IEnt = 1, Counters%nEntity  ! loop over entities
                  if (NPart(IEnt)>0) then
                    NodalDamping(NN,IEnt) = NodalDamping(NN,IEnt) + AverageDamping(IEl,IEnt)
                    NodeTimes(NN,IEnt) = NodeTimes(NN,IEnt) + 1 
                  end if
                end do ! loop over entities
              end do ! loop over element nodes
           
            end if ! only for active elements
           
          end do ! loop over elements
           
          do INode=1,Counters%NodTot ! loop over nodes
            do IEnt = 1, Counters%nEntity ! loop over entities
              
              if (NodeTimes(INode, IEnt)/=0) then
                NodalDamping(INode,IEnt) = NodalDamping (INode, IEnt) / NodeTimes(INode, IEnt)
              end if
              
              do I = 1, NVECTOR
                Ni = ReducedDof(INode) + I   ! global storage coordinate of x-val
                Dampingfactors(Ni,IEnt) = NodalDamping(INode,IEnt)
              end do

            end do ! loop over entities
          end do ! loop over nodes
                
          do I = 1, Counters%N ! loop over degrees-of-freedom
            do J = 1, Counters%nEntity ! loop over entities
              if (TotalVelocitySoil(I,J)/=0.0) then ! needed, beacuse sign() will return (+1) if NodalVelocities(I) = 0, and not zero as expected
                DampedOutOffBalance = - Dampingfactors(I,J) * abs(RateofMomentum(I,J)) * sign(1.d0, TotalVelocitySoil(I,J))
                RateofMomentum(I,J) = RateofMomentum(I,J) + DampedOutOffBalance ! apply damping force
              end if 
            end do ! loop over entities
          end do ! loop over degrees-of-freedom

        end subroutine MapDampingFromParticles
        
         
        subroutine MapWaterDampingFromParticles(RateofMomentumW)
        !**********************************************************************
        !
        ! Function:  To calculate the damping force water 
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: RateofMomentumW
          
          ! Local variables
          integer(INTEGER_TYPE) :: I,J, NElemPart, IPart, ParticleIndex, IEl, INode, NN, EntityID, IEnt
          integer(INTEGER_TYPE), dimension(Counters%nEntity) :: NPart
          integer(INTEGER_TYPE), dimension(Counters%NodTot, Counters%nEntity) :: NodeTimes
          integer(INTEGER_TYPE), dimension(NVECTOR) :: Ni
          
          real(REAL_TYPE) :: DampedOutOffBalance, Damping
          real(REAL_TYPE), dimension(Counters%NEl, Counters%nEntity) :: AverageDamping
          real(REAL_TYPE), dimension(Counters%NodTot, Counters%nEntity) :: NodalDamping
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: Dampingfactors  
          
          NodeTimes = 0
          AverageDamping = 0.0
          NodalDamping = 0.0
          
          do IEl = 1, Counters%NEl ! loop over elements
            
            if (IsActiveElement(IEl)) then ! only for active elements
              NElemPart = NPartEle(IEl)  ! number of material points in element
              NPart = 0
              
              do IPart = 1, NElemPart ! loop over material points
                ParticleIndex = GetParticleIndex(IPart,IEl)
                if((NFORMULATION==1).or. &
                        (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid)) then !NumbOfLayers = 1 or LIQUID Material Point
                  if (.not.CalParams%ApplyContactAlgorithm) then ! NO CONTACT
                    EntityID = 1
                  else ! CONTACT IS APPLIED
                    EntityID = EntityIDArray(ParticleIndex)
                  end if     
                  NPart(EntityID) = NPart (EntityID)+1
                  Damping = Particles(ParticleIndex)%Damping
                  AverageDamping(IEl,EntityID) = AverageDamping(IEl,EntityID) + Damping
                end if !NumbOfLayers = 1 or LIQUID Material Point
              end do ! loop over material points
              
              do IEnt = 1, Counters%nEntity ! loop over entities
                if (NPart(IEnt)/=0) then
                  AverageDamping(IEl,IEnt) = AverageDamping(IEl,IEnt) / NPart(IEnt)
                  end if
              end do ! loop over entities
               
              do INode=1,ELEMENTNODES ! loop over element nodes
                NN=ElementConnectivities(INode,IEl) ! get global node number
                do IEnt = 1, Counters%nEntity ! loop over entities 
                  if (NPart(IEnt)>0) then
                    NodalDamping(NN,IEnt) = NodalDamping (NN,IEnt) + AverageDamping(IEl,IEnt)
                    NodeTimes(NN,IEnt) = NodeTimes(NN,IEnt) + 1 
                  end if
                end do ! loop over entities
              end do ! loop over element nodes
            
            end if ! only for active elements
          end do ! loop over elements
           
          do INode=1,Counters%NodTot ! loop over nodes
            do IEnt = 1, Counters%nEntity ! loop over entities
              
              if (NodeTimes(INode,IEnt)/=0) then
                NodalDamping(INode,IEnt) = NodalDamping(INode,IEnt) / NodeTimes(INode,IEnt)
              end if
             
              do I = 1, NVECTOR
                Ni = ReducedDof(INode) + I   ! global storage coordinate of x-val
                Dampingfactors(Ni,IEnt) = NodalDamping(INode,IEnt)
              end do  

            end do ! loop over entities
          end do ! loop over nodes

          do I = 1, Counters%N ! loop over degrees-of-freedom
            do J = 1, Counters%nEntity ! loop over entities
              if (TotalVelocityWater(I,J)/=0.0) then ! needed, beacuse sign() will return (+1) if NodalVelocities(I) = 0, and not zero as expected
                DampedOutOffBalance = - Dampingfactors(I,J) * abs(RateofMomentumW(I,J)) * sign(1.d0, TotalVelocityWater(I,J))
                RateofMomentumW(I,J) = RateofMomentumW(I,J) + DampedOutOffBalance !apply damping force
              end if
            end do ! loop over entities
          end do ! loop over degrees-of-freedom

        end subroutine MapWaterDampingFromParticles


        subroutine MapGasDampingFromParticles(RateofMomentumG)
        !**********************************************************************
        !
        !   Function:  To calculate the damping force gas 
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N,Counters%nEntity), intent(inout) :: RateofMomentumG
          
          ! Local variables
          integer(INTEGER_TYPE) :: I,J, NElemPart, IPart, ParticleIndex, IEl, INode, NN, EntityID, IEnt
          integer(INTEGER_TYPE), dimension(Counters%nEntity) :: NPart
          integer(INTEGER_TYPE), dimension(Counters%NodTot, Counters%nEntity) :: NodeTimes
          integer(INTEGER_TYPE), dimension(NVECTOR) :: Ni
          
          real(REAL_TYPE) :: DampedOutOffBalance, Damping
          real(REAL_TYPE), dimension(Counters%NEl, Counters%nEntity) :: AverageDamping
          real(REAL_TYPE), dimension(Counters%NodTot, Counters%nEntity) :: NodalDamping
          real(REAL_TYPE), dimension(Counters%N, Counters%nEntity) :: Dampingfactors  
          
          NodeTimes = 0
          AverageDamping = 0.0
          NodalDamping = 0.0
          
          do IEl = 1, Counters%NEl ! loop over elements
            
            if (IsActiveElement(IEl)) then ! only active elements
              NElemPart = NPartEle(IEl)  ! number of material points in element
              NPart = 0
              
              do IPart = 1, NElemPart ! loop over material points
                ParticleIndex = GetParticleIndex(IPart,IEl)
                if (.not.CalParams%ApplyContactAlgorithm) then ! NO CONTACT
                  EntityID = 1
                else ! CONTACT IS APPLIED
                  EntityID = EntityIDArray(ParticleIndex)
                end if ! wheter or not ContactAlgorithm
                NPart(EntityID) = NPart(EntityID) + 1
                Damping = Particles(ParticleIndex)%Damping
                AverageDamping(IEl,EntityID) = AverageDamping(IEl,EntityID) + Damping
              end do ! loop over material points
              
              do IEnt = 1, Counters%nEntity ! loop over entities
                if (NPart(IEnt)/=0) then
                  AverageDamping(IEl,IEnt) = AverageDamping(IEl,IEnt) / NPart(IEnt)
                end if
              end do ! loop over entities
              
              do INode=1,ELEMENTNODES ! loop over element nodes
                NN=ElementConnectivities(INode,IEl) ! get global node number
                do IEnt = 1, Counters%nEntity ! loop over entities
                  if (NPart(IEnt)>0) then
                    NodalDamping(NN,IEnt) = NodalDamping(NN,IEnt) + AverageDamping(IEl,IEnt)
                    NodeTimes(NN,IEnt) = NodeTimes(NN,IEnt) + 1 
                  end if
                end do ! loop over entities
              end do ! loop over element nodes

            end if ! only for active elements
          end do ! loop over elements
           
          do INode=1,Counters%NodTot ! loop over nodes
            do IEnt = 1, Counters%nEntity ! loop over entities
              
              if (NodeTimes(INode,IEnt)/=0) then
                NodalDamping(INode,IEnt) = NodalDamping(INode,IEnt) / NodeTimes(INode,IEnt)
              end if

              do I = 1, NVECTOR
                Ni = ReducedDof(INode) + 1 ! global storage coordinate of x-val
                Dampingfactors(Ni,IEnt) = NodalDamping(INode,IEnt)
              end do   

            end do ! loop over entities
          end do ! loop over nodes
                
          do I = 1, Counters%N ! loop over degrees-of-freedom
            do J = 1, Counters%nEntity ! loop over entities
              if (TotalVelocityGas(I,J)/=0.0) then ! needed, beacuse sign() will return (+1) if NodalVelocities(I) = 0, and not zero as expected
                DampedOutOffBalance = - Dampingfactors(I,J) * abs(RateofMomentumG(I,J)) * sign(1.d0, TotalVelocityGas(I,J))
                RateofMomentumG(I,J) = RateofMomentumG(I,J) + DampedOutOffBalance !apply damping force
              end if 
            end do ! loop over entities
          end do ! loop over degrees-of-freedom
        
        end subroutine MapGasDampingFromParticles         
        
        subroutine ApplyMPPrescribedVelocity(Array)
        !**********************************************************************
        !
        ! Function:  To assign prescribed velocity to structure particles
        !
        ! Note : This subroutine works only if NumbOfLayers = 1
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(inout):: Array
        
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IDir
          integer(INTEGER_TYPE) :: IEntity
          real(REAL_TYPE), dimension(Counters%N) :: NodalPrescribedVelocityFromMP, ISPrescribedDof
          logical :: IsPrescribedVelocity !0=velocity is prescribed, 1=no prescried velocity

          if (.not.(NFORMULATION==1)) RETURN
            
          IsPrescribedVelocity = CalParams%PrescribedVelo%ApplyPrescribedVelo
          if (.not.IsPrescribedVelocity) RETURN
          
            do IParticle = 1, Counters%NParticles ! loop over material points
              ParticleIndex = GetParticleIndexFromList(IParticle)
              do IDir = 1, NVECTOR
                if (Particles(ParticleIndex)%MPPrescribedVeloDir(IDir)) then
                  VelocityArray(ParticleIndex,IDir) = GetParticlePrescribedVelocityI(IDir, ParticleIndex)
                  AccelerationArray(ParticleIndex, IDir) = GetParticlePrescribedAccelerationI(IDir, ParticleIndex) !It is not needed further in the computation
                  if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
                    VelocityWaterArray(ParticleIndex,IDir) = GetParticlePrescribedVelocityI(IDir, ParticleIndex)                    
                  end if
                end if  
              end do
            end do
            
            !Loop over the nodes of the structure to assign the prescribed velocity
            !NB: Prescribed velocity only works together with moving mesh
            if(CalParams%ApplyContactAlgorithm) then
             IEntity = getMaterialEntity(CalParams%MovingMesh%MovingMaterialID)
            else
             IEntity = 1
            end if
            NodalPrescribedVelocityFromMP = 0.0
            call MapMaterialPointPrescribedVelocityToNodes(NodalPrescribedVelocityFromMP,ISPrescribedDof)
            TotalVelocitySoil(:, IEntity) = TotalVelocitySoil(:, IEntity)*ISPrescribedDof(:) + &
                                            NodalPrescribedVelocityFromMP(:)
            AccelerationSoil(:,IEntity) = AccelerationSoil(:,IEntity) * ISPrescribedDof(:)
            if (CalParams%ApplyContactAlgorithm) then
              Array(:) = Array(:)*ISPrescribedDof(:) + NodalPrescribedVelocityFromMP(:)
            end if
            if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
              AccelerationWater(:, IEntity) = AccelerationWater(:, IEntity) * ISPrescribedDof(:)
              TotalVelocityWater(:, IEntity) = TotalVelocityWater(:, IEntity) * ISPrescribedDof(:) + &
                                              NodalPrescribedVelocityFromMP(:)
            end if
        
        end subroutine ApplyMPPrescribedVelocity
        
       subroutine MapMaterialPointPrescribedVelocityToNodes(NodalPrescribedVelocity,ISPrescribedDof)
                              
        !**********************************************************************
        !
        !    Function:  map the prescribed velocity to the nodes
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(inout) :: NodalPrescribedVelocity !value of prescribed velocity at Dof
          real(REAL_TYPE), dimension(Counters%N), intent(inout) :: ISPrescribedDof  !0=velocity is prescribed, 1=no prescried velocity
          ! Local variables
          integer(INTEGER_TYPE) :: I, IEl, IPart, ParticleIndex, J
          integer(INTEGER_TYPE), dimension(NVECTOR, ELEMENTNODES) :: IDof
          real(REAL_TYPE), dimension(NVECTOR) :: ParticleVelocity
          real(REAL_TYPE), dimension(Counters%N) :: PrescribedMomentum, Mass
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
          real(REAL_TYPE) ::  PartilceMass
          
          !!CC - changed to function call - needs to be set to zero here
          PrescribedMomentum = 0.0
          Mass = 0.0
          ISPrescribedDof = 1.0

          do J = 1, GeoParams%PrescribedVeloNElem
            IEL = GeoParams%PrescribedVeloElemID(J)
            do I = 1, NVECTOR
              IDof(I, 1:ELEMENTNODES) = ReducedDof(ElementConnectivities(1:ELEMENTNODES, IEl)) + I
              if (GeoParams%PrescribedVeloElDirection(J,I)==1) then
               ISPrescribedDof(IDof(I, 1:ELEMENTNODES)) = 0.0
              end if 
            end do
            
            do IPart = 1, NPartEle(IEl)                                 ! Loop over all particles in element
              ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the particle ID
              ParticleVelocity = VelocityArray(ParticleIndex,:)
              PartilceMass = MassArray(ParticleIndex)
              ParticleShape = ShapeValuesArray(ParticleIndex,:)
                do I = 1, NVECTOR ! nodal i-momentum
                  if (Particles(ParticleIndex)%MPPrescribedVeloDir(I)) then  
                    PrescribedMomentum(IDof(I,1:ELEMENTNODES)) = PrescribedMomentum(IDof(I,1:ELEMENTNODES)) &
                        + PartilceMass * ParticleShape * ParticleVelocity(I)
                    Mass(IDof(I,1:ELEMENTNODES)) = Mass(IDof(I,1:ELEMENTNODES)) + PartilceMass * ParticleShape
                  end if  !prescribed direction
                end do     
            end do !Loop over particles
          end do !elements
          
          do I=1, Counters%N
            if (Mass(I) /= 0.0) then
              NodalPrescribedVelocity(I) = PrescribedMomentum(I)/Mass(I)
            end if
          end do  
          
         end subroutine MapMaterialPointPrescribedVelocityToNodes
        
        subroutine ApplyNodalPrescribedVelocity(Array)
        !**********************************************************************
        !
        ! Function :  Assign prescribed velocity to nodes
        !
        ! Note : This function works only if NumbOfLayers = 1 and code versions newer than 2018.2
        !
        !**********************************************************************

        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(inout):: Array
        
          ! local variables
          integer(INTEGER_TYPE) :: INode, IDof, I, IEntity, NodeID
         
        if (.not.(NFORMULATION==1)) RETURN
        if (.not.CalParams%PrescribedVelo%ApplyPrescribedVelo) RETURN
        if (CalParams%PrescribedVelo%NNodePrescribedVelo<=0) RETURN !NO prescribed velocity at the nodes, return
          
        do IEntity=1,Counters%NEntity    !loop over entity, prescribed nodal velocity is applied to all entity
          do INode=1, CalParams%PrescribedVelo%NNodePrescribedVelo
              
            NodeID = CalParams%PrescribedVelo%NodePrescribedVelo(INode)
            IDof = ReducedDof(NodeID)
            do I=1,NVECTOR
              if (CalParams%PrescribedVelo%NodalPrescribedVelocityDirection(INode,I)==1) then 
                        if  (CalParams%ApplyContactAlgorithm) then
                            if (Ientity==HARD_ENTITY) then
                                TotalVelocitySoil(IDof + I, IEntity) = GetNodalPrescribedVelocityI(I, INode)
                                AccelerationSoil(IDof + I, IEntity) = GetNodalPrescribedAccelerationI(I, INode) !Not needed further in the computation
                            else !soil entity
                                if (.not.InterfaceNodes(NodeID)) then !Node is not in the contact surface
                                    TotalVelocitySoil(IDof + I, IEntity) = GetNodalPrescribedVelocityI(I, INode)
                                    AccelerationSoil(IDof + I, IEntity) = GetNodalPrescribedAccelerationI(I, INode) !Not needed further in the computation
                                end if
                            end if
                        else
                            TotalVelocitySoil(IDof + I, IEntity) = GetNodalPrescribedVelocityI(I, INode)
                            AccelerationSoil(IDof + I, IEntity) = GetNodalPrescribedAccelerationI(I, INode) !Not needed further in the computation
                        end if
              if (CalParams%ApplyContactAlgorithm) then
                Array(IDof + I) = GetNodalPrescribedVelocityI(I, INode)
              end if  
                
              if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
                  AccelerationWater(IDof + I, IEntity) = GetNodalPrescribedAccelerationI(I, INode)
                  TotalVelocityWater(IDof + I, IEntity) = GetNodalPrescribedVelocityI(I, INode)
              end if
              
              end if
              
            end do
          end do
        end do
            
          
        end subroutine ApplyNodalPrescribedVelocity 
        
      subroutine DetermineContactSurfaceSoilElements()
        !**********************************************************************
        !
        !  Function: Determine the elements which have a node on the contact surface
        !             and are partially filled. It is needed to smoothen the nodal 
        !             acceleration at the contact surface
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
      implicit none
      
        integer(INTEGER_TYPE) :: I, J, IEl, NAdjacentElements, IError
        
        if (allocated(ContactSurfaceSoilElements)) then
          deallocate(ContactSurfaceSoilElements, stat = IError)
        end if
        allocate(ContactSurfaceSoilElements(Counters%NEl), stat = IError)
        ContactSurfaceSoilElements = .false.
        
        if (.not.CalParams%ApplyContactAlgorithm) RETURN
        
        do I = 1, Counters%NodTot
          if (InterfaceNodes(I)) then
            NAdjacentElements = GetNElmOfNode(I)
            do J = 1, NAdjacentElements
              IEl = GetElmIOfNode(I, J)
              if (abs(ElementMaterialID(IEl))==SOFT_ENTITY) then
                ContactSurfaceSoilElements(IEl) = .true.
              end if
            end do
          end if
        end do
        
      end subroutine DetermineContactSurfaceSoilElements

      
      subroutine GetNodalUnitMassGradient(NodalUnitMassGradient, Bcond)
      !**********************************************************************
      !
      !Function: calculate the unit normal mass gradient, which correspond to the normal to the surface
      !
      !**********************************************************************
      implicit none
      
       real(REAL_TYPE), dimension(:,:), intent(inout):: NodalUnitMassGradient
        real(REAL_TYPE), dimension(:), intent(in):: Bcond
       real(REAL_TYPE), dimension(Counters%NodTot,NDIM):: MassGradient
       real(REAL_TYPE), dimension(NDIM,ELEMENTNODES):: B
       real(REAL_TYPE), dimension(NDIM):: LocPos, NodeMassGradient
       real(REAL_TYPE):: Mass, Vlength, DetJac
       integer(INTEGER_TYPE):: IAElem, IElem, IPart, PartIndex, INode, NodeID, IDim,IDof
       
      
      
      MassGradient = 0.0
      do IAElem = 1, Counters%NAEl
          IElem = ActiveElement(IAElem)
          do IPart = 1,NPartEle(IElem)
              PartIndex = GetParticleIndex(IPart,IElem)
              Mass = MassArray(PartIndex)
              LocPos = Particles(PartIndex)%LocPos
              call BMatrix(LocPos, ELEMENTNODES, Counters%NEl, &
                               Counters%NodTot, NDIM, &
                               IElem, ElementConnectivities, &
                               NodalCoordinates, B, DetJac)
              
              do INode = 1, ELEMENTNODES
                  NodeID = ElementConnectivities(INode,IElem)
                  IDof = ReducedDof(NodeID)
                  do IDim = 1, NDIM
                      MassGradient(NodeID,IDim) = MassGradient(NodeID,IDim) + Mass*B(IDim,INode)*Bcond(IDof+IDim)
                  end do                  
              end do
          end do
          
       end do
          
       do NodeID=1,Counters%NodTot   
           NodeMassGradient(:) = MassGradient(NodeID,:)
           Vlength = TINY
           do IDim = 1, NDIM
              Vlength = Vlength + NodeMassGradient(IDim)* NodeMassGradient(IDim)
           end do
           Vlength = sqrt(Vlength)
           NodeMassGradient(:) = NodeMassGradient(:)/Vlength
           NodalUnitMassGradient(NodeID,:)=NodeMassGradient(:)

       end do   
      
      end subroutine GetNodalUnitMassGradient
      
      subroutine RecalculateAccelerationFromCorrectedVelocity(Velocity,VelocityPrevious,Acceleration)
     !*********************************************************************
     ! Function: recalculate acceleration after velocity correction
     !
     !*********************************************************************
                      
        implicit none

        !logical, dimension(:), intent(in) :: IsInfiltrationNode
        real(REAL_TYPE), dimension(:,:), intent(in):: VelocityPrevious
        real(REAL_TYPE), dimension(:,:), intent(inout):: Velocity, Acceleration
        integer(INTEGER_TYPE):: INode, IDof, IDim, IEntity
        

        do INode=1,Counters%NodTot
            do IEntity=1,Counters%NEntity
                IDof = ReducedDof(INode)
                do IDim =1,NVECTOR
                    Acceleration(IDoF + IDim, IEntity) =  (Velocity(IDoF + IDim, IEntity) - VelocityPrevious(IDoF + IDim, IEntity)) /  CalParams%TimeIncrement
                end do
            end do
        end do
    
    end  subroutine RecalculateAccelerationFromCorrectedVelocity
    
               
            
      end module ModLagrangianPhase
