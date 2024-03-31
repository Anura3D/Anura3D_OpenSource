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
 !   and soil�water�structure interaction using the material point method (MPM)
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


module ModMPMDYNBTSigma
   !**********************************************************************
   !
   ! Function:  Contains the routine for integration of the BT*Sigma - assemblage
   !            of the internal force vector For dynamic MPM.
   !
   !            In order to keep the size of this source file reasonably small,
   !            this module only contains routines that are directly related to
   !            the integration of BT*Sigma.
   !
   !     $Revision: 9707 $
   !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
   !
   !**********************************************************************

   use ModCounters
   use ModReadCalculationData
   use ModElementEvaluationTETRA
   use ModElementEvaluationTRI
   use ModElementEvaluationQUAD
   use ModMPMData
   use ModMPMInit
   use ModMeshInfo
   use ModTwoLayerFormulation
   use ModGlobalConstants

   implicit none

contains
   !**********************************************************************
   !
   ! Subroutine: MPMDYNBTSig
   ! Description:  Calculation of the equivalent nodal forces due to
   !               a given stress InternalLD = Integral {BT*Sigma}.
   !
   ! Explaination: Calculates the nodal internal and external forces.
   !               Called during the Lagrangian Phase in the "GetNodalExtAndIntForces()" subroutine.
   !               Intakes External, Internal, and Gravity Forces
   !
   !     @param[inout] ExternalLD : Nodal load array of external load
   !     @param[inout] InternalLD : Nodal load array of internal load
   !     @param[inout] GravityLD :  Nodal load array of gravity load
   !     @param[inout] Reactions : Nodal reactions on a user-defined mesh boundary
   !     @param[inout] ReactionsWater : ??Nodal force from the water??
   !     @param[inout] BulkViscousLD : ??Viscous load term??
   !**********************************************************************
   subroutine MPMDYNBTSig(ExternalLD, InternalLD, GravityLD, Reactions, ReactionsWater, BulkViscousLD)

      implicit none

      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ExternalLD
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: InternalLD
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: GravityLD
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: Reactions
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ReactionsWater
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: BulkViscousLD

      ! Local variables
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
      real(REAL_TYPE), dimension(NTENSOR) :: S
      real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity
      real(REAL_TYPE) :: WtN, Det, WPP, TotIntWeightLiqMP
      real(REAL_TYPE) :: VPressure = 0.0
      real(REAL_TYPE) :: Position
      real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
      integer(INTEGER_TYPE) :: I, J, IntGlo, IEl, Int, NN, IPart, IAEl, NElemPart, iEntityID, iNode, iDof, ParticleIndex, NodeID, GlobDof, iEntity, MaterialIndex,ILoadSystem
      logical :: DoConsiderReactionForces, IsUndrEffectiveStress, IsPrescribedVelocity, IsLoadOnMP

      
      integer(INTEGER_TYPE):: MaxIPart
      
      !reset the needed variables
      IntGlo = 0
      S = 0.0
      ExternalLD = 0.0
      InternalLD = 0.0
      GravityLD = 0.0
      Reactions = 0.0
      ReactionsWater = 0.0
      BulkViscousLD = 0.0
      IsLoadOnMP =(Counters%NLoadedElementSidesSolidMatPoints+Counters%NLoadedElementSidesSolidMatPointsB)>0


      do IAEl = 1, Counters%NAEl ! loop over all active elements
         IEl = ActiveElement(IAEl)
         !========================================================================
         !========================================================================
         if(NFORMULATION==1) then !1 Constituents - Flag for single point model

            !//////////////////// External force calculation ////////////////////////
            do IPart = 1, NPartEle(IEl) ! loop over all material points in element
               ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the material point ID
               PartLoad = 0.0
               PartGravity = 0.0
               if (CalParams%ApplyContactAlgorithm) then
                  iEntity = EntityIDArray(ParticleIndex) ! entity to which material point belongs
               else
                  iEntity = 1
               endif

               do INode = 1, ELEMENTNODES
                  NodeID = ElementConnectivities(INode, IEl)  ! nodal global ID
                  GlobDof = ReducedDof(NodeID)

                  if (IsLoadOnMP) then
                     do ILoadSystem =1, Counters%NSoilLoadSystems
                        !  ILoadSystem=1
                        PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExt(:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
                        do I = 1, NVECTOR
                           ExternalLD(GlobDof + I, iEntity) = ExternalLD(GlobDof + I, iEntity) + PartLoad(I)
                        end do
                     end do
                  end if

                  if (CalParams%ApplySubmergedCalculation) then
                     PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyMixed
                  else
                     PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody
                  end if

                  do I = 1, NVECTOR
                     GravityLD(GlobDof + I, iEntity) = GravityLD(GlobDof + I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                  end do

               end do
            end do

            !//////////////////// Internal force calculation ////////////////////////
            ! Determine number of integration points inside element
            if (IsParticleIntegration(IEl) ) then ! True - material point based integration, false - Gauss point based integration
               NElemPart = NPartEle(IEl)  ! Number of material points in element
            else
               !NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
                
                if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                    ! map the stress using the 4 gauss points 
                    NElemPart = Counters%NGaussPoints !ELEMENTGAUSSPOINTS !this always have to be 4 gauss points
                     else 
                    ! this is left for the triangular elements as when it is 1 for the mixed scheme in general     
                    NElemPart = ELEMENTGAUSSPOINTS
                     end if 
            end if

            

            !------------------------------------ INTEGRATION POINT LOOP --------------------------------
            do Int = 1, NElemPart ! Loop over number of integration points if MIXED Intgration or material points if MP Integration per element IEl

                ! we only need to loop over the gauss points once --> not over the material points
                if (IsParticleIntegration(IEl)) then      
                    ! do nothing
                    MaxIPart = 1!--> effectively no loop as we go through each material point
                    IntGlo = GetParticleIndex(Int, IEl)
                    MaterialIndex = MaterialIDArray(IntGlo)
                else 
                    
                    MaxIPart = SubElementMPOrganization(IEl,Int) ! loop over the material points in the subzone
                    ! Determine global ID of integration point
                    IntGlo = GetParticleIndexInSubElement(IEl, Int, 1) ! pick the first index to find the MaterialIndex
                    MaterialIndex = MaterialIDArray(IntGlo) ! we probably need to be able choose different materials 
                end if
               
               ! if (IsParticleIntegration(IEl)) then
               !     ! recalculating the B matrix for every point in the integration loop
               !     call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN, DShapeValuesArray(IntGlo,:,:)) ! get the B-matrix once per element
               ! else 
               !     ! no need to input the shape function derivatives as we will use the derivatives at the gauss points
               !     call FormB3_GP(Int, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN)
               !end if
                
                if  (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                call FormB3_GP(Int, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get B-matrix
                
                     else 
                         
                         !if (IsParticleIntegration(IEl)) then             
                ! recalculating the B matrix for every point in the integration loop
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN, DShapeValuesArray(IntGlo,:,:)) ! get the B-matrix once per element
               
                end if 

               IsUndrEffectiveStress = &
               !code version 2016 and previous
                  ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
               !code version 2017.1 and following
                  (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))

               ! Set the integration weight
               if (IsParticleIntegration(IEl) ) then
                  if ( ISAXISYMMETRIC ) then
                     ! the integration weight of the MP is not corrected due to the axisymmetry
                     Position = GlobPosArray(IntGlo, 1) ! INDEX_X = 1
                     ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
                  ! use the integration weight of material point if it is partially filled,
                  ! otherwise the weight of gauss point is used
                  WTN = Particles(IntGlo)%IntegrationWeight
               else
                  if ( ISAXISYMMETRIC ) then
                     Position = GPGlobalPositionElement(1, Int, IEl) ! index 1 is radial direction
                     ShapeValues(:) = GPShapeFunction(Int, :)
                     WTN = WTN * Position ! the volume of the element is corrected due to the axisymmetry
                  end if
               end if

               ! Determine stress vector for integration point (material point or Gauss point)
               ! = effective stresses + pore pressure
               if (((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) &
                  .or.IsUndrEffectiveStress) then
                   ! --> we need to use SigmaEffArrayGaussPoints
                  WPP = Particles(IntGlo)%WaterPressure * WTN
               end if
               do J = 1, 3 ! first three stress components both for 2D and 3D
                  if (IsUndrEffectiveStress) then
                      
                     S(J) = SigmaEffArray(IntGlo, J) * WTN + WPP
                  else
                      
                      if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                          S(J) =  SigmaEffArrayGaussPoints(IEl, Int, J) * WTN
                      !if (IsParticleIntegration(IEl) ) then
                       else 
                          S(J) = SigmaEffArray(IntGlo, J) * WTN
                      end if 
                  end if
               end do
               do J = 4, NTENSOR
                   

                   if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                    
                      S(J) =  SigmaEffArrayGaussPoints(IEl, Int, J) * WTN  
                     else      
                         S(J) = SigmaEffArray(IntGlo,J)  * WTN
                   
                     end if 
               end do
               if (CalParams%ApplyBulkViscosityDamping) then
                   if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                    
                       VPressure = SigmaEffArrayGaussPointsBulkViscousPressure(IEl, Int)
                   
                     else 
                         
                  VPressure = Particles(IntGlo)%DBulkViscousPressure * WTN
               
                     end if
               end if
               
               !get material point entity
               if (CalParams%ApplyContactAlgorithm) then ! contact algorithm
                  iEntityID = EntityIDArray(IntGlo) !particles(IntGlo)%EntityID ! get the entity ID for the current particle
               else
                  iEntityID = 1
               end if

               do iNode = 1, ELEMENTNODES ! loop over element nodes

                  nn=ElementConnectivities(iNode,iel) ! get global node number
                  IDof = ReducedDof(nn) ! get dof number

                  if ( NVECTOR == 2 ) then ! 2D
                     ! nodal x-load                                              --> replacing last terms with BTSigma(I) gives error in tests
                     InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4)
                     ! nodal y-load
                     InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2)
                     ! nodal r-load (x-load)
                     if ( ISAXISYMMETRIC ) then
                        ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                        InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + S(3) * ShapeValues(iNode) / Position
                     end if
                  else if ( NVECTOR == 3 ) then ! 3D
                     ! nodal x-load                                              --> replacing last terms with BTSigma(I) gives error in tests
                     InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4) + B(3,iNode)*S(6)
                     ! nodal y-load
                     InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2) + B(3,iNode)*S(5)
                     ! nodal z-load
                     InternalLD(IDof+3, iEntityID) = InternalLD(IDof+3, iEntityID) + B(1,iNode)*S(6) + B(2,iNode)*S(5) + B(3,iNode)*S(3)
                  else
                     call GiveError('Dimension is not correct. It must be 2 or 3.')
                  end if

                  if (CalParams%ApplyBulkViscosityDamping) then
                     do I = 1, NVECTOR
                        BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + B(I, INode) * VPressure
                        if ( ISAXISYMMETRIC ) then
                           ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                           BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + VPressure * ShapeValues(iNode) / Position
                        end if
                     end do
                  end if

                  IsPrescribedVelocity = CalParams%PrescribedVelo%ApplyPrescribedVelo

                  DoConsiderReactionForces =  IsPrescribedVelocity

                  if (DoConsiderReactionForces) then
                     if ( NVECTOR == 3 ) then ! 3D                          ! --> replacing last terms with BTSigma(I) gives error in tests
                        Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(3, INode) * S(6) + B(1, INode) * VPressure
                        Reactions(IDof+2, IEntityID) = Reactions(IDof+2, IEntityID) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(3, INode) * S(5) + B(2, INode) * VPressure
                        Reactions(IDof+3, IEntityID) = Reactions(IDof+3, IEntityID) + B(1, INode) * S(6) + B(2, INode) * S(5) + B(3, INode) * S(3) + B(3, INode) * VPressure
                     else if ( NVECTOR == 2 ) then ! 2D                      ! --> replacing last terms with BTSigma(I) gives error in tests
                        Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(1, INode) * VPressure
                        Reactions(IDof+2, IEntityID) = Reactions(IDof+2, IEntityID) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(2, INode) * VPressure
                        if ( ISAXISYMMETRIC ) then
                           ! Equation below is to be double checked
                           Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + S(3) * ShapeValues(iNode) / Position
                        end if
                     else
                        call GiveError('Dimension is not correct. It must be 2 or 3.')
                     end if

                     if ( ( (CalParams%NumberOfPhases==2) .or. (CalParams%NumberOfPhases==3) ) .or. IsUndrEffectiveStress ) then
                        do I = 1, NVECTOR
                           ReactionsWater(IDof+I, IEntityID) = ReactionsWater(IDof+I, IEntityID) + B(I, INode) * WPP
                        end do
                     end if

                  end if

               end do ! node loop
                end do ! Loop over 1 .. NElemPart
            !end do 
            
            !------------------------------------ END INTEGRATION POINT LOOP -----------------------------
         end if ! 1 Constituent
         !========================================================================
         !========================================================================
         if(NFORMULATION==2) then !2 Constituents

            !//////////////////// External force calculation ////////////////////////
            do IPart = 1, NPartEle(IEl) ! loop over all material points in element
               ParticleIndex = GetParticleIndex(IPart, IEl) ! Get the material point ID

               if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then
                  ! Only SOLID material Point or TwoLayerApplyMixtureLiquidApproach
                  PartLoad = 0.0
                  PartGravity = 0.0
                  iEntity = 1

                  do INode = 1, ELEMENTNODES
                     NodeID = ElementConnectivities(INode, IEl)  ! nodal global ID
                     IDof = ReducedDof(NodeID)

                     if (IsLoadOnMP) then
                        do ILoadSystem=1, Counters%NSoilLoadSystems
                           ! ILoadSystem=1
                           PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExt(:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
                           do I = 1, NVECTOR
                              ExternalLD(IDof+I, iEntity) = ExternalLD(IDof+I, iEntity) + PartLoad(I)
                           end do
                        end do
                     end if

                     PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody
                     do I = 1, NVECTOR
                        GravityLD(IDof+I, iEntity) = GravityLD(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
                     end do

                  end do

               end if ! Only SOLID material Point
            end do

            !//////////////////// Internal force calculation ////////////////////////
            ! Determine number of integration points inside element
            if (IsParticleIntegration(IEl) ) then ! True - material point based integration, false - Gauss point based integration
               NElemPart = NPartEle(IEl)  ! Number of material points in element
            else
               !NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
                if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                    ! map the stress using the 4 gauss points 
                    NElemPart = Counters%NGaussPoints !ELEMENTGAUSSPOINTS !this always have to be 4 gauss points
                     else 
                    ! this is left for the triangular elements as when it is 1 for the mixed scheme in general     
                    NElemPart = ELEMENTGAUSSPOINTS
                     end if 
            end if

            TotIntWeightLiqMP = 0.0
            do J = 1, NPartEle(IEl)
               ParticleIndex = GetParticleIndex(J, IEl)
               if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                  TotIntWeightLiqMP = TotIntWeightLiqMP + Particles(ParticleIndex)%IntegrationWeight
               end if
            end do


            !------------------------------------ INTEGRATION POINT LOOP --------------------------------
            do Int = 1, NElemPart ! Loop over number of integration points per element IEl

                
                !if (IsParticleIntegration(IEl)) then 
                !    ! do nothing
                !    MaxIPart = 1!--> effectively no loop 
                !else 
                !    MaxIPart = SubElementMPOrganization(IEl,Int)
                !end if 
                
                !do IPart = 1, MaxIPart!SubElementMPOrganization(IEl,Int)
                    
                ! Determine global ID of integration point
            if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                IntGlo = GetParticleIndexInSubElement(IEl, Int, 1)     
                
                     else 
                         
                IntGlo = GetParticleIndex(Int, IEl)
                     end if 
                
               ! Determine global ID of integration point
               !IntGlo = GetParticleIndex(Int, IEl)
               

               ! recalculating the B matrix for every point in the integration loop
               !call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN, DShapeValuesArray(IntGlo,:,:) ) ! get the B-matrix once per element

                     if  (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                call FormB3_GP(Int, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get B-matrix
                
                     else 
                         
                         !if (IsParticleIntegration(IEl)) then             
                ! recalculating the B matrix for every point in the integration loop
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN, DShapeValuesArray(IntGlo,:,:)) ! get the B-matrix once per element
               
                end if 

               if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeSolid) then
                  ! Only SOLID material Point or TwoLayerApplyMixtureLiquidApproach

                  ! Set the integration weight
                  if (IsParticleIntegration(IEl) ) then
                     ! use the integration weight of material point if it is partially filled,
                     ! otherwise the weight of gauss point is used
                     WTN = Particles(IntGlo)%IntegrationWeight
                  end if

                  ! Determine stress vector for integration point (material point or Gauss point)
                  
                  if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 

                        
                      do J = 1, NTENSOR
                     S(J) = SigmaEffArrayGaussPoints(IEl,Int,J) * WTN !SigmaEffArray(IntGlo,J)  * WTN
                      end do
                      
                     else 
                         
                  do J = 1, NTENSOR
                     S(J) = SigmaEffArray(IntGlo,J)  * WTN
                  end do
                  
                     end if 
                     
                  if (CalParams%ApplyBulkViscosityDamping) then
                      
                      if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 

                          VPressure =  SigmaEffArrayGaussPointsBulkViscousPressure(IEl, Int) * WTN 
                          
                     else 
                         
                     VPressure = Particles(IntGlo)%DBulkViscousPressure * WTN                     
                     end if
                     
                  end if

                  !get material point entity
                  iEntityID = 1

                  do iNode=1,ELEMENTNODES ! loop over element nodes
                     nn=ElementConnectivities(iNode,iel) ! get global node number
                     IDof = ReducedDof(nn)

                     if ( NVECTOR == 2 ) then ! 2D
                        ! nodal x-load                                           ! --> replacing last terms with BTSigma(I) gives error in tests
                        InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4)
                        ! nodal y-load
                        InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2)
                     else if ( NVECTOR == 3 ) then ! 3D
                        ! nodal x-load                                           ! --> replacing last terms with BTSigma(I) gives error in tests
                        InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4) + B(3,iNode)*S(6)
                        ! nodal y-load
                        InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2) + B(3,iNode)*S(5)
                        ! nodal z-load
                        InternalLD(IDof+3, iEntityID) = InternalLD(IDof+3, iEntityID) + B(1,iNode)*S(6) + B(2,iNode)*S(5) + B(3,iNode)*S(3)
                     else
                        call GiveError('Dimension is not correct. It must be 2 or 3.')
                     end if

                     if (CalParams%ApplyBulkViscosityDamping) then
                        do I = 1, NVECTOR
                           BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + B(I, INode) * VPressure
                        end do
                     end if

                  end do ! node loop
               end if ! Only SOLID material Point

               !--------------------------------------------------------
               if (.not.(CalParams%NumberOfPhases==1)) then !NumberOfPhases > 1
                  if (MaterialPointTypeArray(IntGlo)==MaterialPointTypeLiquid) then
                     ! Only LIQUID material Point if Element = SOLID+LIQUID and NOT TwoLayerApplyMixtureLiquidApproach

                     ! Set the integration weight
                     WTN = Particles(IntGlo)%IntegrationWeight
                     if (TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
                        WTN = ElementSpace(IEl) / TotIntWeightLiqMP * Particles(IntGlo)%IntegrationWeight
                        if (TwoLayerData%Elements(IEl)%IsBoundaryOfLiquidElements) then ! Boundary of Liquid domain
                           if ((TotIntWeightLiqMP/ElementSpace(IEl))<CalParams%RequiredDegreeOfFilling) then
                              WTN = Particles(IntGlo)%IntegrationWeight
                           end if
                        end if
                     end if

                     ! Determine stress vector for integration point (material point or Gauss point)
                     S(:) = SigmaEffArray(IntGlo, :)  * WTN * TwoLayerData%ELEMENTS(IEl)%ConcentrationRatioSolidL

                     !get material point entity
                     iEntityID = 1

                     do iNode = 1, ELEMENTNODES ! loop over element nodes
                        nn=ElementConnectivities(iNode,iel) ! get global node number
                        IDof = ReducedDof(nn)

                        if ( NVECTOR == 2 ) then ! 2D
                           ! nodal x-load                                             ! --> replacing last terms with BTSigma(I) gives error in tests
                           InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4)
                           ! nodal y-load
                           InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2)
                        else if ( NVECTOR == 3 ) then ! 3D
                           ! nodal x-load                                             ! --> replacing last terms with BTSigma(I) gives error in tests
                           InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4) + B(3,iNode)*S(6)
                           ! nodal y-load
                           InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2) + B(3,iNode)*S(5)
                           ! nodal z-load
                           InternalLD(IDof+3, iEntityID) = InternalLD(IDof+3, iEntityID) + B(1,iNode)*S(6) + B(2,iNode)*S(5) + B(3,iNode)*S(3)
                        end if

                     end do ! node loop

                  end if ! Only LIQUID material Point if Element = SOLID+LIQUID
               end if !NumberOfPhases > 1
                !end do ! Loop over 1 .. NElemPart
                
            end do
            
            !------------------------------------ END INTEGRATION POINT LOOP -----------------------------

         end if !2 Constituents
      end do ! Loop over elements

   end subroutine MPMDYNBTSig


   subroutine MPMDYNLoad(ExternalLD, GravityLD)
      !**********************************************************************
      !
      !    Function:  calculate nodal vectors of geostatic loading and external surface loading
      !
      !*********************************************************************

      implicit none

      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: ExternalLD
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: GravityLD

      ! Local variables
      integer(INTEGER_TYPE) :: I, IEl, IPart, IAEl, iNode, ParticleIndex, NodeID, iEntity, IDof, ILoadSystem
      real(REAL_TYPE), dimension(NVECTOR) :: PartLoad, PartGravity
      logical :: IsLoadOnMP

      ExternalLD = 0.0
      GravityLD = 0.0
      IsLoadOnMP =(Counters%NLoadedElementSidesSolidMatPoints+Counters%NLoadedElementSidesSolidMatPointsB)>0

      do IAEl = 1, Counters%NAEl
         IEl = ActiveElement(IAEl)

         do IPart = 1, NPartEle(IEl)
            ParticleIndex = GetParticleIndex(IPart, IEl)
            PartLoad = 0.0
            PartGravity = 0.0
            if (CalParams%ApplyContactAlgorithm) then
               iEntity = EntityIDArray(ParticleIndex)
            else
               iEntity = 1
            endif

            do INode = 1, ELEMENTNODES
               NodeID = ElementConnectivities(INode, IEl)
               IDof = ReducedDof(NodeID)

               if (IsLoadOnMP) then
                  do IloadSystem=1, Counters%NSoilLoadSystems
                     PartLoad = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FExt(:,ILoadSystem) * CalParams%Multipliers%SolidACurrent(ILoadSystem)
                     do I = 1, NVECTOR
                        ExternalLD(IDof+I, iEntity) = ExternalLD(IDof+I, iEntity) + PartLoad(I)
                     end do
                  end do
               end if

               if (CalParams%ApplySubmergedCalculation) then
                  PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBodyMixed
               else
                  PartGravity = ShapeValuesArray(ParticleIndex,INode) * Particles(ParticleIndex)%FBody
               end if

               do I = 1, NVECTOR
                  GravityLD(IDof+I, iEntity) = GravityLD(IDof+I, iEntity) + PartGravity(I) * CalParams%Multipliers%GravityCurrent
               end do

            end do
         end do

      end do

   end subroutine MPMDYNLoad

   subroutine MPMDYNBTSigOnly(InternalLD, Reactions, BulkViscousLD)
      !**********************************************************************
      !
      ! Function:  Calculation of the equivalent nodal forces due to
      !            a given stress InternalLD = Integral {BT*Sigma}.
      !
      !**********************************************************************

      implicit none

      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: InternalLD
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: Reactions
      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(inout) :: BulkViscousLD

      ! Local variables
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
      real(REAL_TYPE), dimension(NTENSOR) :: S
      real(REAL_TYPE) :: Position
      real(REAL_TYPE) :: ShapeValues(ELEMENTNODES)
      integer(INTEGER_TYPE) :: I, IntGlo, IEl, Int, J, NN, IAEl, NElemPart, iEntityID, iNode, IDof, MaterialIndex
      real(REAL_TYPE) :: WtN, Det, WPP, VPressure = 0.0
      logical :: DoConsiderReactionForces, IsUndrEffectiveStress
      
          
      integer(INTEGER_TYPE):: MaxIPart, IPart


      IntGlo = 0
      S = 0.0
      InternalLD = 0.0
      Reactions = 0.0
      BulkViscousLD = 0.0

      do IAEl = 1, Counters%NAEl
         IEl = ActiveElement(IAEl)

         if(NFORMULATION==1) then

            if (IsParticleIntegration(IEl) ) then
               NElemPart = NPartEle(IEl)
            else
                if (ELEMENTTYPE == QUAD4) then 
                    
                    NElemPart = Counters%NGaussPoints
                else 
               
                    NElemPart = ELEMENTGAUSSPOINTS
                end if 
                
            end if


            ! if quadrilateral element, then loop over 4 gaussian points
            do Int = 1, NElemPart

                !if (IsParticleIntegration(IEl)) then 
                !    ! do nothing
                !    MaxIPart = 1!--> effectively no loop 
                !else 
                !    MaxIPart = SubElementMPOrganization(IEl,Int)
                !end if 
                
                !do IPart = 1, MaxIPart!SubElementMPOrganization(IEl,Int)
                    
                if (IsParticleIntegration(IEl)) then 
                    IntGlo = GetParticleIndex(Int, IEl)   ! Determine global ID of integration point 
                else 
                    IntGlo = GetParticleIndexInSubElement(IEl, Int, 1)
                end if 
                !if (.not.IsParticleIntegration(IEl) .and. ELEMENTTYPE == QUAD4) then 
                !IntGlo = GetParticleIndexInSubElement(IEl, Int)
                
                ! I commented these for now....
               !IntGlo = GetParticleIndex(Int, IEl)
               MaterialIndex = MaterialIDArray(IntGlo)
               
               ! get gauss point location 
               
               ! get gauss point shape functions and derivatives 
            
               ! recalculating the B matrix for every point in the integration loop
               if  (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                
                   call FormB3_GP(Int, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get B-matrix
                
                     else 
                         
                         !if (IsParticleIntegration(IEl)) then             
                ! recalculating the B matrix for every point in the integration loop
                call FormB3(1, IEl, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN, DShapeValuesArray(IntGlo,:,:)) ! get the B-matrix once per element
               
                end if 
               


               IsUndrEffectiveStress = &
               !code version 2016 and previous
                  ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
               !code version 2017.1 and following
                  (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))

               if (IsParticleIntegration(IEl) ) then
                  if ( ISAXISYMMETRIC ) then
                     ! the integration weight of the MP is not corrected due to the axisymmetry
                     Position = GlobPosArray(IntGlo, 1) ! INDEX_X = 1
                     ShapeValues(:) = ShapeValuesArray(IntGlo, :)
                  end if
                  WTN = Particles(IntGlo)%IntegrationWeight
               else
                  if ( ISAXISYMMETRIC ) then
                     Position = GPGlobalPositionElement(1, Int, IEl) !(GaussPointGlobalCoord(INDEX_X)
                     ShapeValues(:) = GPShapeFunction(Int, :)
                     WTN = WTN * Position ! the volume of the element is corrected due to the axisymmetry
                  end if
               end if

               if (IsUndrEffectiveStress) then
                  WPP = Particles(IntGlo)%WaterPressure * WTN
               end if
               do J = 1, 3 ! always first three stress components, both for 2D and 3D
                  if (IsUndrEffectiveStress.or.CalParams%ApplyImplicitQuasiStatic) then
                     WPP = Particles(IntGlo)%WaterPressure * WTN
                     S(J) = SigmaEffArray(IntGlo, J) * WTN + WPP
                  else
                      
                      if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                          S(J) =  SigmaEffArrayGaussPoints(IEl, Int, J) * WTN                    
                     else 
                         S(J) = SigmaEffArray(IntGlo, J) * WTN 
                     end if 
                     
                  end if
               end do
               do J = 4, NTENSOR
                   
                   if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                    
                       S(J) =  SigmaEffArrayGaussPoints(IEl, Int, J) * WTN
                       
                     else
                       
                   S(J) = SigmaEffArray(IntGlo,J)  * WTN
                     end if 
               end do

               if (CalParams%ApplyBulkViscosityDamping) then
                   if (ELEMENTTYPE == QUAD4 .and. &
                    (CalParams%ComputationMethod==MPM_MIXED_MG22_INTEGRATION .or. &
                     CalParams%ComputationMethod==MPM_MIXED_MG22_NOINTERPOLATION_INTEGRATION_SPECIFIER)) then 
                       VPressure = SigmaEffArrayGaussPointsBulkViscousPressure(IEl, Int) * WTN !Particles(IntGlo)%DBulkViscousPressure * WTN
                     else 
                         VPressure = Particles(IntGlo)%DBulkViscousPressure * WTN
                     end if
               end if

               if (CalParams%ApplyContactAlgorithm) then
                  iEntityID = EntityIDArray(IntGlo)
               else
                  iEntityID = 1
               end if

               DO iNode=1,ELEMENTNODES
                  nn=ElementConnectivities(iNode,iel)
                  IDof = ReducedDof(nn)

                  if ( NVECTOR == 3 ) then ! 3D
                     InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4) + B(3,iNode)*S(6)
                     InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2) + B(3,iNode)*S(5)
                     InternalLD(IDof+3, iEntityID) = InternalLD(IDof+3, iEntityID) + B(1,iNode)*S(6) + B(2,iNode)*S(5) + B(3,iNode)*S(3)
                  else if ( NVECTOR == 2 ) then ! 2D
                     InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + B(1,iNode)*S(1) + B(2,iNode)*S(4)
                     InternalLD(IDof+2, iEntityID) = InternalLD(IDof+2, iEntityID) + B(1,iNode)*S(4) + B(2,iNode)*S(2)
                     if ( ISAXISYMMETRIC ) then
                        ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                        InternalLD(IDof+1, iEntityID) = InternalLD(IDof+1, iEntityID) + S(3) * ShapeValues(iNode) / Position
                     end if
                  else
                     call GiveError('Dimension is not correct. It must be 2 or 3.')
                  end if

                  if (CalParams%ApplyBulkViscosityDamping) then
                     do I = 1, NVECTOR
                        BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + B(I, INode) * VPressure
                        if ( ISAXISYMMETRIC ) then
                           ! Eq. 27 from Sulsky & Schreyer (1996) : http://www.sciencedirect.com/science/article/pii/S0045782596010912
                           BulkViscousLD(IDof+I, IEntityID) = BulkViscousLD(IDof+I, IEntityID) + VPressure * ShapeValues(iNode) / Position
                        end if
                     end do
                  end if

                  if ( NVECTOR == 3 ) then ! 3D
                     DoConsiderReactionForces = CalParams%PrescribedVelo%ApplyPrescribedVelo
                  else if ( NVECTOR == 2 ) then ! 2D
                     DoConsiderReactionForces =  CalParams%PrescribedVelo%ApplyPrescribedVelo
                  end if

                  if (DoConsiderReactionForces) then
                     if ( NVECTOR == 3 ) then ! 3D
                        Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(3, INode) * S(6) + B(1, INode) * VPressure
                        Reactions(IDof+2, IEntityID) = Reactions(IDof+2, IEntityID) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(3, INode) * S(5) + B(2, INode) * VPressure
                        Reactions(IDof+3, IEntityID) = Reactions(IDof+3, IEntityID) + B(1, INode) * S(6) + B(2, INode) * S(5) + B(3, INode) * S(3) + B(3, INode) * VPressure
                     else if ( NVECTOR == 2 ) then ! 2D
                        Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + B(1, INode) * S(1) + B(2, INode) * S(4) + B(1, INode) * VPressure
                        Reactions(IDof+2, IEntityID) = Reactions(IDof+2, IEntityID) + B(1, INode) * S(4) + B(2, INode) * S(2) + B(2, INode) * VPressure
                        if ( ISAXISYMMETRIC ) then
                           ! Equation below is to be double checked
                           Reactions(IDof+1, IEntityID) = Reactions(IDof+1, IEntityID) + S(3) * ShapeValues(iNode) / Position
                        end if
                     else
                        call GiveError('Dimension is not correct. It must be 2 or 3.')
                     end if
                  end if

               end do
                end do ! loop over points
            !end do
            
         end if
      end do


   end subroutine MPMDYNBTSigOnly

end module ModMPMDYNBTSigma
