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
	
	
      module ModLiquid
      !**********************************************************************
      !
      ! Function: Contains routines for modelling of liquid
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMPMData
      use ModMeshInfo
      use ModTwoLayerFormulation
      use ModMPMStresses
      
      implicit none
      
      contains


        subroutine InitialiseLiquidData()
        !**********************************************************************
        !
        !    Function: initialise 1-phase liquid material points data
        !
        !*********************************************************************  
        
        implicit none
        
          ! local variables
          integer :: IError, NMaterialSets, IMatSet
          
          NMaterialSets = Counters%NLayers
          do IMatSet = 1, NMaterialSets 
              if (NFORMULATION == 1) then ! 1 Constituent
                if (MatParams(IMatSet)%MaterialType == '1-phase-liquid'.or.MatParams(IMatSet)%MaterialPhases == '1-phase-liquid') then
          
                  allocate(NodalDensity(Counters%NodTot), stat = IError)
                  NodalDensity = 0.0

                  allocate(ActiveNodeElementVolume(Counters%NodTot), stat = IError)
                  ActiveNodeElementVolume = 0.0
          
                  ActiveNode = .false.
                  
                  allocate(IsElemWithLiquidFreeSurfMP(Counters%NEl), stat = IError)
                  IsElemWithLiquidFreeSurfMP = .false. ! Initially liquid free surface material points are not detected since all elements are fully filled, and IsElemWithLiquidFreeSurfMP is false                  
              
                end if
              end if
              
              
              if (.not.(NFORMULATION == 1)) then ! 2 constituents
          
                  allocate(NodalDensity(Counters%NodTot), stat = IError)
                  NodalDensity = 0.0

                  allocate(ActiveNodeElementVolume(Counters%NodTot), stat = IError)
                  ActiveNodeElementVolume = 0.0
          
                  ActiveNode = .false.
                  
                  allocate(IsElemWithLiquidFreeSurfMP(Counters%NEl), stat = IError)
                  IsElemWithLiquidFreeSurfMP = .false. ! Initially liquid free surface material points are not detected since all elements are fully filled, and IsElemWithLiquidFreeSurfMP is false
              
              end if
          end do
          
          
        end subroutine InitialiseLiquidData
      

        subroutine DestroyLiquidData()
        !**********************************************************************
        !
        !    Function: destroy 1-phase liquid material points data
        !
        !********************************************************************* 
        
        implicit none
        
          ! Local variables
          integer :: IError, NMaterialSets, IMatSet

          
          if (NFORMULATION == 1) then ! 1 constituent
           NMaterialSets = Counters%NLayers
           do IMatSet = 1, NMaterialSets 
              if (MatParams(IMatSet)%MaterialType == '1-phase-liquid'.or.MatParams(IMatSet)%MaterialPhases == '1-phase-liquid') then
                                  
                  if (allocated(NodalDensity)) then
                    deallocate(NodalDensity, stat = IError)
                  end if  

                  if (allocated(ActiveNodeElementVolume)) then
                    deallocate(ActiveNodeElementVolume, stat = IError)
                  end if
                  
                  if (allocated(IsElemWithLiquidFreeSurfMP) ) then
                     deallocate(IsElemWithLiquidFreeSurfMP, stat = IError)
                  end if
                            
              end if
            end do
           end if
           
              
          if (.not.(NFORMULATION == 1)) then ! 2 constituents
          
                 if (allocated(NodalDensity)) then
                    deallocate(NodalDensity, stat = IError)
                  end if  

                  if (allocated(ActiveNodeElementVolume)) then
                    deallocate(ActiveNodeElementVolume, stat = IError)
                  end if
                  
                  if (allocated(IsElemWithLiquidFreeSurfMP) ) then
                    deallocate(IsElemWithLiquidFreeSurfMP, stat = IError)
                  end if
                        
          end if
              
          
        end subroutine DestroyLiquidData
        
        
        subroutine DetermineLiquidDensityField()
        !**********************************************************************
        !
        !    Function: determine liquid density field
        !
        !********************************************************************* 
        
        implicit none
        
        
          ! Local variables
          integer :: IAEl, IEl, I, J
          integer :: ParticleIndex
          integer :: NAdjacentElements, NMaterialSets, IMatSet
          integer, dimension(ELEMENTNODES) :: NodeIDs
          real(REAL_TYPE) :: LiquidMass
          
          NMaterialSets = Counters%NLayers
          do IMatSet = 1, NMaterialSets 
              if (MatParams(IMatSet)%MaterialType == '1-phase-liquid'.or.MatParams(IMatSet)%MaterialPhases == '1-phase-liquid') then ! only for liquid material
          
                  ! Determine nodal mass from material points of activated elements connected to respective nodes
                  NodalDensity = 0.0
                  ActiveNode = .false.
                  do IAEl = 1, Counters%NAEl 
                    IEl = ActiveElement(IAEl)
                    NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl)
                    ActiveNode(NodeIDs) = .true.
                    do I = 1, NPartEle(IEl)
                      ParticleIndex = GetParticleIndex(I, IEl)
                      if (MaterialIDArray(ParticleIndex) == IMatSet) then 
                        LiquidMass = MassArray(ParticleIndex)
                      else ! not the whole element is filled by liquid, but also solid body is partly in the element
                        ! Replace solid density by fluid density 
                        LiquidMass = MatParams(IMatSet)%FluidThresholdDensity/1000 * Particles(ParticleIndex)%IntegrationWeight
                      end if
                      NodalDensity(NodeIDs) = NodalDensity(NodeIDs) + LiquidMass * ShapeValuesArray(ParticleIndex,:)
                    end do
                  end do
           
                  ActiveNodeElementVolume = 0.0
                  do I = 1, Counters%NodTot
                    if (ActiveNode(I)) then
                      NAdjacentElements = GetNElmOfNode(I)
                      do J = 1, NAdjacentElements
                        IEl = GetElmIOfNode(I, J)
                        if (IsActiveElement(IEl)) then
                          ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) + ElementSpace(IEl)
                        end if
                      end do
                    end if
                  end do
                  do I = 1, Counters%NodTot
                    if (ActiveNode(I)) then
                      ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) / 4.0 ! 4-noded tetrahedral element
                    end if
                  end do
          
                  do I = 1, Counters%NodTot
                    if (ActiveNode(I)) then
                      NodalDensity(I) = NodalDensity(I) / ActiveNodeElementVolume(I)
                    end if
                  end do
                  
              end if
          end do
          
        end subroutine DetermineLiquidDensityField              
      

        subroutine GetStressesLiquidMaterialPoint(ParticleID, IElement)
        !**********************************************************************
        !
        !    Function: get the stresses of 1-phase liquid material points
        !
        !********************************************************************* 
        
        implicit none
        
          integer, intent(in) :: ParticleID
          integer, intent(in) :: IElement
          ! Local variables
          integer :: MaterialIndex
          real(REAL_TYPE), dimension(NTENSOR) :: Sigma
          real(REAL_TYPE), dimension(NTENSOR) :: PrincipalStresses
          real(REAL_TYPE), dimension(NTENSOR) :: StrainRate, Strain, DevStrainRate
          real(REAL_TYPE) :: VolStrainRateComponent
          real(REAL_TYPE) :: VolStrainRate
          real(REAL_TYPE) :: VolStrainIncrement
          real(REAL_TYPE) :: LiquidPressure
          real(REAL_TYPE) :: LiquidPressureIncrement
          real(REAL_TYPE) :: EffectiveViscosity, Viscosity, YieldStress
          real(REAL_TYPE) :: RatioDensity
          logical :: PressureGreaterThreshold, PressureLowerThreshold, FullyFilled, FreeSurface
          logical :: DensityLowerThreshold, ComputeUpdatedStress
          
          
          if (NFORMULATION == 1) then ! 1 constituent

              MaterialIndex = MaterialIDArray(ParticleID)

              Strain = 0.0
              StrainRate = 0.0
              PressureGreaterThreshold = .false.
              DensityLowerThreshold = .false.

              Strain = GetEpsStep(Particles(ParticleID)) ! strain increment in material point
              VolStrainIncrement = Strain(1) + Strain(2) + Strain(3)
              StrainRate = Strain / CalParams%TimeIncrement

              VolStrainRate = StrainRate(1) + StrainRate(2) + StrainRate(3)
              VolStrainRateComponent = VolStrainRate / 3.0
              DevStrainRate(1:3) = StrainRate(1:3) - VolStrainRateComponent
              DevStrainRate(4:NTENSOR) = StrainRate(4:NTENSOR)


              ! Consider stress state of previous time step
              LiquidPressure = (SigmaEffArray(ParticleID,1)  &
                  + SigmaEffArray(ParticleID,2) + SigmaEffArray(ParticleID,3))/3

              if ((CalParams%ApplyBinghamFluid).or.(CalParams%ApplyFrictionalFluid) & !v2016
                  .or. (MatParams(MaterialIndex)%MaterialModel == ESM_FRICTIONAL_LIQUID) & !v2017.1 and following
                  .or. (MatParams(MaterialIndex)%MaterialModel == ESM_BINGHAM_LIQUID))  then ! material is a Bingham fluid

                  call GetBinghamYieldStress (ParticleID, YieldStress)
                  call GetApparentViscosity (StrainRate, YieldStress, MatParams(MaterialIndex)%ViscosityLiquid, Viscosity) !an apparent viscosity is necessary with Bingham fluid
              else
                  Viscosity = MatParams(MaterialIndex)%ViscosityLiquid
              end if

              !Determine interpolated density
              call CalculateRatioDensity(IElement,ParticleID,RatioDensity)

              ComputeUpdatedStress = .false. !initialize

              PressureGreaterThreshold = (LiquidPressure >= CalParams%LiquidPressureCavitationThreshold)
              PressureLowerThreshold = (LiquidPressure < CalParams%LiquidPressureCavitationThreshold)
              FreeSurface = (Particles(ParticleID)%LiquidFreeSurface == 1.0)
              FullyFilled = (RatioDensity > 1.0)

              ComputeUpdatedStress = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
                  ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))


              Sigma(1:3) = Particles(ParticleID)%WaterPressure
              Sigma(4:NTENSOR) = 0.0

              if(ComputeUpdatedStress) then
                  Sigma(1:3) = Particles(ParticleID)%WaterPressure +  &
                      2.0 * Viscosity * DevStrainRate(1:3)
                  Sigma(4:NTENSOR) = Viscosity * DevStrainRate(4:NTENSOR)
              end if ! ComputeUpdatedStress

              SigmaEffArray(ParticleID,:) = Sigma

          end if ! 1 constituent
          
          ! -----------------------------------------------------------------------------
          
          if (NFORMULATION == 2) then ! 2 constituent

              if(.not.(CalParams%NumberOfPhases == 1)) then ! 2 Phases


                  MaterialIndex = MaterialIDArray(ParticleID)

                  Strain = GetEpsStep(Particles(ParticleID))
                  VolStrainIncrement = Strain(1) + Strain(2) + Strain(3)
                  StrainRate = GetEpsStep(Particles(ParticleID)) / CalParams%TimeIncrement
                  VolStrainRate = StrainRate(1) + StrainRate(2) + StrainRate(3)
                  VolStrainRateComponent = VolStrainRate / 3.0
                  DevStrainRate(1:3) = StrainRate(1:3) - VolStrainRateComponent
                  DevStrainRate(4:NTENSOR) = StrainRate(4:NTENSOR)

                  !! Consider stress state of previous time step
                  LiquidPressure = (SigmaEffArray(ParticleID,1)  &
                      + SigmaEffArray(ParticleID,2) + SigmaEffArray(ParticleID,3))/3

                  ComputeUpdatedStress = .false. !initialize
                  PressureGreaterThreshold = (LiquidPressure >= CalParams%LiquidPressureCavitationThreshold)
                  PressureLowerThreshold = (LiquidPressure < CalParams%LiquidPressureCavitationThreshold)
                  FreeSurface = (Particles(ParticleID)%LiquidFreeSurface == 1.0)
                  if(TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL > 0.0) then
                      FullyFilled = (Particles(ParticleID)%FillingRatio > 0.98)
                  else
                      FullyFilled = (Particles(ParticleID)%FillingRatio > 1.0)
                  end if

                  ComputeUpdatedStress = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
                      ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))

                  Sigma(1:3) = Particles(ParticleID)%WaterPressure
                  Sigma(4:NTENSOR) = 0.0
                  if(ComputeUpdatedStress) then
                      if (Particles(ParticleID)%PhaseStatus == PhaseStatusLIQUID) then
                          if(TwoLayerData%Elements(IElement)%ContainedMaterialTypes == ContainedMaterialTypeSOLIDLIQUID) then
                              !Update Liquid Viscosity as a function of the solid volume fraction
                              EffectiveViscosity = MatParams(MaterialIndex)%ViscosityLiquid * ( 1.0 +  &
                                  5.0 / 2.0 * TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL + &
                                  5.2 * TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL *  &
                                  TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL )
                          else
                              EffectiveViscosity = MatParams(MaterialIndex)%ViscosityLiquid
                          end if
                          Sigma(1:3) = Particles(ParticleID)%WaterPressure +  &
                              2.0*EffectiveViscosity * DevStrainRate(1:3)
                          Sigma(4:NTENSOR) = MatParams(MaterialIndex)%ViscosityLiquid * DevStrainRate(4:NTENSOR)
                          !----------------------------------
                          !Correction for Free Water : Always in compression
                          if(CalParams%TwoLayerApplyNoTensStressLiqMPwLiqStatus.and.(Particles(ParticleID)%WaterPressure == 0.0)) then
                              Sigma = 0.0
                          end if
                          !----------------------------------
                      else if(Particles(ParticleID)%PhaseStatus == PhaseStatusSOLID) then
                          Sigma(1:3) = Particles(ParticleID)%WaterPressure
                          Sigma(4:NTENSOR) = 0.0
                      end if
                      !---------------------------------------------------------------------------------
                  end if ! ComputeUpdStress



                  SigmaEffArray(ParticleID,:) = Sigma



              end if ! 2 Phases
              !----------------------------------------------------------

              if((CalParams%NumberOfPhases == 1)) then ! 1 Phases

                  MaterialIndex = MaterialIDArray(ParticleID)


                  Strain = GetEpsStep(Particles(ParticleID))
                  VolStrainIncrement = Strain(1) + Strain(2) + Strain(3)
                  StrainRate = GetEpsStep(Particles(ParticleID)) / CalParams%TimeIncrement


                  VolStrainRate = StrainRate(1) + StrainRate(2) + StrainRate(3)
                  VolStrainRateComponent = VolStrainRate / 3.0


                  ! Consider stress state of previous time step
                  LiquidPressure = (SigmaEffArray(ParticleID,1)  &
                      + SigmaEffArray(ParticleID,2) + SigmaEffArray(ParticleID,3))/3

                  LiquidPressureIncrement = MatParams(MaterialIndex)%BulkModulusLiquid * VolStrainIncrement

                  LiquidPressure = LiquidPressure + LiquidPressureIncrement
                  PressureGreaterThreshold = (LiquidPressure > CalParams%LiquidPressureCavitationThreshold)

                  if (PressureGreaterThreshold) then
                      Sigma(1:3) = CalParams%LiquidPressureCavitationThreshold
                      Sigma(4:NTENSOR) = 0.0
                  else
                      Sigma(1:3) = LiquidPressure + 2.0*MatParams(MaterialIndex)%ViscosityLiquid *  &
                          (StrainRate(1:3) - VolStrainRateComponent)
                      Sigma(4:NTENSOR) = MatParams(MaterialIndex)%ViscosityLiquid * StrainRate(4:NTENSOR)
                  end if


                  SigmaEffArray(ParticleID,:) = Sigma
                  call CalculatePrincipalStresses(ParticleID, Sigma, PrincipalStresses)
                  call SetSigmaPrin(Particles(ParticleID), PrincipalStresses)

                  Particles(ParticleID)%WaterPressure = (Sigma(1) + Sigma(2) + Sigma(3)) / 3.0


              end if ! 1 Phase


          end if ! 2 constituents

                                          
        end subroutine GetStressesLiquidMaterialPoint
        

      subroutine GetPressureLiquidMaterialPoint(ParticleID, IElement)
      !**********************************************************************
      !
      !    Function: get pressure at the 1-phase liquid material points
      !
      !********************************************************************* 
        
        implicit none
        
          integer, intent(in) :: ParticleID
          integer, intent(in) :: IElement
          ! Local variables
          integer :: MaterialIndex
          real(REAL_TYPE), dimension(NTENSOR) :: StrainRate, Strain
          real(REAL_TYPE) :: VolStrainRate
          real(REAL_TYPE) :: VolStrainIncrement
          real(REAL_TYPE) :: LiquidPressure, ViscousDampingPressure
          real(REAL_TYPE) :: LiquidPressureIncrement
          real(REAL_TYPE) :: RatioDensity
          real(REAL_TYPE) :: DilationalWaveSpeed, Density, BulkModulusLiquid
          real(REAL_TYPE) :: ElementLMinLocal
          logical ::  DoComputeUpdatedPressure
            
          if (NFORMULATION == 1) then ! 1 constituent

              MaterialIndex = MaterialIDArray(ParticleID)

              Strain = 0.0
              StrainRate = 0.0

              Strain = GetEpsStep(Particles(ParticleID)) ! strain increment in material point
              VolStrainIncrement = Strain(1) + Strain(2) + Strain(3)
              StrainRate = Strain/CalParams%TimeIncrement

              ! Consider stress state of previous time step
              LiquidPressure = (SigmaEffArray(ParticleID,1)  &
                  + SigmaEffArray(ParticleID,2) + SigmaEffArray(ParticleID,3))/3

              !Determine interpolated density
              call CalculateRatioDensity(IElement,ParticleID,RatioDensity)

              call CheckComputeUpdatedPressure(ParticleID,IElement,LiquidPressure,DoComputeUpdatedPressure)

              if(DoComputeUpdatedPressure) then

                  LiquidPressureIncrement = MatParams(MaterialIndex)%BulkModulusLiquid  * VolStrainIncrement

                  LiquidPressure = LiquidPressure + LiquidPressureIncrement

                  Particles(ParticleID)%WaterPressure = LiquidPressure

                  if(Particles(ParticleID)%WaterPressure > CalParams%LiquidPressureCavitationThreshold) then
                      Particles(ParticleID)%WaterPressure = CalParams%LiquidPressureCavitationThreshold
                  end if

                  !calculate bulk viscosity
                  if (CalParams%ApplyBulkViscosityDamping) then

                      VolStrainRate = StrainRate(1)+StrainRate(2)+StrainRate(3)
                      Density = MatParams(MaterialIndex)%DensityLiquid/1000
                      BulkModulusLiquid = MatParams(MaterialIndex)%BulkModulusLiquid
                      DilationalWaveSpeed = sqrt(BulkModulusLiquid/Density)
                      ElementLMinLocal = ElementLMin(IElement)

                      ViscousDampingPressure = CalParams%BulkViscosityDamping1 *  &
                          Density * DilationalWaveSpeed * ElementLMinLocal * VolStrainRate

                      if ((VolStrainRate < 0.0).and.(CalParams%BulkViscosityDamping2 > 0.0)) then
                          ViscousDampingPressure = ViscousDampingPressure + &
                              Density * (CalParams%BulkViscosityDamping2 * ElementLMinLocal * VolStrainRate)**2
                      end if

                      Particles(ParticleID)%DBulkViscousPressure = ViscousDampingPressure
                  end if
              end if

          end if ! 1 constituent
          
          !! -----------------------------------------------------------------------------
          
          if (NFORMULATION == 2) then ! 2 constituents

              MaterialIndex = MaterialIDArray(ParticleID)

              Strain = 0.0

              Strain = GetEpsStep(Particles(ParticleID)) ! strain increment in material point
              VolStrainIncrement = Strain(1) + Strain(2) + Strain(3)

              ! Consider stress state of previous time step
              LiquidPressure = (SigmaEffArray(ParticleID,1)  &
                  + SigmaEffArray(ParticleID,2) + SigmaEffArray(ParticleID,3))/3

              call CheckComputeUpdatedPressure(ParticleID,IElement,LiquidPressure,DoComputeUpdatedPressure)

              if(DoComputeUpdatedPressure) then

                  LiquidPressureIncrement = MatParams(MaterialIndex)%BulkModulusLiquid  * Particles(ParticleID)%WaterVolumetricStrain

                  LiquidPressure = LiquidPressure + LiquidPressureIncrement

                  Particles(ParticleID)%WaterPressure = LiquidPressure

                  if(Particles(ParticleID)%WaterPressure > CalParams%LiquidPressureCavitationThreshold) then
                      Particles(ParticleID)%WaterPressure = CalParams%LiquidPressureCavitationThreshold
                  end if

                  !--------------------------------------
                  ! Correction for MPL in LIQUID Status in which the (traction) mean pressure is set to zero
                  if(CalParams%TwoLayerApplyNoTensStressLiqMPwLiqStatus) then
                      if((TwoLayerData%Elements(IElement)%ContainedMaterialTypes == ContainedMaterialTypeLIQUID).and. &
                          (Particles(ParticleID)%PhaseStatus == PhaseStatusLIQUID))  then
                          if(Particles(ParticleID)%WaterPressure > 0.0) then
                              Particles(ParticleID)%WaterPressure = 0.0
                          end if
                      end if
                  end if
                  !--------------------------------------


              end if

          end if ! 2 constituent

                                          
      end subroutine GetPressureLiquidMaterialPoint
      
      Subroutine GetBinghamYieldStress (IntGlo, YieldStress)
      !**********************************************************************
      !
      !    Function: get Bingham Yield Stress 
      !
      !********************************************************************* 
         
         implicit none
         integer, intent (in) :: intGlo !particle index
         integer :: MatIndex
         real(REAL_TYPE) :: YieldStress, pressure, phi
         
          Pressure = 0.0
          YieldStress = 0.0
          MatIndex = MaterialIDArray(IntGlo)
          Pressure = (SigmaEffArray(IntGlo,1)+SigmaEffArray(IntGlo,2)+ &
                     SigmaEffArray(IntGlo,3))/3.0
          
          Pressure = min(Pressure, 0.0)
          
          phi = MatParams(MatIndex)%FrictionAngle
          if ((CalParams%ApplyFrictionalFluid).or. &!v2016
           MatParams(MatIndex)%MaterialModel == ESM_FRICTIONAL_LIQUID) then
              YieldStress = -Pressure * SIN(phi*Pi/180.0)
          else
              YieldStress = MatParams(MatIndex)%BinghamYieldStress
          end if    
         
         end subroutine GetBinghamYieldStress
      
      Subroutine GetApparentViscosity (StrainRate, YieldStress, DynamicViscosity, ApparentViscosity)
      !**********************************************************************
      !
      !    Function: get apparent viscosity for Bingham Fluid model
      !
      !********************************************************************* 
      
         implicit none
         real(REAL_TYPE), dimension (6), intent (in):: StrainRate
         real(REAL_TYPE), intent (in):: yieldStress, DynamicViscosity
         real(REAL_TYPE) :: ApparentViscosity, ShearStrainRateNorm
         real(REAL_TYPE), parameter :: tiny=1.0d-17, big = 1.0d+10
         
         !compute second invariant of deviatoric strain tensor
         ShearStrainRateNorm = sqrt(( (StrainRate(1)-StrainRate(2))**2 + &
                                    (StrainRate(1)-StrainRate(3))**2 + &
                                    (StrainRate(2)-StrainRate(3))**2  ) / 3.0 + &
                                  (StrainRate(4)**2 + &
                                   StrainRate(5)**2 + &
                                   StrainRate(6)**2  ) * 2.0)
         
         ShearStrainRateNorm = max(ShearStrainRateNorm,tiny) !shearstrain not less than a minimum
         
         !calculate apparent viscosity
         ApparentViscosity = DynamicViscosity + yieldStress / ShearStrainRateNorm
         
         ApparentViscosity = min(ApparentViscosity,big) !apparent viscosity not larger than a maximum
             
         
         end subroutine GetApparentViscosity
      
      
       subroutine DetectLiquidFreeSurface()
       !**********************************************************************
       !
       !    Function: detect liquid free surface
       !
       !********************************************************************* 
        
        implicit none
        
          ! Local variables
          integer :: MaterialIndex, ParticleID, IParticle, IElement
          integer :: IAEl
          real(REAL_TYPE) :: RatioDensity
         
          if(.not.CalParams%ApplyDetectLiquidFreeSurface) RETURN
          
          do IAEl = 1, Counters%NAEl
              IElement = ActiveElement(IAEl)

              IsElemWithLiquidFreeSurfMP(IElement) = .false. !initialize

              do IParticle = 1, NPartEle(IElement) ! loop over all material points in element
                  ParticleID = GetParticleIndex(IParticle, IElement)

                  if (MaterialPointTypeArray(ParticleID) /= MaterialPointTypeLiquid) CYCLE

                  if (NFORMULATION == 1) then ! 1 constituent

                      MaterialIndex = MaterialIDArray(ParticleID)

                      !Determine interpolated density
                      call CalculateRatioDensity(IElement,ParticleID,RatioDensity)

                      ! Detect Free surface
                      if (RatioDensity < CalParams%FreeSurfaceFactor) then
                          Particles(ParticleID)%LiquidFreeSurface = 1.0
                          Particles(ParticleID)%WaterPressure = 0.0
                          IsElemWithLiquidFreeSurfMP(IElement) = .true. !at least 1 liquid free surface material point is contained in this element
                      end if

                  end if ! 1 constituent

                  !! -----------------------------------------------------------------------------

                  if (NFORMULATION == 2) then ! 2 constituents

                      ! Detect Free surface
                      if (Particles(ParticleID)%FillingRatio < CalParams%FreeSurfaceFactor) then
                          Particles(ParticleID)%LiquidFreeSurface = 1.0
                          Particles(ParticleID)%WaterPressure = 0.0
                          IsElemWithLiquidFreeSurfMP(IElement) = .true. !at least 1 liquid free surface material point is contained in this element
                      end if

                  end if ! 2 constituent

              end do
          end do


                                          
      end subroutine DetectLiquidFreeSurface
      
      subroutine CalculateRatioDensity(IElement,ParticleID,RatioDensity)
       !**********************************************************************
       !
       !    Function: calculate ratio density (liquid material point density/fluid threshold density)
       !
       !********************************************************************* 
      
        implicit none
        
              ! Local variables
              integer(INTEGER_TYPE) :: MaterialIndex, INode, ParticleID, IElement
              real(REAL_TYPE) :: PDensity, PartShape, RatioDensity
              integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs

              MaterialIndex = MaterialIDArray(ParticleID)
              
              !Determine interpolated density
              PDensity = 0.0  !initialize the interpolated density field for the element
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              do INode = 1,ELEMENTNODES
                  PartShape = ShapeValuesArray(ParticleID,INode)
                  PDensity = PDensity + PartShape * NodalDensity(NodeIDs(INode))
              end do
              RatioDensity = PDensity/(MatParams(MaterialIndex)%FluidThresholdDensity/1000)
              
      end subroutine CalculateRatioDensity

      subroutine CheckComputeUpdatedPressure(ParticleID,IElement,LiquidPressure,DoComputeUpdatedPressure)
      !**********************************************************************
      !
      !    Function: Calculates the Ratio Density (rho_{L,MP}/rho_{L,0}) used to determine if the amount an element is filled with a given density of MPs.
      !               The ratio density is compared to the Free Surface Factor (F_{FreeSurf}
      !
      !********************************************************************* 
      
      implicit none
            
          ! input argument 
          integer(INTEGER_TYPE), intent(in) :: ParticleID
          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          real(REAL_TYPE) :: RatioDensity
          logical :: PressureGreaterThreshold, PressureLowerThreshold, FullyFilled, FreeSurface 
          !output arguments
          real(REAL_TYPE) :: LiquidPressure
          logical :: DoComputeUpdatedPressure
          
          call CalculateRatioDensity(IElement,ParticleID,RatioDensity) !ratio density needed for the FullyFilled Calculation
          
               
           PressureGreaterThreshold = (LiquidPressure >= CalParams%LiquidPressureCavitationThreshold)
           PressureLowerThreshold = (LiquidPressure < CalParams%LiquidPressureCavitationThreshold)
           FreeSurface = (Particles(ParticleID)%LiquidFreeSurface == 1.0)
           
           if  (NFORMULATION == 1) then
               FullyFilled = (RatioDensity > 1.0)
           elseif (NFORMULATION == 2) then
               if(TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL > 0.0) then
                   FullyFilled = (Particles(ParticleID)%FillingRatio > 0.98)
               else
                   FullyFilled = (Particles(ParticleID)%FillingRatio > 1.0)
               end if
           end if 
                     
           DoComputeUpdatedPressure = (((.not.FreeSurface).and.PressureLowerThreshold).or. &
              ((PressureGreaterThreshold.or.FreeSurface).and.FullyFilled))  
      
      end subroutine CheckComputeUpdatedPressure
      
     end module ModLiquid
