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
    
    
	  module ModMPMDYNStresses
      !**********************************************************************
      !
      !  Function: Contains the routines related to calculating the stresses at material points
      !            for partially filled elements and Gauss points for fully filled elements.
      !            This module is only used by the dynamic MPM.  
      !
      !            Note: 'Integration point' denotes either a Gauss point or material point. Except for
      !            the data that is used or modified, both types of points are treated identically.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters
      use ModReadCalculationData
      use ModReadMaterialData
      use ModElementEvaluation
      use ModMPMData
      use ModMPMInit
      use ModWriteTestData
      use ModMPMStresses
      use ModStrainSmoothing
      use ModMeshInfo
      use ModLiquid
      use ModTwoLayerFormulation
      use ModExternalSoilModel
      use user32
      use kernel32
      use ModMPMDYN2PhaseSP
       
      contains	  
	  
        

           
        
        
        subroutine GetFluidPressure(IDpt)
        !*************************************************************************************
        !  Function: compute liquid and gas pressures. Solve mass banance equations
        !            Valid for single-point formulation, 2phase and 3phase analysis
        !
        !*************************************************************************************
        
        implicit none
        
        ! local variables
        integer(INTEGER_TYPE), intent(in) :: IDpt ! global integration/material point number

        integer(INTEGER_TYPE) :: I     ! counter
        integer(INTEGER_TYPE) :: IDset ! ID of material parameter set
        integer(INTEGER_TYPE) :: ntens ! Dimension of stress vector to pass to UDSM 

        real(REAL_TYPE), dimension(MatParams(MaterialIDArray_N(IDpt))%UMATDimension) :: StrainIncr ! strain increment in integration/material point
        real(REAL_TYPE), dimension(NTENSOR) :: TempStrainIncr                                    ! incremental strain vector assigned to point 

        real(REAL_TYPE) :: DSigWP    ! Change of water pressure at integration/material point 
        real(REAL_TYPE) :: DSigGP    ! Change of gas pressure at integration/material point 
        real(REAL_TYPE) :: Kf, WaterAdvectiveFlux, WaterDensity, g   !Bulk of water, Advective flow, water density, gravity acceleration
        real(REAL_TYPE) :: N         ! Porosity         
        real(REAL_TYPE) :: DEpsVol   ! Incremental volumetric strain
        real(REAL_TYPE) :: DEpsVolW  ! Incremental volumetric strain (water)
        real(REAL_TYPE) :: dT        ! Change of Temperature at integration point
        real(REAL_TYPE) :: lambda    ! Used to compute water pressure
        real(REAL_TYPE) :: Sr, dSrdp(1) !Degree of saturation, derivative dSr/dp_w
      
       if (.not.(((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)).and.(NFORMULATION==1))) RETURN !Valid for single-point formulation, 2phase and 3phase analysis
    
       IDset = MaterialIDArray_N(IDpt) ! get constitutive model in integration/material point. It is the material number stored in $$MATERIAL_INDEX in the GOM-file
       
       ntens = MatParams(IDset)%UMATDimension     
          
       TempStrainIncr = GetEpsStep(Particles(IDpt)) ! get strain increments in integration/material point. It is the incremental strain vector assigned to point   
        
       StrainIncr = 0.0
       do I=1, NTENSOR
          StrainIncr(I) = StrainIncr(I) + TempStrainIncr(I)
       enddo 
        
       DEpsVol = StrainIncr(1) + StrainIncr(2) + StrainIncr(3) ! volumetric strain, valid for 2D and 3D      
             
       ! initalise water pressure
       DSigWP = 0.0d0
       DSigGP = 0.0d0

       ! for 2-phase analysis
       if ((CalParams%NumberOfPhases==2).and.(NFORMULATION==1)) then
        if (Particles(IDpt)%WaterWeight > 0.0) then
         DEpsVolW = Particles(IDpt)%WaterVolumetricStrain ! Water phase
         N = Particles(IDpt)%Porosity
         Kf = Particles(IDpt)%BulkWater
        
         if (CalParams%ApplyPartialSaturation) then
           Sr = Particles(IDPt)%DegreeSaturation
           WaterAdvectiveFlux = Particles(IDpt)%WaterAdvectiveFlux        
           g = CalParams%GravityData%GAccel                    !Gravity (m/s2)
           WaterDensity = Particles(IDpt)%WaterWeight/g
           call CalculateDerivDegreeSaturation(IDPt,dSrdp,1)
        
           if ((Sr<1).and.(Sr>0)) then              !(consider advective terms)
             lambda = N*Sr*WaterDensity/Kf - N*WaterDensity * dSrdp(1)
             lambda = 1/lambda
             DSigWP = lambda*((Sr*WaterDensity)*DEpsVol + WaterAdvectiveFlux)   
           else
             lambda = Kf/N
             DSigWP = lambda*( N * DEpsVolW + (1.0d0 - N) * DEpsVol)  
           end if
        
         else
           lambda = Kf/N
           DSigWP = lambda*( N * DEpsVolW + (1.0d0 - N) * DEpsVol)
         end if
        
         ! for submerged calculation
         if (CalParams%ApplySubmergedCalculation) then 
           if (CalParams%IStep <= CalParams%NumberSubmergedCalculation) then 
             DSigWP = 0.0 ! excess pore pressure is zero in gravity phase
           end if
         end if ! submerged calculation
              
         else 
          DSigWP = 0.0
         end if
        end if ! 2-phase analysis

       ! for 3-phase analysis (unsaturated soil)
       if (CalParams%NumberOfPhases==3) then ! the incremental water pressure, gas pressure and temperature (for one material point) are obtained
         call SolveBalanceEquations(IDpt, DEpsVol, DSigWP, DSigGP, dT) ! solve mass balances (for unsaturated porous media) and energy balance
       end if ! 3-phase analysis

        ! single phase quasi-static consolidation
        if (CalParams%ApplyImplicitQuasiStatic) then
          DsigWP = Particles(IDpt)%WaterPressure - Particles(IDpt)%WaterPressure0
        end if
    
        call AssignWatandGasPressureToGlobalArray(IDpt, DSigWP, DSigGP)  !Note that the subroutine checks Cavitation Threshold & Gas Pressure
      
        end subroutine GetFluidPressure

    subroutine SmoothenLiquidPressureIncrement(IMatSet)
    !**********************************************************************
    !
    !    Function: Somoothen pressure increment of water phase
    !
    !    IMatSet : The number of the material set
    !
    !**********************************************************************
       
        implicit none
 
          integer(INTEGER_TYPE), intent(in) :: IMatSet
          ! Local variables  
          real(REAL_TYPE), dimension(Counters%NodTot,2) :: dPressureSmooth
          real(REAL_TYPE) :: dWP, Weight, dWPWeighted, EnhancedPressureIncrement
          integer(INTEGER_TYPE) :: IEl, IAEl, I, NElemPart, MaterialID, IParticle, ParticleIndex
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
        
          dPressureSmooth = 0.0

          !Calculate nodal incremental pressure mapping from MP to element nodes    
          do IAEl = 1, Counters%NAEl  ! active element
            IEl = ActiveElement_N(IAEl)
            LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
            if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
              NElemPart = NumberOfIntegrationPoints(IEl)
              do IParticle = 1, NElemPart
                ParticleIndex = GetParticleIndex_N(IParticle, IEl)
                MaterialID = MaterialIDArray_N(ParticleIndex)
                
                if (MaterialID==IMatSet)then ! This particle belongs to the considered entity (Material)
                  ! The volumetric strain in this element for the considered entity (Material)
                  dWP = Particles(ParticleIndex)%WaterPressure - Particles(ParticleIndex)%WaterPressure0
                  Weight = Particles(ParticleIndex)%IntegrationWeight
                  dWPWeighted = dWP * Weight
                  dPressureSmooth(LJ(:), 1) = dPressureSmooth(LJ(:), 1) + dWPWeighted
                  dPressureSmooth(LJ(:), 2) = dPressureSmooth(LJ(:), 2) + Weight
                end if
              end do ! particles

            end if ! entity
          end do ! active element     

          do I = 1, Counters%NodTot
            if(dPressureSmooth(I, 2)>0.0) then
              dPressureSmooth(I, 1) = dPressureSmooth(I, 1) / dPressureSmooth(I, 2) ! Only consider active nodes
            end if
          end do
       
          do IAEl = 1, Counters%NAEl  ! active element 
            IEl = ActiveElement_N(IAEl)
              if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
                LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
                EnhancedPressureIncrement = 0.0
                do I = 1, ELEMENTNODES
                  EnhancedPressureIncrement = EnhancedPressureIncrement + dPressureSmooth(LJ(I), 1) / dble(ELEMENTNODES)
                end do    

              ! Update particle pressure

                NElemPart = NumberOfIntegrationPoints(IEl)
                do IParticle = 1, NElemPart ! Loop over all particles of the element
                  ParticleIndex = GetParticleIndex_N(IParticle, IEl)
                  MaterialID = MaterialIDArray_N(ParticleIndex)
                    if (MaterialID==IMatSet)then ! this particle belongs to material IMatSet. Assign its strains
                      Particles(ParticleIndex)%WaterPressure = Particles(ParticleIndex)%WaterPressure0 + EnhancedPressureIncrement
                    end if
                end do ! particles
              end if ! material
    end  do ! active element
    
    end subroutine SmoothenLiquidPressureIncrement
    
        subroutine CalculateViscousDamping(ParticleID, IEl)
        !**********************************************************************
        !> Computes a pressure term introducing bulk viscosity damping to the equation of motion.
        !>
        !! \param[in] ParticleID ID of considered material point.
        !! \param[in] IEl ID of element of the considered material point.
        !! \param[in] DilationalWaveSpeed Current wave speed computed for the considered material point.
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleID, IEl

          real(REAL_TYPE) :: ViscousDampingPressure = 0.0
          real(REAL_TYPE) :: Density = 0.0
          real(REAL_TYPE) :: ElementLMinLocal = 0.0
          real(REAL_TYPE) :: RateVolStrainLocal = 0.0
          real(REAL_TYPE) :: MaterialIndex = 0.0
          real(REAL_TYPE) :: DilationalWaveSpeed = 0.0
          logical :: IsUndrEffectiveStress

          if (.not.CalParams%ApplyBulkViscosityDamping) return

          MaterialIndex = MaterialIDArray_N(ParticleID)
          
           IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
           
          !if (CalParams%ApplyEffectiveStressAnalysis
           if (IsUndrEffectiveStress &
          .or.((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3))) then
            Density = MatParams(MaterialIndex)%DensityMixture / 1000.0
          else
            Density = (1 - MatParams(MaterialIndex)%InitialPorosity) * MatParams(MaterialIndex)%DensitySolid / 1000.0
          end if

          ElementLMinLocal = ElementLMin(IEl)
          RateVolStrainLocal = RateVolStrain(IEl)

          
          call GetWaveSpeed(ParticleID, DilationalWaveSpeed)

          ViscousDampingPressure = CalParams%BulkViscosityDamping1 *  &
            Density * DilationalWaveSpeed * ElementLMinLocal * RateVolStrainLocal

          if ((RateVolStrainLocal < 0.0).and.(CalParams%BulkViscosityDamping2 > 0.0)) then
            ViscousDampingPressure = ViscousDampingPressure + &
              Density * (CalParams%BulkViscosityDamping2 * ElementLMinLocal * RateVolStrainLocal)**2
          end if

          Particles(ParticleID)%DBulkViscousPressure = ViscousDampingPressure

        end subroutine CalculateViscousDamping
        
        
        end module ModMPMDYNStresses
