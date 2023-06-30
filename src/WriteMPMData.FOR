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
	  
	  
	  module ModWriteMPMData
      !**********************************************************************
      !
      !  Function : Contains routines for writing MPM data to the binary files.
      !
      !             In order to keep the size of this source file reasonably small,
      !             this module only contains routines that write data to the binary files.
      !
      !  Implemented in the frame of the MPM project.
      !
      !     $Revision: 9808 $
      !     $Date: 2022-10-13 16:48:47 +0200 (do, 13 okt 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModMPMData
      use ModGlobalConstants

      implicit none

      contains


        subroutine AppendMPMData(FileUnit)
        !**********************************************************************
        !
        !    Function:  Appends the binary result file by blocks 
        !               which contain the particle data.
        !
        !     FileUnit : Unit assigned to the result file
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Add the individual blocks to the binary result file
          call WriteBlockPARTIC_I(FileUnit) ! Write initial number of particles
          call WriteBlockKERNELFI(FileUnit) ! Write int calculation flags
          call WriteBlockKERNELFR(FileUnit) ! Write double precision calculation flags
          call WriteBlockKERNELVB(FileUnit) ! Write used Kernel version
          call WriteBlockPADWAT_R(FileUnit)
          call WriteBlockACCELW_R(FileUnit) ! Write particle water acceleration data
          call WriteBlockPADUW__R(FileUnit) ! Write water particle incremental displacements
          if (CalParams%NumberOfPhases==3) then ! Write particle data relate to consolidation
            call WriteBlockPADGAS_R(FileUnit)
          end if
          call WriteBlockPADAT__R(FileUnit) ! Write double precision particle data
          call WriteBlockPUDMATPR(FileUnit) ! Write bulk modulus and porosity
          call WriteBlockTWOPNT_R(FileUnit) ! Write initial porosity
          call WriteBlockPADAT__I(FileUnit) ! Write int particle data
          call WriteBlockPARTFACI(FileUnit) ! Write decimal factor
          call WriteBlockNUMPARTI(FileUnit) ! Write number of particles per element
          call WriteBlockELPART_I(FileUnit) ! Write element-particle assignment
          call WriteBlockELPARTHI(FileUnit) ! Write element-particle assignment helper array
          call WriteBlockPARTCOOR(FileUnit) ! Write initial particle coordinates
          call WriteBlockACCEL__R(FileUnit) ! Write particle acceleration data
          call WriteBlockVELOC__R(FileUnit) ! Write particle velocity data
          call WriteBlockPAUTOT_R(FileUnit) ! Write particle total displacements
          call WriteBlockPADU___R(FileUnit) ! Write particle incremental displacements
          call WriteBlockPAPU___R(FileUnit) ! Write particle phase displacements
          call WriteBlockBOUNPARI(FileUnit) ! Write boundary particle status
          call WriteBlockMATPAR_I(FileUnit) ! Write material/virtual particle status
          call WriteBlockMATTYPEI(FileUnit) ! Write particle material type (soil, water ... )
          call WriteBlockPSIGXX_R(FileUnit) ! Write particle stresses xx
          call WriteBlockPSIGYY_R(FileUnit) ! Write particle stresses yy
          call WriteBlockPSIGZZ_R(FileUnit) ! Write particle stresses zz
          call WriteBlockPSIGXY_R(FileUnit) ! Write particle stresses xy
          if (NDIM==3) call WriteBlockPSIGYZ_R(FileUnit) ! Write particle stresses yz
          if (NDIM==3) call WriteBlockPSIGXZ_R(FileUnit) ! Write particle stresses xz
          call WriteBlockWPRESS_R(FileUnit) ! Write particle water pressures
          call WriteBlockGPRESS_R(FileUnit) ! Write particle gas pressures
          call WriteBlockUSEPARTI(FileUnit) ! Write element integration type
          call WriteBlockPARTVARR(FileUnit) ! Write particle arbitrary (testing) data
          call WriteBlockPDEPSXXR(FileUnit) ! Write particle incremental strains xx
          call WriteBlockPDEPSYYR(FileUnit) ! Write particle incremental strains yy
          call WriteBlockPDEPSZZR(FileUnit) ! Write particle incremental strains zz 
          call WriteBlockPDEPSXYR(FileUnit) ! Write particle incremental strains xy
          if (NDIM==3) call WriteBlockPDEPSYZR(FileUnit) ! Write particle incremental strains yz
          if (NDIM==3) call WriteBlockPDEPSXZR(FileUnit) ! Write particle incremental strains xz
          call WriteBlockPEPS_XXR(FileUnit) ! Write particle total strains xx
          call WriteBlockPEPS_YYR(FileUnit) ! Write particle total strains yy
          call WriteBlockPEPS_ZZR(FileUnit) ! Write particle total strains zz
          call WriteBlockPEPS_XYR(FileUnit) ! Write particle total strains xy
          if (NDIM==3) call WriteBlockPEPS_YZR(FileUnit) ! Write particle total strains yz
          if (NDIM==3) call WriteBlockPEPS_XZR(FileUnit) ! Write particle total strains xz
          call WriteBlockPIPL___I(FileUnit) ! Write particle plasticity state
          call WriteBlockDIDL___R(FileUnit) ! damping          
          call WriteBlockEPSP_XXR(FileUnit) ! Write particle total plastic strains xx
          call WriteBlockEPSP_YYR(FileUnit) ! Write particle total plastic strains yy
          call WriteBlockEPSP_ZZR(FileUnit) ! Write particle total plastic strains zz
          call WriteBlockEPSP_XYR(FileUnit) ! Write particle total plastic strains xy
          if (NDIM==3) call WriteBlockEPSP_YZR(FileUnit) ! Write particle total plastic strains yz
          if (NDIM==3) call WriteBlockEPSP_XZR(FileUnit) ! Write particle total plastic strains xz
		  call WriteBlockDPBVD__R(FileUnit)
          call WriteBlockESMstatev(FileUnit)  ! Write user defined material model state variables
          !call WriteBlockESMprops(FileUnit)  ! Write user defined material model properties
          call WriteBlockESM_UnloadingStiffness(FileUnit)  ! Write user defined material model unloading stiffness

        end subroutine AppendMPMData


        subroutine WriteBlockPARTIC_I(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PARTIC_I containing the number of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PARTIC_I is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts
        
          NByts = 4
          write(FileUnit) '$$PARTIC_I$$', NByts 
          write(FileUnit) Counters%NParticles
        
        end subroutine WriteBlockPARTIC_I


        subroutine WriteBlockPARTCOOR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PARTCOOR containing particle coordinates
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PARTCOOR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NDOFL
          NByts = NVals * 8
          write(FileUnit) '$$PARTCOOR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) (GlobPosArray(I,J), J = 1, NDOFL)
          end do

        end subroutine WriteBlockPARTCOOR


        subroutine WriteBlockKERNELFI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block KERNELFI containing calculation flags
        !               of type integer to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block KERNELFI is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals
        
          NVals = 8
          NByts = NVals * 4
          write(FileUnit) '$$KERNELFI$$', NByts
          
          if (IsULFEMComputation()) then
            write(FileUnit) 1
          else
            write(FileUnit) 0
          end if

          write(FileUnit) CalParams%NMaterialPoints

         

          if (CalParams%AssembleSoilLoadVectorFromMP) then
            write(FileUnit) 1
          else
            write(FileUnit) 0
          end if

          write(FileUnit) 1

          if (CalParams%ApplyUpdateWeights) then
            write(FileUnit) 1
          else
            write(FileUnit) 0
          end if

          if (CalParams%ApplyPrescribedDisplacements) then
            write(FileUnit) 1
          else
            write(FileUnit) 0
          end if

          if (CalParams%ApplyMeshSmoothing) then
            write(FileUnit) 1
          else
            write(FileUnit) 0
          end if

          write(FileUnit) CalParams%CalculationType

        end subroutine WriteBlockKERNELFI


        subroutine WriteBlockKERNELFR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block KERNELFR containing calculation flags
        !               of type double to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block KERNELFR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals

          NVals = 3
          NByts = NVals * 8
          write(FileUnit) '$$KERNELFR$$', NByts

          write(FileUnit) CalParams%RequiredDegreeOfFilling
          write(FileUnit) CalParams%FacStiffnessIncrease
          write(FileUnit) 0.0

        end subroutine WriteBlockKERNELFR


        subroutine WriteBlockKERNELVB(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block KERNELVB containing the version of the Kernel used
        !               to calculate the project to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block KERNELVB is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts

          NByts = Len_Trim('1.1')
          write(FileUnit) '$$KERNELVB$$', NByts

          write(FileUnit) '1.1'

        end subroutine WriteBlockKERNELVB


        subroutine WriteBlockPADAT__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PADAT__R containing basic particle information
        !               of type double to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PADAT__R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: FExt
          integer(INTEGER_TYPE) :: NByts, NVals, I, J,K

          NVals = 7 + ( 3 * NVECTOR )+ NVECTOR * MAX_LOAD_SYSTEMS
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$PADAT__R$$', NByts

          do I = 1, Counters%NParticles
            FExt = GetFExt(Particles(I) )

            write(FileUnit) MassArray(I)
            write(FileUnit) Particles(I)%MaterialWeight
            write(FileUnit) Particles(I)%Density
            write(FileUnit) Particles(I)%IntegrationWeight
            write(FileUnit) (GlobPosArray(I,J), J = 1, NVECTOR)
            write(FileUnit) (Particles(I)%LocPos(J), J = 1, NVECTOR)
            write(FileUnit) (Particles(I)%FBody(J), J = 1, NVECTOR)
            do K=1, MAX_LOAD_SYSTEMS
              write(FileUnit) (FExt(J,K), J = 1, NVECTOR)
            end do
            write(FileUnit) Particles(I)%ShearModulus
            write(FileUnit) Particles(I)%CohesionCosPhi
            write(FileUnit) Particles(I)%SFail
          end do
          
        end subroutine WriteBlockPADAT__R


        subroutine WriteBlockPUDMATPR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PUDMATPR containing porosity and bulk modulus of water.
        !
        !     FileUnit : File unit to which the block PUDMATPR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: FileUnit
          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = 2
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$PUDMATPR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%BulkWater
            write(FileUnit) Particles(I)%Porosity
          end do

        end subroutine WriteBlockPUDMATPR


        subroutine WriteBlockTWOPNT_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block TWOPNT_R containing the InitialPorosity
        !
        !     FileUnit : File unit to which the block TWOPNT_R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = 8
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$TWOPNT_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%InitialPorosity
            write(FileUnit) Particles(I)%ConstDensity
            write(FileUnit) Particles(I)%ConcentrationRatioLiquidL
            write(FileUnit) Particles(I)%ConcentrationRatioSolidL
            write(FileUnit) Particles(I)%ConcentrationRatioLiquidS
            write(FileUnit) Particles(I)%EffConcentrationRatioSolid
            write(FileUnit) Particles(I)%EffConcentrationRatioLiquid
            write(FileUnit) Particles(I)%EffPorosity
          end do

		end subroutine WriteBlockTWOPNT_R
		

        subroutine WriteBlockDPBVD__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block DPBVD__R containing BulkViscousPressure 
        !               and RateVolStrain
        !
        !     FileUnit : File unit to which the block DPBVD__R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = 1
          NByts = (Counters%NParticles + Counters%NEl) * NVals * 8
          write(FileUnit) '$$DPBVD__R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%DBulkViscousPressure
          end do
          do I = 1, Counters%NEl
            write(FileUnit) RateVolStrain(I)
          end do

        end subroutine WriteBlockDPBVD__R


        subroutine WriteBlockESMstatev(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block ESMsttvR containing state variables
        !               of the user defined soil modle
        !
        !     FileUnit : File unit to which the block ESMstvR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = 50
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$ESMsttvR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) (ESMstatevArray(I,J), J = 1, NSTATEVAR)
          end do

        end subroutine WriteBlockESMstatev
        
        !subroutine WriteBlockESMprops(FileUnit)
        !!**********************************************************************
        !!
        !!    Function:  Write the block ESMPROPSR containing material properties
        !!               of the user defined soil modle
        !!
        !!     FileUnit : File unit to which the block ESMpropsR is written to
        !!
        !!    Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !implicit none
        !
        !  ! Arguments
        !  integer(INTEGER_TYPE), intent(in) :: FileUnit
        !
        !  ! Local variables
        !  integer(INTEGER_TYPE) :: NByts, NVals, I, J
        !
        !  NVals = NPROPERTIES
        !  NByts = Counters%NParticles * NVals * 8
        !  write(FileUnit) '$$ESMpropsR$$', NByts
        !
        !  do I = 1, Counters%NParticles
        !    write(FileUnit) (ESMpropsArray(I,J), J = 1, NPROPERTIES)
        !  end do
        !
        !end subroutine WriteBlockESMprops

        subroutine WriteBlockESM_UnloadingStiffness(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block ESM_Eunloading containing state variables
        !               of the user defined soil model
        !
        !     FileUnit : File unit to which the variable Unloading Stiffness  is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = 1
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$ESMUnldS$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%ESM_UnloadingStiffness
          end do

        end subroutine WriteBlockESM_UnloadingStiffness
        
        subroutine WriteBlockPADAT__I(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PADAT__I containing basic particle information
        !               of type integer to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PADAT__I is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = 4
          NByts = Counters%NParticles * NVals * 4
          write(FileUnit) '$$PADAT__I$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) IDArray(I)
            write(FileUnit) ElementIDArray(I)
            write(FileUnit) EntityIDArray(I)
            write(FileUnit) MaterialIDArray(I)
          end do

        end subroutine WriteBlockPADAT__I


        subroutine WriteBlockPADWAT_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PADWAT_R containing particle information
        !               (relate to liquid) of type double to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PADWAT_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: PrescrDisp
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: FExtWater
          integer(INTEGER_TYPE) :: NByts, NVals, I, J, MaterialID,K
          logical :: IsUndrEffectiveStress
          
          NVals = 10 + ( 4 * NVECTOR )+ MAX_LOAD_SYSTEMS*NVECTOR
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$PADWAT_R$$', NByts

          do I = 1, Counters%NParticles
            FExtWater = GetFExtWater(Particles(I) )
            PrescrDisp = GetPrescrDisp(Particles(I) )
            MaterialID = MaterialIDArray(I)
            IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialID)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MaterialID)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))

            write(FileUnit) MassWaterArray(I)
            write(FileUnit) Particles(I)%MassMixed
            write(FileUnit) Particles(I)%WaterWeight
            write(FileUnit) Particles(I)%MixedWeight
            write(FileUnit) Particles(I)%Conductivity
            write(FileUnit) Particles(I)%Porosity
            write(FileUnit) Particles(I)%BulkWater
            write(FileUnit) Particles(I)%WaterVolumetricStrain
            write(FileUnit) Particles(I)%WaterPressure
            write(FileUnit) Particles(I)%DegreeSaturation
            write(FileUnit) (Particles(I)%FBodyWater(J), J = 1, NVECTOR)
            write(FileUnit) (Particles(I)%FBodyMixed(J), J = 1, NVECTOR)
            do k=1,MAX_LOAD_SYSTEMS
              write(FileUnit) (FExtWater(J,K), J = 1, NVECTOR)
            END DO
            if (CalParams%ApplyEffectiveStressAnalysis) then 
              write(FileUnit) (UArray(I,J), J = 1, NVECTOR)  !for undrained conditions water displacement equal to soil displacement
              write(FileUnit) (VelocityArray(I, J), J = 1, NVECTOR )   !for undrained conditions water velocity equal to soil velocity
            else
              write(FileUnit) (Particles(I)%UW(J), J = 1, NVECTOR)
              write(FileUnit) (VelocityWaterArray(I, J), J = 1, NVECTOR ) 
            end if
          end do
          
        end subroutine WriteBlockPADWAT_R


        subroutine WriteBlockPADGAS_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PADGAS_R containing particle information
        !               (relate to unsaturated consolidation) of type double to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PADGAS_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: FExtGas
          integer(INTEGER_TYPE) :: NByts, NVals, I, J,K

          NVals = 14 + ( 3 * NVECTOR )+ NVECTOR*MAX_LOAD_SYSTEMS
          NByts = Counters%NParticles * NVals * 8
          write(FileUnit) '$$PADGAS_R$$', NByts

          do I = 1, Counters%NParticles
            FExtGas = GetFExtGas(Particles(I) )

            write(FileUnit) Particles(I)%MassGas
            write(FileUnit) Particles(I)%MassMixed
            write(FileUnit) Particles(I)%GasWeight
            write(FileUnit) Particles(I)%MixedWeight
            write(FileUnit) Particles(I)%ConductivityGas
            write(FileUnit) Particles(I)%Porosity
            write(FileUnit) Particles(I)%DegreeSaturation
            write(FileUnit) Particles(I)%BulkGas
            write(FileUnit) Particles(I)%AirInWaterMassFraction
            write(FileUnit) Particles(I)%VapourInGasMassFraction
            write(FileUnit) Particles(I)%DiffusionAirInWater
            write(FileUnit) Particles(I)%DiffusionVapourInGas
            write(FileUnit) Particles(I)%GasVolumetricStrain
            write(FileUnit) Particles(I)%GasPressure
            write(FileUnit) (Particles(I)%FBodyGas(J), J = 1, NVECTOR)
            do k=1, MAX_LOAD_SYSTEMS
            write(FileUnit) (FExtGas(J,K), J = 1, NVECTOR)
            END DO
            write(FileUnit) (Particles(I)%UG(J), J = 1, NVECTOR)
            write(FileUnit) (VelocityGasArray(I, J), J = 1, NVECTOR) 
          end do

        end subroutine WriteBlockPADGAS_R

        
        subroutine WriteBlockPARTFACI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PARTFACI containing the decimal factor used for
        !               defining the particle - element connectivities in the array EleParticles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PARTFACI is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts
          integer(kind = 8) :: LocalPartFactor
          integer(INTEGER_TYPE) :: Factor1, Factor2

          LocalPartFactor = GetPartFactor()

          if (LocalPartFactor>2147483647) then
            Factor1 = 1000000000
            Factor2 = LocalPartFactor / Factor1
          else
            Factor1 = 1
            Factor2 = LocalPartFactor
          end if

          NByts = 2 * 4
          write(FileUnit) '$$PARTFACI$$', NByts

          write(FileUnit) Factor1, Factor2

        end subroutine WriteBlockPARTFACI


        subroutine WriteBlockNUMPARTI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block NUMPARTI containing the number of particles per element
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block NUMPARTI is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NEl
          NByts = NVals * 4
          write(FileUnit) '$$NUMPARTI$$', NByts

          write(FileUnit) (NPartEle(I), I = 1, NVals)

        end subroutine WriteBlockNUMPARTI


        subroutine WriteBlockELPART_I(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block ELPART_I containing the particle - element connectivities
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block ELPART_I is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8 ! integer (Kind = 8)

          write(FileUnit) '$$ELPART_I$$', NByts

          write(FileUnit) (EleParticles(I),  I = 1, Counters%NParticles)

        end subroutine WriteBlockELPART_I


        subroutine WriteBlockELPARTHI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block ELPARTHI containing the helper array used for
        !               storing the particle - element connectivities
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block ELPARTHI is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NEl
          NByts = NVals *4
          write(FileUnit) '$$ELPARTHI$$', NByts

          write(FileUnit) (EleParticlesHelp(I), I = 1, Counters%NEl)

        end subroutine WriteBlockELPARTHI


        subroutine WriteBlockACCELW_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block ACCELW_R containing the particle accelerations
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block ACCELW_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$ACCELW_R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) Particles(I)%AccelerationWater(NVECTOR)
            end do  
          end do

        end subroutine WriteBlockACCELW_R


        subroutine WriteBlockACCEL__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block ACCEL__R containing the particle accelerations
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block ACCEL__R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$ACCEL__R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR
              write(FileUnit) AccelerationArray(I, J)
            end do  
          end do

        end subroutine WriteBlockACCEL__R


        subroutine WriteBlockVELOC__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block VELOC__R containing the particle velocities
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block VELOC__R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$VELOC__R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) VelocityArray(I, J)
            end do
          end do

        end subroutine WriteBlockVELOC__R


        subroutine WriteBlockPAUTOT_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PAUTOT_R containing the total particle displacements
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PAUTOT_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$PAUTOT_R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) UArray(I,J)
            end do
          end do

        end subroutine WriteBlockPAUTOT_R


        subroutine WriteBlockPADU___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PADU___R containing the incremental particle displacements
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PADU___R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$PADU___R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) UStepArray(I,J)
            end do  
          end do

        end subroutine WriteBlockPADU___R


        subroutine WriteBlockPADUW__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PADUW__R containing the incremental water particle displacements
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PADUW__R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$PADUW__R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) Particles(I)%UStepWater(J)
            end do
          end do

        end subroutine WriteBlockPADUW__R


        subroutine WriteBlockPAPU___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PAPU___R containing the phase particle displacements
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PAPU___R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$PAPU___R$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) VelocityWaterArray(I, J)
            end do
          end do

        end subroutine WriteBlockPAPU___R


        subroutine WriteBlockBOUNPARI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block BOUNPARI containing the status whether a particle
        !               is located on the boundary of a body defined by particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block BOUNPARI is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 4
          write(FileUnit) '$$BOUNPARI$$', NByts

          write(FileUnit) (Particles(I)%IsBoundaryParticle, I = 1, Counters%NParticles)

        end subroutine WriteBlockBOUNPARI


        subroutine WriteBlockMATPAR_I(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block MATPAR_I containing the status whether a particle
        !               is a material or virtual particle to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block MATPAR_I is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 4
          write(FileUnit) '$$MATPAR_I$$', NByts

          write(FileUnit) (Particles(I)%Kind,  I = 1, Counters%NParticles)

        end subroutine WriteBlockMATPAR_I


        subroutine WriteBlockMATTYPEI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block MATTYPEI containing the type of a material point
        !               (mixture, solid, liquid, gas) to the file with unit FileUnit.
        !
        !    FileUnit : File unit to which the block MATTYPEI is written to
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, MatPointType

          NVals = Counters%NParticles
          NByts = NVals * 4
          write(FileUnit) '$$MATTYPEI$$', NByts

          do I=1,Counters%NParticles
              if (MaterialPointTypeArray(I)==MaterialPointTypeMixture) then
                MatPointType = 1
              end if
              if (MaterialPointTypeArray(I)==MaterialPointTypeSolid) then
                MatPointType = 2
              end if
              if (MaterialPointTypeArray(I)==MaterialPointTypeLiquid) then
                MatPointType = 3
              end if

              write(FileUnit) MatPointType
          end do

        end subroutine WriteBlockMATTYPEI


        subroutine WriteBlockPSIGXX_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PSIGXX_R containing the stresses in direction XX of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PSIGXX_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PSIGXX_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) SigmaEffArray(I, 1)
          end do

        end subroutine WriteBlockPSIGXX_R


        subroutine WriteBlockPSIGYY_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PSIGYY_R containing the stresses in direction YY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PSIGYY_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PSIGYY_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) SigmaEffArray(I, 2)
          end do

        end subroutine WriteBlockPSIGYY_R


        subroutine WriteBlockPSIGZZ_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PSIGZZ_R containing the stresses in direction ZZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PSIGZZ_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PSIGZZ_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) SigmaEffArray(I, 3)
          end do

        end subroutine WriteBlockPSIGZZ_R


        subroutine WriteBlockPSIGXY_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PSIGXY_R containing the stresses in direction XY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PSIGXY_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PSIGXY_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) SigmaEffArray(I, 4)
          end do

        end subroutine WriteBlockPSIGXY_R


        subroutine WriteBlockPSIGYZ_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PSIGYZ_R containing the stresses in direction YZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PSIGYZ_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PSIGYZ_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) SigmaEffArray(I, 5)
          end do

        end subroutine WriteBlockPSIGYZ_R


        subroutine WriteBlockPSIGXZ_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PSIGXZ_R containing the stresses in direction XZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PSIGXZ_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PSIGXZ_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) SigmaEffArray(I, 6)
          end do

        end subroutine WriteBlockPSIGXZ_R


        subroutine WriteBlockWPRESS_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block WPRESS_R containing the water pressure of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block WPRESS_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$WPRESS_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%WaterPressure
          end do

        end subroutine WriteBlockWPRESS_R


        subroutine WriteBlockGPRESS_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block GPRESS_R containing the gas pressure of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block GPRESS_R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$GPRESS_R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%GasPressure
          end do

        end subroutine WriteBlockGPRESS_R


        subroutine WriteBlockUSEPARTI(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block USEPARTI containing the flag per element indicating
        !               whether particle based or Gauss Point integration is used
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block USEPARTI is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NEl
          NByts = NVals * 4
          write(FileUnit) '$$USEPARTI$$', NByts
          
          write(FileUnit) (IsParticleIntegration(I), I = 1, Counters%NEl)

        end subroutine WriteBlockUSEPARTI


        subroutine WriteBlockPARTVARR(FileUnit)
        !**********************************************************************
        !
        !  Function:  Write the block PARTVARR containing three values per material pointthat
        !             can be used for displaying any material point information (for testing purposes)
        !             to the file with unit FileUnit.
        !
        !  FileUnit : File unit to which the block PARTVARR is written to
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I, J

          NVals = Counters%NParticles * NVECTOR
          NByts = NVals * 8
          write(FileUnit) '$$PARTVARR$$', NByts

          do I = 1, Counters%NParticles
            do J = 1, NVECTOR  
              write(FileUnit) GetVariableDataI(Particles(I), J)
            end do  
          end do

        end subroutine WriteBlockPARTVARR


        subroutine WriteBlockPDEPSXXR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PDEPSXXR containing the incremental strains in direction XX of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PDEPSXXR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PDEPSXXR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsStepI(Particles(I), 1)
          end do

        end subroutine WriteBlockPDEPSXXR


        subroutine WriteBlockPDEPSYYR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PDEPSYYR containing the incremental strains in direction YY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PDEPSYYR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PDEPSYYR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsStepI(Particles(I), 2)
          end do

        end subroutine WriteBlockPDEPSYYR


        subroutine WriteBlockPDEPSZZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PDEPSZZR containing the incremental strains in direction ZZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PDEPSZZRX is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PDEPSZZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsStepI(Particles(I), 3)
          end do

        end subroutine WriteBlockPDEPSZZR


        subroutine WriteBlockPDEPSXYR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PDEPSXYR containing the incremental strains in direction XY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PDEPSXYR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PDEPSXYR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsStepI(Particles(I), 4)
          end do

        end subroutine WriteBlockPDEPSXYR


        subroutine WriteBlockPDEPSYZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PDEPSYZR containing the incremental strains in direction YZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PDEPSYZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PDEPSYZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsStepI(Particles(I), 5)
          end do

        end subroutine WriteBlockPDEPSYZR


        subroutine WriteBlockPDEPSXZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PDEPSXZR containing the incremental strains in direction XZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PDEPSXZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PDEPSXZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsStepI(Particles(I), 6)
          end do

        end subroutine WriteBlockPDEPSXZR


        subroutine WriteBlockPEPS_XXR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PEPS_XXR containing the total strains in direction XX of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PEPS_XXR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PEPS_XXR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsI(Particles(I), 1)
          end do

        end subroutine WriteBlockPEPS_XXR


        subroutine WriteBlockPEPS_YYR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PEPS_YYR containing the total strains in direction YY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PEPS_YYR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PEPS_YYR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsI(Particles(I), 2)
          end do

        end subroutine WriteBlockPEPS_YYR


        subroutine WriteBlockPEPS_ZZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PEPS_ZZR containing the total strains in direction ZZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PEPS_ZZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PEPS_ZZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsI(Particles(I), 3)
          end do

        end subroutine WriteBlockPEPS_ZZR


        subroutine WriteBlockPEPS_XYR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PEPS_XYR containing the total strains in direction XY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PEPS_XYR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PEPS_XYR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsI(Particles(I), 4)
          end do

        end subroutine WriteBlockPEPS_XYR


        subroutine WriteBlockPEPS_YZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PEPS_YZR containing the total strains in direction YZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PEPS_YZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PEPS_YZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsI(Particles(I), 5)
          end do

        end subroutine WriteBlockPEPS_YZR


        subroutine WriteBlockPEPS_XZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PEPS_XZR containing the total strains in direction XZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PEPS_XZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$PEPS_XZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsI(Particles(I), 6)
          end do

        end subroutine WriteBlockPEPS_XZR


        subroutine WriteBlockEPSP_XXR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EPSP_XXR containing the total plastic strains in direction XX of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block EPSP_XXR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$EPSP_XXR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsPI(Particles(I), 1)
          end do

        end subroutine WriteBlockEPSP_XXR


        subroutine WriteBlockEPSP_YYR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EPSP_YYR containing the total plastic strains in direction YY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block EPSP_YYR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$EPSP_YYR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsPI(Particles(I), 2)
          end do

        end subroutine WriteBlockEPSP_YYR


        subroutine WriteBlockEPSP_ZZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EPSP_ZZR containing the total plastic trains in direction ZZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block EPSP_ZZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$EPSP_ZZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsPI(Particles(I), 3)
          end do

        end subroutine WriteBlockEPSP_ZZR


        subroutine WriteBlockEPSP_XYR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EPSP_XYR containing the total plastic strains in direction XY of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block EPSP_XYR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$EPSP_XYR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsPI(Particles(I), 4)
          end do

        end subroutine WriteBlockEPSP_XYR


        subroutine WriteBlockEPSP_YZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EPSP_YZR containing the total plastic strains in direction YZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block EPSP_YZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$EPSP_YZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsPI(Particles(I), 5)
          end do

        end subroutine WriteBlockEPSP_YZR


        subroutine WriteBlockEPSP_XZR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EPSP_XZR containing the total plastic strains in direction XZ of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block EPSP_XZR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$EPSP_XZR$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) GetEpsPI(Particles(I), 6)
          end do

        end subroutine WriteBlockEPSP_XZR


        subroutine WriteBlockPIPL___I(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block PIPL___I containing the plasticity state of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block PIPL___I is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 4
          write(FileUnit) '$$PIPL___I$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%IPL
          end do

        end subroutine WriteBlockPIPL___I


        subroutine WriteBlockDIDL___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block DIDL___R containing the damping of particles
        !               to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block DIDL___R is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NByts, NVals, I

          NVals = Counters%NParticles
          NByts = NVals * 8
          write(FileUnit) '$$DIDL___R$$', NByts

          do I = 1, Counters%NParticles
            write(FileUnit) Particles(I)%Damping
          end do

        end subroutine WriteBlockDIDL___R

        
      end module ModWriteMPMData
