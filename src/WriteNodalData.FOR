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
	  
	  
	  module ModWriteNodalData
      !**********************************************************************
      !
      !  Function : Contains routines for writing nodal data to the binary files.
      !
      !             In order to keep the size of this source file reasonably small,
      !             this module only contains routines that write data to the binary files.
      !
      !  Implemented in the frame of the MPM project.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModCounters
      use ModMeshInfo
      use ModDynViscousBoudary
      use ModMPMDYN2PhaseSP
      use ModMPMDYN3PhaseSP
      use ModMPMData
      use ModMatrixMath
      
      implicit none

      contains


        subroutine AppendNodalData(FileUnit)
        !**********************************************************************
        !
        !    Function:  Appends the binary result file by blocks 
        !               which contain the nodal data.
        !
        !     FileUnit : Unit assigned to the result file
        !     NodeCoord : Updated nodal coordinates after mesh adjustment
        !     UnrotatedNodalCoordinates : unrotated nodal coordinates after mesh adjustment
        !     PreviousRotationAngle : rotation angle
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          call WriteBlockUTOT___R(FileUnit)
          call WriteBlockVTOT___R(FileUnit)
          
          if (CalParams%ApplyAbsorbingBoundary) then
            call WriteBlockDSPTSD_R(FileUnit)
            call WriteBlockSPSD___R(FileUnit)
            call WriteBlockUTOTVBSR(FileUnit)
            call WriteBlockFVISSLDR(FileUnit)
          end if ! Absorbing boundary

          if ((CalParams%NumberOfPhases==2).or.(CalParams%NumberOfPhases==3)) then
            call WriteBlockUWTOT__R(FileUnit)
            call WriteBlockWTOT___R(FileUnit)
            call WriteBlockEXLDWT_R(FileUnit)

            if (CalParams%ApplyAbsorbingBoundary) then
              call WriteBlockDSPTWT_R(FileUnit)
              call WriteBlockSPWT___R(FileUnit)
              call WriteBlockWTOTVBWR(FileUnit)
              call WriteBlockFVISWATR(FileUnit)
            end if ! Absorbing boundary
          end if ! Consolidation

          if (CalParams%NumberOfPhases==3) then
            call WriteBlockUGTOT__R(FileUnit)
            call WriteBlockGTOT___R(FileUnit)
            call WriteBlockEXLDGS_R(FileUnit)

            if (CalParams%ApplyAbsorbingBoundary) then
              call WriteBlockDSPTGT_R(FileUnit)
              call WriteBlockSPGT___R(FileUnit)
              call WriteBlockWTOTVBGR(FileUnit)
              call WriteBlockFVISGASR(FileUnit)
            end if ! Absorbing boundary
          end if  ! Unsaturated calculation

          if (CalParams%ApplyMeshSmoothing) then
            call WriteBlockNODECOOR(FileUnit) ! Write updated nodal coordinates
          end if

        end subroutine AppendNodalData

    
        subroutine WriteBlockUTOT___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block UTOT___R containing the total displacement 
        !               of the solid phase.
        !
        !     FileUnit : File unit to which the block UTOT___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J


          NValues = Counters%NodTot * NDOFL
          if (Is0Arr(TotalDisplacementSoil, Counters%N)) then
            NValues = 0
          end if
          NBytes = 8 * NValues
          write(FileUnit) '$$UTOT___R$$', NBytes

          if (NBytes>0) then
            do I = 1, Counters%NodTot
              I1 = EdgeNodeTyingsHOE(1, I)
              I3 = EdgeNodeTyingsHOE(2, I)
              write(FileUnit) ( (TotalDisplacementSoil(ReducedDof(I1) + J) + TotalDisplacementSoil(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
            end do
          end if

        end subroutine WriteBlockUTOT___R


        subroutine WriteBlockVTOT___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block VTOT___R containing the velocity 
        !               of the solid phase.
        !
        !     FileUnit : File unit to which the block VTOT___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$VTOT___R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1,I)
            I3 = EdgeNodeTyingsHOE(2,I)
            write(FileUnit) ( (TotalVelocitySoil(ReducedDof(I1) + J, 1) + TotalVelocitySoil(ReducedDof(I3) + J, 1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockVTOT___R

    
        subroutine WriteBlockDSPTSD_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block DSPTSD_R containing the dashpot values 
        !               of the solid phase.
        !
        !     FileUnit : File unit to which the block DSPTSD_R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !*********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$DSPTSD_R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (DashpotSld(ReducedDof(I1) + J) + DashpotSld(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockDSPTSD_R
    
    
        subroutine WriteBlockSPSD___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block SPSD___R containing the spring values 
        !               of the solid phase.
        !
        !     FileUnit : File unit to which the block SPSD___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$SPSD___R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (SpringSld(ReducedDof(I1) + J) + SpringSld(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockSPSD___R


        subroutine WriteBlockUTOTVBSR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block UTOTVBSR containing the displacement 
        !               of the solid phase due to the viscous boundary.
        !
        !     FileUnit : File unit to which the block UTOTVBSR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$UTOTVBSR$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit)( (NodalTotDisplacementVB(ReducedDof(I1) + J, 1) + NodalTotDisplacementVB(ReducedDof(I3) + J, 1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockUTOTVBSR


        subroutine WriteBlockFVISSLDR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block FVISSLDR containing the viscous force 
        !               of the solid phase.
        !
        !     FileUnit : File unit to which the block FVISSLDR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, IDof, IEntity

          NValues = Counters%N*Counters%NEntity
          NBytes = 8 * NValues
          write(FileUnit) '$$FVISSLDR$$', NBytes

          do IDof = 1, Counters%N
            do IEntity = 1, Counters%NEntity
              if (IsMPMComputation()) then
                write(FileUnit) VisDampForceSld (IDof, IEntity) 
              else ! FEM
                write(FileUnit) VisDampForceSld (IDof, IEntity) + SpringSld(IDof) * IncrementalDisplacementSoil(IDof, IEntity) &
                                + DashpotSld(IDof) * AccelerationSoil(IDof, IEntity) * CalParams%TimeIncrement 
              end if
            end do 
          end do 

        end subroutine WriteBlockFVISSLDR

    
        subroutine WriteBlockUWTOT__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block UWTOT__R containing the displacement 
        !               of the water phase.
        !
        !     FileUnit : File unit to which the block UWTOT__R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J


          NValues = Counters%NodTot * NDOFL
          if (Is0Arr(TotalDisplacementWater, Counters%N)) then
            NValues = 0
          end if
          NBytes = 8 * NValues
          write(FileUnit) '$$UWTOT__R$$', NBytes

          if (NBytes>0) then
            do I = 1, Counters%NodTot
              I1 = EdgeNodeTyingsHOE(1, I)
              I3 = EdgeNodeTyingsHOE(2, I)
              write(FileUnit) ((TotalDisplacementWater(ReducedDof(I1)+J)+TotalDisplacementWater(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
            end do
          end if

        end subroutine WriteBlockUWTOT__R


        subroutine WriteBlockWTOT___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block WTOT___R containing the velocity 
        !               of the water phase.
        !
        !     FileUnit : File unit to which the block WTOT___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$WTOT___R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (TotalVelocityWater(ReducedDof(I1) + J, 1) + TotalVelocityWater(ReducedDof(I3) + J, 1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockWTOT___R


        subroutine WriteBlockEXLDWT_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EXLDWT_R containing the external load 
        !               of the water phase.
        !
        !     FileUnit : File unit to which the block EXLDWT_R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$EXLDWT_R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (ExtLoadWaterTotal(ReducedDof(I1) + J, 1,1) + ExtLoadWaterTotal(ReducedDof(I3) + J, 1,1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockEXLDWT_R


        subroutine WriteBlockDSPTWT_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block DSPTWT_R containing the dashpot values 
        !               of the water phase.
        !
        !     FileUnit : File unit to which the block DSPTWT_R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$DSPTWT_R$$', NBytes
          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (DashpotWat(ReducedDof(I1) + J) + DashpotWat(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockDSPTWT_R


        subroutine WriteBlockSPWT___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block SPWT___R containing the spring values 
        !               of the water phase.
        !
        !     FileUnit : File unit to which the block SPWT___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************


        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$SPWT___R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (SpringWat(ReducedDof(I1) + J) + SpringWat(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockSPWT___R


        subroutine WriteBlockWTOTVBWR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block WTOTVBWR containing the displacement 
        !               of the water phase due to the viscous boundary.
        !
        !     FileUnit : File unit to which the block WTOTVBWR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$WTOTVBWR$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (NodalTotDisplacementVBWater(ReducedDof(I1) + J, 1) + NodalTotDisplacementVBWater(ReducedDof(I3) + J, 1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockWTOTVBWR


        subroutine WriteBlockFVISWATR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block FVISWATR containing the viscous force 
        !               of the water phase.
        !
        !     FileUnit : File unit to which the block FVISWATR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************


        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, IDof, IEntity

          NValues = Counters%N*Counters%NEntity
          NBytes = 8 * NValues
          write(FileUnit) '$$FVISWATR$$', NBytes

          do IDof = 1, Counters%N
            do IEntity = 1, Counters%NEntity
              write(FileUnit) VisDampForceWat(IDof, IEntity) + SpringWat(IDof) * IncrementalDisplacementWater(IDof, IEntity) &
                             + DashpotWat(IDof) * AccelerationWater(IDof, IEntity) * CalParams%TimeIncrement
            end do
          end do 

        end subroutine WriteBlockFVISWATR


        subroutine WriteBlockUGTOT__R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block UGTOT__R containing the displacement 
        !               of the gas phase.
        !
        !     FileUnit : File unit to which the block UGTOT__R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J


          NValues = Counters%NodTot * NDOFL
          if (Is0Arr(TotalDisplacementGas, Counters%N)) then
            NValues = 0
          end if
          NBytes = 8 * NValues
          write(FileUnit) '$$UGTOT__R$$', NBytes
          if (NBytes>0) then
            do I = 1, Counters%NodTot
              I1 = EdgeNodeTyingsHOE(1, I)
              I3 = EdgeNodeTyingsHOE(2, I)
              write(FileUnit) ( (TotalDisplacementGas(ReducedDof(I1) + J) + TotalDisplacementGas(ReducedDof(I3) + J)) / 2, J=1,NDOFL)
            end do
          end if

        end subroutine WriteBlockUGTOT__R


        subroutine WriteBlockGTOT___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block GTOT___R containing the velocity 
        !               of the gas phase.
        !
        !     FileUnit : File unit to which the block GTOT___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$GTOT___R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (TotalVelocityGas(ReducedDof(I1) + J, 1) + TotalVelocityGas(ReducedDof(I3) + J, 1)) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockGTOT___R


        subroutine WriteBlockEXLDGS_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block EXLDGS_R containing the external load 
        !               of the gas phase.
        !
        !     FileUnit : File unit to which the block EXLDGS_R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************


        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$EXLDGS_R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (ExtLoadGasTotal(ReducedDof(I1) + J, 1,1) + ExtLoadGasTotal(ReducedDof(I3) + J, 1,1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockEXLDGS_R


        subroutine WriteBlockDSPTGT_R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block DSPTGT_R containing the dashpot values 
        !               of the gas phase.
        !
        !     FileUnit : File unit to which the block DSPTGT_R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************


        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$DSPTGT_R$$', NBytes
          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (DashpotGas(ReducedDof(I1) + J) + DashpotGas(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockDSPTGT_R


        subroutine WriteBlockSPGT___R(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block SPGT___R containing the spring values 
        !               of the gas phase.
        !
        !     FileUnit : File unit to which the block SPGT___R is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************


        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$SPGT___R$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (SpringGas(ReducedDof(I1) + J) + SpringGas(ReducedDof(I3) + J) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockSPGT___R


        subroutine WriteBlockWTOTVBGR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block WTOTVBGR containing the displacement 
        !               of the gas phase due to the viscous boundary.
        !
        !     FileUnit : File unit to which the block WTOTVBGR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, I1, I3, J

          NValues = Counters%NodTot * NDOFL
          NBytes = 8 * NValues
          write(FileUnit) '$$WTOTVBGR$$', NBytes

          do I = 1, Counters%NodTot
            I1 = EdgeNodeTyingsHOE(1, I)
            I3 = EdgeNodeTyingsHOE(2, I)
            write(FileUnit) ( (NodalTotDisplacementVBGas(ReducedDof(I1) + J, 1) + NodalTotDisplacementVBGas(ReducedDof(I3) + J, 1) ) / 2, J=1,NDOFL)
          end do

        end subroutine WriteBlockWTOTVBGR


        subroutine WriteBlockFVISGASR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block FVISGASR containing the viscous force 
        !               of the gas phase.
        !
        !     FileUnit : File unit to which the block FVISGASR is written to
        !
        !    Implemented in the frame of the MPM project.
        !
        !**********************************************************************


        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, IDof, IEntity

          NValues = Counters%N*Counters%NEntity
          NBytes = 8 * NValues
          write(FileUnit) '$$FVISGASR$$', NBytes

          do IDof = 1, Counters%N
            do IEntity = 1, Counters%NEntity
               write(FileUnit) VisDampForceGas (IDof, IEntity) + SpringGas(IDof) * IncrementalDisplacementGas(IDof, IEntity) &
                              + DashpotGas(IDof) * AccelerationGas(IDof, IEntity) * CalParams%TimeIncrement
            end do
          end do  

        end subroutine WriteBlockFVISGASR

        
        subroutine WriteBlockNODECOOR(FileUnit)
        !**********************************************************************
        !
        !    Function:  Write the block NODECOOR containing the nodal coordinates x, y, z
        !               after mesh adjustment has been performed to the file with unit FileUnit.
        !
        !     FileUnit : File unit to which the block NODECOOR is written to
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Arguments
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! Local variables
          integer(INTEGER_TYPE) :: NBytes, NValues, I, J

          NValues = Counters%NodTot * NDOFL
          NBytes = NValues * 8
          write(FileUnit) '$$NODECOOR$$', NBytes
          
          if (NBytes>0) then
            do I = 1, Counters%NodTot
              write(FileUnit) (NodalCoordinates(I, J), J = 1, NDOFL)
            end do
          end if

        end subroutine WriteBlockNODECOOR

      end module