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


      module ModReadMPMData
      !**********************************************************************
      !
      !    Function : Contains routines for reading MPM data from the binary result files,
      !               the PPD and RON file (in case of prescribed displacements).
      !            
      !               This module should only contain routines that read data that define 
      !               the material point discretisation.
      !
      !     $Revision: 9808 $
      !     $Date: 2022-10-13 16:48:47 +0200 (do, 13 okt 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModReadCalculationData
      use ModMPMData
      use ModParticle
      use ModMeshInfo
      use ModMPMDynContact
      use ModDynViscousBoudary
      use ModMPMDYN2PhaseSP
      use ModMPMDYN3PhaseSP
      use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
      use ModFileIO
      
      implicit none

      contains
      
        !subroutine ReadPreviousRotationAngle(BRFfilename, PreviousRotationAngle)
        !!**********************************************************************
        !!
        !!  Function : Reads the rotation angle from the BRF file
        !!             of the previous phase if the block exists.
        !!             The block exists, if mesh adjustment has also
        !!             been applied in the previous phase.
        !!
        !!  I    BRFfilename : Name of the BRF file of the previous load phase
        !! 
        !!  O    PreviousRotationAngle : previous rotation angle
        !!
        !!**********************************************************************
        !implicit none
        !
        !  ! Arguments
        !  character(len = *), intent(in) :: BRFfilename
        !  real(REAL_TYPE), intent(out) :: PreviousRotationAngle
        !  
        !  ! Local variables
        !  integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 1
        !  integer(INTEGER_TYPE) :: BlockCheck, NBytes
        !  character(len = 80) :: Header
        !  character(len = 12) :: BlockName
        !
        !  call FileOpenAction(BRFunit, BRFfilename, 'R')
        !  
        !  BlockCheck = 0
        !  if (FExist(BRFfilename) ) then
        !    read(BRFunit) Header
        !
        !    do
        !      read(BRFunit) BlockName, NBytes
        !
        !      if ( (BlockName=='$$ENDOFFII$$').or. &
        !           (BlockCheck==RELEVANTBLOCKS) ) then
        !        EXIT
        !      else
        !        if (NBytes>0) then
        !        
        !          if (BlockName=='$$PREVRA_R$$') then
        !            call ReadBlockPREVRA_R(BRFunit, NBytes, BlockCheck, PreviousRotationAngle)
        !          else
        !            call Skip(BRFunit, NBytes)
        !          end if  
        !        
        !        end if
        !      end if
        !    end do
        !    
        !    close(BRFunit)
        !  end if         
        !
        !end subroutine ReadPreviousRotationAngle
        !
        !
        !subroutine ReinitialiseUnrotatedNodalCoordinates(BRFfilename, NodTot, IDim, UnrotatedNodalCoordinates)
        !!**********************************************************************
        !!
        !!  Function : Reads the unrotated nodal coordinates from the BRF 
        !!             file of the previous phase if the block exists.
        !!             The block exists, if mesh adjustment has also been
        !!             applied in the previous phase.
        !!
        !! I    BRFfilename : Name of the BRF file of the previous load phase
        !! I    NodTot : Total number of nodes
        !! I    IDim : Number of dimensions
        !!
        !! O    UnrotatedNodalCoordinates : Unrotated nodal coordinates
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  ! Arguments
        !  integer(INTEGER_TYPE), intent(in) :: NodTot, IDim
        !  character(len = *), intent(in) :: BRFfilename
        !  real(REAL_TYPE), dimension(NodTot, IDim), intent(inout) :: UnrotatedNodalCoordinates
        !  
        !  ! Local variables
        !  integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 1
        !  integer(INTEGER_TYPE) :: BlockCheck, NBytes
        !  character(len = 80) :: Header
        !  character(len = 12) :: BlockName
        !
        !  call FileOpenAction(BRFunit, BRFfilename, 'R')
        !  
        !  BlockCheck = 0
        !  if (FExist(BRFfilename) ) then
        !    read(BRFunit) Header
        !
        !    do
        !      read(BRFunit) BlockName, NBytes
        !
        !      if ( (BlockName=='$$ENDOFFII$$').or. &
        !           (BlockCheck==RELEVANTBLOCKS) ) then
        !        EXIT
        !      else
        !        if (NBytes>0) then
        !        
        !          if (BlockName=='$$UNROTC_R$$') then
        !            call ReadBlockNODECOORtobechanged(BRFunit, NBytes, BlockCheck,  &
        !                                   NodTot, IDim,  &
        !                                   UnrotatedNodalCoordinates)
        !          else
        !            call Skip(BRFunit, NBytes)
        !          end if  
        !        
        !        end if
        !      end if
        !    end do
        !    
        !    close(BRFunit)
        !  end if
        !
        !end subroutine ReinitialiseUnrotatedNodalCoordinates
      

        subroutine ReadNodalCoordinates(BRFfilename, NodTot, IDim, InitialNodalCoord)
        !**********************************************************************
        !
        !  Function : Reads the initial nodal coordinates from the BRF
        !             file of the previous phase if the block exists.
        !             The block exists, if mesh adjustment has also
        !             been applied in the previous phase.
        !
        !     BRFfilename : Name of the BRF file of the previous load phase
        !     NodTot : Total number of nodes
        !     IDim : Number of dimensions
        !
        ! O   InitialNodalCoord : Initial nodal coordinates
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim
          character(len = *), intent(in) :: BRFfilename
          real(REAL_TYPE), dimension(NodTot, IDim), intent(inout) :: InitialNodalCoord
          
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 1
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')
          
          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$NODECOOR$$') then ! Read nodal coordinates
                    call ReadBlockNODECOORtobechanged(BRFunit, NBytes, BlockCheck, NodTot, IDim, InitialNodalCoord)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if  
                
                end if
              end if
            end do
            
            close(BRFunit)
          end if

        end subroutine ReadNodalCoordinates
        
        
        subroutine ReadBlockNODECOORtobechanged(FileUnit, NBytes, BlockCheck, &
                                                NodTot, IDim, NodeCoord)
        !**********************************************************************
        !
        !    Function:  Reads data from the block NODECOOR containing nodal coordinates
        !               updated through mesh adjustment from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block NODECOOR is read
        !   NBytes : Number of bytes inside the block
        !   NodTot : Total number of nodes
        !   IDim : Number of dimensions
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        ! O NodeCoord : Nodal coordinates
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          real(REAL_TYPE), dimension(NodTot, IDim),  &
            intent(inout) :: NodeCoord
          ! Local variables
          integer(INTEGER_TYPE) :: I, J

          if ( (NBytes / 8 / IDim)==NodTot) then

            do I = 1, NodTot
              read(FileUnit) (NodeCoord(I, J), J = 1, IDim)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
         
        end subroutine ReadBlockNODECOORtobechanged


        subroutine ReadHouseKeepingCountersFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Reads the total number of particles Counters%NParticles from 
        !             the BRF file of the previous phase. 
        !             The value is needed to reinitialise the particle house-keeping 
        !             arrays from the previous phase.
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************

        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 1
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')
          
          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header
            
            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$PARTIC_I$$') then ! Read number of particles
                    call ReadBlockPARTIC_I(BRFunit, NBytes, BlockCheck)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if  
                
                end if
              end if
            end do
            
            close(BRFunit)
          end if

          call Assert(BlockCheck == RELEVANTBLOCKS, 'Error in ReadHouseKeepingCountersFromFile')

        end subroutine ReadHouseKeepingCountersFromFile


        subroutine ReadHouseKeepingDataFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Fills the house-keeping data structures with data
        !             read from the BRF file of the previous phase.
        !             The particle data itself (coordinates, stresses, ... )
        !             is not read - only house-keeping data.
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************
        
        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 5
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')
          
          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$PARTFACI$$') then ! Read decimal factor
                    call ReadBlockPARTFACI(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$ELPART_I$$') then ! Read element-particle assignment
                    call ReadBlockELPART_I(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$ELPARTHI$$') then ! Read element-particle assignment helper array
                    call ReadBlockELPARTHI(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$USEPARTI$$') then ! Read element integration type
                    call ReadBlockUSEPARTI(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$NUMPARTI$$') then ! Read number of particles per element
                    call ReadBlockNUMPARTI(BRFunit, NBytes, BlockCheck)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if  
                
                end if
              end if
            end do
            
            close(BRFunit)
          end if
          
          call Assert(BlockCheck == RELEVANTBLOCKS, 'error in ReadHouseKeepingDataFromFile')

        end subroutine ReadHouseKeepingDataFromFile


        subroutine ReadParticleDataFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Fills the particle data with data read from the BRF 
        !             file of the previous phase (material ID, coordinates, ...). 
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************
        
        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 5
          integer(INTEGER_TYPE) :: BlockCheck, NBytes, BlockUnchecked
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')
         
          BlockCheck = 0 ! No relevant blocks found yet
          BlockUnchecked = 0
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$PADAT__R$$') then ! Read real(REAL_TYPE) particle data
                    call ReadBlockPADAT__R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PUDMATPR$$') then
                    call ReadBlockPUDMATPR(BRFunit, NBytes, BlockUnchecked)
                  else if (BlockName=='$$TWOPNT_R$$') then ! initial porosity
                    call ReadBlockTWOPNT_R(BRFunit, NBytes, BlockUnchecked)
                  else if (BlockName=='$$SSOFTMCR$$') then ! state parameters Strain Softening MC model
                    call ReadBlockSSOFTMCR(BRFunit, NBytes, BlockUnChecked)
                  else if (BlockName=='$$MCCMatPR$$') then ! state parameters of MCC model
                    call ReadBlockMCCMATPR(BRFunit, NBytes, BlockUnchecked)
                  else if (BlockName=='$$ESMsttvR$$') then ! state parameters of external model
                    call ReadBlockESMstatev(BRFunit, NBytes, BlockUnchecked)
                  !else if (BlockName=='$$ESMpropsR$$') then ! material parameters of external model
                  !  call ReadBlockESMprops(BRFunit, NBytes, BlockUnchecked)
                  else if (BlockName=='$$ESMUnldS$$') then ! Eunloading
                    call ReadBlockESMUnload(BRFunit, NBytes, BlockUnchecked)
                  else if (BlockName=='$$PICHPSPR$$') then ! state parameters of ICH and HP models
                    call ReadBlockPICHPSPR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$DPBVD__R$$') then
                    call ReadBlockDPBVD__R(BRFunit, NBytes)
                  else if (BlockName=='$$PADAT__I$$') then ! Read integer(INTEGER_TYPE) particle data
                    call ReadBlockPADAT__I(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$BOUNPARI$$') then ! Read boundary particle status
                    call ReadBlockBOUNPARI(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$MATPAR_I$$') then ! Read virtual/material particle status
                    call ReadBlockMATPAR_I(BRFunit, NBytes, BlockUnchecked)
                  else if (BlockName=='$$MATTYPEI$$') then
                    call ReadBlockMATTYPEI(BRFunit, NBytes, BlockUnchecked)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if  
                
                end if
              end if
            end do
            
            close(BRFunit)
          end if


        end subroutine ReadParticleDataFromFile


        subroutine ReadParticleStateParametersFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Fills the particle data with data read from the BRF
        !             file of the previous phase (displacements, stresses, strains, ...). 
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************
        
        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 24
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')

          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$PAUTOT_R$$') then ! Read total particle displacements
                    call ReadBlockPAUTOT_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PSIGXX_R$$') then ! Read particle stresses xx
                    call ReadBlockPSIGXX_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PSIGYY_R$$') then ! Read particle stresses yy
                    call ReadBlockPSIGYY_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PSIGZZ_R$$') then ! Read particle stresses zz
                    call ReadBlockPSIGZZ_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PSIGXY_R$$') then ! Read particle stresses xy
                    call ReadBlockPSIGXY_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PSIGYZ_R$$') then ! Read particle stresses yz
                    call ReadBlockPSIGYZ_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PSIGXZ_R$$') then ! Read particle stresses xz
                    call ReadBlockPSIGXZ_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$WPRESS_R$$') then ! Read particle water pressures
                    call ReadBlockWPRESS_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$GPRESS_R$$') then ! Read particle gas pressures
                    call ReadBlockGPRESS_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PEPS_XXR$$') then ! Read particle total strains xx
                    call ReadBlockPEPS_XXR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PEPS_YYR$$') then ! Read particle total strains yy
                    call ReadBlockPEPS_YYR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PEPS_ZZR$$') then ! Read particle total strains zz
                    call ReadBlockPEPS_ZZR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PEPS_XYR$$') then ! Read particle total strains xy
                    call ReadBlockPEPS_XYR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PEPS_YZR$$') then ! Read particle total strains yz
                    call ReadBlockPEPS_YZR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PEPS_XZR$$') then ! Read particle total strains xz
                    call ReadBlockPEPS_XZR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$EPSP_XXR$$') then ! Read particle total strains xx
                    call ReadBlockEPSP_XXR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$EPSP_YYR$$') then ! Read particle total strains yy
                    call ReadBlockEPSP_YYR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$EPSP_ZZR$$') then ! Read particle total strains zz
                    call ReadBlockEPSP_ZZR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$EPSP_XYR$$') then ! Read particle total strains xy
                    call ReadBlockEPSP_XYR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$EPSP_YZR$$') then ! Read particle total strains yz
                    call ReadBlockEPSP_YZR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$EPSP_XZR$$') then ! Read particle total strains xz
                    call ReadBlockEPSP_XZR(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PIPL___I$$') then ! Read particle plasticity state
                    call ReadBlockPIPL___I(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$DIDL___R$$') then ! Read Damping
                    call ReadBlockDIDL___R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PADU___R$$') then ! Read particle incremental displacements
                    call ReadBlockPADU___R(BRFunit, NBytes, BlockCheck)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if
                 
                end if
              end if
            end do
          
            close(BRFunit)
          end if

        end subroutine ReadParticleStateParametersFromFile


        subroutine ReadParticleDynamicDataFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Fills the dynamic related data with data read from 
        !             the BRF file of the previous load phase
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************
        
        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 3
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')

          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$VELOC__R$$') then ! Read particle velocity data
                    call ReadBlockVELOC__R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$ACCEL__R$$') then ! Read particle acceleration data
                    call ReadBlockACCEL__R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PAPU___R$$') then ! Read particle phase displacements
                    call ReadBlockPAPU___R(BRFunit, NBytes, BlockCheck)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if
                 
                end if
              end if
            end do
          
            close(BRFunit)
          end if
        
        end subroutine ReadParticleDynamicDataFromFile


        subroutine ReadConsolidationDataFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Fills the consolidation particle data with data read
        !             from the BRF file of the previous load phase 
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************
        
        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 3
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')

          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$PADWAT_R$$') then ! Read consolidation data 
                    call ReadBlockPADWAT_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$ACCELW_R$$') then ! Read particle acceleration data
                    call ReadBlockACCELW_R(BRFunit, NBytes, BlockCheck)
                  else if (BlockName=='$$PADUW__R$$') then ! Read particle incremental displacements
                    call ReadBlockPADUW__R(BRFunit, NBytes, BlockCheck)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if
                 
                end if
              end if
            end do
          
            close(BRFunit)
          end if
        
        end subroutine ReadConsolidationDataFromFile


        subroutine ReadConsolidationDataGasFromFile(BRFfilename)
        !**********************************************************************
        !
        !  Function : Fills the unsaturated (gas) particle data with 
        !             data from the BFR file of the previous load phase
        !
        !  BRFfilename : Name of the BRF file of the previous load phase
        !
        !**********************************************************************
        
        implicit none
        
          character(len = *), intent(in) :: BRFfilename
          ! Local variables
          integer(INTEGER_TYPE), parameter :: RELEVANTBLOCKS = 1
          integer(INTEGER_TYPE) :: BlockCheck, NBytes
          character(len = 80) :: Header
          character(len = 12) :: BlockName

          call FileOpenAction(BRFunit, BRFfilename, 'R')

          BlockCheck = 0 ! No relevant blocks found yet
          if (FExist(BRFfilename) ) then
            read(BRFunit) Header

            do
              read(BRFunit) BlockName, NBytes

              if ( (BlockName=='$$ENDOFFII$$').or.(BlockCheck==RELEVANTBLOCKS) ) then ! Close file
                EXIT
              else
                if (NBytes>0) then ! The block is not empty
                
                  if (BlockName=='$$PADGAS_R$$') then ! Read consolidation unsaturated data (gas data) 
                    call ReadBlockPADGAS_R(BRFunit, NBytes, BlockCheck)
                  else ! Skip bytes of block
                    call Skip(BRFunit, NBytes)
                  end if
                 
                end if
              end if
            end do
          
            close(BRFunit)
          end if
        
        end subroutine ReadConsolidationDataGasFromFile

        subroutine ReadBlockPARTIC_I(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PARTIC_I containing the number of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PARTIC_I is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
        
          if (NBytes==4) then
          
            read(FileUnit) Counters%NParticles
        
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockPARTIC_I

        
        subroutine ReadBlockPREVRA_R(FileUnit, NBytes, BlockCheck, PreviousRotationAngle)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PREVRA_R containing the previous 
        !               rotation angle of type double from the file with unit FileUnit.
        !
        ! I  FileUnit : File unit from which the block PREVRA_R is read
        ! I  NBytes : Number of bytes inside the block
        !
        ! O  BlockCheck : Increase by 1, if the block was successfully read
        ! O  PreviousRotationAngle : previous rotation angle
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          real(REAL_TYPE), intent(out) :: PreviousRotationAngle
        
          if (NBytes==8) then
          
            read(FileUnit) PreviousRotationAngle
            
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockPREVRA_R
                

        subroutine ReadBlockPADAT__R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PADAT__R containing basic particle information
        !               of type double from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PADAT__R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, MAX_LOAD_SYSTEMS) :: FExt
          integer(INTEGER_TYPE) :: I, J, NVals, K

          NVals = 7 + ( 3 * NVECTOR ) + MAX_LOAD_SYSTEMS * NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              read(FileUnit) MassArray(I)
              read(FileUnit) Particles(I)%MaterialWeight
              read(FileUnit) Particles(I)%Density
              read(FileUnit) Particles(I)%IntegrationWeight
              read(FileUnit) (GlobPosArray(I,J), J = 1, NVECTOR)
              read(FileUnit) (Particles(I)%LocPos(J), J = 1, NVECTOR)
              read(FileUnit) (Particles(I)%FBody(J), J = 1, NVECTOR)
              do K=1,MAX_LOAD_SYSTEMS
                read(FileUnit) (FExt(J,K), J = 1, NVECTOR)
              end do
              call SetFExt(Particles(I), FExt)
              read(FileUnit) Particles(I)%ShearModulus
              read(FileUnit) Particles(I)%CohesionCosPhi
              read(FileUnit) Particles(I)%SFail
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockPADAT__R
        
        subroutine ReadBlockPUDMATPR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PUDMATPR containing water properties.
        !
        !   FileUnit : File unit from which the block PUDMATPR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        !**********************************************************************
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals

          NVals = 2

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              read(FileUnit) Particles(I)%BulkWater
              read(FileUnit) Particles(I)%Porosity
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockPUDMATPR


        subroutine ReadBlockTWOPNT_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !   Function:  Reads data from the block TWOPNT_R containing initial porosity
        !              of type double from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block TWOPNT_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals

          NVals = 8

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              read(FileUnit) Particles(I)%InitialPorosity
              read(FileUnit) Particles(I)%ConstDensity
              read(FileUnit) Particles(I)%ConcentrationRatioLiquidL
              read(FileUnit) Particles(I)%ConcentrationRatioSolidL
              read(FileUnit) Particles(I)%ConcentrationRatioLiquidS
              read(FileUnit) Particles(I)%EffConcentrationRatioSolid
              read(FileUnit) Particles(I)%EffConcentrationRatioLiquid
              read(FileUnit) Particles(I)%EffPorosity
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
        
        end subroutine ReadBlockTWOPNT_R
       
        
        subroutine ReadBlockSSOFTMCR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block SSOFTMCR containing basic particle information
        !               of type double from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block SSOFTMCR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals

          NVals = 3

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              read(FileUnit) Particles(I)%CohesionStSoft
              read(FileUnit) Particles(I)%PhiStSoft
              read(FileUnit) Particles(I)%PsiStSoft
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockSSOFTMCR
        
        subroutine ReadBlockMCCMATPR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block MccMatPR containing the state
        !               variables of the MCC model
        !
        !   FileUnit : File unit from which the block MCCMATPR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals

          NVals = 2

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              read(FileUnit) Particles(I)%init_void_ratio
              read(FileUnit) Particles(I)%pp
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockMCCMATPR
        
        subroutine ReadBlockESMstatev(FileUnit, NBytes, BlockCheck)
          implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, NVals
        
          NVals = 50
          if ( (NBytes / NVals / 8)==Counters%NParticles) then
            do I = 1, Counters%NParticles
              read(FileUnit) (ESMstatevArray(I,J), J=1, NSTATEVAR)
            end do
            BlockCheck = BlockCheck + 1
          end if
        
        end subroutine ReadBlockESMstatev
        
        !subroutine ReadBlockESMprops(FileUnit, NBytes, BlockCheck)
        !  implicit none
        !
        !  integer(INTEGER_TYPE), intent(in) :: FileUnit
        !  integer(INTEGER_TYPE), intent(in) :: NBytes
        !  integer(INTEGER_TYPE), intent(inout) :: BlockCheck
        !  ! Local variables
        !  integer(INTEGER_TYPE) :: I, J, NVals
        !
        !  NVals = NPROPERTIES
        !  if ( (NBytes / NVals / 8)==Counters%NParticles) then
        !    do I = 1, Counters%NParticles
        !      read(FileUnit) (ESMpropsArray(I,J), J=1, NPROPERTIES)
        !    end do
        !    BlockCheck = BlockCheck + 1
        !  end if
        !
        !end subroutine ReadBlockESMprops

        subroutine ReadBlockESMUnload(FileUnit, NBytes, BlockCheck)
          implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals
        
          NVals = 1

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              read(FileUnit) Particles(I)%ESM_UnloadingStiffness
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
        
        end subroutine ReadBlockESMUnload
        
        subroutine ReadBlockPICHPSPR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PICHPSPR containing the state
        !               variables of the ICH and HP models.
        !
        !   FileUnit : File unit from which the block PICHPSPR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE) :: HP(2), HPIG(7), Eoed, Alpha, Beta
          integer(INTEGER_TYPE) :: I, J, NVals

          NVals = 16

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
              
              read(FileUnit) (HP(J), J = 1, 2)
              call SetHPStateVariables(Particles(I), HP)
              read(FileUnit) (HP(J), J = 1, 2)
              call SetModifiedHPStateVariables(Particles(I), HP)
              read(FileUnit) (HPIG(J), J = 1, 7)
              call SetHPIGStateVariables(Particles(I), HPIG)
              read(FileUnit) Eoed
              Particles(I)%EoedHP = Eoed
              read(FileUnit) Alpha
              Particles(I)%Alpha = Alpha
              read(FileUnit) Beta
              Particles(I)%Beta = Beta
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockPICHPSPR

        subroutine ReadBlockPADAT__I(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PADAT__I containing basic particle information
        !               of type integer(INTEGER_TYPE) from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PADAT__I is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals
        
          NVals = 4
          
          if ( (NBytes / NVals / 4)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) IDArray(I)
              read(FileUnit) ElementIDArray(I)
              read(FileUnit) EntityIDArray(I)
              read(FileUnit) MaterialIDArray(I)
            end do

            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPADAT__I


        subroutine ReadBlockDPBVD__R(FileUnit, NBytes)
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes

          integer(INTEGER_TYPE) :: I, NVals

          NVals = 1

          if ( (NBytes / NVals / 8)==(Counters%NParticles + Counters%NEl)) then
            do I = 1, Counters%NParticles
              read(FileUnit) Particles(I)%DBulkViscousPressure
            end do
            do I = 1, Counters%NEl
              read(FileUnit) RateVolStrain(I)
            end do
          end if
 
        end subroutine ReadBlockDPBVD__R
        
        subroutine ReadBlockPARTFACI(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PARTFACI containing the decimal factor used for
        !               defining the particle - element connectivities in the array EleParticles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PARTFACI is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: Factor1, Factor2
          integer(INTEGER_TYPE) :: LocalFactor1, LocalFactor2
          integer(kind = 8) :: NewPartFactor
        
          if (NBytes==8) then
          
            read(FileUnit) Factor1, Factor2
          
            LocalFactor1 = Factor1
            LocalFactor2 = Factor2
          
            NewPartFactor = LocalFactor1 * LocalFactor2          
            call SetPartFactor(NewPartFactor)

            BlockCheck = BlockCheck + 1
            
          end if
          
        end subroutine ReadBlockPARTFACI

        subroutine ReadBlockNUMPARTI(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block NUMPARTI containing the number of particles per element
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block NUMPARTI is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
        
          if ( (NBytes / 4)==Counters%NEl) then
            
            read(FileUnit) (NPartEle(I), I = 1, Counters%NEl)
            
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockNUMPARTI

        subroutine ReadBlockELPART_I(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block ELPART_I containing the particle - element connectivities
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block ELPART_I is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I

          if ( (NBytes / 8)==(Counters%NParticles) ) then

            read(FileUnit) (EleParticles(I), I = 1, Counters%NParticles)

            BlockCheck = BlockCheck + 1
            
          end if
                  
        end subroutine ReadBlockELPART_I

        subroutine ReadBlockELPARTHI(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block ELPARTHI containing the helper array used for
        !               storing the particle - element connectivities
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block ELPARTHI is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
        
          if ( (NBytes / 4)==Counters%NEl) then
            
            read(FileUnit) (EleParticlesHelp(I), I = 1, Counters%NEl)
            
            BlockCheck = BlockCheck + 1
            
          end if
         
        end subroutine ReadBlockELPARTHI

        subroutine ReadBlockVELOC__R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block VELOC__R containing the
        !               particle velocity.
        !
        !   FileUnit : File unit from which the block VELOC__R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: dum
          integer(INTEGER_TYPE) :: I, NVals, id

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
                read(FileUnit) (dum(id), id = 1, NVECTOR)
                VelocityArray(I,:) = dum
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockVELOC__R

        subroutine ReadBlockACCELW_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block ACCELW_R containing the
        !               particle acceleration.
        !
        !   FileUnit : File unit from which the block ACCELW_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: dum
          integer(INTEGER_TYPE) :: I, NVals, id

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
                read(FileUnit) (dum(id), id = 1, NVECTOR)
                Particles(I)%AccelerationWater(1:NVECTOR) = dum(1:NVECTOR)
            end do

            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockACCELW_R

        subroutine ReadBlockACCEL__R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block ACCEL__R containing the
        !               particle acceleration.
        !
        !   FileUnit : File unit from which the block ACCEL__R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR) :: dum
          integer(INTEGER_TYPE) :: I, NVals, id

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
                read(FileUnit) (dum(id), id = 1, NVECTOR)
                AccelerationArray(I,:) = dum
            end do

            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockACCEL__R

        subroutine ReadBlockPAPU___R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PAPU___R containing the
        !               particle phase dispalcement.
        !
        !   FileUnit : File unit from which the block PAPU___R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals, id

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles
                read(FileUnit) (UPhaseArray(I, id), id = 1, NVECTOR)
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockPAPU___R

        subroutine ReadBlockPADWAT_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PADWAT_R containing the
        !               consolidation particle information.
        !
        !   FileUnit : File unit from which the block PADWAT_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE) :: WatPressure
          real(REAL_TYPE), dimension(NVECTOR) :: PartVelocityWater
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: FExtWater
          integer(INTEGER_TYPE) :: I, J, NVals, K

          NVals = 10 + ( 4 * NVECTOR ) + NVECTOR*MAX_LOAD_SYSTEMS

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles

              read(FileUnit) MassWaterArray(I)
              read(FileUnit) Particles(I)%MassMixed
              read(FileUnit) Particles(I)%WaterWeight
              read(FileUnit) Particles(I)%MixedWeight
              read(FileUnit) Particles(I)%Conductivity
              read(FileUnit) Particles(I)%Porosity
              read(FileUnit) Particles(I)%BulkWater
              read(FileUnit) Particles(I)%WaterVolumetricStrain
              read(FileUnit) WatPressure
              Particles(I)%WaterPressure = WatPressure
              Particles(I)%WaterPressure0 = WatPressure
              read(FileUnit) Particles(I)%DegreeSaturation
              read(FileUnit) (Particles(I)%FBodyWater(J), J = 1, NVECTOR)
              read(FileUnit) (Particles(I)%FBodyMixed(J), J = 1, NVECTOR)
              do K=1,MAX_LOAD_SYSTEMS
                read(FileUnit) (FExtWater(J,K), J = 1, NVECTOR)
              end do
              call SetFExtWater(Particles(I), FExtWater)
              read(FileUnit) (Particles(I)%UW(J), J = 1, NVECTOR)
              read(FileUnit) (PartVelocityWater(J), J = 1, NVECTOR)
              VelocityWaterArray(I,:) = PartVelocityWater
             
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockPADWAT_R

        subroutine ReadBlockPADGAS_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PADGAS_R containing the
        !               unsaturated consolidation particle information.
        !
        !   FileUnit : File unit from which the block PADGAS_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          real(REAL_TYPE) :: GasPressure
          real(REAL_TYPE), dimension(NVECTOR) :: PartVelocityGas
          real(REAL_TYPE), dimension(NVECTOR,MAX_LOAD_SYSTEMS) :: FExtGas
          integer(INTEGER_TYPE) :: I, J, NVals,K

          NVals = 14 + ( 3 * NVECTOR )+ NVECTOR*MAX_LOAD_SYSTEMS

          if ( (NBytes / NVals / 8)==Counters%NParticles) then

            do I = 1, Counters%NParticles

              read(FileUnit) Particles(I)%MassGas
              read(FileUnit) Particles(I)%MassMixed
              read(FileUnit) Particles(I)%GasWeight
              read(FileUnit) Particles(I)%MixedWeight
              read(FileUnit) Particles(I)%ConductivityGas
              read(FileUnit) Particles(I)%Porosity
              read(FileUnit) Particles(I)%DegreeSaturation
              read(FileUnit) Particles(I)%BulkGas
              read(FileUnit) Particles(I)%AirInWaterMassFraction
              read(FileUnit) Particles(I)%VapourInGasMassFraction
              read(FileUnit) Particles(I)%DiffusionAirInWater
              read(FileUnit) Particles(I)%DiffusionVapourInGas
              read(FileUnit) Particles(I)%GasVolumetricStrain
              read(FileUnit) GasPressure
              Particles(I)%GasPressure = GasPressure
              Particles(I)%GasPressure0 = GasPressure
              read(FileUnit) (Particles(I)%FBodyGas(J), J = 1, NVECTOR)
              do K=1, MAX_LOAD_SYSTEMS
                  read(FileUnit) (FExtGas(J,K), J = 1, NVECTOR)
               end do
              call SetFExtGas(Particles(I), FExtGas)
              
              read(FileUnit) (Particles(I)%UG(J), J = 1, NVECTOR)
              read(FileUnit) (PartVelocityGas(J), J = 1, NVECTOR)
              VelocityGasArray(I,:) = PartVelocityGas
             
            end do
 
            BlockCheck = BlockCheck + 1
          
          end if
 
        end subroutine ReadBlockPADGAS_R
        

        subroutine ReadBlockPAUTOT_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PAUTOT_R containing the total particle displacements
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PAUTOT_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals, J

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) (UArray(I, J), J = 1, NVECTOR)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPAUTOT_R
        
        subroutine ReadBlockPADU___R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PADU___R containing the incremental particle displacements
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PADU___R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals, J

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) (UStepArray(I, J), J = 1, NVECTOR)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPADU___R

        subroutine ReadBlockPADUW__R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PADUW__R containing the incremental water particle displacements
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PADUW__R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, NVals, J

          NVals = NVECTOR

          if ( (NBytes / NVals / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) (Particles(I)%UStepWater(J), J = 1, NVECTOR)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPADUW__R

        subroutine ReadBlockBOUNPARI(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block BOUNPARI containing the status whether a particle
        !               is located on the boundary of a body defined by particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block BOUNPARI is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          
          if ( (NBytes / 4)==Counters%NParticles) then
          
            read(FileUnit) (Particles(I)%IsBoundaryParticle, I = 1, Counters%NParticles)
            
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockBOUNPARI

        subroutine ReadBlockMATPAR_I(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block MATPAR_I containing the status whether a particle
        !               is a virtual or material particle from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block MATPAR_I is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, IParticleCounter
          
          if ( (NBytes / 4)==Counters%NParticles) then
            IParticleCounter = 1
            do I = 1, Counters%NParticles
              read(FileUnit) Particles(IParticleCounter)%Kind
              IParticleCounter = IParticleCounter + 1
            
            end do
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockMATPAR_I

        subroutine ReadBlockMATTYPEI(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !   Function:  Read the block MATTYPEI containing the type of a material point
        !              (mixture, solid, liquid, gas) to the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block MATTYPEI is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        !**********************************************************************
        use ModGlobalConstants
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I, IParticleCounter, MatPointType
          
          if ( (NBytes / 4)==Counters%NParticles) then
            IParticleCounter = 1
            do I = 1, Counters%NParticles
              read(FileUnit) MatPointType
              
              if (MatPointType==1) then
                MaterialPointTypeArray(I) = MaterialPointTypeMixture
              end if
              if (MatPointType==2) then
                MaterialPointTypeArray(I) = MaterialPointTypeSolid
              end if
              if (MatPointType==3) then
                MaterialPointTypeArray(I) = MaterialPointTypeLiquid
              end if
              
              IParticleCounter = IParticleCounter + 1
            
            end do
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockMATTYPEI
        
        subroutine ReadBlockPSIGXX_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PSIGXX_R containing the stresses in direction XX of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PSIGXX_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StressValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StressValue
              SigmaEff0Array(I, 1) = StressValue
              SigmaEffArray(I, 1)  = StressValue
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
         
        end subroutine ReadBlockPSIGXX_R

        subroutine ReadBlockPSIGYY_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PSIGYY_R containing the stresses in direction YY of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PSIGYY_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StressValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StressValue
              SigmaEff0Array(I, 2) = StressValue
              SigmaEffArray(I, 2)  = StressValue
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPSIGYY_R

        subroutine ReadBlockPSIGZZ_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PSIGZZ_R containing the stresses in direction ZZ of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PSIGZZ_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StressValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StressValue
              SigmaEff0Array(I, 3) = StressValue
              SigmaEffArray(I, 3)  = StressValue
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPSIGZZ_R

        subroutine ReadBlockPSIGXY_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PSIGXY_R containing the stresses in direction XY of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PSIGXY_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StressValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StressValue
              SigmaEff0Array(I, 4) = StressValue
              SigmaEffArray(I, 4)  = StressValue
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPSIGXY_R

        subroutine ReadBlockPSIGYZ_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PSIGYZ_R containing the stresses in direction YZ of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PSIGYZ_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StressValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StressValue
              SigmaEff0Array(I, 5) = StressValue
              SigmaEffArray(I, 5)  = StressValue
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPSIGYZ_R

        subroutine ReadBlockPSIGXZ_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PSIGXZ_R containing the stresses in direction XZ of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PSIGXZ_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StressValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StressValue
              SigmaEff0Array(I, 6) = StressValue
              SigmaEffArray(I, 6)  = StressValue
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPSIGXZ_R

        subroutine ReadBlockWPRESS_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block WPRESS_R containing the water pressure of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block WPRESS_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I

          if ( (NBytes / 8)==Counters%NParticles) then
          
            read(FileUnit) (Particles(I)%WaterPressure0, I = 1, Counters%NParticles)
          
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockWPRESS_R

        subroutine ReadBlockGPRESS_R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block GPRESS_R containing the gas pressure of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block GPRESS_R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I

          if ( (NBytes / 8)==Counters%NParticles) then
          
            read(FileUnit) (Particles(I)%GasPressure0, I = 1, Counters%NParticles)
          
            BlockCheck = BlockCheck + 1
            
          end if
        
        end subroutine ReadBlockGPRESS_R

        subroutine ReadBlockUSEPARTI(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block USEPARTI containing the flag per element indicating
        !               whether particle based or Gauss Point integration is used
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block USEPARTI is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
        
          if ( (NBytes / 4)==Counters%NEl) then
            
            read(FileUnit) (IsParticleIntegration(I), I = 1, Counters%NEl)
            
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockUSEPARTI

        subroutine ReadBlockPEPS_XXR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PEPS_XXR containing the total strains 
        !               in direction XX of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PEPS_XXR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StrainValue
              call SetEpsI(Particles(I), 1, StrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPEPS_XXR

        subroutine ReadBlockPEPS_YYR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PEPS_YYR containing the total strains 
        !               in direction YY of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PEPS_YYR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StrainValue
              call SetEpsI(Particles(I), 2, StrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPEPS_YYR

        subroutine ReadBlockPEPS_ZZR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PEPS_ZZR containing the total strains 
        !               in direction ZZ of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PEPS_ZZR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StrainValue
              call SetEpsI(Particles(I), 3, StrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPEPS_ZZR

        subroutine ReadBlockPEPS_XYR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PEPS_XYR containing the total 
        !               strains in direction XY of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PEPS_XYR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StrainValue
              call SetEpsI(Particles(I), 4, StrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPEPS_XYR

        subroutine ReadBlockPEPS_YZR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PEPS_YZR containing the total 
        !               strains in direction YZ of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PEPS_YZR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StrainValue
              call SetEpsI(Particles(I), 5, StrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPEPS_YZR

        subroutine ReadBlockPEPS_XZR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PEPS_XZR containing the total strains 
        !               in direction XZ of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PEPS_XZR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: StrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) StrainValue
              call SetEpsI(Particles(I), 6, StrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockPEPS_XZR

        subroutine ReadBlockEPSP_XXR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block EPSP_XXR containing the total plastic strains 
        !               in direction XX of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block EPSP_XXR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: PlasticStrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) PlasticStrainValue
              call SetEpsPI(Particles(I), 1, PlasticStrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockEPSP_XXR

        subroutine ReadBlockEPSP_YYR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block EPSP_YYR containing the total plastic strains 
        !               in direction YY of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block EPSP_YYR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: PlasticStrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) PlasticStrainValue
              call SetEpsPI(Particles(I), 2, PlasticStrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockEPSP_YYR

        subroutine ReadBlockEPSP_ZZR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block EPSP_ZZR containing the total plastic strains 
        !               in direction ZZ of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block EPSP_ZZR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: PlasticStrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) PlasticStrainValue
              call SetEpsPI(Particles(I), 3, PlasticStrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockEPSP_ZZR

        subroutine ReadBlockEPSP_XYR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block EPSP_XYR containing the total plastic 
        !               strains in direction XY of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block EPSP_XYR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: PlasticStrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) PlasticStrainValue
              call SetEpsPI(Particles(I), 4, PlasticStrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockEPSP_XYR

        subroutine ReadBlockEPSP_YZR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block EPSP_YZR containing the total plastic
        !               strains in direction YZ of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block EPSP_YZR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: PlasticStrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) PlasticStrainValue
              call SetEpsPI(Particles(I), 5, PlasticStrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockEPSP_YZR

        subroutine ReadBlockEPSP_XZR(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block EPSP_XZR containing the total plastic strains 
        !               in direction XZ of particles from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block EPSP_XZR is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: PlasticStrainValue

          if ( (NBytes / 8)==Counters%NParticles) then
          
            do I = 1, Counters%NParticles
              read(FileUnit) PlasticStrainValue
              call SetEpsPI(Particles(I), 6, PlasticStrainValue)
            end do
          
            BlockCheck = BlockCheck + 1
            
          end if
 
        end subroutine ReadBlockEPSP_XZR

        subroutine ReadBlockPIPL___I(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !
        !    Function:  Reads data from the block PIPL___I containing the plasticity state of particles
        !               from the file with unit FileUnit.
        !
        !   FileUnit : File unit from which the block PIPL___I is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I

          if ( (NBytes / 4)==Counters%NParticles) then
          
            read(FileUnit) (Particles(I)%IPL, I = 1, Counters%NParticles)
          
            BlockCheck = BlockCheck + 1
            
          end if
         
        end subroutine ReadBlockPIPL___I

        
        subroutine ReadBlockDIDL___R(FileUnit, NBytes, BlockCheck)
        !**********************************************************************
        !        
        !    Function:  Reads data from the block DIDL___R containing .
        !
        !   FileUnit : File unit from which the block DIDL___R is read
        !   NBytes : Number of bytes inside the block
        !
        ! O BlockCheck : Increase by 1, if the block was successfully read
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE), intent(inout) :: BlockCheck
          ! Local variables
          integer(INTEGER_TYPE) :: I

          if ( (NBytes / 8)==Counters%NParticles) then
          
            read(FileUnit) (Particles(I)%Damping, I = 1, Counters%NParticles)
          
            BlockCheck = BlockCheck + 1
            
          end if
         
        end subroutine ReadBlockDIDL___R
        
            
        logical function ReadPPD()
        !**********************************************************************
        !
        !    Function:  Reads the defined prescribed particle displacements from the
        !               PPD file. If this file exists .true. is returned, else .false. .
        !
        ! O   ReadPPD : Returns .true. if the PPD file exists for the considered project
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE), parameter :: PRESCRPARTICLES = 1
          integer(INTEGER_TYPE), parameter :: PRESCRMATERIALS = 2

          character(len = 255) :: CompleteFileName
          integer(INTEGER_TYPE) :: DefinitionType

          ReadPPD = .false.
        
          CompleteFileName = trim(CalParams%FileNames%ProjectName)//'.PPD'

          if (FExist(CompleteFileName) ) then
            ReadPPD = .true.
            call FileOpen(1, CompleteFileName)
            
            read(1, *) DefinitionType
            select case (DefinitionType)
              case(PRESCRPARTICLES)
                call ReadPrescrDispParticles(1)
              case(PRESCRMATERIALS)
                call ReadPrescrDispMaterial(1)
            end select
            
            close(1)
          end if
        
        end function ReadPPD

        subroutine ReadPrescrDispParticles(FileUnit)
        !**********************************************************************
        !
        !    Function:  Reads the defined prescribed particle displacements from the
        !               PPD file. Prescribed displacements are defined for particles.
        !
        !     FileUnit : Unit of the opened PPD file
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, ParticleID, NPrescribedDisplacementParticles
          real(REAL_TYPE), dimension(NVECTOR) :: PrescrDisp

          read(FileUnit, *) NPrescribedDisplacementParticles
          
          do I = 1, NPrescribedDisplacementParticles
            read(FileUnit, *) ParticleID, PrescrDisp(1:NVECTOR)
            
            if ( (ParticleID>0) .and. (ParticleID<=GetIDCounter() ) ) then ! The particle exists
              do J = 1, NVECTOR
                  call SetPrescrDispI(Particles(I), J, PrescrDisp(J))
              end do
            end if
          end do
            
        end subroutine ReadPrescrDispParticles

        subroutine ReadPrescrDispMaterial(FileUnit)
        !**********************************************************************
        !
        !    Function:  Reads the defined prescribed particle displacements from the
        !               PPD file. Prescribed displacements are defined for particles
        !               of the same material.
        !
        !     FileUnit : Unit of the opened PPD file
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: FileUnit
          ! Local variables
          integer(INTEGER_TYPE) :: I, J, MaterialID
          real(REAL_TYPE), dimension(NVECTOR) :: PrescrDisp
          
          
          read(FileUnit, *) MaterialID, PrescrDisp(1:NVECTOR)
          
          do I = 1, Counters%NParticles
            if (MaterialIDArray(I)==MaterialID) then
                do J = 1, NVECTOR
                    call SetPrescrDispI(Particles(I), J, PrescrDisp(J))
                end do
            end if
          end do
            
        end subroutine ReadPrescrDispMaterial
         
              
        subroutine ReadNodalDataFromFile()
        !**********************************************************************
        !
        !  Function : Fills the nodal data with data read from previous load step 
        !
        !**********************************************************************
        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: I, J, NBytes
          character :: BName*12, Header*80
          character(len = MAX_FILENAME_LENGTH) :: BRFFileName

          if (.not.IsFollowUpPhase()) RETURN 
          
          ! open BRF file for reading
          BRFFileName = trim(CalParams%FileNames%ProjectName)//BRF_FILE_EXTENSION//trim(CalParams%FileNames%PreviousStepExt)
          call FileOpenAction(BRFunit, BRFFileName, 'R')

          ! read total displacement and total velocity from BRF file
          if (FExist(BRFFileName) ) then
            read(BRFunit) Header
            do
              read(BRFunit) BName, NBytes
              if (trim(BName)=='$$ENDOFFII$$') then ! close file
                EXIT
              elseif (NBytes>0) then ! The block is not empty
                if (trim(BName)=='$$UTOT___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) TotalDisplacementSoil(I) 
                  End Do
                else if (trim(BName)=='$$UWTOT__R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) TotalDisplacementWater(I) 
                  End Do
                else if (trim(BName)=='$$UGTOT__R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) TotalDisplacementGas(I) 
                  End Do
                else if (trim(BName)=='$$VTOT___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) TotalVelocitySoil(I, 1) 
                  End Do
                else if (trim(BName)=='$$WTOT___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) TotalVelocityWater(I, 1) ! Global coordinate system
                  End Do
                else if (trim(BName)=='$$GTOT___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) TotalVelocityGas(I, 1) 
                  End Do
                else if (trim(BName)=='$$DSPTSD_R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) DashpotSld(I) 
                  End Do
                else if (trim(BName)=='$$SPSD___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) SpringSld(I) 
                  End Do
                else if (trim(BName)=='$$UTOTVBSR$$') then
                  do I = 1, Counters%N
                    read(BRFunit) NodalTotDisplacementVB(I,1) 
                  End Do
                else if (trim(BName)=='$$DSPTWT_R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) DashpotWat(I)
                  End Do
                else if (trim(BName)=='$$SPWT___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) SpringWat(I)
                  End Do
                else if (trim(BName)=='$$WTOTVBWR$$') then
                  do I = 1, Counters%N
                    read(BRFunit) NodalTotDisplacementVBWater(I,1) 
                  End Do
                else if (trim(BName)=='$$DSPTGT_R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) DashpotGas(I)
                  End Do
                else if (trim(BName)=='$$SPGT___R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) SpringGas(I)
                  End Do
                else if (trim(BName)=='$$WTOTVBGR$$') then
                  do I = 1, Counters%N
                    read(BRFunit) NodalTotDisplacementVBGas(I,1) 
                  End Do
                else if (trim(BName)=='$$FVISSLDR$$') then
                  do I = 1, Counters%N
                    do J = 1, Counters%NEntity 
                      read(BRFunit) VisDampForceSld(I,J)
                    end do  
                  End Do
                else if (trim(BName)=='$$FVISWATR$$') then
                  do I = 1, Counters%N
                    do J = 1, Counters%NEntity 
                      read(BRFunit) VisDampForceWat(I,J) 
                    end do  
                    end do
                else if (trim(BName)=='$$FVISGASR$$') then
                  do I = 1, Counters%N
                    do J = 1, Counters%NEntity 
                      read(BRFunit) VisDampForceGas(I,J) 
                    end do  
                End Do
                else if (trim(BName)=='$$EXLDWT_R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) ExtLoadWaterTotal(I, 1,1)
                  End Do
                else if (trim(BName)=='$$EXLDGS_R$$') then
                  do I = 1, Counters%N
                    read(BRFunit) ExtLoadGasTotal(I, 1,1)
                  End Do
                else ! Skip bytes of block
                  call Skip(BRFunit, NBytes)
                end if
              end if
            end do
          end if
         
          close(BRFunit)
          
        end subroutine ReadNodalDataFromFile
               


      end module ModReadMPMData
