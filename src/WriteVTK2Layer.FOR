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


module ModWriteVTKTwoLayer
    !**********************************************************************
    !
    !    Function: Contains pointer routines for writing output data in VTK format
    !    for 2 layer formulation
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************
    
    use ModGlobalConstants
    use ModMPMData
    use ModWriteVTKOutput
    use ModFileIO
    
    contains
    
    subroutine WriteVTKOutput_2LayForm_Solid()
    !**********************************************************************
    !
    ! Layer formulation :: Solid Material Point info
    ! Function:  Write the VTK output file
    !
    !**********************************************************************

    implicit none

    ! local variables
    integer(INTEGER_TYPE) :: ElementID, SolidMaterialPointID, MaterialPointID, DataSetID,  &
        NumberElements, NumberMaterialPoints, DataSetSize, NumberSolidMaterialPoints, k
    character(len=50), dimension(:), allocatable :: DataSetName
    real(REAL_TYPE), dimension(:,:), allocatable :: DataSetScalar
    real(REAL_TYPE), dimension(:,:,:), allocatable :: DataSetVector
    real(REAL_TYPE), dimension(:,:,:,:), allocatable :: DataSetTensor
    real(REAL_TYPE):: EpsD
    real(REAL_TYPE), dimension(NVECTOR) :: VectorValue

    !********** START: writing scalar material point data ****************************************
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    DataSetSize = 16 ! total number of datasets for material points to be written
    allocate(DataSetName(DataSetSize))

    NumberSolidMaterialPoints = 0

    !Count Number of Solid MatPoint, allocate and initialize
    do MaterialPointID = 1, NumberMaterialPoints
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            NumberSolidMaterialPoints =  NumberSolidMaterialPoints + 1 ! total number of solid material points
        end if
    end do
    allocate(DataSetScalar(NumberSolidMaterialPoints,DataSetSize))
    DataSetScalar = 0.0

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'mass_solid' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = MassArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 2 ! index of dataset
    DataSetName(DataSetID) = 'weight_solid' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%MaterialWeight ! assign output data to material points
        end if
    end do

    DataSetID = 3 ! index of dataset
    DataSetName(DataSetID) = 'initial_porosity' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%InitialPorosity ! assign output data to material points
        end if
    end do

    DataSetID = 4 ! index of dataset
    DataSetName(DataSetID) = 'integration_weight' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%IntegrationWeight ! assign output data to material points
        end if
    end do

    DataSetID = 5 ! index of dataset
    DataSetName(DataSetID) = 'material_point_id' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) =IDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 6 ! index of dataset
    DataSetName(DataSetID) = 'element_id' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = ElementIDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 7 ! index of dataset
    DataSetName(DataSetID) = 'entity_id' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = EntityIDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 8 ! index of dataset
    DataSetName(DataSetID) = 'material_id' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = MaterialIDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 9 ! index of dataset
    DataSetName(DataSetID) = 'density' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%Density
        end if
    end do

    DataSetID = 10 ! index of dataset
    DataSetName(DataSetID) = 'damping' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%damping ! assign output data to material points
        end if
    end do

    DataSetID = 11 ! index of dataset
    DataSetName(DataSetID) = 'mean_eff_stress' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = (SigmaEffArray(MaterialPointID,1) +  &
                SigmaEffArray(MaterialPointID,2) + SigmaEffArray(MaterialPointID,3) ) / 3.0 ! valid for 2D and 3D
        end if
    end do

    DataSetID = 12 ! index of dataset
    DataSetName(DataSetID) = 'Phase_Status' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            if(Particles(MaterialPointID)%PhaseStatus==PhaseStatusSOLID) then
                DataSetScalar(SolidMaterialPointID,DataSetID) = 1.0
            else if(Particles(MaterialPointID)%PhaseStatus==PhaseStatusLIQUID) then
                DataSetScalar(SolidMaterialPointID,DataSetID) = 2.0
            end if
        end if
    end do

    DataSetID = 13 ! index of dataset
    DataSetName(DataSetID) = 'porosity_Eff' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%EffPorosity ! assign output data to material points
        end if
    end do

    DataSetID = 14 ! index of dataset
    DataSetName(DataSetID) = 'ConcRatioLiquid_SolElm' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%ConcentrationRatioLiquidS ! assign output data to material points
        end if
    end do

    DataSetID = 15 ! index of dataset
    DataSetName(DataSetID) = 'ConcRatioSolid_Eff' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetScalar(SolidMaterialPointID,DataSetID) = Particles(MaterialPointID)%EffConcentrationRatioSolid ! assign output data to material points
        end if
    end do
    
    DataSetID = 16 ! index of dataset
    DataSetName(DataSetID) = 'deviatoric_strain_solid' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            EpsD = getEpsq( GetEps(Particles(MaterialPointID)) ) ! deviatoric strain
            DataSetScalar(SolidMaterialPointID,DataSetID) = EpsD ! assign output data to material points
        end if  
    end do

    call WriteVTKMaterialPointScalar_2LayForm_Solid(DataSetName,DataSetScalar) ! write all scalar datasets for material points

    deallocate(DataSetName,DataSetScalar)
    !********** END: writing scalar material point data ****************************************

    !********** START: writing vector material point data ****************************************
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    DataSetSize = 7 ! total number of datasets for material points to be written
    allocate(DataSetName(DataSetSize))

    NumberSolidMaterialPoints = 0

    !Count Number of Solid MatPoint, allocate and initialize
    do MaterialPointID = 1, NumberMaterialPoints
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            NumberSolidMaterialPoints =  NumberSolidMaterialPoints + 1 ! total number of solid material points
        end if
    end do
    allocate(DataSetVector(NumberSolidMaterialPoints,NVECTOR,DataSetSize))
    DataSetVector = 0.0

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'acceleration' ! assign name to dataset
    SolidMaterialPointID = 0
    DataSetVector(:,:,DataSetID) = 0.0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetVector(SolidMaterialPointID,:,DataSetID) = AccelerationArray(MaterialPointID,:)
        end if
    end do

    DataSetID = 2 ! index of dataset
    DataSetName(DataSetID) = 'velocity_solid' ! assign name to dataset
    SolidMaterialPointID = 0
    DataSetVector(:,:,DataSetID) = 0.0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetVector(SolidMaterialPointID,:,DataSetID) = VelocityArray(MaterialPointID,:)
        end if
    end do

    DataSetID = 3 ! index of dataset
    DataSetName(DataSetID) = 'displacement_solid' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetVector(SolidMaterialPointID,:,DataSetID) = UArray(MaterialPointID,:)
        end if
    end do

    DataSetID = 4 ! index of dataset
    DataSetName(DataSetID) = 'global_position' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetVector(SolidMaterialPointID,:,DataSetID) = GlobPosArray(MaterialPointID,:)
        end if
    end do

    DataSetID = 5 ! index of dataset
    DataSetName(DataSetID) = 'local_position' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetVector(SolidMaterialPointID,:,DataSetID) = Particles(MaterialPointID)%LocPos(:)
        end if
    end do

    DataSetID = 6 ! index of dataset
    DataSetName(DataSetID) = 'bodyforce' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            DataSetVector(SolidMaterialPointID,:,DataSetID) = Particles(MaterialPointID)%FBody(:)
        end if
    end do

    DataSetID = 7 ! index of dataset
    DataSetName(DataSetID) = 'externalforce' ! assign name to dataset
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            VectorValue = 0.0
            do k = 1,Counters%NSoilLoadSystems
             VectorValue(:) = VectorValue (:) + Particles(MaterialPointID)%FExt(:,k)*CalParams%Multipliers%SolidACurrent(k)
            end do
            DataSetVector(SolidMaterialPointID,:,DataSetID) = VectorValue(:)
        end if
    end do

    call WriteVTKMaterialPointVector_2LayForm_Solid(DataSetName,DataSetVector) ! write all vector datasets for material points

    deallocate(DataSetName,DataSetVector)
    !********** END: writing vector material point data ****************************************

    !********** START: writing tensor material point data ****************************************
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    DataSetSize = 2 ! total number of datasets for material points to be written
    allocate(DataSetName(DataSetSize))

    NumberSolidMaterialPoints = 0

    !Count Number of Solid MatPoint, allocate and initialize
    do MaterialPointID = 1, NumberMaterialPoints
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            NumberSolidMaterialPoints =  NumberSolidMaterialPoints + 1 ! total number of solid material points
        end if
    end do
    allocate(DataSetTensor(NumberSolidMaterialPoints,NPRINCIPAL,NPRINCIPAL,DataSetSize))
    DataSetTensor = 0.0

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'effective_stress' ! assign name to dataset
    DataSetTensor(:,:,:,DataSetID) = 0.0
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            if (NDIM == 3) then ! 3D case
                DataSetTensor(SolidMaterialPointID,1,1,DataSetID) = SigmaEffArray(MaterialPointID,1)
                DataSetTensor(SolidMaterialPointID,1,2,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(SolidMaterialPointID,1,3,DataSetID) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(SolidMaterialPointID,2,1,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(SolidMaterialPointID,2,2,DataSetID) = SigmaEffArray(MaterialPointID,2)
                DataSetTensor(SolidMaterialPointID,2,3,DataSetID) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(SolidMaterialPointID,3,1,DataSetID) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(SolidMaterialPointID,3,2,DataSetID) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(SolidMaterialPointID,3,3,DataSetID) = SigmaEffArray(MaterialPointID,3)
            elseif (NDIM == 2) then ! 2D case
                DataSetTensor(SolidMaterialPointID,1,1,DataSetID) = SigmaEffArray(MaterialPointID,1)
                DataSetTensor(SolidMaterialPointID,1,2,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(SolidMaterialPointID,2,1,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(SolidMaterialPointID,2,2,DataSetID) = SigmaEffArray(MaterialPointID,2)
                DataSetTensor(SolidMaterialPointID,3,3,DataSetID) = SigmaEffArray(MaterialPointID,3)
            end if
        end if
    end do

    DataSetID = 2 ! index of dataset
    DataSetName(DataSetID) = 'strain' ! assign name to dataset
    DataSetTensor(:,:,:,DataSetID) = 0.0
    SolidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            SolidMaterialPointID = SolidMaterialPointID + 1
            if (NDIM == 3) then ! 3D case
                DataSetTensor(SolidMaterialPointID,1,1,DataSetID) = Particles(MaterialPointID)%Eps(1)
                DataSetTensor(SolidMaterialPointID,1,2,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(SolidMaterialPointID,1,3,DataSetID) = Particles(MaterialPointID)%Eps(6)/2.0
                DataSetTensor(SolidMaterialPointID,2,1,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(SolidMaterialPointID,2,2,DataSetID) = Particles(MaterialPointID)%Eps(2)
                DataSetTensor(SolidMaterialPointID,2,3,DataSetID) = Particles(MaterialPointID)%Eps(5)/2.0
                DataSetTensor(SolidMaterialPointID,3,1,DataSetID) = Particles(MaterialPointID)%Eps(6)/2.0
                DataSetTensor(SolidMaterialPointID,3,2,DataSetID) = Particles(MaterialPointID)%Eps(5)/2.0
                DataSetTensor(SolidMaterialPointID,3,3,DataSetID) = Particles(MaterialPointID)%Eps(3)
            elseif (NDIM == 2) then ! 2D case
                DataSetTensor(SolidMaterialPointID,1,1,DataSetID) = Particles(MaterialPointID)%Eps(1)
                DataSetTensor(SolidMaterialPointID,1,2,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(SolidMaterialPointID,2,1,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(SolidMaterialPointID,2,2,DataSetID) = Particles(MaterialPointID)%Eps(2)
                DataSetTensor(SolidMaterialPointID,3,3,DataSetID) = Particles(MaterialPointID)%Eps(3)
            end if
        end if
    end do

    call WriteVTKMaterialPointTensor_2LayForm_Solid(DataSetName,DataSetTensor) ! write all tensor datasets for material points

    deallocate(DataSetName,DataSetTensor)
    !********** END: writing tensor material point data ****************************************

    !********** START: writing scalar mesh data ****************************************
    NumberElements = Counters%NEl ! total number of elements
    DataSetSize = 1 ! total number of datasets for elements to be written
    allocate(DataSetName(DataSetSize))
    allocate(DataSetScalar(NumberElements,DataSetSize))

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'active_elements' ! assign name to dataset
    do ElementID = 1, NumberElements ! loop over elements
        if (IsActiveElement(ElementID)) then
            DataSetScalar(ElementID,DataSetID) = 1.0 ! assign output data to elements
        else
            DataSetScalar(ElementID,DataSetID) = 0.0
        end if
    end do

    call WriteVTKMeshData(DataSetName,DataSetScalar) ! write all datasets for elements

    deallocate(DataSetName,DataSetScalar)
    !********** END: writing scalar mesh data ****************************************

    end subroutine WriteVTKOutput_2LayForm_Solid
    
    subroutine WriteVTKMaterialPointScalar_2LayForm_Solid(DataSetName,DataSetScalar)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************
    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:), intent(in) :: DataSetScalar ! dimension(NumberMaterialPoints,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID, DataSetID, NumberSolidMaterialPoints, NumberMaterialPoints, DataSetSize
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    NumberSolidMaterialPoints = size(DataSetScalar,1) ! total number of SOLID material points
    DataSetSize = size(DataSetScalar,2) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPScalarSOLID'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberSolidMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            NoCo(1) = GlobPosArray(MaterialPointID, 1)
            NoCo(2) = GlobPosArray(MaterialPointID, 2)
            if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,*) NoCo
        end if
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberSolidMaterialPoints,2*NumberSolidMaterialPoints
    do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberSolidMaterialPoints
    do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberSolidMaterialPoints
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'SCALARS ',trim(DataSetName(DataSetID)),' float 1' ! write dataset name
        write(VTKUnit,'(A)')'LOOKUP_TABLE default'
        do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
            write(VTKUnit, '(f14.8)')DataSetScalar(MaterialPointID,DataSetID) ! write data of material point
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMaterialPointScalar_2LayForm_Solid


    subroutine WriteVTKMaterialPointVector_2LayForm_Solid(DataSetName,DataSetVector)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************
    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:,:), intent(in) :: DataSetVector ! dimension(NumberMaterialPoints,3,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID, DataSetID, NumberMaterialPoints, DataSetSize, NumberSolidMaterialPoints
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    NumberSolidMaterialPoints = size(DataSetVector,1) ! total number of SOLID material points
    DataSetSize = size(DataSetVector,3) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPVectorSOLID'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberSolidMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            NoCo(1) = GlobPosArray(MaterialPointID, 1)
            NoCo(2) = GlobPosArray(MaterialPointID, 2)
            if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,*) NoCo
        end if
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberSolidMaterialPoints,2*NumberSolidMaterialPoints
    do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberSolidMaterialPoints
    do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberSolidMaterialPoints
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'VECTORS ',trim(DataSetName(DataSetID)),' float' ! write dataset name
        do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
            NoCo(1) = DataSetVector(MaterialPointID, 1, DataSetID)
            NoCo(2) = DataSetVector(MaterialPointID, 2, DataSetID)
            if (NDIM==3) NoCo(3) = DataSetVector(MaterialPointID, 3, DataSetID) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,'(f14.8)') NoCo
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMaterialPointVector_2LayForm_Solid


    subroutine WriteVTKMaterialPointTensor_2LayForm_Solid(DataSetName,DataSetTensor)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************
    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:,:,:), intent(in) :: DataSetTensor ! dimension(NumberMaterialPoints,NVECTOR,NVECTOR,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: I, MaterialPointID, DataSetID, NumberMaterialPoints, DataSetSize, NumberSolidMaterialPoints
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    real(REAL_TYPE), dimension(9) :: Tensor = 0.0 ! dimension 3x3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    NumberSolidMaterialPoints = size(DataSetTensor,1) ! total number of material points
    DataSetSize = size(DataSetTensor,4) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPTensorSOLID'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

   
    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberSolidMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid) then
            NoCo(1) = GlobPosArray(MaterialPointID, 1)
            NoCo(2) = GlobPosArray(MaterialPointID, 2)
            if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,*) NoCo
        end if
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberSolidMaterialPoints,2*NumberSolidMaterialPoints
    do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberSolidMaterialPoints
    do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberSolidMaterialPoints
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'TENSORS ',trim(DataSetName(DataSetID)),' float' ! write dataset name
        do MaterialPointID = 1, NumberSolidMaterialPoints ! loop over material points
                      
          Tensor(1) = DataSetTensor(MaterialPointID,1,1,DataSetID)
          Tensor(2) = DataSetTensor(MaterialPointID,1,2,DataSetID)
          Tensor(4) = DataSetTensor(MaterialPointID,2,1,DataSetID)
          Tensor(5) = DataSetTensor(MaterialPointID,2,2,DataSetID)
          Tensor(9) = DataSetTensor(MaterialPointID,3,3,DataSetID)

          if (NDIM==3) then
            Tensor(3) = DataSetTensor(MaterialPointID,1,3,DataSetID)
            Tensor(6) = DataSetTensor(MaterialPointID,2,3,DataSetID)
            Tensor(7) = DataSetTensor(MaterialPointID,3,1,DataSetID)
            Tensor(8) = DataSetTensor(MaterialPointID,3,2,DataSetID)
          end if
          
          do I = 1, 3
            write(VTKUnit,'(f14.8)') Tensor(I + 2 * (I - 1):I + 2 * I)
          end do
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMaterialPointTensor_2LayForm_Solid


    subroutine WriteVTKOutput_2LayForm_Liquid()
    !**********************************************************************
    !
    ! Layer formulation :: Liquid Material Point info
    ! Function:  Write the VTK output file
    !
    !**********************************************************************
    implicit none

    ! local variables
    integer(INTEGER_TYPE) :: LiquidMaterialPointID, MaterialPointID, DataSetID,  &
        NumberMaterialPoints, DataSetSize, NumberLiquidMaterialPoints, k
    character(len=50), dimension(:), allocatable :: DataSetName
    real(REAL_TYPE), dimension(:,:), allocatable :: DataSetScalar
    real(REAL_TYPE), dimension(:,:,:), allocatable :: DataSetVector
    real(REAL_TYPE), dimension(:,:,:,:), allocatable :: DataSetTensor
    real(REAL_TYPE), dimension(NVECTOR) :: VectorValue

    !********** START: writing scalar material point data ****************************************
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    DataSetSize = 19 ! total number of datasets for material points to be written
    allocate(DataSetName(DataSetSize))

    NumberLiquidMaterialPoints = 0

    !Count Number of Liquid MatPoint, allocate and initialize
    do MaterialPointID = 1, NumberMaterialPoints
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            NumberLiquidMaterialPoints =  NumberLiquidMaterialPoints + 1 ! total number of Liquid material points
        end if
    end do
    allocate(DataSetScalar(NumberLiquidMaterialPoints,DataSetSize))
    DataSetScalar = 0.0

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'mass_Liquid' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = MassArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 2 ! index of dataset
    DataSetName(DataSetID) = 'weight_Liquid' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%MaterialWeight ! assign output data to material points
        end if
    end do

    DataSetID = 3 ! index of dataset
    DataSetName(DataSetID) = 'porosity_LiqElm' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%PorosityL ! assign output data to material points
        end if
    end do

    DataSetID = 4 ! index of dataset
    DataSetName(DataSetID) = 'initial_porosity' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%InitialPorosity ! assign output data to material points
        end if
    end do

    DataSetID = 5 ! index of dataset
    DataSetName(DataSetID) = 'integration_weight' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%IntegrationWeight ! assign output data to material points
        end if
    end do

    DataSetID = 6 ! index of dataset
    DataSetName(DataSetID) = 'material_point_id' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = IDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 7 ! index of dataset
    DataSetName(DataSetID) = 'element_id' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = ElementIDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 8 ! index of dataset
    DataSetName(DataSetID) = 'entity_id' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = EntityIDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 9 ! index of dataset
    DataSetName(DataSetID) = 'material_id' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = MaterialIDArray(MaterialPointID) ! assign output data to material points
        end if
    end do

    DataSetID = 10 ! index of dataset
    DataSetName(DataSetID) = 'density' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%Density ! assign output data to material points
        end if
    end do

    DataSetID = 11 ! index of dataset
    DataSetName(DataSetID) = 'damping' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%damping ! assign output data to material points
        end if
    end do

    DataSetID = 12 ! index of dataset
    DataSetName(DataSetID) = 'ConcRatioLiquid_LiqElm' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%ConcentrationRatioLiquidL ! assign output data to material points
        end if
    end do

    DataSetID = 13 ! index of dataset
    DataSetName(DataSetID) = 'ConcRatioLiquid_Eff' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%EffConcentrationRatioLiquid ! assign output data to material points
        end if
    end do

    DataSetID = 14 ! index of dataset
    DataSetName(DataSetID) = 'ConcRatioSolid_LiqElm' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%ConcentrationRatioSolidL ! assign output data to material points
        end if
    end do

    DataSetID = 15 ! index of dataset
    DataSetName(DataSetID) = 'liquid_pres' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = (SigmaEffArray(MaterialPointID,1) +  &
                SigmaEffArray(MaterialPointID,2) + SigmaEffArray(MaterialPointID,3) ) / 3.0 ! valid for 2D and 3D
        end if
    end do

    DataSetID = 16 ! index of dataset
    DataSetName(DataSetID) = 'Phase_Status' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            if(Particles(MaterialPointID)%PhaseStatus==PhaseStatusSOLID) then
                DataSetScalar(LiquidMaterialPointID,DataSetID) = 1.0
            else if(Particles(MaterialPointID)%PhaseStatus==PhaseStatusLIQUID) then
                DataSetScalar(LiquidMaterialPointID,DataSetID) = 2.0
            end if
        end if
    end do

    DataSetID = 17 ! index of dataset
    DataSetName(DataSetID) = 'Liquid_free_surface' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%LiquidFreeSurface ! assign output data to material points
        end if
    end do

    DataSetID = 18 ! index of dataset
    DataSetName(DataSetID) = 'Volumetric_Strain_Liquid' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%WaterVolumetricStrain ! assign output data to material points
        end if
    end do

    DataSetID = 19 ! index of dataset
    DataSetName(DataSetID) = 'FillingRatio' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetScalar(LiquidMaterialPointID,DataSetID) = Particles(MaterialPointID)%FillingRatio ! assign output data to material points
        end if
    end do

    ! write all scalar datasets for material points
    call WriteVTKMaterialPointScalar_2LayForm_Liquid(DataSetName,DataSetScalar)

    deallocate(DataSetName,DataSetScalar)
    !********** END: writing scalar material point data ****************************************

    !********** START: writing vector material point data ****************************************
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    DataSetSize = 7 ! total number of datasets for material points to be written
    allocate(DataSetName(DataSetSize))

    NumberLiquidMaterialPoints = 0

    !Count Number of Liquid MatPoint, allocate and initialize
    do MaterialPointID = 1, NumberMaterialPoints
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            NumberLiquidMaterialPoints =  NumberLiquidMaterialPoints + 1 ! total number of Liquid material points
        end if
    end do
    allocate(DataSetVector(NumberLiquidMaterialPoints,NVECTOR,DataSetSize))
    DataSetVector = 0.0

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'acceleration' ! assign name to dataset
    LiquidMaterialPointID = 0
    DataSetVector(:,:,DataSetID) = 0.0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = AccelerationArray(MaterialPointID,:) ! assign output data
        end if
    end do

    DataSetID = 2 ! index of dataset
    DataSetName(DataSetID) = 'velocity_Liquid' ! assign name to dataset
    LiquidMaterialPointID = 0
    DataSetVector(:,:,DataSetID) = 0.0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = VelocityArray(MaterialPointID,:) ! assign output data
        end if
    end do

    DataSetID = 3 ! index of dataset
    DataSetName(DataSetID) = 'displacement_Liquid' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = UArray(MaterialPointID,:) ! assign output data
        end if
    end do

    DataSetID = 4 ! index of dataset
    DataSetName(DataSetID) = 'global_position' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = GlobPosArray(MaterialPointID,:) ! assign output data
        end if
    end do

    DataSetID = 5 ! index of dataset
    DataSetName(DataSetID) = 'locbal_position' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = Particles(MaterialPointID)%LocPos(:)
        end if
    end do

    DataSetID = 6 ! index of dataset
    DataSetName(DataSetID) = 'bodyforce' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = Particles(MaterialPointID)%FBody(:)
        end if
    end do

    DataSetID = 7 ! index of dataset
    DataSetName(DataSetID) = 'externalforce' ! assign name to dataset
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            VectorValue = 0.0
            do k = 1,Counters%NSoilLoadSystems
             VectorValue = VectorValue + Particles(MaterialPointID)%FExt(:,k)*CalParams%Multipliers%SolidACurrent(k)
            end do            
            DataSetVector(LiquidMaterialPointID,:,DataSetID) = VectorValue
        end if
    end do

    ! write all vector datasets for material points
    call WriteVTKMaterialPointVector_2LayForm_Liquid(DataSetName,DataSetVector)

    deallocate(DataSetName,DataSetVector)
    !********** END: writing vector material point data ****************************************

    !********** START: writing tensor material point data ****************************************
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    DataSetSize = 2 ! total number of datasets for material points to be written
    allocate(DataSetName(DataSetSize))

    NumberLiquidMaterialPoints = 0

    !Count Number of Liquid MatPoint, allocate and initialize
    do MaterialPointID = 1, NumberMaterialPoints
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            NumberLiquidMaterialPoints =  NumberLiquidMaterialPoints + 1 ! total number of Liquid material points
        end if
    end do
    allocate(DataSetTensor(NumberLiquidMaterialPoints,NPRINCIPAL,NPRINCIPAL,DataSetSize))
    DataSetTensor = 0.0

    DataSetID = 1 ! index of dataset
    DataSetName(DataSetID) = 'liquid_stress' ! assign name to dataset
    DataSetTensor(:,:,:,DataSetID) = 0.0
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            if (NDIM == 3) then ! 3D case
                DataSetTensor(LiquidMaterialPointID,1,1,DataSetID) = SigmaEffArray(MaterialPointID,1)
                DataSetTensor(LiquidMaterialPointID,1,2,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(LiquidMaterialPointID,1,3,DataSetID) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(LiquidMaterialPointID,2,1,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(LiquidMaterialPointID,2,2,DataSetID) = SigmaEffArray(MaterialPointID,2)
                DataSetTensor(LiquidMaterialPointID,2,3,DataSetID) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(LiquidMaterialPointID,3,1,DataSetID) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(LiquidMaterialPointID,3,2,DataSetID) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(LiquidMaterialPointID,3,3,DataSetID) = SigmaEffArray(MaterialPointID,3)
            elseif (NDIM == 2) then ! 2D case
                DataSetTensor(LiquidMaterialPointID,1,1,DataSetID) = SigmaEffArray(MaterialPointID,1)
                DataSetTensor(LiquidMaterialPointID,1,2,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(LiquidMaterialPointID,2,1,DataSetID) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(LiquidMaterialPointID,2,2,DataSetID) = SigmaEffArray(MaterialPointID,2)
                DataSetTensor(LiquidMaterialPointID,3,3,DataSetID) = SigmaEffArray(MaterialPointID,3)
            end if
        end if
    end do

    DataSetID = 2 ! index of dataset
    DataSetName(DataSetID) = 'strain' ! assign name to dataset
    DataSetTensor(:,:,:,DataSetID) = 0.0
    LiquidMaterialPointID = 0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            LiquidMaterialPointID = LiquidMaterialPointID + 1
            if (NDIM == 3) then ! 3D case
                DataSetTensor(LiquidMaterialPointID,1,1,DataSetID) = Particles(MaterialPointID)%Eps(1)
                DataSetTensor(LiquidMaterialPointID,1,2,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(LiquidMaterialPointID,1,3,DataSetID) = Particles(MaterialPointID)%Eps(6)/2.0
                DataSetTensor(LiquidMaterialPointID,2,1,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(LiquidMaterialPointID,2,2,DataSetID) = Particles(MaterialPointID)%Eps(2)
                DataSetTensor(LiquidMaterialPointID,2,3,DataSetID) = Particles(MaterialPointID)%Eps(5)/2.0
                DataSetTensor(LiquidMaterialPointID,3,1,DataSetID) = Particles(MaterialPointID)%Eps(6)/2.0
                DataSetTensor(LiquidMaterialPointID,3,2,DataSetID) = Particles(MaterialPointID)%Eps(5)/2.0
                DataSetTensor(LiquidMaterialPointID,3,3,DataSetID) = Particles(MaterialPointID)%Eps(3)
            elseif (NDIM == 2) then ! 2D case
                DataSetTensor(LiquidMaterialPointID,1,1,DataSetID) = Particles(MaterialPointID)%Eps(1)
                DataSetTensor(LiquidMaterialPointID,1,2,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(LiquidMaterialPointID,2,1,DataSetID) = Particles(MaterialPointID)%Eps(4)/2.0
                DataSetTensor(LiquidMaterialPointID,2,2,DataSetID) = Particles(MaterialPointID)%Eps(2)
                DataSetTensor(LiquidMaterialPointID,3,3,DataSetID) = Particles(MaterialPointID)%Eps(3)
            end if
        end if
    end do

    ! write all tensor datasets for material points
    call WriteVTKMaterialPointTensor_2LayForm_Liquid(DataSetName,DataSetTensor)

    deallocate(DataSetName,DataSetTensor)
    !********** END: writing tensor material point data ****************************************

    end subroutine WriteVTKOutput_2LayForm_Liquid


    subroutine WriteVTKMaterialPointScalar_2LayForm_Liquid(DataSetName,DataSetScalar)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************
    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:), intent(in) :: DataSetScalar ! dimension(NumberMaterialPoints,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID, DataSetID, NumberLiquidMaterialPoints, NumberMaterialPoints, DataSetSize
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    NumberLiquidMaterialPoints = size(DataSetScalar,1) ! total number of Liquid material points
    DataSetSize = size(DataSetScalar,2) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPScalarLIQUID'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberLiquidMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        if (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            NoCo(1) = GlobPosArray(MaterialPointID, 1)
            NoCo(2) = GlobPosArray(MaterialPointID, 2)
            if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,*) NoCo
        end if
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberLiquidMaterialPoints,2*NumberLiquidMaterialPoints
    do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberLiquidMaterialPoints
    do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberLiquidMaterialPoints
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'SCALARS ',trim(DataSetName(DataSetID)),' float 1' ! write dataset name
        write(VTKUnit,'(A)')'LOOKUP_TABLE default'
        do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
            write(VTKUnit,'(f14.8)')DataSetScalar(MaterialPointID,DataSetID) ! write data of material point
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMaterialPointScalar_2LayForm_Liquid


    subroutine WriteVTKMaterialPointVector_2LayForm_Liquid(DataSetName,DataSetVector)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************
    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:,:), intent(in) :: DataSetVector ! dimension(NumberMaterialPoints,NVECTOR,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID, DataSetID, NumberMaterialPoints, DataSetSize, NumberLiquidMaterialPoints
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    NumberLiquidMaterialPoints = size(DataSetVector,1) ! total number of Liquid material points
    DataSetSize = size(DataSetVector,3) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPVectorLIQUID'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberLiquidMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            NoCo(1) = GlobPosArray(MaterialPointID, 1)
            NoCo(2) = GlobPosArray(MaterialPointID, 2)
            if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,*) NoCo
        end if
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberLiquidMaterialPoints,2*NumberLiquidMaterialPoints
    do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberLiquidMaterialPoints
    do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberLiquidMaterialPoints
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'VECTORS ',trim(DataSetName(DataSetID)),' float' ! write dataset name
        do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
            NoCo(1) = DataSetVector(MaterialPointID, 1, DataSetID)
            NoCo(2) = DataSetVector(MaterialPointID, 2, DataSetID)
            if (NDIM==3) NoCo(3) = DataSetVector(MaterialPointID, 3, DataSetID) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,'(f14.8)') NoCo
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMaterialPointVector_2LayForm_Liquid


    subroutine WriteVTKMaterialPointTensor_2LayForm_Liquid(DataSetName,DataSetTensor)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************
    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:,:,:), intent(in) :: DataSetTensor ! dimension(NumberMaterialPoints,NVECTOR,NVECTOR,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: I, MaterialPointID, DataSetID, NumberMaterialPoints, DataSetSize, NumberLiquidMaterialPoints
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    real(REAL_TYPE), dimension(9) :: Tensor = 0.0 ! dimension 3x3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    NumberLiquidMaterialPoints = size(DataSetTensor,1) ! total number of material points
    DataSetSize = size(DataSetTensor,4) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPTensorLIQUID'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberLiquidMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        if(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid) then
            NoCo(1) = GlobPosArray(MaterialPointID, 1)
            NoCo(2) = GlobPosArray(MaterialPointID, 2)
            if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
            write(VTKUnit,*) NoCo
        end if
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberLiquidMaterialPoints,2*NumberLiquidMaterialPoints
    do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberLiquidMaterialPoints
    do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberLiquidMaterialPoints
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'TENSORS ',trim(DataSetName(DataSetID)),' float' ! write dataset name

        do MaterialPointID = 1, NumberLiquidMaterialPoints ! loop over material points
                      
          Tensor(1) = DataSetTensor(MaterialPointID,1,1,DataSetID)
          Tensor(2) = DataSetTensor(MaterialPointID,1,2,DataSetID)
          Tensor(4) = DataSetTensor(MaterialPointID,2,1,DataSetID)
          Tensor(5) = DataSetTensor(MaterialPointID,2,2,DataSetID)
          Tensor(9) = DataSetTensor(MaterialPointID,3,3,DataSetID)

          if (NDIM==3) then
            Tensor(3) = DataSetTensor(MaterialPointID,1,3,DataSetID)
            Tensor(6) = DataSetTensor(MaterialPointID,2,3,DataSetID)
            Tensor(7) = DataSetTensor(MaterialPointID,3,1,DataSetID)
            Tensor(8) = DataSetTensor(MaterialPointID,3,2,DataSetID)
          end if
          
          do I = 1, 3
            write(VTKUnit,'(f14.8)') Tensor(I + 2 * (I - 1):I + 2 * I)
          end do
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMaterialPointTensor_2LayForm_Liquid

    end module ModWriteVTKTwoLayer
