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


module ModWriteVTKOutput
    !**********************************************************************
    !
    ! Function: Contains routines for writing output data in VTK format
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************

    use ModReadCalculationData
    use ModReadMaterialData
    use ModGlobalConstants
    use ModCounters
    use ModMeshInfo
    use ModParticle
    use ModMPMData
    use ModMeshInfo
    use ModWriteVTKASCII
    use ModWriteVTKBinary
    use ModFileIO

    implicit none

    contains
       
    subroutine SetVTKPointers()
    !**********************************************************************
    !
    !    Function: set VTK pointers
    !
    !**********************************************************************
        implicit none
        
        if (CalParams%isVTKBinary) then ! Set vtk writing pointers    
            InitialiseVTKFilePointer => InitialiseVTKFileBinary
            WriteVTKFloatScalarDataPointer => WriteVTKFloatScalarDataBinary
            WriteVTKIntegerScalarDataPointer => WriteVTKIntegerScalarDataBinary
            WriteVTKFloatVectorDataPointer => WriteVTKFloatVectorDataBinary
            WriteVTKIntegerVectorDataPointer => WriteVTKIntegerVectorDataBinary
            WriteVTKFloatTensorDataPointer => WriteVTKFloatTensorDataBinary
            WriteVTKIntegerTensorDataPointer => WriteVTKIntegerTensorDataBinary         
        else ! Set vtk writing pointers  
            InitialiseVTKFilePointer => InitialiseVTKFileASCII
            WriteVTKFloatScalarDataPointer => WriteVTKFloatScalarDataASCII
            WriteVTKIntegerScalarDataPointer => WriteVTKIntegerScalarDataASCII
            WriteVTKFloatVectorDataPointer => WriteVTKFloatVectorDataASCII
            WriteVTKIntegerVectorDataPointer => WriteVTKIntegerVectorDataASCII
            WriteVTKFloatTensorDataPointer => WriteVTKFloatTensorDataASCII
            WriteVTKIntegerTensorDataPointer => WriteVTKIntegerTensorDataASCII
        endif
    end subroutine SetVTKPointers
    
    subroutine WriteVTKOutput()
    !**********************************************************************
    !
    ! Function:  Write VTK output files
    !
    !**********************************************************************
    implicit none

        call CollectAndWriteDataSetScalar()
        call CollectAndWriteDataSetVector()
        call CollectAndWriteDataSetTensor()
        call CollectAndWriteMeshData()
        
    end subroutine WriteVTKOutput

    subroutine WriteVTKNodeScalar(DataSetNameNode,DataSetScalarNode)
    !**********************************************************************
    !
    ! Function:  Write Node Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetNameNode ! dimension(DataSetSizeNode)
    real(REAL_TYPE), dimension(:,:), intent(in) :: DataSetScalarNode ! dimension(NumberNodes,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: NodeID, DataSetID, NumberNodes, DataSetSizeNode
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberNodes = size(DataSetScalarNode,1) ! total number of material points
    DataSetSizeNode = size(DataSetScalarNode,2) ! total number of datasets

    ! open the VTK file for material point data
    VTKFileName = trim(CalParams%FileNames%ProjectName)//'_NodeScalar'//'_'// &
        trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberNodes,' float'
    do NodeID = 1, NumberNodes ! loop over Nodes
        NoCo(1) = NodalCoordinates(NodeID, 1)
        NoCo(2) = NodalCoordinates(NodeID, 2)
        if (NDIM==3) NoCo(3) = NodalCoordinates(NodeID, 3) ! for 2D third coordinate is set equal to 0.0
        write(VTKUnit,*) NoCo
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberNodes,2*NumberNodes
    do NodeID = 1, NumberNodes ! loop over Nodes
        write(VTKUnit,'(A,I8)')'1',NodeID-1 ! write "connectivities" of Nodes (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberNodes
    do NodeID = 1, NumberNodes ! loop over Nodes
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberNodes
    do DataSetID = 1, DataSetSizeNode ! loop over datasets
        write(VTKUnit,'(3A)')'SCALARS ',trim(DataSetNameNode(DataSetID)),' float 1' ! write dataset name
        write(VTKUnit,'(A)')'LOOKUP_TABLE default'
        do NodeID = 1, NumberNodes ! loop over Nodes
            write(VTKUnit,*)DataSetScalarNode(NodeID,DataSetID) ! write data of Node
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKNodeScalar


    subroutine WriteVTKMeshData(DataSetName,DataSetScalar)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), dimension(:), intent(in) :: DataSetName ! dimension(DataSetSize)
    real(REAL_TYPE), dimension(:,:), intent(in) :: DataSetScalar ! dimension(NumberElements,DataSetSize)

    ! local variables
    integer(INTEGER_TYPE) :: I, ElementID, NodeID, DataSetID, NumberElements, NumberNodes, DataSetSize
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName

    NumberNodes = Counters%NodTot ! total number of nodes
    NumberElements = size(DataSetScalar,1) ! total number of elements
    DataSetSize = size(DataSetScalar,2) ! total number of datasets

    ! open the VTK file for mesh data
    if (.not.CalParams%ApplyImplicitQuasiStatic) then
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MeshData'//'_'// &
            trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'
    else
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MeshData'//'_'// &
            trim(CalParams%FileNames%TimeStepExt) &
            //'.vtk'
    end if

    call FileOpen(VTKUnit, VTKFileName)

   
    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)')'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)')'POINTS ',NumberNodes,' float'
    do NodeID = 1, NumberNodes ! loop over number of nodes
        NoCo(1) = NodalCoordinates(NodeID, 1)
        NoCo(2) = NodalCoordinates(NodeID, 2)
        if (NDIM==3) NoCo(3) = NodalCoordinates(NodeID, 3) ! for 2D third coordinate is set equal to 0.0
        write(VTKUnit,*) NoCo
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberElements,(ELEMENTNODES+1)*NumberElements
    do ElementID = 1, NumberElements ! loop over elements
        write(VTKUnit,*) ELEMENTNODES, (((ElementConnectivities(I,ElementID)-1)), I = 1, size(ElementConnectivities,1)) ! write element connectivities (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberElements
    do ElementID = 1, NumberElements ! loop over elements
        write(VTKUnit,*) VTK_CELL ! specification for element, i.e. CELL_TYPE_VTK
    end do

    ! writing dataset attributes for elements (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'CELL_DATA ',NumberElements
    do DataSetID = 1, DataSetSize ! loop over datasets
        write(VTKUnit,'(3A)')'SCALARS ',trim(DataSetName(DataSetID)),' float 1' ! write dataset name
        write(VTKUnit,'(A)')'LOOKUP_TABLE default'
        do ElementID = 1, NumberElements ! loop over elements
            write(VTKUnit,*)DataSetScalar(ElementID,DataSetID) ! write data of elements
        end do
    end do

    close(VTKUnit) ! close the VTK file

    end subroutine WriteVTKMeshData

    subroutine CollectAndWriteDataSetScalar()
    !**********************************************************************
    !
    ! Function:  Write Scalar data VTK output files
    !
    !**********************************************************************

    implicit none

    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID, NumberMaterialPoints, StId
    real(REAL_TYPE), dimension(:), allocatable :: DataSetScalar
    real(REAL_TYPE) :: EpsD, SigD
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName
    character(len=20) :: strst
    character(len=80) :: format_string

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    allocate(DataSetScalar(NumberMaterialPoints))

    ! open the VTK file for material point data
    if (.not.CalParams%ApplyImplicitQuasiStatic) then
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPScalar'//'_'// &
            trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'
    else
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPScalar'//'_'// &
            trim(CalParams%FileNames%TimeStepExt) &
            //'.vtk'
    end if

    call InitialiseVTKFilePointer(NumberMaterialPoints, VTKFileName)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
        DataSetScalar(MaterialPointID) = (SigmaEffArray(MaterialPointID,1) + & ! assign output data to material points
            SigmaEffArray(MaterialPointID,2) + &
            SigmaEffArray(MaterialPointID,3) ) / 3.0 ! valid for 2D and 3D
        else
            DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%WaterPressure
        end if
    end do

    call WriteVTKScalarDataPointer('pressure_liquid',DataSetScalar)
    
    if ((CalParams%PrescribedHead%HydraulicHead == .true.) .and. (CalParams%Multipliers%HydraulicHeadType == LOAD_TYPE_FILE)) then
     do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
         Particles(MaterialPointID)%HydraulicHead =  -Particles(MaterialPointID)%WaterPressure/CalParams%GammaWater + GlobPosArray(MaterialPointID,2)
         DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%HydraulicHead
    end do
    end if
    
    call WriteVTKScalarDataPointer('hydraulic_head',DataSetScalar)
    
    call WriteVTKScalarDataPointer('virtual_material_point',Particles(:)%Kind)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetScalar(MaterialPointID) = GetEpsI(Particles(MaterialPointID),1) + &
            GetEpsI(Particles(MaterialPointID),2) + &
            GetEpsI(Particles(MaterialPointID),3) ! volumetric strain ! valid for 2D and 3D
    end do
    call WriteVTKScalarDataPointer('volumetric_strain_solid',DataSetScalar)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        EpsD = getEpsq( GetEps(Particles(MaterialPointID)) ) ! deviatoric strain
        DataSetScalar(MaterialPointID) = EpsD ! assign output data to material points
    end do
    call WriteVTKScalarDataPointer('deviatoric_strain_solid',DataSetScalar)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetScalar(MaterialPointID) =(  SigmaEffArray(MaterialPointID,1) &
            + SigmaEffArray(MaterialPointID,2) &
            + SigmaEffArray(MaterialPointID,3) ) / 3.0 ! mean effective stress ! valid for 2D and 3D
    end do
    call WriteVTKScalarDataPointer('mean_effective_stress',DataSetScalar)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        SigD = getq( GetSigmaPrin(Particles(MaterialPointID)) ) ! deviatoric stress
        DataSetScalar(MaterialPointID) = SigD ! assign output data to material points
    end do
    call WriteVTKScalarDataPointer('deviatoric_stress',DataSetScalar)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetScalar(MaterialPointID) = GetEpsStepI(Particles(MaterialPointID),1) + &
            GetEpsStepI(Particles(MaterialPointID),2) + &
            GetEpsStepI(Particles(MaterialPointID),3) ! valid for 2D and 3D
    end do
    call WriteVTKScalarDataPointer('incremental_volumetric_strain_solid',DataSetScalar)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        EpsD = getEpsq( GetEpsStep(Particles(MaterialPointID)) ) ! deviatoric strain
        DataSetScalar(MaterialPointID) = EpsD ! assign output data to material points
    end do
    call WriteVTKScalarDataPointer('incremental_deviatoric_strain_solid',DataSetScalar)


    if (.not.CalParams%OutputBasicData) then
        ! these results are written only if the user wants to write all results
        call WriteVTKScalarDataPointer('mass_solid',MassArray(:))

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetScalar(MaterialPointID) = MassArray(MaterialPointID) ! assign output data to material points
            else
                DataSetScalar(MaterialPointID) = MassWaterArray(MaterialPointID)
            end if
        end do

        call WriteVTKScalarDataPointer('mass_liquid', DataSetScalar)

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetScalar(MaterialPointID) = 0.0d0 ! assign output data to material points
            else
                DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%MaterialWeight
            end if
        end do
        call WriteVTKScalarDataPointer('weight_solid',DataSetScalar)

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%MaterialWeight ! assign output data to material points
            else
                DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%WaterWeight
            end if
        end do
        call WriteVTKScalarDataPointer('weight_liquid',DataSetScalar)
        call WriteVTKScalarDataPointer('porosity',                  Particles(:)%Porosity)
        call WriteVTKScalarDataPointer('integration_weight',        Particles(:)%IntegrationWeight)
        call WriteVTKScalarDataPointer('material_point_id',         IDArray(:))
        call WriteVTKScalarDataPointer('element_id',                ElementIDArray(:))
        call WriteVTKScalarDataPointer('entity_id',                 EntityIDArray(:))
        call WriteVTKScalarDataPointer('material_id',               MaterialIDArray(:))
        call WriteVTKScalarDataPointer('damping',                   Particles(:)%Damping)
        call WriteVTKScalarDataPointer('liquid_free_surface',       Particles(:)%LiquidFreeSurface)
        call WriteVTKScalarDataPointer('liquid_free_surface_cumul', Particles(:)%LiquidFreeSurfaceCumul)
        call WriteVTKScalarDataPointer('material_point_kind', Particles(:)%Kind)
        
        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
          if (Particles(MaterialPointID)%IsBoundaryParticle) then
            DataSetScalar(MaterialPointID) = 1
          else
            DataSetScalar(MaterialPointID) = 0  
          end if
        end do
        call WriteVTKScalarDataPointer('boundary_material_point',DataSetScalar)
        
        !Plot the state variables values
        do StId= 1, NSTATEVAR
        if (StId < 10) then
            format_string = "(A4,I1)"
        else
            format_string = "(A4,I2)"
        endif
            write(strst, format_string) "stv_", StId
            call WriteVTKScalarDataPointer(trim(strst),   ESMstatevArray(:,StId))
        end do

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetScalar(MaterialPointID) = GetEpsI(Particles(MaterialPointID),1) + & ! assign output data to material points
                GetEpsI(Particles(MaterialPointID),2) + &
                GetEpsI(Particles(MaterialPointID),3) ! valid for 2D and 3D
            else
                DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%watervolumetricstrain
            end if
        end do
        call WriteVTKScalarDataPointer('volumetric_strain_liquid',DataSetScalar)

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetScalar(MaterialPointID) = Particles(MaterialPointID)%Density ! assign output data to material points
            else
                DataSetScalar(MaterialPointID) = 0.0d0
            end if
        end do
        call WriteVTKScalarDataPointer('density_liq',DataSetScalar)
    endif

    if (isTwoPhaseCalculation() .or. isThreePhaseCalculation()) then
        if (.not.CalParams%OutputBasicData) then ! these results are written only user wants to write all results
            call WriteVTKScalarDataPointer('mass_gas',              Particles(:)%MassGas)
            call WriteVTKScalarDataPointer('mass_mixture',          Particles(:)%MassMixed)
            call WriteVTKScalarDataPointer('pressure_gas',          Particles(:)%GasPressure)
            call WriteVTKScalarDataPointer('weight_gas',            Particles(:)%GasWeight)
            call WriteVTKScalarDataPointer('weight_mixture',        Particles(:)%MixedWeight)
            call WriteVTKScalarDataPointer('initial_porosity',      Particles(:)%InitialPorosity)
            call WriteVTKScalarDataPointer('volumetric_strain_gas', Particles(:)%GasVolumetricStrain)
            call WriteVTKScalarDataPointer('degree_saturation_liquid', Particles(:)%DegreeSaturation)

        end if
    endif

    close (VTKUnit)
    deallocate(DataSetScalar)

    end subroutine CollectAndWriteDataSetScalar
    
    

    subroutine CollectAndWriteDataSetVector()
    !**********************************************************************
    !
    ! Function:  Write vectorial data VTK output files
    !
    !**********************************************************************
    implicit none
    ! local variables
    integer(INTEGER_TYPE) :: N, K, MaterialPointID, NumberMaterialPoints
    real(REAL_TYPE), dimension(:,:), allocatable :: DataSetVector
    real(REAL_TYPE), dimension(NVECTOR) :: VectorValue
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName
    logical :: HasValue

    ! open the VTK file for material point data
    if (.not.CalParams%ApplyImplicitQuasiStatic) then
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPVector'//'_'// trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'
    else
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPVector'//'_'// trim(CalParams%FileNames%TimeStepExt)//'.vtk'
    end if

    NumberMaterialPoints = Counters%NParticles ! total number of material points
    allocate(DataSetVector(NumberMaterialPoints,NVECTOR))

    call InitialiseVTKFilePointer(NumberMaterialPoints, VTKFileName)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
        DataSetVector(MaterialPointID,:) = 0.0
        else
            DataSetVector(MaterialPointID,:) = VelocityArray(MaterialPointID,:)
        end if
    end do
    call WriteVTKVectorDataPointer('velocity_solid',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid)))then
        DataSetVector(MaterialPointID,:) = VelocityArray(MaterialPointID,:)
        else
            DataSetVector(MaterialPointID,:) = VelocityWaterArray(MaterialPointID,:)
        end if
    end do
    call WriteVTKVectorDataPointer('velocity_liquid',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
        DataSetVector(MaterialPointID,:) = 0.0
        else
            DataSetVector(MaterialPointID,:) = UArray(MaterialPointID,:)
        end if
    end do
    call WriteVTKVectorDataPointer('displacement_solid',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
        DataSetVector(MaterialPointID,:) = UArray(MaterialPointID,:)
        else
            DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%UW(:)
        end if
    end do
    call WriteVTKVectorDataPointer('displacement_liquid',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetVector(MaterialPointID,:) = AccelerationArray(MaterialPointID,:)
    end do
    call WriteVTKVectorDataPointer('acceleration_solid',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetVector(MaterialPointID,:) = GlobPosArray(MaterialPointID,:)
    end do
    call WriteVTKVectorDataPointer('global_position',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%LocPos(:)
    end do
    call WriteVTKVectorDataPointer('local_position',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%FBody(:)
    end do
    call WriteVTKVectorDataPointer('bodyforce',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
        DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%FBody(:)
        else
            DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%FBodyWater(:)
        end if
    end do
    call WriteVTKVectorDataPointer('bodyforce_liquid',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%FBodyMixed(:)
    end do
    call WriteVTKVectorDataPointer('bodyforce_mixture',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        HasValue = .false.
        do N = 1, NVECTOR
          do K=1,MAX_LOAD_SYSTEMS  
            HasValue = HasValue .or. (Particles(MaterialPointID)%FExt(N,K) /= 0.0)
          end do
        end do
        if (HasValue) then
            VectorValue = 0.0
          do K=1,MAX_LOAD_SYSTEMS    
            VectorValue(:) = VectorValue(:)+Particles(MaterialPointID)%FExt(:,K)*CalParams%Multipliers%SolidACurrent(k)
          end do
            DataSetVector(MaterialPointID,:) = VectorValue(:)
        else
            DataSetVector(MaterialPointID,:) = 0.0
        end if
    end do
    call WriteVTKVectorDataPointer('externalforce',DataSetVector)

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        HasValue = .false.
        do N = 1, NVECTOR
          do K=1,MAX_LOAD_SYSTEMS  
            HasValue = HasValue .or. (Particles(MaterialPointID)%FExtWater(N,K) /= 0.0)
          end do
        end do
        if (HasValue) then
            VectorValue = 0.0
          do K=1,MAX_LOAD_SYSTEMS    
            VectorValue(:) = VectorValue(:)+Particles(MaterialPointID)%FExtWater(:,K)*CalParams%Multipliers%WaterACurrent(k)
          end do
            DataSetVector(MaterialPointID,:) = VectorValue(:)
        else
            DataSetVector(MaterialPointID,:) = 0.0
        end if
    end do
    call WriteVTKVectorDataPointer('externalforce_liquid',DataSetVector)

    if (isThreePhaseCalculation()) then
        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            DataSetVector(MaterialPointID,:) = VelocityGasArray(MaterialPointID,:)
        end do
        call WriteVTKVectorDataPointer('velocity_gas',DataSetVector)

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%UG(:)
        end do
        call WriteVTKVectorDataPointer('displacement_gas',DataSetVector)

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            DataSetVector(MaterialPointID,:) = Particles(MaterialPointID)%FBodyGas(:)
        end do
        call WriteVTKVectorDataPointer('bodyforce_gas',DataSetVector)

        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
          HasValue = .false.
          do N = 1, NVECTOR
            do K=1,MAX_LOAD_SYSTEMS  
              HasValue = HasValue .or. (Particles(MaterialPointID)%FExtGas(N,K) /= 0.0)
            end do
          end do
          if (HasValue) then
            VectorValue = 0.0
            do K=1,MAX_LOAD_SYSTEMS    
              VectorValue(:) = VectorValue(:)+Particles(MaterialPointID)%FExtGas(:,K)*CalParams%Multipliers%GasACurrent(k)
            end do
            DataSetVector(MaterialPointID,:) = VectorValue(:)
          else
            DataSetVector(MaterialPointID,:) = 0.0
          end if
        end do
        call WriteVTKVectorDataPointer('externalforce_gas',DataSetVector)
    end if

    close(VTKunit)
    deallocate(DataSetVector)

    end subroutine CollectAndWriteDataSetVector


    subroutine CollectAndWriteDataSetTensor()
    !**********************************************************************
    !
    ! Function:  Write tensorial data VTK output files
    !
    !**********************************************************************
    implicit none
    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID, NumberMaterialPoints
    real(REAL_TYPE), dimension(:,:,:), allocatable :: DataSetTensor
    real(REAL_TYPE)::Wp, Sr
    character(len=MAX_FILENAME_LENGTH) :: VTKFileName
    logical :: IsUndrEffectiveStress
    integer(INTEGER_TYPE) :: IDset ! ID of material parameter set
    
    
    NumberMaterialPoints = Counters%NParticles ! total number of material points
    allocate(DataSetTensor(NumberMaterialPoints,NPRINCIPAL,NPRINCIPAL))

    ! open the VTK file for material point data
    if (.not.CalParams%ApplyImplicitQuasiStatic) then
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPTensor'//'_'// trim(CalParams%FileNames%LoadStepExt)//trim(CalParams%FileNames%TimeStepExt)//'.vtk'
    else
        VTKFileName = trim(CalParams%FileNames%ProjectName)//'_MPTensor'//'_'// trim(CalParams%FileNames%TimeStepExt) //'.vtk'
    end if

    call InitialiseVTKFilePointer(NumberMaterialPoints, VTKFileName)

    DataSetTensor = 0.0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if (NDIM == 3) then ! 3D case
            DataSetTensor(MaterialPointID,1,1) = Particles(MaterialPointID)%Eps(1)
            DataSetTensor(MaterialPointID,1,2) = Particles(MaterialPointID)%Eps(4)
            DataSetTensor(MaterialPointID,1,3) = Particles(MaterialPointID)%Eps(5)
            DataSetTensor(MaterialPointID,2,1) = Particles(MaterialPointID)%Eps(4)
            DataSetTensor(MaterialPointID,2,2) = Particles(MaterialPointID)%Eps(2)
            DataSetTensor(MaterialPointID,2,3) = Particles(MaterialPointID)%Eps(6)
            DataSetTensor(MaterialPointID,3,1) = Particles(MaterialPointID)%Eps(5)
            DataSetTensor(MaterialPointID,3,2) = Particles(MaterialPointID)%Eps(6)
            DataSetTensor(MaterialPointID,3,3) = Particles(MaterialPointID)%Eps(3)
        else if (NDIM == 2) then ! 2D case
            DataSetTensor(MaterialPointID,1,1) = Particles(MaterialPointID)%Eps(1)
            DataSetTensor(MaterialPointID,1,2) = Particles(MaterialPointID)%Eps(4)
            DataSetTensor(MaterialPointID,2,1) = Particles(MaterialPointID)%Eps(4)
            DataSetTensor(MaterialPointID,2,2) = Particles(MaterialPointID)%Eps(2)
            DataSetTensor(MaterialPointID,3,3) = Particles(MaterialPointID)%Eps(3)
        end if
    end do
    call WriteVTKTensorDataPointer('strain',DataSetTensor)

    DataSetTensor = 0.0
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid).or. &
            (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeMixture)) then
            if (NDIM == 3) then ! 3D case
                DataSetTensor(MaterialPointID,1,1) = SigmaEffArray(MaterialPointID,1)
                DataSetTensor(MaterialPointID,1,2) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,1,3) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(MaterialPointID,2,1) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,2) = SigmaEffArray(MaterialPointID,2)
                DataSetTensor(MaterialPointID,2,3) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(MaterialPointID,3,1) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(MaterialPointID,3,2) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(MaterialPointID,3,3) = SigmaEffArray(MaterialPointID,3)
            elseif (NDIM == 2) then ! 2D case
                DataSetTensor(MaterialPointID,1,1) = SigmaEffArray(MaterialPointID,1)
                DataSetTensor(MaterialPointID,1,2) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,1) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,2) = SigmaEffArray(MaterialPointID,2)
                DataSetTensor(MaterialPointID,3,3) = SigmaEffArray(MaterialPointID,3)
            end if
        end if
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetTensor(MaterialPointID,:,:) = 0.0
        end if
    end do
    call WriteVTKTensorDataPointer('eff_stress_solid',DataSetTensor)
    
    DataSetTensor = 0.0
    Wp=0
    Sr=1
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        if ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid).or. &
            (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeMixture)) then
            WP = Particles(MaterialPointID)%WaterPressure
            Sr = Particles(MaterialPointID)%DegreeSaturation
            if (NDIM == 3) then ! 3D case
                DataSetTensor(MaterialPointID,1,1) = SigmaEffArray(MaterialPointID,1)+WP*Sr
                DataSetTensor(MaterialPointID,1,2) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,1,3) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(MaterialPointID,2,1) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,2) = SigmaEffArray(MaterialPointID,2)+WP*Sr
                DataSetTensor(MaterialPointID,2,3) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(MaterialPointID,3,1) = SigmaEffArray(MaterialPointID,6)
                DataSetTensor(MaterialPointID,3,2) = SigmaEffArray(MaterialPointID,5)
                DataSetTensor(MaterialPointID,3,3) = SigmaEffArray(MaterialPointID,3)+WP*Sr
            elseif (NDIM == 2) then !2D case
                DataSetTensor(MaterialPointID,1,1) = SigmaEffArray(MaterialPointID,1)+WP*Sr
                DataSetTensor(MaterialPointID,1,2) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,1) = SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,2) = SigmaEffArray(MaterialPointID,2)+WP*Sr
                DataSetTensor(MaterialPointID,3,3) = SigmaEffArray(MaterialPointID,3)+WP*Sr
            end if
        end if
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetTensor(MaterialPointID,:,:) = 0.0
        end if
    end do
    call WriteVTKTensorDataPointer('total_stress_solid',DataSetTensor)

    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    DataSetTensor = 0.0
    Wp=0
    Sr=1
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        !if ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid).or. &
        !    (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeMixture)) then
        IDset = MaterialIDArray(MaterialPointID) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file

        IsUndrEffectiveStress = &
        !code version 2016 and previous
        ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(IDSet)%MaterialType)=='2-phase')) .or. &
        !code version 2017.1 and following
        (trim(MatParams(IDSet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
        
        if (IsUndrEffectiveStress) then 
            WP = Particles(MaterialPointID)%WaterPressure
            Sr = Particles(MaterialPointID)%DegreeSaturation
            if (NDIM == 3) then ! 3D case
                DataSetTensor(MaterialPointID,1,1) = WP*Sr!SigmaEffArray(MaterialPointID,1)+WP*Sr
                DataSetTensor(MaterialPointID,1,2) = 0.0!SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,1,3) = 0.0!SigmaEffArray(MaterialPointID,6)
                DataSetTensor(MaterialPointID,2,1) = 0.0!SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,2) = WP*Sr!SigmaEffArray(MaterialPointID,2)+WP*Sr
                DataSetTensor(MaterialPointID,2,3) = 0.0!SigmaEffArray(MaterialPointID,5)
                DataSetTensor(MaterialPointID,3,1) = 0.0!SigmaEffArray(MaterialPointID,6)
                DataSetTensor(MaterialPointID,3,2) = 0.0!SigmaEffArray(MaterialPointID,5)
                DataSetTensor(MaterialPointID,3,3) = WP*Sr!SigmaEffArray(MaterialPointID,3)+WP*Sr
            elseif (NDIM == 2) then !2D case
                DataSetTensor(MaterialPointID,1,1) = WP*Sr!SigmaEffArray(MaterialPointID,1)+WP*Sr
                DataSetTensor(MaterialPointID,1,2) = 0.0!SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,1) = 0.0!SigmaEffArray(MaterialPointID,4)
                DataSetTensor(MaterialPointID,2,2) = WP*Sr!SigmaEffArray(MaterialPointID,2)+WP*Sr
                DataSetTensor(MaterialPointID,3,3) = WP*Sr!SigmaEffArray(MaterialPointID,3)+WP*Sr
            end if
        end if
        if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
            ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
            DataSetTensor(MaterialPointID,:,:) = 0.0
        end if
    end do
    
    
    
    call WriteVTKTensorDataPointer('stress_liquid',DataSetTensor)
    
    
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
    
    
    if (.not.CalParams%OutputBasicData) then
        DataSetTensor = 0.0
        do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
            if ((MatParams(MaterialIDArray(MaterialPointID))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(MaterialPointID))%MaterialPhases=='1-phase-liquid'.or. &
                ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeLiquid))) then
                if (NDIM == 3) then ! 3D case
                    DataSetTensor(MaterialPointID,1,1) = SigmaEffArray(MaterialPointID,1)
                    DataSetTensor(MaterialPointID,1,2) = SigmaEffArray(MaterialPointID,4)
                    DataSetTensor(MaterialPointID,1,3) = SigmaEffArray(MaterialPointID,6)
                    DataSetTensor(MaterialPointID,2,1) = SigmaEffArray(MaterialPointID,4)
                    DataSetTensor(MaterialPointID,2,2) = SigmaEffArray(MaterialPointID,2)
                    DataSetTensor(MaterialPointID,2,3) = SigmaEffArray(MaterialPointID,5)
                    DataSetTensor(MaterialPointID,3,1) = SigmaEffArray(MaterialPointID,6)
                    DataSetTensor(MaterialPointID,3,2) = SigmaEffArray(MaterialPointID,5)
                    DataSetTensor(MaterialPointID,3,3) = SigmaEffArray(MaterialPointID,3)
                elseif (NDIM == 2) then !2D case
                    DataSetTensor(MaterialPointID,1,1) = SigmaEffArray(MaterialPointID,1)
                    DataSetTensor(MaterialPointID,1,2) = SigmaEffArray(MaterialPointID,4)
                    DataSetTensor(MaterialPointID,2,1) = SigmaEffArray(MaterialPointID,4)
                    DataSetTensor(MaterialPointID,2,2) = SigmaEffArray(MaterialPointID,2)
                    DataSetTensor(MaterialPointID,3,3) = SigmaEffArray(MaterialPointID,3)
                end if
            end if
            if ((MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeSolid).or. &
                (MaterialPointTypeArray(MaterialPointID)==MaterialPointTypeMixture)) then
                DataSetTensor(MaterialPointID,:,:) = 0.0
            end if
        end do
        call WriteVTKTensorDataPointer('stress_liquid',DataSetTensor)
    endif

    close(VTKunit)
    deallocate(DataSetTensor)

    end subroutine CollectAndWriteDataSetTensor


    subroutine CollectAndWriteMeshData()
    !**********************************************************************
    !
    ! Function:  Write mesh data VTK output files
    !
    !**********************************************************************
    implicit none
    ! local variables
    integer(INTEGER_TYPE) :: ElementID, DataSetID, &
        NumberElements, DataSetSize
    character(len=50), dimension(:), allocatable :: DataSetName
    real(REAL_TYPE), dimension(:,:), allocatable :: DataSetScalar
    
    !********** START: writing scalar mesh data ****************************************
    
    NumberElements = Counters%NEl ! total number of elements
    DataSetSize = 2 ! total number of datasets for elements to be written
    allocate(DataSetName(DataSetSize))
    allocate(DataSetScalar(NumberElements,DataSetSize))

    !DataSetID = 1 ! index of dataset
   
    do DataSetID = 1, 2
        if (DataSetID == 1) then
            DataSetName(DataSetID) = 'active_elements' ! assign name to dataset
            do ElementID = 1, NumberElements ! loop over elements
                if (IsActiveElement(ElementID)) then
                    DataSetScalar(ElementID,DataSetID) = 1.0 ! assign output data to elements
                else
                    DataSetScalar(ElementID,DataSetID) = 0.0
                end if
            end do 
            call WriteVTKMeshData(DataSetName,DataSetScalar) ! write all datasets for elements
        else
            if (allocated(HydraulicHeadLoadedElemID)) then
                DataSetName(DataSetID) = 'border_elements' ! assign name to dataset
                do ElementID = 1, NumberElements ! loop over elements
                    if (HydraulicHeadLoadedElemID(ElementID)) then
                        DataSetScalar(ElementID,DataSetID) = 1.0 ! assign output data to elements
                    else
                        DataSetScalar(ElementID,DataSetID) = 0.0
                    end if
                end do
                call WriteVTKMeshData(DataSetName,DataSetScalar)! write all datasets for elements
            end if
        end if
    end do
     
    !********** END: writing scalar mesh data ****************************************

    deallocate(DataSetName,DataSetScalar)

    end subroutine CollectAndWriteMeshData

    subroutine CollectAndWriteNodeScalar()
    !**********************************************************************
    !
    ! Function:  Write node data VTK output files
    !
    !**********************************************************************
    implicit none
    ! local variables
    integer(INTEGER_TYPE) :: DataSetID, NodeID, &
        NumberNodes, NumberMaterialPoints, DataSetSize
    character(len=50), dimension(:), allocatable :: DataSetNameNode
    real(REAL_TYPE), dimension(:,:), allocatable :: DataSetScalarNode

    NumberMaterialPoints = Counters%NParticles ! total number of material points

    !********** START: writing scalar node data ****************************************
    NumberNodes = Counters%NodTot ! total number of nodes
    DataSetSize = 1 ! total number of datasets for material points to be written
    allocate(DataSetNameNode(DataSetSize))
    allocate(DataSetScalarNode(NumberNodes,DataSetSize))

    DataSetScalarNode = 0.0

    DataSetID = 1 ! index of dataset
    DataSetNameNode(DataSetID) = 'Nodal_Density' ! assign name to dataset
    if (allocated(NodalDensity)) then
        do NodeID = 1, NumberNodes ! loop over material points
            DataSetScalarNode(NodeID,DataSetID) = NodalDensity(NodeID) ! assign output data to node
        end do
    end if

    !********** END: writing scalar node data ****************************************
    end subroutine CollectAndWriteNodeScalar


    end module ModWriteVTKOutput