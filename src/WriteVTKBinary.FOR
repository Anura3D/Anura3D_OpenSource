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


module ModWriteVTKBinary
    !**********************************************************************
    !
    !    Function: Contains pointer routines for writing binary output data in VTK format
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************
    
    use ModGlobalConstants
    use ModReadCalculationData
    use ModMPMData
    use ModFileIO
    
    character(1), parameter:: end_rec = char(10)    ! end-character for binary-record finalize
    integer(INTEGER_TYPE), parameter:: maxlen = 500         ! max number of characters os static string
    character(len=maxlen) ::  cbuffer ! string buffer for binary data
    
    contains
    
    subroutine InitialiseVTKFileBinary(NumberMaterialPoints, VTKFileName)
    !**********************************************************************
    !
    !    Function: initialise VTK files
    !
    ! I NumberMaterialPoints: Number of material  points
    ! O VTKFileName: VTK output file name
    !
    !**********************************************************************    
    implicit none

    integer(INTEGER_TYPE), intent(in) :: NumberMaterialPoints
    character*(*), intent(in) :: VTKFileName
    ! local variables
    integer(INTEGER_TYPE) :: MaterialPointID
    real(REAL_TYPE), dimension(3) :: NoCo = 0.0 ! dimension 3 as ParaView requires 3D data structure

    ! open the VTK file for material point data

    call FileOpenWriteBinary(VTKUnit, VTKFileName)

    write(VTKUnit) trim(VTK_VERSION)//end_rec ! writing file version and identifier
    write(VTKUnit) trim(CalParams%FileNames%ProjectName)//end_rec ! writing header
    write(VTKUnit) trim(BINARY_FORMAT)//end_rec ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit) 'DATASET UNSTRUCTURED_GRID'//end_rec

    write(cbuffer, '(A,I8,A)') 'POINTS ',NumberMaterialPoints,' DOUBLE'
    write(VTKUnit) trim(cbuffer)//end_rec

    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        NoCo(1) = GlobPosArray(MaterialPointID, 1)
        NoCo(2) = GlobPosArray(MaterialPointID, 2)
        if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
        write(VTKUnit) NoCo
    end do

    write(cbuffer,'(A,2I8)')'CELLS ',NumberMaterialPoints,2*NumberMaterialPoints
    write(VTKUnit) trim(cbuffer)//end_rec

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        write(VTKUnit) 1, MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(cbuffer,'(A,I8)')'CELL_TYPES ',NumberMaterialPoints
    write(VTKUnit) trim(cbuffer)//end_rec

    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        write(VTKUnit) 1 ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(cbuffer,'(A,I8)')'POINT_DATA ',NumberMaterialPoints
    write(VTKUnit) trim(cbuffer)//end_rec

    end subroutine InitialiseVTKFileBinary
    
    subroutine WriteVTKFloatScalarDataBinary(DataSetName,DataSetScalar)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    real(REAL_TYPE), dimension(:), intent(in) :: DataSetScalar

    ! local variables
    integer(INTEGER_TYPE) :: i 

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit)'SCALARS '//trim(DataSetName)//' DOUBLE 1'//end_rec ! write dataset name
    write(VTKUnit)'LOOKUP_TABLE default'//end_rec

    do i = 1, size(DataSetScalar) ! loop over material points
        write(VTKUnit) DataSetScalar(i)
    end do

    end subroutine WriteVTKFloatScalarDataBinary
    
    subroutine WriteVTKIntegerScalarDataBinary(DataSetName,DataSetScalar)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    integer(INTEGER_TYPE), dimension(:), intent(in) :: DataSetScalar

    ! local variables
    integer(INTEGER_TYPE) :: i

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit)'SCALARS '//trim(DataSetName)//' int 1'//end_rec ! write dataset name
    write(VTKUnit)'LOOKUP_TABLE default'//end_rec
    do i = 1, size(DataSetScalar) ! loop over material points
        write(VTKUnit) DataSetScalar(i) ! write data of material point
    end do

    end subroutine WriteVTKIntegerScalarDataBinary
    
    subroutine WriteVTKFloatVectorDataBinary(DataSetName,DataSet)
    !**********************************************************************

    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    real(REAL_TYPE), dimension(:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i
    real(REAL_TYPE), dimension(3) :: Value ! dimension is always 3 as this is required by ParaView, for 2D Value(3) is equal to 0.0

    Value = 0.0
    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit) 'VECTORS '//trim(DataSetName)//' DOUBLE'//end_rec ! write dataset name
    do i = 1, size(DataSet,1) ! loop over material points
        Value(1:NVECTOR) = DataSet(i, 1:NVECTOR) ! use original data
        write(VTKUnit) Value
    end do

    end subroutine WriteVTKFloatVectorDataBinary
    
    subroutine WriteVTKIntegerVectorDataBinary(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    integer(INTEGER_TYPE), dimension(:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i
    integer(INTEGER_TYPE), dimension(3) :: Value ! dimension is always 3 as this is required by ParaView, for 2D Value(3) is equal to 0.0

    Value = 0
    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit) 'VECTORS '//trim(DataSetName)//' int'//end_rec ! write dataset name

    do i = 1, size(DataSet,1) ! loop over material points
        Value(1:NVECTOR) = DataSet(i, 1:NVECTOR)
        write(VTKUnit) Value
    end do

    end subroutine WriteVTKIntegerVectorDataBinary
    
    subroutine WriteVTKFloatTensorDataBinary(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    real(REAL_TYPE), dimension(:,:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit)'TENSORS '//trim(DataSetName)//' DOUBLE'//end_rec ! write dataset name
    do i = 1, size(DataSet,1) ! loop over material points
        write(VTKUnit) DataSet(i,1,1), &
            DataSet(i,1,2), &
            DataSet(i,1,3), &
            DataSet(i,2,1), &
            DataSet(i,2,2), &
            DataSet(i,2,3), &
            DataSet(i,3,1), &
            DataSet(i,3,2), &
            DataSet(i,3,3)
    end do

    end subroutine WriteVTKFloatTensorDataBinary
    
    subroutine WriteVTKIntegerTensorDataBinary(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    integer(INTEGER_TYPE), dimension(:,:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit)'TENSORS '//trim(DataSetName)//' int'//end_rec ! write dataset name
    do i = 1, size(DataSet,1) ! loop over material points
        write(VTKUnit) DataSet(i,1,1), &
            DataSet(i,1,2), &
            DataSet(i,1,3), &
            DataSet(i,2,1), &
            DataSet(i,2,2), &
            DataSet(i,2,3), &
            DataSet(i,3,1), &
            DataSet(i,3,2), &
            DataSet(i,3,3)
    end do

    end subroutine WriteVTKIntegerTensorDataBinary
    end module 