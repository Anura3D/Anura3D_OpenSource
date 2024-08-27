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


module ModWriteVTKASCII
    !**********************************************************************
    !
    !    Function: Contains pointer routines for writing ASCII output data in VTK format
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************
    
    use ModGlobalConstants
    use ModReadCalculationData
    use ModMPMData
    use ModFileIO
    
    contains
    
    subroutine InitialiseVTKFileASCII(NumberMaterialPoints, VTKFileName)
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

    call FileOpen(VTKUnit, VTKFileName)

    write(VTKUnit,'(A)') trim(VTK_VERSION) ! writing file version and identifier
    write(VTKUnit,'(A)') trim(CalParams%FileNames%ProjectName) ! writing header
    write(VTKUnit,'(A)') trim(ASCII_FORMAT) ! writing file format

    ! writing dataset structure (defines geometry and topology of dataset)
    write(VTKUnit,'(A)') 'DATASET UNSTRUCTURED_GRID'

    write(VTKUnit,'(A,I8,A)') 'POINTS ',NumberMaterialPoints,' float'
    do MaterialPointID = 1, NumberMaterialPoints ! loop over number of material points
        NoCo(1) = GlobPosArray(MaterialPointID, 1)
        NoCo(2) = GlobPosArray(MaterialPointID, 2)
        if (NDIM==3) NoCo(3) = GlobPosArray(MaterialPointID, 3) ! for 2D third coordinate is set equal to 0.0
        write(VTKUnit,*) NoCo
    end do

    write(VTKUnit,'(A,2I8)')'CELLS ',NumberMaterialPoints,2*NumberMaterialPoints
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        write(VTKUnit,'(A,I8)')'1',MaterialPointID-1 ! write "connectivities" of material points (note that the index has to be decreased by 1)
    end do

    write(VTKUnit,'(A,I8)')'CELL_TYPES ',NumberMaterialPoints
    do MaterialPointID = 1, NumberMaterialPoints ! loop over material points
        write(VTKUnit,*)'1' ! i.e. specification for single point, VTK_VERTEX=1
    end do

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(A,I8)')'POINT_DATA ',NumberMaterialPoints

    end subroutine InitialiseVTKFileASCII
    
    subroutine WriteVTKFloatScalarDataASCII(DataSetName,DataSetScalar)
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
    real(REAL_TYPE) :: Value

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(3A)')'SCALARS ',trim(DataSetName),' float 1' ! write dataset name
    write(VTKUnit,'(A)')'LOOKUP_TABLE default'
    do i = 1, size(DataSetScalar) ! loop over material points
        if ( DataSetScalar(i) < SMALL .and. DataSetScalar(i) > -SMALL ) then
            Value = 0.0 ! set value to zero
        else
            Value = DataSetScalar(i) ! use original data
        end if
        write(VTKUnit,*) Value ! write data of material point
    end do

    end subroutine WriteVTKFloatScalarDataASCII
    
    subroutine WriteVTKIntegerScalarDataASCII(DataSetName,DataSetScalar)
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
    write(VTKUnit,'(3A)')'SCALARS ',trim(DataSetName),' int 1' ! write dataset name
    write(VTKUnit,'(A)')'LOOKUP_TABLE default'
    do i = 1, size(DataSetScalar) ! loop over material points
        write(VTKUnit,*) DataSetScalar(i) ! write data of material point
    end do

    end subroutine WriteVTKIntegerScalarDataASCII
    
    subroutine WriteVTKFloatVectorDataASCII(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    real(REAL_TYPE), dimension(:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i, j
    real(REAL_TYPE), dimension(3) :: Value ! dimension is always 3 as this is required by ParaView, for 2D Value(3) is equal to 0.0

    Value = 0.0
    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(3A)') 'VECTORS ',trim(DataSetName),' float' ! write dataset name

    do i = 1, size(DataSet,1) ! loop over material points

        do j = 1, NVECTOR ! loop over dimension
            if ( DataSet(i,j) < SMALL .and. DataSet(i,j) > -SMALL ) then
                Value(j) = 0.0 ! set value to zero
            else
                Value(j) = DataSet(i,j) ! use original data
            end if
        end do

        write(VTKUnit,*) Value

    end do

    end subroutine WriteVTKFloatVectorDataASCII
    
    subroutine WriteVTKIntegerVectorDataASCII(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    integer(INTEGER_TYPE), dimension(:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i, j
    integer(INTEGER_TYPE), dimension(3) :: Value ! dimension is always 3 as this is required by ParaView, for 2D Value(3) is equal to 0.0

    Value = 0
    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(3A)') 'VECTORS ',trim(DataSetName),' int' ! write dataset name

    do i = 1, size(DataSet,1) ! loop over material points

        do j = 1, NVECTOR ! loop over dimension
            Value(j) = DataSet(i, j)
        end do

        write(VTKUnit,*) Value

    end do

    end subroutine WriteVTKIntegerVectorDataASCII
    
    subroutine WriteVTKFloatTensorDataASCII(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    real(REAL_TYPE), dimension(:,:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i, j, k, n, m
    real(REAL_TYPE), dimension(size(DataSet,2),size(DataSet,3)) :: Value

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(3A)')'TENSORS ',trim(DataSetName),' float' ! write dataset name
    do i = 1, size(DataSet,1) ! loop over material points

        do j = 1, size(Value,1)
            do k = 1, size(Value,2)
                if ( DataSet(i,j,k) < SMALL .and. DataSet(i,j,k) > -SMALL ) then
                    Value(j,k) = 0.0 ! set value to zero
                else
                    Value(j,k) = DataSet(i,j,k) ! use original data
                end if
            end do
        end do
        do N = 1, size(Value,1)
            write(VTKUnit,*) (Value(N, M), M = 1, size(Value,2))
        end do
    end do

    end subroutine WriteVTKFloatTensorDataASCII
    
    subroutine WriteVTKIntegerTensorDataASCII(DataSetName,DataSet)
    !**********************************************************************
    !
    ! Function:  Write Material Point Data in VTK format
    !
    !**********************************************************************

    implicit none

    character(len=*), intent(in) :: DataSetName
    integer(INTEGER_TYPE), dimension(:,:,:), intent(in) :: DataSet

    ! local variables
    integer(INTEGER_TYPE) :: i, j, k

    ! writing dataset attributes for material points (defines the actual dataset values)
    write(VTKUnit,'(3A)')'TENSORS ',trim(DataSetName),' int' ! write dataset name
    do i = 1, size(DataSet,1) ! loop over material points
        do j = 1, size(DataSet,2)
            write(VTKUnit,*) (DataSet(i,j,k), k = 1, size(DataSet,3))
        end do
    end do

    end subroutine WriteVTKIntegerTensorDataASCII
    
    end module