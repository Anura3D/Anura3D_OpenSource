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


    module ModInitialiseElementType
    !**********************************************************************
    !
    ! Function : Contains initialisation routines of kernel and project
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************
      
    use ModGlobalConstants
    use ModElementEvaluationTETRA
    use ModElementEvaluationTRI
    use ModMeshAdjacencies
      
    contains 
      
      subroutine InitialiseElementType()
      !**********************************************************************
      !
      ! Function : Initialises the element type in global variables
      !            Note: Only the following element types are available, 
      !                  triangular 3-noded (TRI3)
      !                  tetrahedral_old (TETRAOLD)
      !
      !**********************************************************************
 
      implicit none
 
      select case (ELEMENTTYPE)
            
          case(TRI3) ! 'triangular_3-noded' 
              CheckForGlobPosPointer => CheckTRIForGlobPos
              Gauss_Q1Pointer => GaussTRI_Q1         
              InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsLINE2
              InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsTRI3
              IsInsideElementLocPosPointer => IsInsideElementLocPosTRI 
              GetMinAltitudePointer => GetMinAltitudeTri
              InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesTRI
              ShapeLocPosPointer => ShapeLocPosTRI3
              RearrangeConnectivitiesPointer => RearrangeConnectivitiesLINE2
          case(TETRAOLD) ! 'tetrahedral_old' 
              CheckForGlobPosPointer => CheckTetrahedronForGlobPos
              Gauss_Q1Pointer => GaussTETRA_Q1
              InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsTRI3
              InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsTETRA4
              IsInsideElementLocPosPointer => IsInsideElementLocPosTETRA
              GetMinAltitudePointer => GetMinAltitudeTetra
              InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesTETRA
              RearrangeConnectivitiesPointer => RearrangeConnectivitiesTRI6
         
          !**************************NOT AVAILABLE***************************
          !case(TRI6) ! 'triangular_6-noded'
          !    CheckForGlobPosPointer => CheckTRIForGlobPos
          !    Gauss_Q1Pointer => GaussTRI_Q1      
          !    !InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsLINE3  --> to be added
          !    !InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsTRI6  --> to be added
          !    IsInsideElementLocPosPointer => IsInsideElementLocPosTRI 
          !    GetMinAltitudePointer => GetMinAltitudeTri
          !    InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesTRI
          !    ShapeLocPosPointer => ShapeLocPosTRI6
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesLINE3
          !    
          !case(QUAD4) ! 'quadrilateral_4-noded' 
          !    CheckForGlobPosPointer => CheckQUADForGlobPos
          !    Gauss_Q1Pointer => GaussQUAD_Q1 
          !    InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsLINE2
          !    InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsQUAD4
          !    IsInsideElementLocPosPointer => IsInsideElementLocPosQUAD 
          !    GetMinAltitudePointer => GetMinAltitudeQUAD
          !    InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesQUAD
          !    ShapeLocPosPointer => ShapeLocPosQUAD4
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesLINE2
          !    
          !case(QUAD8) ! 'quadrilateral_8-noded'    
          !    CheckForGlobPosPointer => CheckQUADForGlobPos
          !    Gauss_Q1Pointer => GaussQUAD_Q1 
          !    InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsLINE3  --> to be added
          !    InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsQUAD8  --> to be added
          !    IsInsideElementLocPosPointer => IsInsideElementLocPosQUAD 
          !    GetMinAltitudePointer => GetMinAltitudeQUAD
          !    InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesQUAD
          !    ShapeLocPosPointer => ShapeLocPosQUAD4
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesLINE2
          !                 
          !case(TETRA4) ! 'tetrahedral_4-noded'       
          !    CheckForGlobPosPointer => CheckTetrahedronForGlobPos
          !    Gauss_Q1Pointer => GaussTETRA_Q1
          !    InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsTRI3
          !    InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsTETRA4
          !    IsInsideElementLocPosPointer => IsInsideElementLocPosTETRA
          !    GetMinAltitudePointer => GetMinAltitudeTetra
          !    InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesTETRA
          !    ShapeLocPosPointer => ShapeLocPosTETRA4
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesTRI3
          !  
          !case(TETRA10) ! 'tetrahedral_10-noded'
          !    CheckForGlobPosPointer => CheckTetrahedronForGlobPos
          !    Gauss_Q1Pointer => GaussTETRA_Q1
          !    InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsTRI6  --> to be added
          !    InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsTETRA10   --> to be added
          !    IsInsideElementLocPosPointer => IsInsideElementLocPosTETRA
          !    GetMinAltitudePointer => GetMinAltitudeTetra
          !    InitialLocalMaterialPointCoordinatesPointer => InitialLocalMaterialPointCoordinatesTETRA
          !    ShapeLocPosPointer => ShapeLocPosTETRA10
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesTRI6
          !  
          !case(HEXA8) ! 'hexahedral_8-noded'
          !    InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsQUAD4   --> to be added
          !    InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsHEXA8  --> to be added 
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesQUAD4
          !
          !case(HEXA20) ! 'hexahedral_20-noded'
          !    InitialiseShapeFunctionsBoundaryPointer => InitialiseShapeFunctionsQUAD8  --> to be added
          !    InitialiseShapeFunctionsPointer => InitialiseShapeFunctionsHEXA20  --> to be added
          !    RearrangeConnectivitiesPointer => RearrangeConnectivitiesQUAD8
          !******************************************************************
             
          case default  ! not defined
            call GiveError('Element type not defined. [subroutine SetElementType()].')
            
      end select
          
      end subroutine InitialiseElementType
 end module ModInitialiseElementType