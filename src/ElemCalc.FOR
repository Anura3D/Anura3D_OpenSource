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
    !	Copyright (C) 2022  Members of the Anura3D MPM Research Community 
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


      module ModElementEvaluation
      !**********************************************************************
      !
      !     ModElementEvaluation:  This module contains routines related to Finite Element
      !                            evaluation independent of the type of element
      !                            (Shape functions, strain interpolation matrix,
      !                            Jacobian matrix).
      !
      !                            Whenever, element specific data is needed it is referred
      !                            to the corresponding source files:
      !                             - ModElementEvaluationTETRA
      !                             - ModElementEvaluationTRI
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************
          
      use ModElementEvaluationTETRA
      use ModElementEvaluationTRI
      use ModElementEvaluationQUAD
      use ModString
      use ModFeedback
      use ModGlobalConstants
      use ModCounters
      use ModGeometryMath
      use ModNURBS
      
      implicit none
      
        real(REAL_TYPE), dimension(:), allocatable :: GPWeight 
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative
        
        
        ! Originaly in Anura3D, GPShapeFunctionDerivative is the same for all elements. 
        ! So we do not need to calculate it for each element. But in NURBS, we have a super element
        ! so we need to know for what element we are doing this for. We can add a dimension with
        ! element number
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeight_AllElements  
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunction_AllElements
        real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivative_AllElements

        ! NURBS xi direction  
        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary ! for each boundary (1D) there are a number of weights 
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary ! similarly here for shape function evaluation (1D) 
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary ! similarly here for shape function derivative (1D) 
        
        !NURBS eta direction 
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeight_ETA 
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction_ETA
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative_ETA

        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_ETA
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_ETA
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ETA
        
        !NURBS zeta direction 
        !real(REAL_TYPE), dimension(:), allocatable :: GPWeight_ZETA 
        !real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunction_ZETA
        !real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivative_ZETA

        real(REAL_TYPE), dimension(:), allocatable :: GPWeightBoundary_ZETA
        real(REAL_TYPE), dimension(:,:), allocatable :: GPShapeFunctionBoundary_ZETA
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ZETA
        
        

        ! ------------------ all elements needed for traction implementation -------------------------------------------------

        
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_XI_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_XI_AllElements
        real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_XI_AllElements
        
        
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_ETA_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_ETA_AllElements
        real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ETA_AllElements
        
        
        real(REAL_TYPE), dimension(:,:), allocatable :: GPWeightBoundary_ZETA_AllElements
        real(REAL_TYPE), dimension(:,:,:), allocatable :: GPShapeFunctionBoundary_ZETA_AllElements
        real(REAL_TYPE), dimension(:,:,:,:), allocatable :: GPShapeFunctionDerivativeBoundary_ZETA_AllElements

        
        ! ------------------ all elements needed for traction implementation -------------------------------------------------


      contains ! Routines of this module

      
      subroutine GaussPointLocalCoordinates(IGaussPoint, WeiGP, PosGP, ELEMENTGAUSSPOINTS)
      !**********************************************************************
      !
      !>   GaussPointLocalCoordinates:  Determines the local coordinates and integration weight assigned to Gauss point
      !                                 IGaussPoint which are returned through WeiGP and PosGP.
      !
      !> IN:
      !> IGaussPoint : Number of the Gauss point
      !
      !> OUT:
      !> WeiGP : Initial weight assigned to Gauss point IGaussPoint
      !> PosGP : Initial local position of Gauss point IGaussPoint
      !
      !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          
          integer(INTEGER_TYPE), intent(in) :: ELEMENTGAUSSPOINTS

          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
        
          integer(INTEGER_TYPE) :: ID
          
          WeiGP = 0.0 ! why is this zero?  --> should be 4
          PosGP = 0.0
    
          !call Gauss_Q1Pointer(IGaussPoint, PosGP, WeiGP) ! --> it is assigned as 4 in this subroutine     
        
          !case (HEXA)
          
          ID = IGaussPoint
            
          select case(ELEMENTGAUSSPOINTS)
            case (1) ! 1 material point per element
              call InitialHEXA_MP1(ID, PosGP, WeiGP)
            case (8) ! 8 material point per element
              call InitialHEXA_MP8(ID, PosGP, WeiGP)
            case (27) ! 27 material point per element
              call InitialHEXA_MP27(ID, PosGP, WeiGP)
            case (64) ! 1 material point per element
              call InitialHEXA_MP64(ID, PosGP, WeiGP)

            case (125) ! 1 material point per element
              call InitialHEXA_MP125(ID, PosGP, WeiGP)
            
            case (216) ! 1 material point per element
              call InitialHEXA_MP216(ID, PosGP, WeiGP)
              
              
              case (343) ! 1 material point per element
              call InitialHEXA_MP343(ID, PosGP, WeiGP)
          
          end select
          
        end subroutine GaussPointLocalCoordinates
        

        integer(INTEGER_TYPE) function GetNSideNodes(IElTyp)
        !**********************************************************************
        !
        !> GetNSideNodes:  Returns the number of nodes on each element side
        !!
        !> Implemented in the frame of the MPM project.
        !!
        !> IN:
        !!
        !>    IElTyp : Number of node connectivities
        !! - = 4 : 4-noded tetrahedral element
        !! - = 10 : 10-noded tetrahedral element
        !!
        !> OUT:
        !!
        !> GetNSideNodes : Number of nodes per element side
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp !< Number of node connectivities

          GetNSideNodes = 0 
          
          select case(IElTyp)
            case(10) ! 10-noded tetrahedral element  ! for 3D only  
              GetNSideNodes = 6
            case(4) ! 4-noded tetrahedral element  ! for 3D only
              GetNSideNodes = 3
          end select
        
        end function GetNSideNodes 

        
        subroutine ShapeXiEtaT(NSideNodes, Xi, Eta, HS, DHS)
        !**********************************************************************
        !
        !    Function:  Determines the nodal shape function values HS and derivatives DHS for
        !               (Xi, Eta) for either a 3-noded or 6-noded tetrahedral element - depending
        !               on NSideNodes (either 3 or 6).
        !
        !     Xi, Eta : Local coordinates
        !
        ! O   HS : Shape function values at Xi, Eta
        ! O   DHS : Shape function derivatives at Xi, Eta
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NSideNodes
          real(REAL_TYPE), intent(in) :: Xi, Eta
          real(REAL_TYPE), dimension(NSideNodes), intent(out) :: HS
          real(REAL_TYPE), dimension(NSideNodes, 2), intent(out) :: DHS
        
          select case(NSideNodes)
            case(6) ! 6-noded triangular element
              call ShapeXiEtaTetrahedronHOE(Xi, Eta, HS, DHS)
            case(3) ! 3-noded triangular element
              call ShapeXiEtaTetrahedronLOE(Xi, Eta, HS, DHS)
          end select
        
        end subroutine ShapeXiEtaT

        
        !subroutine ShapeFunctionData(LocPos, IElTyp, ShapeValues, DShapeValues, IElement, &
        !                         ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
        !                         ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta)
        !!**********************************************************************
        !!
        !!    Function:  Evaluates the shape functions and shape function derivatives
        !!               at a local coordinate LocPos inside an element.
        !!
        !!     LocPos : Point inside an element in the local coordinate system
        !!     IElTyp : Number of nodes per element
        !!     IDim : Number of dimensions
        !!
        !! O   ShapeValues : Shape function values at LocPos
        !! O   DShapeValues : Shape function derivatives at LocPos
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  real(REAL_TYPE), dimension(:), intent(in) :: LocPos !local position is an input... how is this calculated?
        !  integer(INTEGER_TYPE), intent(in) :: IElTyp
        !  real(REAL_TYPE), dimension(:), intent(inout) :: ShapeValues
        !  real(REAL_TYPE), dimension(:, :), intent(inout) :: DShapeValues
        !  integer(INTEGER_TYPE), intent(in) :: IElement
        !  integer(INTEGER_TYPE) :: ni, nj
        !  real(REAL_TYPE), dimension(NXiGaussPoints) :: Xi_ParametricDomain
        !  real(REAL_TYPE), dimension(NEtaGaussPoints) :: Eta_ParametricDomain
        !
        !  ! NURBS
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS 
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt !classic inout parameters
        !                                    
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Xi
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Xi
        !
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta          
        !  
        !  ! NURBS Xi output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Xi
        !                         
        !  ! NURBS Eta output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Eta
        !  
        !  !ShapeValuesArray_Eta 
        !  !DShapeValuesArray_Eta
        !  
        !  
        !  
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: ShapeValues_Xi 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: DShapeValues_Xi 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: Weights_Xi                 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: ShapeValues_Eta 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: DShapeValues_Eta 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: Weights_Eta
        !  
        !  
        !  
        !  
        !  
        !  
        !  ! NURBS debugging 
        !  allocate( Wt(ELEMENTGAUSSPOINTS) )
        !  allocate( HS(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        !  allocate( dHS(ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
        !
        !  ! Xi direction 
        !  allocate( Wt_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI) )
        !  allocate( HS_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI) )
        !  allocate( dHS_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI, NDOFL-1) )
        !  
        !  allocate(ShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
        !  
        !  
        !  ! Eta direction 
        !  allocate( Wt_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA) )
        !  allocate( HS_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA) )
        !  allocate( dHS_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA, NDOFL-1) )
        !
        !  allocate(ShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
        !  
        !  
        !  
        !  
        !  ! Xi direction 
        !  
        !  
        !  
        !  
        !  
        !  ! Eta direction 
        !  
        !  
        !  
        !  
        !  
        !  
        !!Temp Solution for TetraOld           
        !  select case(ELEMENTTYPE) ! -> this is just for backward compatibility
        !  case(TETRAOLD)
        !      select case(IElTyp)
        !      case(10) ! 10-noded tetrahedral element  ! for 3D only
        !        ShapeLocPosPointer => ShapeLocPosTETRA10
        !      case(4) ! 4-noded tetrahedral element  ! for 3D only
        !        ShapeLocPosPointer => ShapeLocPosTETRA4
        !      end select
        !  end select
        ! 
        ! !Triangle and Quadilateral use this with LocPos as an input 
        ! !call ShapeLocPosPointer(LocPos, ShapeValues, DShapeValues)  ! -> NURBS implementation 
        !
        ! ni = INN(IEN(nen_NURBS,IElement),1)    
        ! nj = INN(IEN(nen_NURBS,IElement),2)
        ! 	
        ! if (NDIM == 3) then  
        !     nk = INN(IEN(nen_NURBS,IElement),3)
        ! end if 
        ! 
        !  != LocPos(1)
        !  != LocPos(2)
        ! 
        ! Xi_ParametricDomain =  ( (XiKnotEntries(ni+1) - XiKnotEntries(ni) ) * LocPos(1) &
        !                            + (XiKnotEntries(ni+1) + XiKnotEntries(ni)) ) * 0.5; ! this should be a scalar always
        !    
        ! Eta_ParametricDomain =  ( (EtaKnotEntries(nj+1) - EtaKnotEntries(nj) ) * LocPos(2) &
        !                            + (EtaKnotEntries(nj+1) + EtaKnotEntries(nj)) ) * 0.5; ! this should be a scalar always
        ! 
        ! !if (Eta_ParametricDomain(NXiGaussPoints) < EtaKnotEntries(nj)) then 
        ! !    Eta_ParametricDomain(NXiGaussPoints) = Eta_ParametricDomain(NXiGaussPoints) + EtaKnotEntries(nj)
        ! !end if 
        ! !
        ! !if (Xi_ParametricDomain(NXiGaussPoints) < XiKnotEntries(ni)) then 
        ! !    Xi_ParametricDomain(NXiGaussPoints) = Xi_ParametricDomain(NXiGaussPoints) + XiKnotEntries(ni)
        ! !end if 
        ! 
        !!InitialiseShapeFunctionsQUAD4_NURBS -> NURBS implementation 
        !call InitialiseShapeFunctionsQUAD4_NURBS(HS, dHS, Wt, & !classic inout parameters
        !                                         HS_Xi, dHS_Xi, Wt_Xi, &
        !                                         HS_Eta, dHS_Eta, Wt_Eta, &
        !                                         XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
        !                                         EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
        !                                         ni, nj) !NURBS related inputs in the eta direction 
        !! need to normalize the gradients by dividing by the element width....
        !! In this case, the element width is 0.5 in every direction. 
        !! So I will hardcode this in here by just multiplying by 0.5 
        !
        !
        !!Let's store HS_Xi and HS_Eta and see how these look like for the MPs as they cross in the column
        !! the second index needs to be particle number 
        !
        !!ShapeValuesArray_Xi(:,  IElement) = HS_Xi
        !!DShapeValuesArray_Xi(:, IElement) = dHS_Xi
        !!WeightArray_Xi(:,  IElement) = Wt_Xi
        !
        !!
        !!ShapeValuesArray_Eta = HS_Eta
        !!DShapeValuesArray_Eta = dHS_Eta
        !!Weight_EtaArray = Wt_Eta
        !
        !
        !
        !
        !  !real(REAL_TYPE), dimension(:,:), intent(out) ::  
        !  !real(REAL_TYPE), dimension(:,:), intent(out) :: 
        !  !                       
        !  !! NURBS Eta output for debugging purposes 
        !  !real(REAL_TYPE), dimension(:,:), intent(out) ::  
        !  !real(REAL_TYPE), dimension(:,:), intent(out) :: 
        !
        !
        !
        !! -----------------------------------------------------------------------------------
        !
        !      
        !! initialize 
        !!
        !!NN_IncludesZeroValues_Print = 0.0
        !!dN_dxi_IncludesZeroValues_Print = 0.0
        !!      
        !!
        !!MM_IncludesZeroValues_Print = 0.0
        !!dM_deta_IncludesZeroValues_Print = 0.0
        !!      
        !!      
        !!! write debug parameters 
        !!NN_IncludesZeroValues_Print = NN_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !!dN_dxi_IncludesZeroValues_Print = dN_dxi_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !!      
        !!      
        !!MM_IncludesZeroValues_Print = MM_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !!dM_deta_IncludesZeroValues_Print = dM_deta_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !      
        !      
        !! -----------------------------------------------------------------------------------
        !
        !
        !
        !
        !
        !! Xi 
        !ShapeValuesArray_Temp_Xi = HS_Xi(1,:)
        !DShapeValuesArray_Temp_Xi = dHS_Eta(1,:,1)
        !  
        !  
        !! Eta 
        !ShapeValuesArray_Temp_Eta = HS_Eta(1,:)
        !DShapeValuesArray_Temp_Eta = dHS_Eta(1,:,1)
        !
        !
        !ShapeValues = HS(1,:) !* 0.333 !hardcoded
        !DShapeValues = dHS(1, :, :)  !* 0.333 !hardcoded...
        !
        !                         end subroutine ShapeFunctionData
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine ShapeFunctionData(LocPos, IElTyp, ShapeValues, DShapeValues, IElement, IPatch) !, &
                                 !IElTyp, 
                                 !ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
                                 !ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta, &
                                 !ShapeValuesArray_Temp_Zeta, DShapeValuesArray_Temp_Zeta)
        !**********************************************************************
        !
        !    Function:  Evaluates the shape functions and shape function derivatives
        !               at a local coordinate LocPos inside an element.
        !
        !     LocPos : Point inside an element in the local coordinate system
        !     IElTyp : Number of nodes per element
        !     IDim : Number of dimensions
        !
        ! O   ShapeValues : Shape function values at LocPos
        ! O   DShapeValues : Shape function derivatives at LocPos
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos !local position is an input... how is this calculated?
          integer(INTEGER_TYPE), intent(in) :: IElTyp
          real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: ShapeValues
          real(REAL_TYPE), allocatable, dimension(:, :), intent(inout) :: DShapeValues
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE) :: ni, nj, nk
          real(REAL_TYPE), dimension(1) :: Xi_ParametricDomain !NXiGaussPoints
          real(REAL_TYPE), dimension(1) :: Eta_ParametricDomain !NEtaGaussPoints
          real(REAL_TYPE), dimension(1) :: Zeta_ParametricDomain !NZetaGaussPoints


          ! NURBS
          real(REAL_TYPE), allocatable, dimension(:,:) :: HS
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS 
          real(REAL_TYPE), allocatable, dimension(:) :: Wt !classic inout parameters
                                            
          real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Xi 
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Xi
          real(REAL_TYPE), allocatable, dimension(:) :: Wt_Xi

          real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
          real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta          
          
          real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Zeta 
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Zeta
          real(REAL_TYPE), allocatable, dimension(:) :: Wt_Zeta    
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary = 1
          integer(INTEGER_TYPE), intent(in) :: IPatch
          
          
          ! Temporary debugging variables --> commented because no longer needed but kept just in case.
          ! NURBS Xi output for debugging purposes 
          !real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Xi !, intent(out)
          !real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Xi !, intent(out)
          !                       
          !! NURBS Eta output for debugging purposes 
          !real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Eta !, intent(out)
          !real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Eta !, intent(out)
          !
          !! NURBS Zeta output for debugging purposes 
          !real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Zeta !, intent(out)
          !real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Zeta !, intent(out) 
          
          
          ! NURBS debugging 
          allocate( Wt(1) ) !ELEMENTGAUSSPOINTS --> This is for 1 material point
          allocate( HS(1, ELEMENTNODES) ) !ELEMENTGAUSSPOINTS --> This is for 1 material point
          allocate( dHS(1, ELEMENTNODES, NDIM) ) !ELEMENTGAUSSPOINTS, NDOFL --> This is for 1 material point
        
          ! Xi direction 
          allocate( Wt_Xi(1) ) !ELEMENTBOUNDARYGAUSSPOINTS_XI --> This is for 1 material point
          allocate( HS_Xi(1, ELEMENTBOUNDARYNODES_XI) ) !ELEMENTBOUNDARYGAUSSPOINTS_XI --> This is for 1 material point
          allocate( dHS_Xi(1, ELEMENTBOUNDARYNODES_XI, 1)) !NDOFL-1) ) !ELEMENTBOUNDARYGAUSSPOINTS_XI --> This is for 1 material point

          ! Temporary debugging variables --> commented because no longer needed but kept just in case.
          !allocate(ShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
          !allocate(DShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
          
          
          ! Eta direction 
          allocate( Wt_Eta(1) ) !ELEMENTBOUNDARYGAUSSPOINTS_ETA
          allocate( HS_Eta(1, ELEMENTBOUNDARYNODES_ETA) ) !ELEMENTBOUNDARYGAUSSPOINTS_ETA
          allocate( dHS_Eta(1, ELEMENTBOUNDARYNODES_ETA, 1)) !NDOFL-1) ) !ELEMENTBOUNDARYGAUSSPOINTS_ETA
        

          ! Temporary debugging variables --> commented because no longer needed but kept just in case.
          !allocate(ShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
          !allocate(DShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
          
          
          ! Zeta direction
          allocate( Wt_Zeta(1) ) !ELEMENTBOUNDARYGAUSSPOINTS_ZETA
          allocate( HS_Zeta(1, ELEMENTBOUNDARYNODES_ZETA) ) 
          allocate( dHS_Zeta(1, ELEMENTBOUNDARYNODES_ZETA, 1)) !NDOFL-1) ) !ELEMENTBOUNDARYGAUSSPOINTS_ZETA 
        

          ! Temporary debugging variables --> commented because no longer needed but kept just in case.
          !allocate(ShapeValuesArray_Temp_Zeta(NZetaKnotOrder+1) )
          !allocate(DShapeValuesArray_Temp_Zeta(NZetaKnotOrder+1) )
          
          !NNodes = size(ICon,1)
          !if (allocated(ShapeValues)) then
          !    deallocate(ShapeValues)
          !    
          !end if 
          !    
          !
          !allocate( ShapeValues(IElTyp) )
          !ShapeValues = 0
          !
          !
          !
          !if (allocated(DShapeValues)) then
          !    deallocate(DShapeValues)
          !    
          !end if 
          !    
          !
          !allocate( DShapeValues(IElTyp, NDIM) )
          !DShapeValues = 0
          
          
          
          
          
          
          
        !Temp Solution for TetraOld           
          select case(ELEMENTTYPE) ! -> this is just for backward compatibility
          case(TETRAOLD)
              select case(IElTyp)
              case(10) ! 10-noded tetrahedral element  ! for 3D only
                ShapeLocPosPointer => ShapeLocPosTETRA10
              case(4) ! 4-noded tetrahedral element  ! for 3D only
                ShapeLocPosPointer => ShapeLocPosTETRA4
              end select
          end select
         
         !Triangle and Quadilateral use this with LocPos as an input 
         !call ShapeLocPosPointer(LocPos, ShapeValues, DShapeValues)  ! -> NURBS implementation 
         ni = INN(IEN(1,IElement,IPatch),1,IPatch)    
         nj = INN(IEN(1,IElement,IPatch),2,IPatch)
         if (NDIM == 3) then 
             nk = INN(IEN(1,IElement,IPatch),3,IPatch)
         end if 
         
         ! calculate parametric domain values 
         Xi_ParametricDomain =  ( (XiKnotEntries(ni+1,IPatch) - XiKnotEntries(ni,IPatch) ) * LocPos(1) &
                                    + (XiKnotEntries(ni+1,IPatch) + XiKnotEntries(ni,IPatch)) ) * 0.5 ! this should be a scalar always
         
         Eta_ParametricDomain =  ( (EtaKnotEntries(nj+1,IPatch) - EtaKnotEntries(nj,IPatch) ) * LocPos(2) &
                                    + (EtaKnotEntries(nj+1,IPatch) + EtaKnotEntries(nj,IPatch)) ) * 0.5 ! this should be a scalar always
         if (NDIM == 3) then 
         Zeta_ParametricDomain =  ( (ZetaKnotEntries(nk+1,IPatch) - ZetaKnotEntries(nk,IPatch) ) * LocPos(3) &
                                    + (ZetaKnotEntries(nk+1,IPatch) + ZetaKnotEntries(nk,IPatch)) ) * 0.5 ! this should be a scalar always
         end if
        
         if (NDIM == 2) then 
         
         call InitialiseShapeFunctionsQUAD4_NURBS(HS, dHS, Wt, & !classic inout parameters
                                                 HS_Xi, dHS_Xi, Wt_Xi, &
                                                 HS_Eta, dHS_Eta, Wt_Eta, &
                                                 XiKnotEntries(:,IPatch), NXiKnotEntries(IPatch), Xi_ParametricDomain, NXiKnotOrder(IPatch), & !NURBS related inputs in the xi direction 
                                                 EtaKnotEntries(:,IPatch), NEtaKnotEntries(IPatch), Eta_ParametricDomain, NEtaKnotOrder(IPatch), &
                                                 ni, nj, &
                                                 IPatch) !NURBS related inputs in the eta direction 
         
         elseif (NDIM == 3) then 
             
         call InitialiseShapeFunctionsHEXA_NURBS(HS, dHS, Wt, & !classic inout parameters
                                                 HS_Xi, dHS_Xi, Wt_Xi, &
                                                 HS_Eta, dHS_Eta, Wt_Eta, &
                                                 HS_Zeta, dHS_Zeta, Wt_Zeta, &
                                                 XiKnotEntries(:,IPatch), NXiKnotEntries(IPatch), Xi_ParametricDomain, NXiKnotOrder(IPatch), & !NURBS related inputs in the xi direction 
                                                 EtaKnotEntries(:,IPatch), NEtaKnotEntries(IPatch), Eta_ParametricDomain, NEtaKnotOrder(IPatch), &
                                                 ZetaKnotEntries(:,IPatch), NZetaKnotEntries(IPatch), Zeta_ParametricDomain, NZetaKnotOrder(IPatch), &
                                                 ni, nj, nk, &
                                                 IElement, &
                                                 IPatch) !NURBS related inputs in the eta direction     
             
         end if 
         

        
        ! -----------------------------------------------------------------------------------
        ! Abandoned debugging variables 
              
        ! initialize 
        !
        !NN_IncludesZeroValues_Print = 0.0
        !dN_dxi_IncludesZeroValues_Print = 0.0
        !      
        !
        !MM_IncludesZeroValues_Print = 0.0
        !dM_deta_IncludesZeroValues_Print = 0.0
        !      
        !      
        !! write debug parameters 
        !NN_IncludesZeroValues_Print = NN_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !dN_dxi_IncludesZeroValues_Print = dN_dxi_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !      
        !      
        !MM_IncludesZeroValues_Print = MM_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !dM_deta_IncludesZeroValues_Print = dM_deta_IncludesZeroValues(1,:, NEtaKnotOrder+1)
              
        ! Xi 
        !ShapeValuesArray_Temp_Xi = HS_Xi(1,:)
        !DShapeValuesArray_Temp_Xi = dHS_Eta(1,:,1)
        !  
        !  
        !! Eta 
        !ShapeValuesArray_Temp_Eta = HS_Eta(1,:)
        !DShapeValuesArray_Temp_Eta = dHS_Eta(1,:,1)
        !
        !! Zeta 
        !ShapeValuesArray_Temp_Zeta = HS_Zeta(1,:)
        !DShapeValuesArray_Temp_Zeta = dHS_Zeta(1,:,1)      
        ! -----------------------------------------------------------------------------------
        
         
         ! assign outputs 
         ShapeValues = HS(1,:) !--> RR
         DShapeValues = dHS(1, :, :)  ! --> dR/dxi
        
                                 
                                 
        end subroutine ShapeFunctionData
                                 
                                 
                                 
                                 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                                 
                                 
        ! ---------------------------------------------------------------------------------------------
        ! SINGLE PARTICLE 
        ! subroutine ShapeFunctionData_SingleParticle(LocPos, IElTyp, ShapeValues, DShapeValues, IElement, &
        !                         ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
        !                         ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta)
        !!**********************************************************************
        !!
        !!    Function:  Evaluates the shape functions and shape function derivatives
        !!               at a local coordinate LocPos inside an element.
        !!
        !!     LocPos : Point inside an element in the local coordinate system
        !!     IElTyp : Number of nodes per element
        !!     IDim : Number of dimensions
        !!
        !! O   ShapeValues : Shape function values at LocPos
        !! O   DShapeValues : Shape function derivatives at LocPos
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  real(REAL_TYPE), dimension(:), intent(in) :: LocPos !local position is an input... how is this calculated?
        !  integer(INTEGER_TYPE), intent(in) :: IElTyp
        !  real(REAL_TYPE), dimension(:), intent(inout) :: ShapeValues
        !  real(REAL_TYPE), dimension(:, :), intent(inout) :: DShapeValues
        !  integer(INTEGER_TYPE), intent(in) :: IElement
        !  integer(INTEGER_TYPE) :: ni, nj
        !  real(REAL_TYPE), dimension(1) :: Xi_ParametricDomain
        !  real(REAL_TYPE), dimension(1) :: Eta_ParametricDomain
        !
        !  ! NURBS
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS 
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt !classic inout parameters
        !                                    
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Xi
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Xi
        !
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta          
        !  
        !  ! NURBS Xi output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Xi
        !                         
        !  ! NURBS Eta output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Eta
        !  
        !  !ShapeValuesArray_Eta 
        !  !DShapeValuesArray_Eta
        !  
        !  
        !  
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: ShapeValues_Xi 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: DShapeValues_Xi 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: Weights_Xi                 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: ShapeValues_Eta 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: DShapeValues_Eta 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: Weights_Eta
        !  
        !  
        !  
        !  
        !  
        !  
        !  ! NURBS debugging 
        !  allocate( Wt(ELEMENTGAUSSPOINTS) )
        !  allocate( HS(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        !  allocate( dHS(ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
        !
        !  ! Xi direction 
        !  allocate( Wt_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI) )
        !  allocate( HS_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI) )
        !  allocate( dHS_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI, NDOFL-1) )
        !  
        !  allocate(ShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
        !  
        !  
        !  ! Eta direction 
        !  allocate( Wt_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA) )
        !  allocate( HS_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA) )
        !  allocate( dHS_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA, NDOFL-1) )
        !
        !  allocate(ShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
        !  
        !  
        !  
        !  
        !  ! Xi direction 
        !  
        !  
        !  
        !  
        !  
        !  ! Eta direction 
        !  
        !  
        !  
        !  
        !  
        !  
        !!Temp Solution for TetraOld           
        !  select case(ELEMENTTYPE) ! -> this is just for backward compatibility
        !  case(TETRAOLD)
        !      select case(IElTyp)
        !      case(10) ! 10-noded tetrahedral element  ! for 3D only
        !        ShapeLocPosPointer => ShapeLocPosTETRA10
        !      case(4) ! 4-noded tetrahedral element  ! for 3D only
        !        ShapeLocPosPointer => ShapeLocPosTETRA4
        !      end select
        !  end select
        ! 
        ! !Triangle and Quadilateral use this with LocPos as an input 
        ! !call ShapeLocPosPointer(LocPos, ShapeValues, DShapeValues)  ! -> NURBS implementation 
        !
        !ni = INN(IEN(nen_NURBS,IElement),1)    
        !nj = INN(IEN(nen_NURBS,IElement),2)
        !
        !	
        !if (NDIM == 3) then 
        !  nk = INN(IEN(nen_NURBS,IElement),3)
        !end if 
        ! 
        !  != LocPos(1)
        !  != LocPos(2)
        ! 
        ! Xi_ParametricDomain =  ( (XiKnotEntries(ni+1) - XiKnotEntries(ni) ) * LocPos(1) &
        !                            + (XiKnotEntries(ni+1) + XiKnotEntries(ni)) ) * 0.5; ! this should be a scalar always
        !    
        ! Eta_ParametricDomain =  ( (EtaKnotEntries(nj+1) - EtaKnotEntries(nj) ) * LocPos(2) &
        !                            + (EtaKnotEntries(nj+1) + EtaKnotEntries(nj)) ) * 0.5; ! this should be a scalar always
        ! 
        ! !if (Eta_ParametricDomain(NXiGaussPoints) < EtaKnotEntries(nj)) then 
        ! !    Eta_ParametricDomain(NXiGaussPoints) = Eta_ParametricDomain(NXiGaussPoints) + EtaKnotEntries(nj)
        ! !end if 
        ! !
        ! !if (Xi_ParametricDomain(NXiGaussPoints) < XiKnotEntries(ni)) then 
        ! !    Xi_ParametricDomain(NXiGaussPoints) = Xi_ParametricDomain(NXiGaussPoints) + XiKnotEntries(ni)
        ! !end if 
        ! 
        !!InitialiseShapeFunctionsQUAD4_NURBS -> NURBS implementation 
        !call InitialiseShapeFunctionsQUAD4_NURBS_SINGLEPARTICLE(HS, dHS, Wt, & !classic inout parameters
        !                                         HS_Xi, dHS_Xi, Wt_Xi, &
        !                                         HS_Eta, dHS_Eta, Wt_Eta, &
        !                                         XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
        !                                         EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
        !                                         ni, nj) !NURBS related inputs in the eta direction 
        !! need to normalize the gradients by dividing by the element width....
        !! In this case, the element width is 0.5 in every direction. 
        !! So I will hardcode this in here by just multiplying by 0.5 
        !
        !
        !!Let's store HS_Xi and HS_Eta and see how these look like for the MPs as they cross in the column
        !! the second index needs to be particle number 
        !
        !!ShapeValuesArray_Xi(:,  IElement) = HS_Xi
        !!DShapeValuesArray_Xi(:, IElement) = dHS_Xi
        !!WeightArray_Xi(:,  IElement) = Wt_Xi
        !
        !!
        !!ShapeValuesArray_Eta = HS_Eta
        !!DShapeValuesArray_Eta = dHS_Eta
        !!Weight_EtaArray = Wt_Eta
        !
        !
        !
        !
        !  !real(REAL_TYPE), dimension(:,:), intent(out) ::  
        !  !real(REAL_TYPE), dimension(:,:), intent(out) :: 
        !  !                       
        !  !! NURBS Eta output for debugging purposes 
        !  !real(REAL_TYPE), dimension(:,:), intent(out) ::  
        !  !real(REAL_TYPE), dimension(:,:), intent(out) :: 
        !
        !
        !
        !! -----------------------------------------------------------------------------------
        !
        !      
        !! initialize 
        !!
        !!NN_IncludesZeroValues_Print = 0.0
        !!dN_dxi_IncludesZeroValues_Print = 0.0
        !!      
        !!
        !!MM_IncludesZeroValues_Print = 0.0
        !!dM_deta_IncludesZeroValues_Print = 0.0
        !!      
        !!      
        !!! write debug parameters 
        !!NN_IncludesZeroValues_Print = NN_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !!dN_dxi_IncludesZeroValues_Print = dN_dxi_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !!      
        !!      
        !!MM_IncludesZeroValues_Print = MM_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !!dM_deta_IncludesZeroValues_Print = dM_deta_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !      
        !      
        !! -----------------------------------------------------------------------------------
        !
        !
        !
        !
        !
        !! Xi 
        !ShapeValuesArray_Temp_Xi = HS_Xi(1,:)
        !DShapeValuesArray_Temp_Xi = dHS_Eta(1,:,1)
        !  
        !  
        !! Eta 
        !ShapeValuesArray_Temp_Eta = HS_Eta(1,:)
        !DShapeValuesArray_Temp_Eta = dHS_Eta(1,:,1)
        !
        !
        !ShapeValues = HS(1,:) !* 0.333 !hardcoded
        !DShapeValues = dHS(1, :, :)  !* 0.333 !hardcoded...
        !
        !end subroutine ShapeFunctionData_SingleParticle
                                 
                
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
        !    subroutine ShapeFunctionData_SingleParticle(LocPos, IElTyp, ShapeValues, DShapeValues, IElement, &
        !                         ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
        !                         ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta, &
        !                         ShapeValuesArray_Temp_Zeta, DShapeValuesArray_Temp_Zeta)
        !!**********************************************************************
        !!
        !!    Function:  Evaluates the shape functions and shape function derivatives
        !!               at a local coordinate LocPos inside an element.
        !!
        !!     LocPos : Point inside an element in the local coordinate system
        !!     IElTyp : Number of nodes per element
        !!     IDim : Number of dimensions
        !!
        !! O   ShapeValues : Shape function values at LocPos
        !! O   DShapeValues : Shape function derivatives at LocPos
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  real(REAL_TYPE), dimension(:), intent(in) :: LocPos !local position is an input... how is this calculated?
        !  integer(INTEGER_TYPE), intent(in) :: IElTyp
        !  real(REAL_TYPE), dimension(:), intent(inout) :: ShapeValues
        !  real(REAL_TYPE), dimension(:, :), intent(inout) :: DShapeValues
        !  integer(INTEGER_TYPE), intent(in) :: IElement
        !  integer(INTEGER_TYPE) :: ni, nj, nk
        !  
        !  
        !  real(REAL_TYPE), dimension(1) :: Xi_ParametricDomain
        !  real(REAL_TYPE), dimension(1) :: Eta_ParametricDomain
        !  real(REAL_TYPE), dimension(1) :: Zeta_ParametricDomain
        !  
        !
        !  ! NURBS
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS 
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt !classic inout parameters
        !                                    
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Xi
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Xi
        !
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta          
        !  
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Zeta 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Zeta
        !  real(REAL_TYPE), allocatable, dimension(:) :: Wt_Zeta
        !  !
        !  !real(REAL_TYPE), allocatable, dimension(:,:) :: HS_Eta 
        !  !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dHS_Eta
        !  !real(REAL_TYPE), allocatable, dimension(:) :: Wt_Eta         
        !  
        !  
        !  
        !  ! NURBS Xi output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Xi 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Xi
        !                         
        !  ! NURBS Eta output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Eta 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Eta
        !  
        !  ! NURBS Eta output for debugging purposes 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: ShapeValuesArray_Temp_Zeta 
        !  real(REAL_TYPE), allocatable, dimension(:), intent(out) :: DShapeValuesArray_Temp_Zeta
        !  
        !  
        !  !ShapeValuesArray_Eta 
        !  !DShapeValuesArray_Eta
        !  
        !  
        !  
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: ShapeValues_Xi 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: DShapeValues_Xi 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: Weights_Xi                 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: ShapeValues_Eta 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: DShapeValues_Eta 
        !  !real(REAL_TYPE), dimension(:,:), intent(in) :: Weights_Eta
        !  
        !  
        !  
        !  
        !  
        !  
        !  ! NURBS debugging 
        !  allocate( Wt(ELEMENTGAUSSPOINTS) )
        !  allocate( HS(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        !  allocate( dHS(ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
        !
        !  ! Xi direction 
        !  allocate( Wt_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI) )
        !  allocate( HS_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI) )
        !  allocate( dHS_Xi(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI, 1))!NDOFL-1) )
        !  
        !  allocate(ShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Xi(NXiKnotOrder+1) )
        !  
        !  
        !  ! Eta direction 
        !  allocate( Wt_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA) )
        !  allocate( HS_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA) )
        !  allocate( dHS_Eta(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA, 1))!NDOFL-1) )
        !
        !  allocate(ShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Eta(NEtaKnotOrder+1) )
        !  
        !  
        !  ! Zeta direction 
        !  allocate( Wt_Zeta(ELEMENTBOUNDARYGAUSSPOINTS_ZETA) )
        !  allocate( HS_Zeta(ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA) )
        !  allocate( dHS_Zeta(ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA, 1))!NDOFL-1) )
        !
        !  allocate(ShapeValuesArray_Temp_Zeta(NZetaKnotOrder+1) )
        !  allocate(DShapeValuesArray_Temp_Zeta(NZetaKnotOrder+1) )
        !  
        !  
        !  ! Xi direction 
        !  
        !  
        !  
        !  
        !  
        !  ! Eta direction 
        !  
        !  
        !  
        !  
        !  
        !  
        !!Temp Solution for TetraOld           
        !  select case(ELEMENTTYPE) ! -> this is just for backward compatibility
        !  case(TETRAOLD)
        !      select case(IElTyp)
        !      case(10) ! 10-noded tetrahedral element  ! for 3D only
        !        ShapeLocPosPointer => ShapeLocPosTETRA10
        !      case(4) ! 4-noded tetrahedral element  ! for 3D only
        !        ShapeLocPosPointer => ShapeLocPosTETRA4
        !      end select
        !  end select
        ! 
        ! !Triangle and Quadilateral use this with LocPos as an input 
        ! !call ShapeLocPosPointer(LocPos, ShapeValues, DShapeValues)  ! -> NURBS implementation 
        !
        ! ni = INN(IEN(nen_NURBS,IElement),1)    
        ! nj = INN(IEN(nen_NURBS,IElement),2)
        ! 
        ! if (NDIM == 3) then 
        !     nk = INN(IEN(1,IElement),3)
        ! end if 
        !  != LocPos(1)
        !  != LocPos(2)
        ! 
        ! Xi_ParametricDomain =  ( (XiKnotEntries(ni+1) - XiKnotEntries(ni) ) * LocPos(1) &
        !                            + (XiKnotEntries(ni+1) + XiKnotEntries(ni)) ) * 0.5; ! this should be a scalar always
        !    
        ! Eta_ParametricDomain =  ( (EtaKnotEntries(nj+1) - EtaKnotEntries(nj) ) * LocPos(2) &
        !                            + (EtaKnotEntries(nj+1) + EtaKnotEntries(nj)) ) * 0.5; ! this should be a scalar always
        ! 
        ! if (NDIM == 3) then 
        ! 
        ! Zeta_ParametricDomain =  ( (ZetaKnotEntries(nk+1) - ZetaKnotEntries(nk) ) * LocPos(3) &
        !                            + (ZetaKnotEntries(nk+1) + ZetaKnotEntries(nk)) ) * 0.5; ! this should be a scalar always
        !
        ! end if 
        ! 
        ! !if (Eta_ParametricDomain(NXiGaussPoints) < EtaKnotEntries(nj)) then 
        ! !    Eta_ParametricDomain(NXiGaussPoints) = Eta_ParametricDomain(NXiGaussPoints) + EtaKnotEntries(nj)
        ! !end if 
        ! !
        ! !if (Xi_ParametricDomain(NXiGaussPoints) < XiKnotEntries(ni)) then 
        ! !    Xi_ParametricDomain(NXiGaussPoints) = Xi_ParametricDomain(NXiGaussPoints) + XiKnotEntries(ni)
        ! !end if 
        ! 
        !!InitialiseShapeFunctionsQUAD4_NURBS -> NURBS implementation 
        !
        ! ! remember that for each element we have a list of connectivity that comes from the IEN array
        ! 
        ! 
        ! 
        ! 
        ! if (NDIM == 2) then 
        ! 
        !     call InitialiseShapeFunctionsQUAD4_NURBS_SINGLEPARTICLE(HS, dHS, Wt, & !classic inout parameters
        !                                         HS_Xi, dHS_Xi, Wt_Xi, &
        !                                         HS_Eta, dHS_Eta, Wt_Eta, &
        !                                         XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
        !                                         EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
        !                                         ni, nj) !NURBS related inputs in the eta direction 
        !
        ! elseif (NDIM == 3) then 
        ! 
        !     call InitialiseShapeFunctionsHEXA_NURBS_SINGLEPARTICLE(HS, dHS, Wt, & !classic inout parameters
        !                                         HS_Xi, dHS_Xi, Wt_Xi, &
        !                                         HS_Eta, dHS_Eta, Wt_Eta, &
        !                                         HS_Zeta, dHS_Zeta, Wt_Zeta, &
        !                                         XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
        !                                         EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
        !                                         ZetaKnotEntries, NZetaKnotEntries, Zeta_ParametricDomain, NZetaKnotOrder, &
        !                                         ni, nj, nk, &  
        !                                         IElement) !NURBS related inputs in the eta direction     
        !     
        ! end if 
        ! 
        ! ! need to normalize the gradients by dividing by the element width....
        !! In this case, the element width is 0.5 in every direction. 
        !! So I will hardcode this in here by just multiplying by 0.5 
        !
        !
        !!Let's store HS_Xi and HS_Eta and see how these look like for the MPs as they cross in the column
        !! the second index needs to be particle number 
        !
        !!ShapeValuesArray_Xi(:,  IElement) = HS_Xi
        !!DShapeValuesArray_Xi(:, IElement) = dHS_Xi
        !!WeightArray_Xi(:,  IElement) = Wt_Xi
        !
        !!
        !!ShapeValuesArray_Eta = HS_Eta
        !!DShapeValuesArray_Eta = dHS_Eta
        !!Weight_EtaArray = Wt_Eta
        !
        !
        !
        !
        !  !real(REAL_TYPE), dimension(:,:), intent(out) ::  
        !  !real(REAL_TYPE), dimension(:,:), intent(out) :: 
        !  !                       
        !  !! NURBS Eta output for debugging purposes 
        !  !real(REAL_TYPE), dimension(:,:), intent(out) ::  
        !  !real(REAL_TYPE), dimension(:,:), intent(out) :: 
        !
        !
        !
        !! -----------------------------------------------------------------------------------
        !
        !      
        !! initialize 
        !!
        !!NN_IncludesZeroValues_Print = 0.0
        !!dN_dxi_IncludesZeroValues_Print = 0.0
        !!      
        !!
        !!MM_IncludesZeroValues_Print = 0.0
        !!dM_deta_IncludesZeroValues_Print = 0.0
        !!      
        !!      
        !!! write debug parameters 
        !!NN_IncludesZeroValues_Print = NN_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !!dN_dxi_IncludesZeroValues_Print = dN_dxi_IncludesZeroValues(1,:, NXiKnotOrder+1)
        !!      
        !!      
        !!MM_IncludesZeroValues_Print = MM_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !!dM_deta_IncludesZeroValues_Print = dM_deta_IncludesZeroValues(1,:, NEtaKnotOrder+1)
        !      
        !      
        !! -----------------------------------------------------------------------------------
        !
        !
        !
        !
        !
        !! Xi 
        !ShapeValuesArray_Temp_Xi = HS_Xi(1,:)
        !DShapeValuesArray_Temp_Xi = dHS_Eta(1,:,1)
        !  
        !  
        !! Eta 
        !ShapeValuesArray_Temp_Eta = HS_Eta(1,:)
        !DShapeValuesArray_Temp_Eta = dHS_Eta(1,:,1)
        !
        !
        !! Zeta 
        !ShapeValuesArray_Temp_Zeta = HS_Zeta(1,:)
        !DShapeValuesArray_Temp_Zeta = dHS_Zeta(1,:,1)
        !
        !ShapeValues = HS(1,:) !* 0.333 !hardcoded
        !DShapeValues = dHS(1, :, :)  !* 0.333 !hardcoded...
        !
        !end subroutine ShapeFunctionData_SingleParticle
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
        ! SINGLE PARTICLE                          
        ! ---------------------------------------------------------------------------------------------------------                         
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 

        
        subroutine DetJacob(HS, dHS, LocPos, &
                            NEl, NodTot, IDim, IElement, &
                            ICon, Co, &
                            dx_dxi, &
                            dx_dtildexi, dtildexi_dx, Det_dx_dtildexi, & !, &
                            dxi_dx, dtildexi_dxi, &
                            IPatch)
        !**********************************************************************
        !
        !  Function : Determination of the Jacobian matrix and the determinant of
        !             the Jacobian matrix for the location LocPos inside an element IElement
        !             whose nodal connectivities are defined by ICon and Co depending on
        !             the type of element.
        !
        !  I  LocPos : Local coordinates of the considered point inside an element
        !  I  NEl : Total number of elements
        !  I  NodTot : Total number of nodes
        !  I  IDim : Number of dimensions
        !  I  IElement : ID of the element
        !  I  ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !  I  Co : Global nodal coordinates Co(I, J): j-coordinate of node I
        !
        !  O  RJac : Jacobian matrix
        !  O  InvRJac : Inverse of the Jacobian matrix
        !  O  DetJac : Determinate of the Jacobian matrix
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: NEl, NodTot, IDim, IElement
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:, :), intent(in) :: Co
          real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dx_dtildexi !RJac
          real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dtildexi_dx !InvRJac
          real(REAL_TYPE), intent(out) :: Det_dx_dtildexi
          
          real(REAL_TYPE), dimension(:), intent(in) :: HS
          real(REAL_TYPE), dimension(:,:), intent(in) :: dHS 
          
          ! local variables
          integer(INTEGER_TYPE) :: I, J, INode, NodeID, NNodes
          !real(REAL_TYPE), dimension(:), allocatable :: HS
          !real(REAL_TYPE), dimension(:, :), allocatable :: dHS
          real(REAL_TYPE) :: Det1
          
          ! Hardcoded : nen_NURBS should be patch specific but here we hardcode to 1 
          real(REAL_TYPE), dimension(nen_NURBS((IPatch)),NDIM) :: dR_dx ! IPatch_Temporary
          real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
          real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dtildexi_dxi

          real(REAL_TYPE), dimension(NDIM,NDIM) :: J_mat
          
          real(REAL_TYPE) :: Det_dxi_dtildexi, Det_dtildexi_dxi
          
          integer(INTEGER_TYPE) :: aa, bb, cc, loc_num, ii, jj
          integer(INTEGER_TYPE) :: ni, nj, nk
                    
          !real(REAL_TYPE), dimension(NDIM,NDIM) :: dx_dtildexi
          !real(REAL_TYPE), dimension(NDIM,NDIM) :: dtildexi_dx
          !real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dx_dxi
          real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dxi_dx
          real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dx_dxi
          
          real(REAL_TYPE) :: det_dx_dxi, det_dxi_dx, det_dtildexi_dx ! det_dx_dtildexi, 
          
          ! NURBS debugging
          !real(REAL_TYPE), dimension(:), allocatable :: ShapeValuesArray_Temp_Xi
          !real(REAL_TYPE), dimension(:), allocatable :: DShapeValuesArray_Temp_Xi
          !
          !real(REAL_TYPE), dimension(:), allocatable :: ShapeValuesArray_Temp_Eta
          !real(REAL_TYPE), dimension(:), allocatable :: DShapeValuesArray_Temp_Eta
          !
          !
          !real(REAL_TYPE), dimension(:), allocatable :: ShapeValuesArray_Temp_Zeta
          !real(REAL_TYPE), dimension(:), allocatable :: DShapeValuesArray_Temp_Zeta
          
          !aa, bb, loc_num, ii, jj
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          !integer(INTEGER_TYPE) :: aa

          
          NNodes = size(ICon,1)
          !allocate( HS(NNodes), dHS(NNodes, IDim) )

          
          !this subroutine asks for the element ID to calcuate the Jacobian. This might need to be cast in NURBS parametric domain using ni, nj, nk...
          
          !HS = 0.0
          !dHS = 0.0
          !
          !! get dR_dxi 
          !! Determine the shape functions HS and shape function derivatives dHS for the element for LocPos.
          !call ShapeFunctionData(LocPos, NNodes, HS, dHS, IElement, &
          !                       ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
          !                       ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta, &
          !                       ShapeValuesArray_Temp_Zeta, DShapeValuesArray_Temp_Zeta) ! -> dHS here is dR_dxi

          ! we would need to divide by denominators to complete the definitions of functions and derivatives w.r.t. parametric coodinates 
          
          !call InitialiseShapeFunctionsQUAD4_NURBS
          
          ! finding the NURBS coordinates for the element
          ! Think about this to minimize the times we have to evaluate ni, nj, and nk
          ni = INN(IEN(1,IElement,IPatch),1,IPatch)    
          nj = INN(IEN(1,IElement,IPatch),2,IPatch)
          if (NDIM == 3) then 
              nk = INN(IEN(1,IElement,IPatch),3,IPatch)
          end if 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dx/dxi
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call Evaluate_dx_dxi(dHS, Co, ICon, dx_dxi, IElement)
          ! Determine the Jacobian matrix RJac
          !RJac = 0.0
          !dx_dxi = 0.0
          !dxi_dx = 0.0
          !do INode = 1, NNodes ! loop nodes of each element
          !  
          !    !ni = INN(IEN(1,IElement),1)    
          !    !nj = INN(IEN(2,IElement),2)
          !    
          !    !XiKnotEntries(ni) 
          !    !note that the coordinates of the control points are the ones we need here in 'Co'
          !    NodeID = ICon(INode, IElement)
          !    ! This is gradient mapping from the parameter space to the physical space 
          !  do I = 1, IDim
          !    do J = 1, IDim
          !      dx_dxi(I, J) = dx_dxi(I, J) + dHS(INode, I) * Co(NodeID, J) !calculation of the jacobian matrix from the coordinates of one element ! -> this is checked and now it is working
          !    end do
          !  end do
          !end do
          
        !  loc_num = 0
        !  
        !  !do kk = 0, rr ! loop accross the order of the z dimension 
        !      do jj = 0, NEtaKnotOrder ! loop accross the order of the y dimension 
        !          do ii = 0, NXiKnotOrder ! loop accross the order of x dimension 
        !              loc_num = loc_num + 1 
        !              
        !              do aa = 1, NDIM 
        !                  do bb = 1, NDIM 
        !                      
        !                      RJac(aa, bb) = RJac(aa, bb) + Co(ni - ii, nj - jj) * dHS(loc_num, bb)
        !                  end do 
        !              end do 
        !          end do 
        !      end do 
        !    !end do 
        !!end do 
        
        
                      
                      
                      
                      !XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
                      !                           EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
                      !                           ni, nj) 
                      
                      
                      
                      
                      
              
          
          !if (RJac(1,1) == 0.0) then 
          !
          !    return 
          !
          !end if 
            
          ! RJac -> dx_dxi in Cottrell et al. (2011)
           
          !if (RJac(1,1) == 0.0) then 
          !
          !    return 
          !
          !end if 
          !
          !RJac(1,1) = 1.0
          !RJac(1,2) = 0.0
          !RJac(2,1) = 0.0
          !RJac(2,2) = 1.0
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dxi/dx
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculate inverse and determinant of Jacobian matrix
          
          ! debugger
          if (dx_dtildexi(3,3) == 0.0) then 
              aa=1
          end if

          
          call RJacInv(IDim, dx_dxi, dxi_dx, Det_dx_dxi, Det_dxi_dx)

          
         call Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries(IPatch), XiKnotEntries(:,IPatch), ni, &
                                                                        NEtaKnotEntries(IPatch), EtaKnotEntries(:,IPatch), nj, &  
                                                                        NZetaKnotEntries(IPatch), ZetaKnotEntries(:,IPatch), nk, &       
                                                                        dxi_dtildexi)
                  
         call Get_dx_dtildexi_MappingFromParentElementToParameterSpace(dx_dxi, dxi_dtildexi, dx_dtildexi) ! -> RJac = dx/dxi * dxi/dtildexi
         
         ! debugger
      if (dxi_dtildexi(3,3) == 0.0) then 
          aa=1
      end if

         
         call RJacInv(IDim, dxi_dtildexi, dtildexi_dxi, Det_dxi_dtildexi, Det_dtildexi_dxi)
         
         
         if (dx_dtildexi(3,3) == 0.0) then 
          aa=1
      end if
         call RJacInv(NDIM, dx_dtildexi, dtildexi_dx, Det_dx_dtildexi, Det_dtildexi_dx)
         
         
         ! assign the determinant to dx_dtildexi
         !RJac = dx_dtildexi
         !InvRJac = dtildexi_dx
         !DetJac = Det_dx_dtildexi
         !Det1 = Det_dtildexi_dx
          
          !_______________________________
          
          
          
          
          
          ! InvRJac -> dxi_dx 
          !
          !
          ! EVERYTHING BELOW HERE SHOULD GO IN THE PARENT SUBROUTINE 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dR/dx
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !Calculate derivative of basis functions with respect to physical coordinates
          !call Get_dR_dx_MappingShapeFunctionDerivativeWRTPhysicalCoordinates(InvRJac, dHS, dR_dx) ! dR_dx = dxi_dx * dR_dxi 
          
          !dHS = dR_dx
          
          ! -> We would need to store with dR_dx in the shape function array 
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Calculating dxi/dtildexi
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !call Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries, XiKnotEntries, ni, &
          !                                                               NEtaKnotEntries, EtaKnotEntries, nj, & 
          !                                                               NZetaKnotEntries, ZetaKnotEntries, nk, &   
          !                                                               dxi_dtildexi)
          !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Calculating dx/dtildexi
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !call Get_dx_dtildexi_MappingFromParentElementToParameterSpace(RJac, dxi_dtildexi, dx_dtildexi)
       
          
          
          
        
          
          ! reassign the jacobian here with respect to the parent coordinates 
          !RJac = dx_dxi 
          
  
          !call RJacInv(IDim, RJac, InvRJac, DetJac, Det1)
          
          
          
          
          if ( (Det_dx_dtildexi < 0.0) .or. (Det_dtildexi_dx < 0.0) ) then
            call WriteInLogFile('Negative determinant '// trim(String(Det_dx_dtildexi)) // trim(String(Det_dtildexi_dx)) // ' element ' // trim(String(IElement)))
          end if

        end subroutine DetJacob
        
        
        
        !subroutine DetJacob_B(LocPos, NEl, NodTot, IDim, IElement, ICon, Co, RJac, InvRJac, DetJac)
        !!**********************************************************************
        !!
        !!  Function : Determination of the Jacobian matrix and the determinant of
        !!             the Jacobian matrix for the location LocPos inside an element IElement
        !!             whose nodal connectivities are defined by ICon and Co depending on
        !!             the type of element.
        !!
        !!  I  LocPos : Local coordinates of the considered point inside an element
        !!  I  NEl : Total number of elements
        !!  I  NodTot : Total number of nodes
        !!  I  IDim : Number of dimensions
        !!  I  IElement : ID of the element
        !!  I  ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !!  I  Co : Global nodal coordinates Co(I, J): j-coordinate of node I
        !!
        !!  O  RJac : Jacobian matrix
        !!  O  InvRJac : Inverse of the Jacobian matrix
        !!  O  DetJac : Determinate of the Jacobian matrix
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  integer(INTEGER_TYPE), intent(in) :: NEl, NodTot, IDim, IElement
        !  real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        !  integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        !  real(REAL_TYPE), dimension(:, :), intent(in) :: Co
        !  real(REAL_TYPE), dimension(:, :), intent(inout) :: RJac
        !  real(REAL_TYPE), dimension(:, :), intent(inout) :: InvRJac
        !  real(REAL_TYPE), intent(inout) :: DetJac
        !  
        !  ! local variables
        !  integer(INTEGER_TYPE) :: I, J, INode, NodeID, NNodes
        !  real(REAL_TYPE), dimension(:), allocatable :: HS
        !  real(REAL_TYPE), dimension(:, :), allocatable :: dHS
        !  real(REAL_TYPE) :: Det1
        !  
        !  real(REAL_TYPE), dimension(nen_NURBS,NDIM) :: dR_dx 
        !  real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
        !  real(REAL_TYPE), dimension(NDIM,NDIM) :: J_mat
        !  integer(INTEGER_TYPE) :: aa, bb, cc
        !  integer(INTEGER_TYPE) :: ni, nj, nk
        !  
        !   
        !  real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Xi
        !  real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Xi
        !  
        !  real(REAL_TYPE), allocatable, dimension(:) :: ShapeValuesArray_Temp_Eta
        !  real(REAL_TYPE), allocatable, dimension(:) :: DShapeValuesArray_Temp_Eta
        !  
        !  real(REAL_TYPE), dimension(:), allocatable :: ShapeValuesArray_Temp_Zeta
        !  real(REAL_TYPE), dimension(:), allocatable :: DShapeValuesArray_Temp_Zeta
        !  
        !  !real(REAL_TYPE) :: Det ! Temporary debug 
        !  !integer(INTEGER_TYPE) :: counter
        !  NNodes = size(ICon,1)
        !  allocate( HS(NNodes), dHS(NNodes, IDim) )
        !  
        !  !allocate( HS(NNodes), dHS(NNodes, IDim) )
        !  !allocate( HS(NNodes), dHS(NNodes, IDim) )
        !  !allocate( HS(NNodes), dHS(NNodes, IDim) )
        !  
        !  !this subroutine asks for the element ID to calcuate the Jacobian. This might need to be cast in NURBS parametric domain using ni, nj, nk...
        !  
        !  HS = 0.0
        !  dHS = 0.0
        !  !counter = 1
        !  
        !  ! get dR_dxi 
        !  ! Determine the shape functions HS and shape function derivatives dHS for the element for LocPos.
        !  call ShapeFunctionData(LocPos, NNodes, HS, dHS, IElement, &
        !                         ShapeValuesArray_Temp_Xi, DShapeValuesArray_Temp_Xi, &
        !                         ShapeValuesArray_Temp_Eta, DShapeValuesArray_Temp_Eta, & 
        !                         ShapeValuesArray_Temp_Zeta, DShapeValuesArray_Temp_Zeta) ! -> dHS here is dR_dxi
        !
        !  ! we would need to divide by denominators to complete the definitions of functions and derivatives w.r.t. parametric coodinates 
        !  
        !  !call InitialiseShapeFunctionsQUAD4_NURBS
        !  
        !  ! finding the NURBS coordinates for the element
        !  ni = INN(IEN(nen_NURBS,IElement),1)    
        !  nj = INN(IEN(nen_NURBS,IElement),2)
        !  
        !  if (NDIM == 3) then 
        !      nk = INN(IEN(nen_NURBS,IElement),3)
        !  end if 
        !  
        !  ! Determine the Jacobian matrix RJac
        !  RJac = 0.0
        !  do INode = 1, NNodes ! loop nodes of each element
        !    
        !      !ni = INN(IEN(1,IElement),1)    
        !      !nj = INN(IEN(2,IElement),2)
        !      
        !      !XiKnotEntries(ni) 
        !      !note that the coordinates of the control points are the ones we need here in 'Co'
        !      NodeID = ICon(INode, IElement)
        !      ! This is gradient mapping from the parameter space to the physical space 
        !    do I = 1, IDim
        !      do J = 1, IDim
        !        RJac(I, J) = RJac(I, J) + dHS(INode, I) * Co(NodeID, J) !calculation of the jacobian matrix from the coordinates of one element
        !      end do
        !    end do
        !  end do
        !  
        !  !Det = RJac(1,1) * RJac(2,2) - RJac(1,2) * RJac(2,1) ! Temporary debug !determinant of the jacobian matrix for 2D case
        !  
        !  ! RJac -> dx_dxi in Cottrell et al. (2011)
        !  
        !  !if (RJac(1,1) == 0.0) then 
        !  !
        !  !    call ShapeFunctionData(LocPos, NNodes, HS, dHS, IElement)
        !  !    return 
        !  !
        !  !end if 
        !  ! Calculate inverse and determinant of Jacobian matrix
        !  
        !  !counter = counter + 1
        !  
        !            
        !  !if (RJac(1,1) == 0.0) then 
        !  !
        !  !    return 
        !  !
        !  !end if 
        !  
        !  call RJacInv(IDim, RJac, InvRJac, DetJac, Det1)
        !  
        !  
        !  
        !  !_______________________________
        !  
        !  
        !  
        !  
        !  
        !  ! InvRJac -> dxi_dx 
        !  !
        !  !
        !  ! EVERYTHING BELOW HERE SHOULD GO IN THE PARENT SUBROUTINE 
        !  
        !  
        !   !Calculate derivative of basis functions with respect to physical coordinates
        !  call Get_dR_dx_MappingShapeFunctionDerivativeWRTPhysicalCoordinates(InvRJac, dHS, dR_dx) ! dR_dx = dxi_dx * dR_dxi 
        !  
        !  !dR_dx = dHS
        !  
        !  ! -> We would need to store with dR_dx in the shape function array 
        !  
        !  
        !  !call Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries, XiKnotEntries, ni, &
        !  !                                                               NEtaKnotEntries, EtaKnotEntries, nj, &   
        !  !                                                               dxi_dtildexi)
        !  !
        !  !          
        !  !call Get_dx_dtildexi_MappingFromParentElementToParameterSpace(RJac, dxi_dtildexi, dx_dtildexi)
        !  !
        !  !
        !  !
        !  !
        !  !
        !  !
        !  !! reassign the jacobian here with respect to the parent coordinates 
        !  !RJac_NURBS = J_mat 
        !  !
        !  !! recalculate jacobian determinant -> recalculate the deteminant and its inverse 
        !  !call RJacInv(IDim, RJac_NURBS, InvRJac_NURBS, DetJac, Det1)
        !  
        !  
        !  
        !  
        !  ! InvRJac -> dxi_dx 
        !  !
        !  !
        !  !! EVERYTHING BELOW HERE SHOULD GO IN THE PARENT SUBROUTINE 
        !  !
        !  !
        !  ! Calculate derivative of basis functions with respect to physical coordinates
        !  call Get_dR_dx_MappingShapeFunctionDerivativeWRTPhysicalCoordinates(InvRJac, dHS, dR_dx) ! dR_dx = dxi_dx * dR_dxi 
        !  !
        !  !!dR_dx = dHS
        !  !
        !  !! -> We would need to store with dR_dx in the shape function array 
        !  !
        !  !
        !  !call Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries, XiKnotEntries, ni, &
        !  !                                                               NEtaKnotEntries, EtaKnotEntries, nj, &   
        !  !                                                               dxi_dtildexi)
        !  !
        !  !          
        !  !
        !  !
        !  !
        !  !J_mat = 0.0
        !  !
        !  !
        !  !! calculation of J_matrix ! -> this is dx/dtildexi -> tilde symbolizes the parent element  
        !  !do aa = 1, NDIM
        !  !    do bb = 1, NDIM
        !  !        do cc = 1, NDIM
        !  !            J_mat(aa,bb) = J_mat(aa,bb) &
        !  !                           + RJac(aa,cc) * dxi_dtildexi(cc,bb)
        !  !        end do 
        !  !    end do 
        !  !end do 
        !  !
        !  !! reassign the jacobian here with respect to the parent coordinates 
        !  !RJac = J_mat 
        !  !
        !  !! recalculate jacobian determinant -> recalculate the deteminant and its inverse 
        !  !call RJacInv(IDim, RJac, InvRJac, DetJac, Det1)
        !  !
        !  
        !  if ( (DetJac < 0.0) .or. (Det1 < 0.0) ) then
        !    call WriteInLogFile('Negative determinant '// trim(String(DetJac)) // trim(String(Det1)) // ' element ' // trim(String(IElement)))
        !  end if
        !
        !end subroutine DetJacob_B
        
        
        
        
        
        subroutine Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries, XiKnotEntries, ni, &
                                                                             NEtaKnotEntries, EtaKnotEntries, nj, &  
                                                                             NZetaKnotEntries, ZetaKnotEntries, nk, &
                                                                             dxi_dtildexi)
        !------------------------------------------------------------------------------------------------------------
        ! Subroutine information: 
        !
        !
        !
        
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
        integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
        integer(INTEGER_TYPE), intent(in) :: NZetaKnotEntries

        
        integer(INTEGER_TYPE), intent(in) :: ni
        integer(INTEGER_TYPE), intent(in) :: nj
        integer(INTEGER_TYPE), intent(in) :: nk

        
        real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
        real(REAL_TYPE), dimension(NEtaKnotEntries), intent(in) :: EtaKnotEntries
        real(REAL_TYPE), dimension(NEtaKnotEntries), intent(in) :: ZetaKnotEntries


        real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi 


        dxi_dtildexi = 0.0
        ! gradient of mapping from parent element to parameter space 
        dxi_dtildexi(1,1) = ( XiKnotEntries(ni+1) - XiKnotEntries(ni) )/2
        dxi_dtildexi(2,2) = ( EtaKnotEntries(nj+1) - EtaKnotEntries(nj) )/2
        if (NDIM ==3) then 
        dxi_dtildexi(3,3) = ( ZetaKnotEntries(nk+1) - ZetaKnotEntries(nk) )/2
        end if 
        
        end subroutine Get_dxi_dtildexi_MappingFromParentElementToParameterSpace
        
        
        
        
        
        subroutine Get_dR_dx_MappingShapeFunctionDerivativeWRTPhysicalCoordinates(InvRJac, dHS, dR_dx, IPatch)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        implicit none 
        
        ! Hardcoded index to 1
        real(REAL_TYPE), dimension(nen_NURBS(1),NDIM), intent(out) :: dR_dx ! IPatch_Temporary
        
        real(REAL_TYPE), dimension(:, :), intent(in) :: dHS
        real(REAL_TYPE), dimension(:, :), intent(in) :: InvRJac
        
        integer(INTEGER_TYPE) :: loc_num, aa, bb
        
        integer(INTEGER_TYPE), intent(in) :: IPatch
        !integer(INTEGER_TYPE) :: IPatch_Temporary
        
        !do INode = 1, NNodes
        !    do I = 1, IDim            
        !        do J = 1, IDim
        !            dHSMultipledByInvRJac(I,J) = dHSMultipledByInvRJac(I,J) + ( dHS(INode, I) * InvRJac ) 
        !
        !        end do 
        !    end do 
        !end do 
        
        !dx_dxi = InvRJac
        !dx_dxi = InvRJac
        
        dR_dx = 0.0
        
        do loc_num = 1, nen_NURBS(IPatch)
            do aa = 1, NDIM 
                do bb = 1, NDIM 
                    dR_dx(loc_num, aa) = dR_dx(loc_num, aa) + &
                                                         (dHS(loc_num, bb) * InvRJac(bb, aa) )
                end do
            end do
        end do
        
        
        end subroutine Get_dR_dx_MappingShapeFunctionDerivativeWRTPhysicalCoordinates

        
        
        
        
        
        
        
        
        
        
        
        
        subroutine BMatrix(LocPos, IElTyp, NEl, NodTot, IDim, IElement, ICon, NodeCoord, B, IPatch)!, DetJac) ! -> DetJac here is | dx/dtilexi |
        !**********************************************************************
        !
        !    Function:  Determination of the strain interpolation matrix (B matrix)
        !               for the location LocPos inside an element IElement whose nodal
        !               connectivities are defined by ICon and NodeCoord.
        !
        !     LocPos : Position inside an element in local coordinates
        !     IElTyp : Number of nodes per element
        !     NEl : Total number of elements
        !     NodTot : Total number of nodes
        !     IDim : Number of dimensions
        !     IElement : ID of the element
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NodeCoord : Global nodal coordinates Co(I, J): J-coordinate of node I
        !
        ! O   B : IDim x IElTyp matrix containing the strain interpolation terms
        ! O   DetJac : Determinante of Jacobian
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          real(REAL_TYPE), dimension(IDim, IElTyp), intent(out) :: B
          real(REAL_TYPE) :: DetJac !, intent(out)
          real(REAL_TYPE) :: Det1

          ! Local variables
          integer(INTEGER_TYPE) :: I, J, K
          real(REAL_TYPE), allocatable, dimension(:) :: HS ! Shape functions
          real(REAL_TYPE), allocatable, dimension(:, :) :: dHS ! Derivatives of shape functions
          real(REAL_TYPE), dimension(IDim, IDim) :: RJac!, RJacInv ! Jacobian matrix, inverse of Jacobian matrix
          real(REAL_TYPE), dimension(IDim, IDim) :: dx_dxi
          real(REAL_TYPE), dimension(IDim, IDim) :: dxi_dx

          ! Hardcoded index to 1
          real(REAL_TYPE), dimension(nen_NURBS(1),NDIM) :: dR_dx ! IPatch_Temporary
          real(REAL_TYPE), dimension(NDIM,NDIM) :: J_mat
          
          real(REAL_TYPE) :: Det_dx_dxi
          real(REAL_TYPE) :: Det_dxi_dx
          
          real(REAL_TYPE) :: Det_dtildexi_dx
          real(REAL_TYPE) :: Det_dx_dtildexi
          !real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
                  
          
          integer(INTEGER_TYPE) :: ni
          integer(INTEGER_TYPE) :: nj
          integer(INTEGER_TYPE) :: nk
          
         
          
          ! allocation of variables
          real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
          real(REAL_TYPE), dimension(NDIM,NDIM) :: dx_dtildexi
          real(REAL_TYPE), dimension(NDIM,NDIM) :: dtildexi_dx
          
          ! Multipatch variables 
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          integer(INTEGER_TYPE) :: aa
          
          ni = INN(IEN(1,IElement,IPatch),1,IPatch)    
          nj = INN(IEN(1,IElement,IPatch),2,IPatch)
          if (NDIM == 3) then 
              nk = INN(IEN(1,IElement,IPatch),3,IPatch)
          end if 



          ! Determine the shape functions HS and shape function derivatives dHS for LocPos.
          call ShapeFunctionData(LocPos, IElTyp, HS, dHS, IElement, IPatch)

          ! Evaluate jacobian between physical domain and parameter domain
          call Evaluate_dx_dxi(dHS, NodeCoord, ICon, dx_dxi, IElement)
          
          ! debugger
      if (dx_dxi(3,3) == 0.0) then 
          aa=1
      end if

          
          ! Evaluate the inverse of jacobian. (dx_dxi)^-1 = dxi_dx
          call RJacInv(NDIM, dx_dxi, dxi_dx, Det_dx_dxi, Det_dxi_dx)
          
          ! Check for negative jacobian matrix. If it is negative, you need to fix the connectivity. 
          if ( (Det_dx_dxi < 0.0) .or. (Det_dxi_dx < 0.0) ) then
            call WriteInLogFile('Negative determinant '// trim(String(Det_dx_dxi)) // trim(String(Det_dxi_dx)) // ' element ' // trim(String(IElement)))
          end if
          
         ! Assemble B matrix (cartesian derivatives)
         B = 0.0
         do J = 1, IDim
           do I = 1, IElTyp
             do K = 1, IDim
               ! dR/dx =          dxi/dx       * dR/dxi
               B(K, I) = B(K,I) + dxi_dx(K, J) * dHS(I, J)
               
               ! dHS here needs to be dR_dxi (physical coordinates)
               ! RJacInv here is dxi_dx 
               ! B is a matrix collecting all dR_dx values
                                 
               ! This is how we usually write B in books on computational mechanics. 
               ! 2D----------------------------
               ! B = [dR_dx     0       
               !      0        dR_dy    
               !      dR_dy    dR_dx]
               ! 3D----------------------------
               ! B = [dR_dx     0       0
               !      0        dR_dy    0
               !      0         0       dR_dz
               !      dR_dy    dR_dx    0
               !      0        dR_dz    dR_dy
               !      dR_dz     0       dR_dx]
               ! However, in this case, this is how it is arranged.
               ! B = [dRi_dx     ---       dRIElTyp_dx          --> x global direction  
               !      dRi_dy     ---       dRIElTyp_dx          --> y global direction
               !      dRi_dz     ---       dRIElTyp_dx]         --> z global direction
               !      
               
             end do
           end do
         end do       
        end subroutine BMatrix
        
        
        
        
        
        
        
        
        
        
        subroutine Evaluate_dx_dxi(dHS, NodeCoord, ICon, dx_dxi, IElement)
          ! this subroutine is part of the BMatrix evaluation. 
          !  
          ! 
          
          implicit none 
          

          integer(INTEGER_TYPE) :: I, J, INode, NodeID, NNodes
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          
          real(REAL_TYPE), dimension(:, :), intent(in) :: dHS
          integer(INTEGER_TYPE), intent(in) :: IElement


          real(REAL_TYPE), dimension(NDIM, NDIM), intent(out) :: dx_dxi


          NNodes = size(ICon,1)
          
          dx_dxi = 0.0
          do INode = 1, NNodes ! loop nodes of each element
            
              !note that the coordinates of the control points are the ones we need here in 'Co'
              NodeID = ICon(INode, IElement)
              ! This is gradient mapping from the parameter space to the physical space 
            do I = 1, NDIM 
              do J = 1, NDIM
                dx_dxi(I, J) = dx_dxi(I, J) + (dHS(INode, I) * NodeCoord(NodeID, J)) !calculation of the jacobian matrix from the coordinates of one element ! -> this is checked and now it is working
              end do
            end do
          end do
          
          
          
          
          
          
          end subroutine Evaluate_dx_dxi
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        subroutine Get_dx_dtildexi_MappingFromParentElementToParameterSpace(RJac, dxi_dtildexi, dx_dtildexi)
          
        ! not J_mat is same as dx_dtildexi
        real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: dx_dtildexi
        real(REAL_TYPE), dimension(NDIM,NDIM), intent(in) :: dxi_dtildexi
        real(REAL_TYPE), dimension(NDIM,NDIM), intent(in) :: RJac
        
        integer(INTEGER_TYPE) :: aa, bb, cc


          dx_dtildexi = 0.0
          
          
          ! calculation of J_matrix ! -> this is dx/dtildexi -> tilde symbolizes the parent element  
          do aa = 1, NDIM
              do bb = 1, NDIM
                  do cc = 1, NDIM
                      dx_dtildexi(aa,bb) = dx_dtildexi(aa,bb) &
                                     + RJac(aa,cc) * dxi_dtildexi(cc,bb)
                  end do 
              end do 
          end do 
         
          
        end subroutine Get_dx_dtildexi_MappingFromParentElementToParameterSpace
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STRUCTURED MESH - GetLocalCoordinates - NURBS implementation 
        
        
        
        
    !    subroutine GetLocalCoordinates_StructuredMesh(GlobPos, IElement, IElTyp, NEl, NodTot, IDim, NodeCoord, ICon, LocPos, OutsideElement, Success)
    !    !**********************************************************************
    !    !
    !    !    Function:  Determination of local coordinates from global coordinates,
    !    !               assuming that the point with global coordinates lies inside 
    !    !               the element within IElement.
    !    !               Success returns .false. if the local position could not be found
    !    !               inside the element within 10 iterations with sufficient accuracy.
    !    !
    !    !     GlobPos : Global coordinates of a point inside IElement
    !    !     IElement : ID of the element
    !    !     IElTyp : Number of node connectivities of IElement
    !    !     NEl : Number of elements
    !    !     NodTot : Total number of nodes
    !    !     IDim : Dimension of the mesh
    !    !     NodeCoord : Global nodal coordinates
    !    !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
    !    !
    !    ! O   LocPos : Local coordinates of the considered point inside IElement
    !    ! O   OutsideElement : True, if the local coordinate lie outside IElement
    !    ! O   Success : True, if the local position could be found
    !    !
    !    ! Implemented in the frame of the MPM project.
    !    !
    !    !**********************************************************************
    !  
    !    implicit none
    !
    !      integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
    !      real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
    !      integer(INTEGER_TYPE), intent(in) :: IElement
    !      real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
    !      integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
    !      real(REAL_TYPE), dimension(IDim), intent(out) :: LocPos !local position of point is an output 
    !      logical, intent(out) :: OutsideElement, Success
    !      
    !      ! Local variables
    !      real(REAL_TYPE), dimension(IDim, IDim) :: RJac, RJacInv
    !      real(REAL_TYPE), dimension(IDim) :: GlobPosIteration
    !      real(REAL_TYPE), dimension(IDim) :: DeltaGlobPos
    !      integer(INTEGER_TYPE) :: I, J, Iteration
    !      real(REAL_TYPE) :: Tolerance, Difference, DetJac
    !      
    !      real(REAL_TYPE) :: KnotSpanCoord1X, KnotSpanCoord1Y
    !      real(REAL_TYPE) :: KnotSpanCoord2X, KnotSpanCoord2Y
    !      real(REAL_TYPE) :: KnotSpanCoord3X, KnotSpanCoord3Y
    !      real(REAL_TYPE) :: KnotSpanCoord4X, KnotSpanCoord4Y
    !      
    !      real(REAL_TYPE) :: LowerBoundXCoordinate, UpperBoundXCoordinate
    !      real(REAL_TYPE) :: LowerBoundYCoordinate, UpperBoundYCoordinate
    !      
    !      real(REAL_TYPE) :: DenominatorLocPosRange_1, DenominatorLocPosRange_2
    !            
    !      if (IDim==2) then
    !        Tolerance = 1d-15 * (dabs(GlobPos(1)) + dabs(GlobPos(2)) ) + 1d-15 !1d-10
    !      else
    !        Tolerance = 1d-15 * (dabs(GlobPos(1) ) + dabs(GlobPos(2) ) + dabs(GlobPos(3) ) ) + 1d-10  
    !      end if    
    !      
    !      Success = .true.
    !      OutsideElement = .false.
    !      
    !      
    !      !ICon = KnotBezierMeshConnectivity(IElement, :) ! This variable needs to be 1 row 4 columns 
    !      ! for element 1, this should be ICon = [1, 2, 3, 4]
    !
    !      !Based on this identify the coordinates of the knot spans
    !      KnotSpanCoord1X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      KnotSpanCoord1Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      
    !      KnotSpanCoord2X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      KnotSpanCoord2Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      
    !      KnotSpanCoord3X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      KnotSpanCoord3Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      
    !      KnotSpanCoord4X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      KnotSpanCoord4Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
    !      
    !      
    !      !We can use the knot coordinates for "1" and "4". That will give you the extent of the element 
    !      
    !      LowerBoundXCoordinate = KnotSpanCoord1X !should be 0.0
    !      UpperBoundXCoordinate = KnotSpanCoord4X !should be 1.0
    !      
    !      LowerBoundYCoordinate = KnotSpanCoord1Y !should be 0.0
    !      UpperBoundYCoordinate = KnotSpanCoord4Y !should be 0.5
    !      
    !      
    !      
    !      !Based on this, we can interpolate to find the LocPos in each direction
    !      
    !      ! LocPos(1)
    !      DenominatorLocPosRange_1 = UpperBoundXCoordinate - LowerBoundXCoordinate
    !      
    !      LocPos(1) = ( (  (GlobPos(1)-LowerBoundXCoordinate)/DenominatorLocPosRange_1 ) * (1--1) ) + -1
    !      
    !      
    !      ! LocPos(2)
    !      DenominatorLocPosRange_2 = UpperBoundYCoordinate - LowerBoundYCoordinate 
    !      
    !      LocPos(2) = ( ( (GlobPos(2)-LowerBoundYCoordinate) /DenominatorLocPosRange_2 ) * (1--1) ) + -1
    !      
    !      
    !      !! X-axis loc position
    !      !do I = 1, IElTyp
    !      !  NodeID = ICon(I, IElement)
    !      !  do J = 1, IDim
    !      !      NodeCoord(NodeID, J)
    !      !  end do
    !      !end do
    !      
    !      
    ! 
    !      ! Initial guess at local coordinates
    !      !LocPos = (/0.2, 0.211111, 0.188888/)
    !      !LocPos = (/0.0, 0.0, 0.0/)
    !      ! Global coordinates for guess
    !      !call GetGlobalCoordinates(LocPos,  &
    !      !                          IElTyp, NEl, NodTot, IDim, &
    !      !                          IElement, ICon, NodeCoord, &
    !      !                          GlobPosIteration)
    !
    !      !Iteration = 1
    !      !do ! Iterate towards equality of local and global coordinates
    !        ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix RJacInv
    !        
    !        !  if (LocPos(1) < -1) then 
    !        !
    !        !  return 
    !        !
    !        !
    !        !end if 
    !          
    !          
    !          !
    !          !call DetJacob(LocPos, NEl, NodTot, IDim, &
    !          !            IElement, ICon, NodeCoord, &
    !          !            RJac, RJacInv, DetJac) !update Jacob
    !
    !        
    !      !  
    !      !  ! Determine new LocPos
    !      !  DeltaGlobPos = GlobPos - GlobPosIteration
    !      !
    !      !  do I = 1, IDim
    !      !    do J = 1, IDim
    !      !      LocPos(I) = LocPos(I) + RJacInv(J, I) * DeltaGlobPos(J) !update Jacob
    !      !    end do
    !      !  end do
    !      !
    !      !  ! Global coordinates for guess
    !      !  call GetGlobalCoordinates(LocPos,  &
    !      !                            IElTyp, NEl, NodTot, IDim, &
    !      !                            IElement, ICon, NodeCoord, &
    !      !                            GlobPosIteration)
    !      !
    !      !  ! Check whether loop can be aborted
    !      !  Difference = 0.0
    !      !  do I = 1, IDim
    !      !    Difference = Difference + dabs(GlobPos(I) - GlobPosIteration(I))
    !      !  end do
    !      !
    !      !  if (Difference < Tolerance) then ! Found local coordinates
    !      !    Success = .true.
    !      !    EXIT
    !      !  else ! Difference greater or equal Tolerance
    !      !    if (Iteration >= 1000) then ! Too many iterations needed, something went wrong
    !      !      Success = .false.
    !      !      EXIT
    !      !    else ! Continue iteration
    !      !      Iteration = Iteration + 1
    !      !    end if
    !      !  end if
    !      !end do ! Iteration loop
    !
    !      if (.not.IsInsideElementLocPos(LocPos) ) then
    !        OutsideElement = .true.
    !        Success = .true.
    !      end if
    !
    !      if (.not.Success) then
    !        call GiveError('Did not find local coordinates of particle in '// &
    !                       trim(String(IElement)) // &
    !                       ' within limit number of iterations.')
    !      end if
    !
    !end subroutine GetLocalCoordinates_StructuredMesh
    !    
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            
    
    
           
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !subroutine GetLocalCoordinates_StructuredMesh_3D(GlobPos, IElement, IElTyp, NEl, NodTot, IDim, NodeCoord, ICon, LocPos, OutsideElement, Success)
        !!**********************************************************************
        !!
        !!    Function:  Determination of local coordinates from global coordinates,
        !!               assuming that the point with global coordinates lies inside 
        !!               the element within IElement.
        !!               Success returns .false. if the local position could not be found
        !!               inside the element within 10 iterations with sufficient accuracy.
        !!
        !!     GlobPos : Global coordinates of a point inside IElement
        !!     IElement : ID of the element
        !!     IElTyp : Number of node connectivities of IElement
        !!     NEl : Number of elements
        !!     NodTot : Total number of nodes
        !!     IDim : Dimension of the mesh
        !!     NodeCoord : Global nodal coordinates
        !!     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !!
        !! O   LocPos : Local coordinates of the considered point inside IElement
        !! O   OutsideElement : True, if the local coordinate lie outside IElement
        !! O   Success : True, if the local position could be found
        !!
        !! Implemented in the frame of the MPM project.
        !!
        !!**********************************************************************
        !
        !implicit none
        !
        !  integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
        !  real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
        !  integer(INTEGER_TYPE), intent(in) :: IElement
        !  real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
        !  integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        !  real(REAL_TYPE), dimension(IDim), intent(out) :: LocPos !local position of point is an output 
        !  logical, intent(out) :: OutsideElement, Success
        !  
        !  ! Local variables
        !  real(REAL_TYPE), dimension(IDim, IDim) :: RJac, RJacInv
        !  real(REAL_TYPE), dimension(IDim) :: GlobPosIteration
        !  real(REAL_TYPE), dimension(IDim) :: DeltaGlobPos
        !  integer(INTEGER_TYPE) :: I, J, Iteration
        !  real(REAL_TYPE) :: Tolerance, Difference, DetJac
        !  
        !  real(REAL_TYPE) :: KnotSpanCoord1X, KnotSpanCoord1Y, KnotSpanCoord1Z
        !  real(REAL_TYPE) :: KnotSpanCoord2X, KnotSpanCoord2Y, KnotSpanCoord2Z
        !  real(REAL_TYPE) :: KnotSpanCoord3X, KnotSpanCoord3Y, KnotSpanCoord3Z
        !  real(REAL_TYPE) :: KnotSpanCoord4X, KnotSpanCoord4Y, KnotSpanCoord4Z
        !  real(REAL_TYPE) :: KnotSpanCoord5X, KnotSpanCoord5Y, KnotSpanCoord5Z
        !  real(REAL_TYPE) :: KnotSpanCoord6X, KnotSpanCoord6Y, KnotSpanCoord6Z
        !  real(REAL_TYPE) :: KnotSpanCoord7X, KnotSpanCoord7Y, KnotSpanCoord7Z
        !  real(REAL_TYPE) :: KnotSpanCoord8X, KnotSpanCoord8Y, KnotSpanCoord8Z
        !  
        !  
        !  real(REAL_TYPE) :: LowerBoundXCoordinate, UpperBoundXCoordinate
        !  real(REAL_TYPE) :: LowerBoundYCoordinate, UpperBoundYCoordinate
        !  real(REAL_TYPE) :: LowerBoundZCoordinate, UpperBoundZCoordinate
        !  
        !  
        !  real(REAL_TYPE) :: DenominatorLocPosRange_1, DenominatorLocPosRange_2, DenominatorLocPosRange_3
        !        
        !  if (IDim==2) then
        !    Tolerance = 1d-15 * (dabs(GlobPos(1)) + dabs(GlobPos(2)) ) + 1d-15 !1d-10
        !  else
        !    Tolerance = 1d-15 * (dabs(GlobPos(1) ) + dabs(GlobPos(2) ) + dabs(GlobPos(3) ) ) + 1d-10  
        !  end if    
        !  
        !  Success = .true.
        !  OutsideElement = .false.
        !  
        !  
        !  !ICon = KnotBezierMeshConnectivity(IElement, :) ! This variable needs to be 1 row 4 columns 
        !  ! for element 1, this should be ICon = [1, 2, 3, 4]
        !
        !  !Based on this identify the coordinates of the knot spans
        !  !KnotSpanCoord1X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord1Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord1Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord2X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord2Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord2Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord3X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord3Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord3Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord4X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord4Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord4Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord5X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 5),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord5Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 5),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord5Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 5),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord6X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 6),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord6Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 6),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord6Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 6),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord7X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 7),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord7Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 7),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord7Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 7),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !
        !  !KnotSpanCoord8X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 8),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord8Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 8),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  !KnotSpanCoord8Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 8),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  KnotSpanCoord1X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord1Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord1Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 1),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  KnotSpanCoord2X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord2Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord2Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 2),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  KnotSpanCoord3X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord3Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord3Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 3),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  KnotSpanCoord4X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord4Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord4Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 4),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  KnotSpanCoord5X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 5),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord5Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 5),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord5Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 5),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  KnotSpanCoord6X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 6),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord6Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 6),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord6Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 6),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  KnotSpanCoord7X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 7),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord7Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 7),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord7Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 7),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !
        !  KnotSpanCoord8X = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 8),1) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord8Y = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 8),2) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  KnotSpanCoord8Z = KnotBezierMesh(KnotBezierMeshConnectivity(IElement, 8),3) ! first entry is from the connectivity, second entry is from the IDim (i.e., x or y)
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  
        !  !We can use the knot coordinates for "1" and "4". That will give you the extent of the element 
        !  ! Are we keeping this as 1, 4, ..., ???
        !  LowerBoundXCoordinate = KnotSpanCoord1X !should be 0.0
        !  UpperBoundXCoordinate = KnotSpanCoord8X !should be 1.0
        !  
        !  LowerBoundYCoordinate = KnotSpanCoord1Y !should be 0.0
        !  UpperBoundYCoordinate = KnotSpanCoord8Y !should be 0.5
        !  
        !  LowerBoundZCoordinate = KnotSpanCoord1Z !should be 0.0
        !  UpperBoundZCoordinate = KnotSpanCoord8Z !should be 0.5
        !  
        !  
        !  !Based on this, we can interpolate to find the LocPos in each direction
        !  
        !  ! LocPos(1)
        !  DenominatorLocPosRange_1 = UpperBoundXCoordinate - LowerBoundXCoordinate
        !  
        !  LocPos(1) = ( (  (GlobPos(1)-LowerBoundXCoordinate)/DenominatorLocPosRange_1 ) * (1--1) ) + -1
        !  
        !  
        !  ! LocPos(2)
        !  DenominatorLocPosRange_2 = UpperBoundYCoordinate - LowerBoundYCoordinate 
        !  
        !  LocPos(2) = ( ( (GlobPos(2)-LowerBoundYCoordinate) /DenominatorLocPosRange_2 ) * (1--1) ) + -1
        !  
        !  
        !  ! LocPos(3)
        !  DenominatorLocPosRange_3 = UpperBoundZCoordinate - LowerBoundZCoordinate 
        !  
        !  LocPos(3) = ( ( (GlobPos(3)-LowerBoundZCoordinate) /DenominatorLocPosRange_3 ) * (1--1) ) + -1
        !  
        !  !! X-axis loc position
        !  !do I = 1, IElTyp
        !  !  NodeID = ICon(I, IElement)
        !  !  do J = 1, IDim
        !  !      NodeCoord(NodeID, J)
        !  !  end do
        !  !end do
        !  
        !  
        !
        !  ! Initial guess at local coordinates
        !  !LocPos = (/0.2, 0.211111, 0.188888/)
        !  !LocPos = (/0.0, 0.0, 0.0/)
        !  ! Global coordinates for guess
        !  !call GetGlobalCoordinates(LocPos,  &
        !  !                          IElTyp, NEl, NodTot, IDim, &
        !  !                          IElement, ICon, NodeCoord, &
        !  !                          GlobPosIteration)
        !
        !  !Iteration = 1
        !  !do ! Iterate towards equality of local and global coordinates
        !    ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix RJacInv
        !    
        !    !  if (LocPos(1) < -1) then 
        !    !
        !    !  return 
        !    !
        !    !
        !    !end if 
        !      
        !      
        !      !
        !      !call DetJacob(LocPos, NEl, NodTot, IDim, &
        !      !            IElement, ICon, NodeCoord, &
        !      !            RJac, RJacInv, DetJac) !update Jacob
        !
        !    
        !  !  
        !  !  ! Determine new LocPos
        !  !  DeltaGlobPos = GlobPos - GlobPosIteration
        !  !
        !  !  do I = 1, IDim
        !  !    do J = 1, IDim
        !  !      LocPos(I) = LocPos(I) + RJacInv(J, I) * DeltaGlobPos(J) !update Jacob
        !  !    end do
        !  !  end do
        !  !
        !  !  ! Global coordinates for guess
        !  !  call GetGlobalCoordinates(LocPos,  &
        !  !                            IElTyp, NEl, NodTot, IDim, &
        !  !                            IElement, ICon, NodeCoord, &
        !  !                            GlobPosIteration)
        !  !
        !  !  ! Check whether loop can be aborted
        !  !  Difference = 0.0
        !  !  do I = 1, IDim
        !  !    Difference = Difference + dabs(GlobPos(I) - GlobPosIteration(I))
        !  !  end do
        !  !
        !  !  if (Difference < Tolerance) then ! Found local coordinates
        !  !    Success = .true.
        !  !    EXIT
        !  !  else ! Difference greater or equal Tolerance
        !  !    if (Iteration >= 1000) then ! Too many iterations needed, something went wrong
        !  !      Success = .false.
        !  !      EXIT
        !  !    else ! Continue iteration
        !  !      Iteration = Iteration + 1
        !  !    end if
        !  !  end if
        !  !end do ! Iteration loop
        !
        !  if (.not.IsInsideElementLocPos(LocPos) ) then
        !    OutsideElement = .true.
        !    Success = .true.
        !  end if
        !
        !  if (.not.Success) then
        !    call GiveError('Did not find local coordinates of particle in '// &
        !                   trim(String(IElement)) // &
        !                   ' within limit number of iterations.')
        !  end if
        !
        !end subroutine GetLocalCoordinates_StructuredMesh_3D
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
        ! STRUCTURED MESH - GetLocalCoordinates - NURBS implementation 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        subroutine GetLocalCoordinates(GlobPos, UStepArray, IElement, IElTyp, NEl, NodTot, IDim, NodeCoord, ICon, LocPos, OutsideElement, Success, ElementAdjacencies, IPatch)!, IncrementalDisplacementSoil)
        !**********************************************************************
        !
        !    Function:  Determination of local coordinates from global coordinates,
        !               assuming that the point with global coordinates lies inside 
        !               the element within IElement.
        !               Success returns .false. if the local position could not be found
        !               inside the element within 10 iterations with sufficient accuracy.
        !
        !     GlobPos : Global coordinates of a point inside IElement
        !     IElement : ID of the element
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   LocPos : Local coordinates of the considered point inside IElement
        ! O   OutsideElement : True, if the local coordinate lie outside IElement
        ! O   Success : True, if the local position could be found
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(:), intent(in) :: UStepArray
          integer(INTEGER_TYPE), intent(in) :: IElement
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          !real(REAL_TYPE), dimension(:, :), intent(in) :: IncrementalDisplacementSoil
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim), intent(inout) :: LocPos !local position of point is an output 
          logical, intent(out) :: OutsideElement, Success
          integer(INTEGER_TYPE), dimension(:,:), intent(in) :: ElementAdjacencies
          real(REAL_TYPE), dimension(NDIM, NDIM) :: dtildexi_dx
          ! Local variables
          real(REAL_TYPE), dimension(IDim, IDim) :: RJac, RJacInv
          real(REAL_TYPE), dimension(IDim, IDim) :: dx_dxi
          real(REAL_TYPE), dimension(IDim) :: GlobPosIteration
          real(REAL_TYPE), dimension(IDim) :: DeltaGlobPos
          integer(INTEGER_TYPE) :: I, J, Iteration, k
          real(REAL_TYPE) :: Tolerance, Difference, DetJac
          
          real(REAL_TYPE), dimension(IDim) :: ParametricDomain
          real(REAL_TYPE), dimension(IDim, IDim) :: dxi_dx
          real(REAL_TYPE), dimension(IDim, IDim) :: dtildexi_dxi
          
          integer(INTEGER_TYPE) :: ni, nj, nk, test
          
          real(REAL_TYPE), allocatable, dimension(:) :: HS
          real(REAL_TYPE), allocatable, dimension(:,:) :: dHS 
          
          ! NURBS local variables 
          real(REAL_TYPE), dimension(IDim) :: LocPos_Xi
          real(REAL_TYPE), dimension(IDim) :: LocPos_XiTilde      
          real(REAL_TYPE), dimension(IDim, IDim) :: dx_dtildexi!, dtildexi_dx
          real(REAL_TYPE) :: Det_dx_dtildexi
          real(REAL_TYPE), dimension(IDim) :: tolerance_Xi
          integer(INTEGER_TYPE) :: NewElementID, CrossedSide
          integer(INTEGER_TYPE) :: counter_xi_outsideelement, counter_eta_outsideelement, counter_zeta_outsideelement, counter_outside
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch
          integer(INTEGER_TYPE) :: aa ! debugging place holder

          
          counter_xi_outsideelement = 0
          counter_eta_outsideelement = 0
          counter_zeta_outsideelement = 0
          counter_outside = 0
          

          !integer(INTEGER_TYPE), intent(in) :: GetAdjacentElement
          ! evaluate tolerances based on input "true" global positon 
          if (IDim==2) then
            Tolerance = 1d-15 * (dabs(GlobPos(1)) + dabs(GlobPos(2)) ) + 1d-15 !1d-10
          else
            Tolerance = 1d-12 * (dabs(GlobPos(1) ) + dabs(GlobPos(2) ) + dabs(GlobPos(3) ) ) + 1d-12
          end if    
          
          ! Set the tolerance value to 1%
            tolerance_Xi(1) = 1e-4!Tolerance!!0.02*(XiKnotEntries(ni+1)-XiKnotEntries(ni)) ! Adjust this value as needed
            tolerance_Xi(2) = 1e-4!Tolerance!0.02*(EtaKnotEntries(nj+1)-EtaKnotEntries(nj)) ! Adjust this value as needed
            tolerance_Xi(3) = 1e-4!Tolerance!0.02*(ZetaKnotEntries(nk+1)-ZetaKnotEntries(nk)) ! Adjust this value as needed
          
          ! local positions can be expressed in either xi or xitilde format
          
          Success = .true.
          OutsideElement = .false.
          
          !if (.not.IsInsideElementLocPos(LocPos) ) then
          !  LocPos_Xi = (/0.0, 0.0, 0.0/)  
          !end if
     
          ! Initial guess at local coordinates
          !LocPos = (/0.1, 0.3, 0.5/)
                
          !LocPos = 0! (/-0.1, 0.211111, -0.1/)
          !ParametricDomain = (/0.5, 0.5, 0.5/)
          !LocPos = (/0.2, 0.211111, 0/)
          !LocPos = 0.0
          ! Global coordinates for guess
          
          
         ni = INN(IEN(1,IElement,IPatch),1,IPatch)    
         nj = INN(IEN(1,IElement,IPatch),2,IPatch)
         
         if (NDIM == 3) then 
         nk = INN(IEN(1,IElement,IPatch),3,IPatch)
         end if 
         ! != LocPos(1)
         ! != LocPos(2)
         !
         !ParametricDomain(1) =  (( (XiKnotEntries(ni+1) - XiKnotEntries(ni) ) * LocPos(1)) &
         !                           + (XiKnotEntries(ni+1) + XiKnotEntries(ni)) ) * 0.5; ! this should be a scalar always
         !   
         !ParametricDomain(2) =  ( ((EtaKnotEntries(nj+1) - EtaKnotEntries(nj) ) * LocPos(2)) &
         !                           + (EtaKnotEntries(nj+1) + EtaKnotEntries(nj)) ) * 0.5; ! this should be a scalar always
         !
         !if (NDIM == 3) then 
         !ParametricDomain(3) =  ( ((ZetaKnotEntries(nk+1) - ZetaKnotEntries(nk) ) * LocPos(3)) &
         !                           + (ZetaKnotEntries(nk+1) + ZetaKnotEntries(nk)) ) * 0.5; ! this should be a scalar always
         !end if 
          
          ! Local position in the xi parametric domain 
          LocPos_Xi = (/0.511111111*(XiKnotEntries(ni,IPatch)+XiKnotEntries(ni+1,IPatch)), 0.522222222*(EtaKnotEntries(nj,IPatch)+EtaKnotEntries(nj+1,IPatch)), 0.511111111*(ZetaKnotEntries(nk,IPatch)+ZetaKnotEntries(nk+1,IPatch))/)  
          
          LocPos_XiTilde = 0
          
          ! Corresponding local position in the xi tilde parent domain
          LocPos_XiTilde(1) = ( (2*LocPos_Xi(1)) - (XiKnotEntries(ni+1,IPatch)  +  XiKnotEntries(ni,IPatch)) )    /(XiKnotEntries(ni+1,IPatch)-XiKnotEntries(ni,IPatch))
          LocPos_XiTilde(2) = ( (2*LocPos_Xi(2)) - (EtaKnotEntries(nj+1,IPatch) + EtaKnotEntries(nj,IPatch)) )    /(EtaKnotEntries(nj+1,IPatch)-EtaKnotEntries(nj,IPatch))
          LocPos_XiTilde(3) = ( (2*LocPos_Xi(3)) - (ZetaKnotEntries(nk+1,IPatch)+ZetaKnotEntries(nk,IPatch)) )    /(ZetaKnotEntries(nk+1,IPatch)-ZetaKnotEntries(nk,IPatch))
          
          ! Calculate global position based on the local position
          call GetGlobalCoordinates(LocPos_XiTilde,  &
                                    IElTyp, NEl, NodTot, IDim, &
                                    IElement, ICon, NodeCoord, &
                                    GlobPosIteration, &
                                    IPatch)!, &
                                    !IncrementalDisplacementSoil)

          Iteration = 1
          do ! Iterate towards equality of local and global coordinates
            ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix RJacInv
            
              !evaluate shape functions here and input into the jacobian evaluation 
               call ShapeFunctionData(LocPos_XiTilde, ELEMENTNODES, HS, dHS, IElement, IPatch)
              
              ! Jacobian calculation
              call DetJacob(HS, dHS, & !-> dHS = dR/dx
                          LocPos_XiTilde, NEl, NodTot, IDim, &
                          IElement, ICon, NodeCoord, & 
                          dx_dxi, &  
                          dx_dtildexi, dtildexi_dx, Det_dx_dtildexi, &
                          dxi_dx, dtildexi_dxi, &
                  IPatch) !update Jacob

                            
              !end if 

            
            
            ! Determine new LocPos
            DeltaGlobPos = GlobPos - GlobPosIteration

            !do I = 1, IDim
            !  do J = 1, IDim
            !    LocPos(I) = LocPos(I) + RJacInv(J, I) * DeltaGlobPos(J) !update Jacob
            !  end do
            !end do
        
            !do i = 1, NDIM
            !    do j = 1, NDIM
            !        dtildexi_dx(i, j) = 0
            !        do k = 1, NDIM
            !            dtildexi_dx(i, j) = dtildexi_dx(i, j) + dxi_dx(i, k) * dtildexi_dxi(k, j)
            !        end do
            !    end do
            !end do
            ni = INN(IEN(1,IElement,IPatch),1,IPatch)    
            nj = INN(IEN(1,IElement,IPatch),2,IPatch)
         
            if (NDIM == 3) then 
                nk = INN(IEN(1,IElement,IPatch),3,IPatch)
            end if 
            
            do I = 1, IDim
              do J = 1, IDim
                !LocPos(I) = LocPos(I) + (dxi_dx(J, I) * dtildexi_dxi(J, I) * DeltaGlobPos(J))!update Jacob
                !ParametricDomain(I) = ParametricDomain(I) + (dxi_dx(I,J) * DeltaGlobPos(J))!update Jacob

                  LocPos_Xi(I) = LocPos_Xi(I) + (dxi_dx(J,I) * DeltaGlobPos(J))!update Jacob

              end do
            end do
            
            
            ! Corresponding local position in the xi tilde parent domain
            LocPos_XiTilde(1) = ( (2*LocPos_Xi(1)) - (XiKnotEntries(ni+1,IPatch)  +  XiKnotEntries(ni,IPatch)) )    /(XiKnotEntries(ni+1,IPatch)-XiKnotEntries(ni,IPatch))
            LocPos_XiTilde(2) = ( (2*LocPos_Xi(2)) - (EtaKnotEntries(nj+1,IPatch) + EtaKnotEntries(nj,IPatch)) )    /(EtaKnotEntries(nj+1,IPatch)-EtaKnotEntries(nj,IPatch))
            LocPos_XiTilde(3) = ( (2*LocPos_Xi(3)) - (ZetaKnotEntries(nk+1,IPatch)+ZetaKnotEntries(nk,IPatch)) )    /(ZetaKnotEntries(nk+1,IPatch)-ZetaKnotEntries(nk,IPatch))
            
            if (.not.IsInsideElementLocPos(LocPos_XiTilde) ) then
                    ! based on which side you crossed you can knot which element you are in
                    if (LocPos_XiTilde(1) <-1 ) then 
                        CrossedSide = 3!side 3
                    elseif (LocPos_XiTilde(1) >1 ) then
                        CrossedSide = 2!side 2
                    elseif (LocPos_XiTilde(2) <-1 ) then 
                        CrossedSide = 6!side 6
                    elseif (LocPos_XiTilde(2) >1 ) then
                        CrossedSide = 5!side 5
                    elseif (LocPos_XiTilde(3) <-1 ) then
                        CrossedSide = 1!side 1
                    elseif (LocPos_XiTilde(3) >1 ) then
                        CrossedSide = 4!side 4
                    end if 
                    
                    NewElementID = ElementAdjacencies(IElement, CrossedSide)
                    
                    counter_outside = counter_outside +1
                    
                    if ((NewElementID .ne. 0) .and. counter_outside>3) then 
            !.and. ((counter_xi_outsideelement>2) .or. (counter_eta_outsideelement>2) .or. (counter_zeta_outsideelement>2))) then
            OutsideElement = .true.
            Success = .false.
            LocPos = LocPos_XiTilde
            !exit
            return
                    else 
                        ! Dimension 1 (Xi)
                    if (LocPos_Xi(1) < XiKnotEntries(ni,IPatch) + tolerance_Xi(1)) then
                        LocPos_Xi(1) = XiKnotEntries(ni,IPatch) + tolerance_Xi(1)
                        counter_xi_outsideelement = counter_xi_outsideelement + 1
                    !test = ElementAdjacencies(1,1)
                    elseif (LocPos_Xi(1) > XiKnotEntries(ni+1,IPatch) - tolerance_Xi(1)) then
                        LocPos_Xi(1) = XiKnotEntries(ni+1,IPatch) - tolerance_Xi(1)
                        counter_xi_outsideelement = counter_xi_outsideelement + 1
                    endif

                    ! Dimension 2 (Eta)
                    if (LocPos_Xi(2) < EtaKnotEntries(nj,IPatch) + tolerance_Xi(2)) then
                        LocPos_Xi(2) = EtaKnotEntries(nj,IPatch) + tolerance_Xi(2)
                        counter_eta_outsideelement = counter_eta_outsideelement + 1
                    elseif (LocPos_Xi(2) > EtaKnotEntries(nj+1,IPatch) - tolerance_Xi(2)) then
                        LocPos_Xi(2) = EtaKnotEntries(nj+1,IPatch) - tolerance_Xi(2)
                        counter_eta_outsideelement = counter_eta_outsideelement + 1
                    endif

                    ! Dimension 3 (Zeta)
                    if (LocPos_Xi(3) < ZetaKnotEntries(nk,IPatch) + tolerance_Xi(3)) then
                        LocPos_Xi(3) = ZetaKnotEntries(nk,IPatch) + tolerance_Xi(3)
                        counter_zeta_outsideelement = counter_zeta_outsideelement + 1
                    elseif (LocPos_Xi(3) > ZetaKnotEntries(nk+1,IPatch) - tolerance_Xi(3)) then
                        LocPos_Xi(3) = ZetaKnotEntries(nk+1,IPatch) - tolerance_Xi(3)
                        counter_zeta_outsideelement = counter_zeta_outsideelement + 1
                    endif
                    
                    
                    
                    
                    !  ! Corresponding local position in the xi tilde parent domain
            LocPos_XiTilde(1) = ( (2*LocPos_Xi(1)) - (XiKnotEntries(ni+1,IPatch)  +  XiKnotEntries(ni,IPatch)) )    /(XiKnotEntries(ni+1,IPatch)-XiKnotEntries(ni,IPatch))
            LocPos_XiTilde(2) = ( (2*LocPos_Xi(2)) - (EtaKnotEntries(nj+1,IPatch) + EtaKnotEntries(nj,IPatch)) )    /(EtaKnotEntries(nj+1,IPatch)-EtaKnotEntries(nj,IPatch))
            LocPos_XiTilde(3) = ( (2*LocPos_Xi(3)) - (ZetaKnotEntries(nk+1,IPatch)+ZetaKnotEntries(nk,IPatch)) )    /(ZetaKnotEntries(nk+1,IPatch)-ZetaKnotEntries(nk,IPatch))
            
            
                    
                    end if
                    
            
                
                
                    
            endif
            
              
                    ! assign new element
                    
                    
                    
                    
            
            !if (.not.IsInsideElementLocPos(LocPos_XiTilde) .and. (NewElementID .ne. 0)) then 
            !!.and. ((counter_xi_outsideelement>2) .or. (counter_eta_outsideelement>2) .or. (counter_zeta_outsideelement>2))) then
            !OutsideElement = .true.
            !Success = .false.
            !LocPos = LocPos_XiTilde
            !!exit
            !return
            !end if
            

            


            ! Check and adjust LocPos_Xi for each dimension
! if outside knot span 
              !if (.not.IsInsideElementLocPos(LocPos_XiTilde) ) then
              !      ! based on which side you crossed you can knot which element you are in
              !      if (LocPos_XiTilde(1) <-1 ) then 
              !          CrossedSide = 3!side 3
              !      elseif (LocPos_XiTilde(1) >1 ) then
              !          CrossedSide = 2!side 2
              !      elseif (LocPos_XiTilde(2) <-1 ) then 
              !          CrossedSide = 6!side 6
              !      elseif (LocPos_XiTilde(2) >1 ) then
              !          CrossedSide = 5!side 5
              !      elseif (LocPos_XiTilde(3) <-1 ) then
              !          CrossedSide = 1!side 1
              !      elseif (LocPos_XiTilde(3) >1 ) then
              !          CrossedSide = 4!side 4
              !      end if 
              !
              !      ! assign new element
              !      NewElementID = ElementAdjacencies(IElement, CrossedSide)
              !      
              !      
              !      
              !      if (NewElementID == 0) then 
                    
                    
                    
                    !end if
            
            !if ( LocPos_Xi(1)<0 ) then  !(iteration < 3) .and. 
            !    LocPos_Xi(1) = XiKnotEntries(ni)    
            !end if 
            !if ( LocPos_Xi(2)<0 ) then !(iteration < 3) .and. 
            !    LocPos_Xi(2) = EtaKnotEntries(nj)   
            !end if 
            !if ( LocPos_Xi(3)<0 ) then !(iteration < 3) .and. 
            !    LocPos_Xi(3) = ZetaKnotEntries(nk)    
            !end if 
            !if ( LocPos_Xi(1)>1 ) then !(iteration < 3) .and. 
            !    LocPos_Xi(1) = XiKnotEntries(ni+1)      
            !end if 
            !if ( LocPos_Xi(2)>1 ) then !(iteration < 3) .and. 
            !    LocPos_Xi(2) = EtaKnotEntries(nj+1) 
            !end if 
            !if ( LocPos_Xi(3)>1 ) then  !(iteration < 3) .and. 
            !    LocPos_Xi(3) = ZetaKnotEntries(nk+1)    
            !end if 
                    
                    
                    
                    
                    
                    
                    
                    
                    !! automatically assign success 
                    !Success = .true.
                    !
                    !! reset local position in the new element
                    !NewLocPos = 0.0
              !end if ! if outside knot span 
            
            
            
            
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Correct if local particle position is outside element due to rounding
            !if (LocPos(1)<-1) then
            !    LocPos(1) = -1
            !    
            !    !if (CalParams%OutputDebugData) then
            !    !  call WriteInLogFile('Particle neg loc pos (x) '// trim(String(IParticle))    //' '// &
            !    !                                                trim(String(NewElementID)) //' '// &
            !    !                                                trim(String(LocPos(1))))
            !    !end if
            !end if
            !
            !if (LocPos(1)>1) then
            !    LocPos(1) = 1
            !    !
            !    !if (CalParams%OutputDebugData) then
            !    !  call WriteInLogFile('Particle neg loc pos (x) '// trim(String(IParticle))    //' '// &
            !    !                                                trim(String(NewElementID)) //' '// &
            !    !                                                trim(String(LocPos(1))))
            !    !end if
            !end if
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !if (LocPos(2)<-1) then
            !    LocPos(2) = -1
            !    
            !    !if (CalParams%OutputDebugData) then
            !    !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
            !    !                                                trim(String(NewElementID)) //' '// &
            !    !                                                trim(String(LocPos(2))))
            !    !end if
            !end if
            !
            !
            !if (LocPos(2)>1) then
            !    LocPos(2) = 1
            !    
            !    !if (CalParams%OutputDebugData) then
            !    !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
            !    !                                                trim(String(NewElementID)) //' '// &
            !    !                                                trim(String(LocPos(2))))
            !    !end if
            !end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !if (LocPos(3)<-1) then
            !    LocPos(3) = -1
            !    
            !    !if (CalParams%OutputDebugData) then
            !    !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
            !    !                                                trim(String(NewElementID)) //' '// &
            !    !                                                trim(String(LocPos(2))))
            !    !end if
            !end if
            !
            !
            !if (LocPos(3)>1) then
            !    LocPos(3) = 1
            !    
            !    !if (CalParams%OutputDebugData) then
            !    !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
            !    !                                                trim(String(NewElementID)) //' '// &
            !    !                                                trim(String(LocPos(2))))
            !    !end if
            !end if
              
            
            ! Global coordinates for guess
            call GetGlobalCoordinates(LocPos_XiTilde,  &
                                      IElTyp, NEl, NodTot, IDim, &
                                      IElement, ICon, NodeCoord, &
                                      GlobPosIteration, &
                                      IPatch)!, &
            !IncrementalDisplacementSoil)
        
            ! Check whether loop can be aborted
            Difference = 0.0
            do I = 1, IDim
              Difference = Difference + dabs(GlobPos(I) - GlobPosIteration(I))
            end do
            
            
            
     
            if ((Difference < Tolerance) ) then ! Found local coordinates
              Success = .true.
              LocPos = LocPos_XiTilde
              EXIT
            else ! Difference greater or equal Tolerance
              if (Iteration >= 30) then ! Too many iterations needed, something went wrong
                Success = .false.
                EXIT
              else ! Continue iteration
                Iteration = Iteration + 1
              end if
            end if
          end do ! Iteration loop

          

          if (.not.Success) then
            !call GiveError('Did not find local coordinates of particle in element'// &
            !               trim(String(IElement)) // &
            !               ' within limit number of iterations.')
              
              aa = 1 ! place holder for debugging
              
              OutsideElement = .true.
            Success = .false.
              
          end if

        end subroutine GetLocalCoordinates
!KnotBezierMesh

        subroutine GetLocalCoordinates3(GlobPos, IDim, IElement, IElTyp, NEl, ICon, LocPos, OutsideElement, MInv, MIX1, CrossedSide)
        ! This subroutine is never called...
        !**********************************************************************
        !
        !    Function:  Determination of local coordinates from global coordinates,
        !               assuming that the point with global coordinates lies inside 
        !               the element within IElement. Only for 3D
        !               Success returns .false. if the local position could not be found
        !               inside the element within 10 iterations with sufficient accuracy.
        !
        !     GlobPos : Global coordinates of a point inside IElement
        !     IDim : dimension
        !     IElement : ID of the element
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   LocPos : Local coordinates of the considered point inside IElement
        ! O   OutsideElement : True, if the local coordinate lie outside IElement
        ! O   Success : True, if the local position could be found
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IDim
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(IDim), intent(out) :: LocPos
          logical, intent(out) :: OutsideElement
          ! Local variables
          integer(INTEGER_TYPE) , dimension(IElTyp) :: NodeID
          integer(INTEGER_TYPE) :: J, INode
          real(REAL_TYPE), dimension(IDIm, IDim) :: MInv
          real(REAL_TYPE), dimension(IDim) :: MIX1
          integer(INTEGER_TYPE) :: CrossedSide
          
          OutsideElement = .true.

          do INode = 1, IElTyp
            NodeID(INode) = ICon(INode, IElement)
          end do

          
          LocPos = 0.0
          do J = 1, 3
            LocPos(3) = LocPos(3) + MInv(1, J) *  GlobPos(J)
            LocPos(1) = LocPos(1) + MInv(2, J) *  GlobPos(J)
            LocPos(2) = LocPos(2) + MInv(3, J) *  GlobPos(J)
          end do
         
          LocPos(3) = LocPos(3) - MIX1(1)
          LocPos(1) = LocPos(1) - MIX1(2)
          LocPos(2) = LocPos(2) - MIX1(3)

         CrossedSide = -1
        if     (LocPos(1) < 0) then
         CrossedSide = 2
        elseif (LocPos(2) < 0) then
         CrossedSide = 1
        elseif (LocPos(3) < 0) then
         CrossedSide = 3
        elseif ((LocPos(1) + LocPos(2) + LocPos(3)) > 1) then
         CrossedSide = 4
        else
        OutsideElement = .false.
        end if
        
        end subroutine GetLocalCoordinates3

        
        subroutine GetGlobalCoordinates(LocPos_XiTilde, IElTyp, NEl, NodTot, IDim, IElement, ICon, NodeCoord, GlobPos, IPatch)!, IncrementalDisplacementSoil)
        !**********************************************************************
        !
        !    Function:  Determination of global coordinates from global coordinates.
        !
        !     LocPos : Local coordinates of the considered point inside an element
        !     IElTyp : Number of node connectivities of IElement
        !     NEl : Number of elements
        !     NodTot : Total number of nodes
        !     IDim : Dimension of the mesh
        !     IElement : ID of the element that the point is located in
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !     NodeCoord : Global nodal coordinates Co(I, J): j-coordinate of node I
        !
        ! O   GlobPos : Global coordinates of the point inside IElement
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NEl, NodTot, IDim
          real(REAL_TYPE), dimension(IDim), intent(inout) :: LocPos_XiTilde
          integer(INTEGER_TYPE), intent(in) :: IElement
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          !real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: IncrementalDisplacementSoil
          real(REAL_TYPE), dimension(IDim), intent(out) :: GlobPos
          ! Local variables
          real(REAL_TYPE), allocatable, dimension(:) :: ShapeValues ! Shape functions
          real(REAL_TYPE), allocatable, dimension(:, :) :: DShapeValues ! Derivatives of shape functions
          !real(REAL_TYPE), dimension(IElTyp) :: ShapeValues ! Shape functions
          !real(REAL_TYPE), dimension(IElTyp, IDim) :: DShapeValues ! Derivatives of shape functions
          integer(INTEGER_TYPE) :: I, J, NodeID
          
          ! Multipatch variables
          integer(INTEGER_TYPE) :: IPatch_Temporary
          integer(INTEGER_TYPE), intent(in) :: IPatch

          
          !do II = 1,NVECTOR 
          !    
          !    if (LosPos(II) < -1) then 
          !        LocPos
          !    
          !end do 
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  ! Correct if local particle position is outside element due to rounding
          !  if (LocPos(1)<-0.95) then
          !      LocPos(1) = -0.95
          !      
          !      !if (CalParams%OutputDebugData) then
          !      !  call WriteInLogFile('Particle neg loc pos (x) '// trim(String(IParticle))    //' '// &
          !      !                                                trim(String(NewElementID)) //' '// &
          !      !                                                trim(String(LocPos(1))))
          !      !end if
          !  end if
          !  
          !  if (LocPos(1)>0.95) then
          !      LocPos(1) = 0.95
          !      
          !      !if (CalParams%OutputDebugData) then
          !      !  call WriteInLogFile('Particle neg loc pos (x) '// trim(String(IParticle))    //' '// &
          !      !                                                trim(String(NewElementID)) //' '// &
          !      !                                                trim(String(LocPos(1))))
          !      !end if
          !  end if
          !  
          !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !
          !  if (LocPos(2)<-1) then
          !      LocPos(2) = -0.95
          !      
          !      !if (CalParams%OutputDebugData) then
          !      !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
          !      !                                                trim(String(NewElementID)) //' '// &
          !      !                                                trim(String(LocPos(2))))
          !      !end if
          !  end if
          !  
          !  
          !  if (LocPos(2)>1) then
          !      LocPos(2) = 0.95
          !      
          !      !if (CalParams%OutputDebugData) then
          !      !  call WriteInLogFile('Particle neg loc pos (y) '// trim(String(IParticle))    //' '// &
          !      !                                                trim(String(NewElementID)) //' '// &
          !      !                                                trim(String(LocPos(2))))
          !      !end if
          !  end if
          !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          
          
          
          ! Determine the shape functions HS and shape function derivatives dHS for LocPos.
          call ShapeFunctionData(LocPos_XiTilde, ELEMENTNODES, ShapeValues, DShapeValues, IElement, IPatch)

          GlobPos = 0.0
          do I = 1, IElTyp
            NodeID = ICon(I, IElement)
            do J = 1, IDim
              GlobPos(J) = GlobPos(J) + (ShapeValues(I) * NodeCoord(NodeID, J))
              !GlobPos(J) = GlobPos(J) + (ShapeValues(I) * (IncrementalDisplacementSoil(NodeID, J)))!)
            end do  
          end do

        end subroutine GetGlobalCoordinates

        
        logical function IsInsideElementLocPos(LocPos)
        !**********************************************************************
        !
        !    Function:  Returns .true. if LocPos (local coordinates) lies inside the 
        !               volume of the considered element.
        !
        !     IElTyp : Number of node connectivities of IElement
        !
        ! O   IsInsideElementLocPos : True, if the point lies inside the element.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          
        IsInsideElementLocPos = IsInsideElementLocPosPointer(LocPos)
     
        end function IsInsideElementLocPos

        
        logical function IsInsideElementGlobPos(GlobPos, ElementID, NodTot, IDim, IElTyp, NEl, NodeCoord, ICon)
        !**********************************************************************
        !
        !    Function:  Returns .true. if GlobPos (global coordinates) lies inside the 
        !               volume of the considered element.
        !
        !     GlobPos : Global coordinates of the considered point inside an element
        !     ElementID : ID of the considered element
        !     NodTot : Total number of nodes
        !     IDim : Number of dimensions
        !     IElTyp : Number of nodes per element
        !     NEl : Total number of elements
        !     NodeCoord : Nodal coordinates
        !     ICon : Element connectivities
        !
        ! O   IsInsideElementGlobPos : True, if the point lies inside the element.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NodTot, IDim, IElTyp, NEl
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
        
          select case(IElTyp) 
            case(10) ! 10-noded tetrahedral element
              IsInsideElementGlobPos =  IsInsideElementGlobPosTETRA(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
            case(4) ! 4-noded tetrahedral element
              IsInsideElementGlobPos =  IsInsideElementGlobPosTETRA(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
          end select

        end function IsInsideElementGlobPos

        
        logical function IsCornerNode(INode, IElTyp)
        !**********************************************************************
        !
        !    Function:  Returns .true. if INode is a corner node.
        !
        !     INode : Local number of a node (1 .. IElTyp)
        !     IElTyp : Number of node connectivities
        !
        ! O   IsCornerNode : True, if INode is a corner node
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: INode
          integer(INTEGER_TYPE), intent(in) :: IElTyp
          
          IsCornerNode = .false.
          
          select case(ELEMENTTYPE)
              
            case(TRI3) ! 'triangular_3-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
              
            case(TRI6) ! 'triangular_6-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) ) IsCornerNode = .true. ! first three nodes are corner nodes
              
            case(QUAD4) ! 'quadrilateral_4-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
              
            case(QUAD8) ! 'quadrilateral_8-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) ) IsCornerNode = .true. ! first four nodes are corner nodes 
            
            !NURBS - it depends on the order --> need to generalize this here   
            case(QUAD4_NURBS)             
                IsCornerNode = .true. ! all nodes are corner nodes -> this is not correct and need a better way to link this to the order of the shape functions

              
            case(TETRA4) ! 'tetrahedral_4-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
                         
            case(TETRA10) ! 'tetrahedral_10-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) ) IsCornerNode = .true. ! first four nodes are corner nodes
                          
            case(HEXA8) ! 'hexahedral_8-noded'
              IsCornerNode = .true. ! all nodes are corner nodes
              
            case(HEXA20) ! 'hexahedral_20-noded'
              if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) .or. &
                   (INode == 5) .or. (INode == 6) .or. (INode == 7) .or. (INode == 8) ) IsCornerNode = .true. ! first eight nodes are corner nodes 
         
              
              
              
              
              
              
              
  
            case(HEXA_NURBS)
                !IsCornerNode = .true. ! all nodes are corner nodes --> this is not true and depends on the element order
                ! it seems that the first 8 nodes are corner nodes and the rest are inner nodes
                if ( (INode == 1) .or. (INode == 2) .or. (INode == 3) .or. (INode == 4) .or. &
                    (INode == 5) .or. (INode == 6) .or. (INode == 7) .or. (INode == 8) ) IsCornerNode = .true. ! first eight nodes are corner nodes 
              
              
              
              
              
              
              
              
              
              
            case(TETRAOLD)
              select case(IElTyp) 
                case(10) ! 10-noded tetrahedral element
                  if ( (INode == 1).or. (INode == 2).or. (INode == 3).or. (INode == 4) ) IsCornerNode = .true.
                case(4) ! 4-noded tetrahedral element
                  if ( (INode==1).or. (INode==2).or. (INode==3).or. (INode==4) ) IsCornerNode = .true.
              end select
              
          case default
            call GiveError('Element type not defined in function IsCornerNode().')
            
          end select  
        
        end function IsCornerNode

        
        subroutine RearrangeConnectivitiesTRI3(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 3-noded triangular element.
        !               For this element type no re-ordering is required as only corner nodes exist.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(2)
          IConLocal(3) = IConGlobal(3)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(2, :)
          ValuesLocal(3, :) = ValuesGlobal(3, :)

        end subroutine RearrangeConnectivitiesTRI3
        
        
        subroutine RearrangeConnectivitiesTRI6(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 6-noded triangular element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I    IConGlobal : Node connectivities
        ! I    ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO   IConLocal : Rearranged node connectivities
        ! IO   ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(5)
          IConLocal(4) = IConGlobal(2)
          IConLocal(5) = IConGlobal(4)
          IConLocal(6) = IConGlobal(6)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(5, :)
          ValuesLocal(4, :) = ValuesGlobal(2, :)
          ValuesLocal(5, :) = ValuesGlobal(4, :)
          ValuesLocal(6, :) = ValuesGlobal(6, :)

        end subroutine RearrangeConnectivitiesTRI6

        
        subroutine RearrangeConnectivitiesLINE2(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 2-noded line element.
        !               For this element type no re-ordering is required as only corner nodes exist.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(2)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(2, :)

        end subroutine RearrangeConnectivitiesLINE2
        
        
        subroutine RearrangeConnectivitiesLINE3(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 3-noded line element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(2)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(2, :)

        end subroutine RearrangeConnectivitiesLINE3
        
        
        subroutine RearrangeConnectivitiesQUAD4(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 4-noded quadrilateral element.
        !               For this element type no re-ordering is required as only corner nodes exist.
        !
        ! I   IConGlobal : Node connectivities
        ! I   ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO  IConLocal : Rearranged node connectivities
        ! IO  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(2)
          IConLocal(3) = IConGlobal(3)
          IConLocal(4) = IConGlobal(4)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(2, :)
          ValuesLocal(3, :) = ValuesGlobal(3, :)
          ValuesLocal(4, :) = ValuesGlobal(4, :)

        end subroutine RearrangeConnectivitiesQUAD4 
        
        
        subroutine RearrangeConnectivitiesQUAD8(IConGlobal, ValuesGlobal, IConLocal, ValuesLocal)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal for a 8-noded quadrilateral element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I    IConGlobal : Node connectivities
        ! I    ValuesGlobal : Values belonging to the connectivities of IConGlobal
        !
        ! IO   IConLocal : Rearranged node connectivities
        ! IO   ValuesLocal : Rearranged node values
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), dimension(:), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(:, :), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(:), intent(inout) :: IConLocal
          real(REAL_TYPE), dimension(:, :), intent(inout) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(5)
          IConLocal(4) = IConGlobal(7)
          IConLocal(5) = IConGlobal(2)
          IConLocal(6) = IConGlobal(4)
          IConLocal(7) = IConGlobal(6)
          IConLocal(7) = IConGlobal(8)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(5, :)
          ValuesLocal(4, :) = ValuesGlobal(7, :)
          ValuesLocal(5, :) = ValuesGlobal(2, :)
          ValuesLocal(6, :) = ValuesGlobal(4, :)
          ValuesLocal(7, :) = ValuesGlobal(6, :)
          ValuesLocal(7, :) = ValuesGlobal(8, :)

        end subroutine RearrangeConnectivitiesQUAD8
        
        
        function IntegrateVectorSurface(NSurfaceNodes, NGP, NodTot, NodeCoord, IConSurface, NodeValues, IElements, NumberOfTractionElements, NURBSTractionNodes)
        !**********************************************************************
        !
        !    Function:  Returns the integral of the distributed values defined by NodeValues
        !               for a 6-noded triangular surface (3D) and 2-noded line (2D) by Gauss point integration.
        !
        !     NSurfaceNodes : Number of nodes of the surface (=6 for 6-noded triangular element, =2 for 2-noded linear element)
        !     NGP : Number of Gauss points
        !     NodTot : Total number of nodes
        !     NodeCoord : Nodal coordinates
        !     IConSurface : Surface connectivities
        !     NodeValues : Node values integrated over the surface
        !
        ! O   IntegrateVectorSurface : Integral over surface of NodeValues
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
          implicit none
        
          integer(INTEGER_TYPE), intent(in) :: NSurfaceNodes, NGP, NodTot
          real(REAL_TYPE), dimension(NodTot, 3), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConSurface
          real(REAL_TYPE), dimension(NSurfaceNodes, NDOFL), intent(in) :: NodeValues
          !real(REAL_TYPE), allocatable, dimension(:, :) :: NodeValues_AllNodes
          
          real(REAL_TYPE), dimension(NDOFL) :: NodeValues_Debug ! hardcoded
          
          real(REAL_TYPE), dimension(NDOFL,4) :: IntegrateVectorSurface
          ! Local variables
          real(REAL_TYPE), dimension(NDOFL, 4) :: GPValue 

          real(REAL_TYPE), dimension(NDOFL) :: VectorN
          integer(INTEGER_TYPE) :: IGP, INode, IDim, nNode
          !real(REAL_TYPE) :: DetJ
          integer(INTEGER_TYPE), intent(in) :: NumberOfTractionElements
          integer(INTEGER_TYPE), dimension(:), intent(in) :: NURBSTractionNodes
          
          !integer(INTEGER_TYPE), allocatable, dimension(:) :: NURBSTractionNodes !need to update
          integer(INTEGER_TYPE), allocatable, dimension(:,:) :: IConLElement !need to update
          integer(INTEGER_TYPE) :: ii 
          integer(INTEGER_TYPE), dimension(:) :: IElements
          
          
          real(REAL_TYPE), dimension(4) :: DetJ ! this represents 
          
          IntegrateVectorSurface = 0.0
          
          NodeValues_Debug(1) = NodeValues(1,1) 
          NodeValues_Debug(2) = NodeValues(1,2)

          ! hardcoded
          !allocate(NURBSTractionNodes(IConSurface(2)))
          !NURBSTractionNodes(1) = 1
          !!NURBSTractionNodes(2) = 2
          !!NURBSTractionNodes(3) = 3
          !NumberOfTractionElements = 1
          !
          !!IConSurface(1) -> start 
          !!IConSurface(2) -> end 
          !
          !do ii = 2,IConSurface(2)
          !    NURBSTractionNodes(ii) = ii     
          !end do 
              
              
              
          
          
          
          
          !do IGP = 1, NGP

              
              ! All these GPShapeFunctionDerivativeBoundary need to change to DShapeValuesArray otherwise they will stay the same throughout the calculation
              
              
              
            ! Determine determinante of the Jacobian
            select case(NSurfaceNodes)
              case(6) ! 6-noded triangular boundary element, only for 3D
    !            nNode = 3
    !            !call Normal_T3(IGP, NodeCoord, IConSurface, nNode, GPShapeFunctionDerivativeBoundary, VectorN, DetJ)
				!DetJ = DetJ * 0.5
			  case(2) ! 2-noded linear element, for 2D
                nNode = 4
                call NormalOnLine(IGP, NodeCoord, IConSurface, GPShapeFunctionDerivativeBoundary_XI_AllElements, VectorN, DetJ, NumberOfTractionElements, NURBSTractionNodes, IElements)!, IConLElement)
            end select
          
            ! Determine distributed value at Gauss point
            GPValue = 0.0
            
            do IGP = 1, NGP
            do INode = 1, nNode 
              do IDim = 1, NDOFL
                GPValue(IDim, INode) = GPValue(IDim, INode) + GPShapeFunctionBoundary(IGP, INode) * NodeValues_Debug(IDim)* DetJ(INode)*2 !NodeValues_Debug
                !GPValue(IDim) = GPValue(IDim) + GPShapeFunctionBoundary(IGP, INode) * NodeValues(INode, IDim) !NodeValues_Debug
                ! hardcoding GPWeightBoundary to 2
              end do
            end do
            end do
          
            !! Integrate over the element surface
            !do IDim = 1, NDOFL
            !  IntegrateVectorSurface(IDim) = IntegrateVectorSurface(IDim) + GPValue(IDim) * GPWeightBoundary(IGP) !* DetJ(IDim) ! Hardcoded: IDim here needs to be fixed
            !end do
            IntegrateVectorSurface = GPValue

        
        end function IntegrateVectorSurface
        
        subroutine InitialiseShapeFunctions()
        !*************************************************************************************   
        !    FUNCTION:     Initialise shape functions
        ! 
        !    DESCRIPTION:        
        !>   Initialises the shape functions for all element types.
        !
        !>   @note: 
        !
        !>   @param[in] 
        !
        !>   @return
        !
        !*************************************************************************************
        implicit none
        
        real(REAL_TYPE), dimension(NXiGaussPoints) :: Xi_ParametricDomain, Eta_ParametricDomain, Zeta_ParametricDomain 
        real(REAL_TYPE), dimension(NXiGaussPoints) :: xi_tilde 
        real(REAL_TYPE), dimension(NEtaGaussPoints) :: eta_tilde
        real(REAL_TYPE), dimension(NEtaGaussPoints) :: zeta_tilde
        integer(INTEGER_TYPE) :: ni, nj, nk, ee
        
        integer(INTEGER_TYPE) :: gauss_point_counter, gp_zeta, gp_eta, gp_xi
        
        ! Multipatch variables
        integer(INTEGER_TYPE) :: IPatch_Temporary =1 
        integer(INTEGER_TYPE) :: IPatch
        
        allocate( GPWeight(ELEMENTGAUSSPOINTS) )
        allocate( GPShapeFunction(ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        allocate( GPShapeFunctionDerivative(ELEMENTGAUSSPOINTS, ELEMENTNODES, 1) )!NDOFL) ) ! --> 1 corresponds to one gauss point
        
        ! In NURBS we need these values for all elements
        allocate( GPWeight_AllElements(nel_NURBS(IPatch_Temporary), ELEMENTGAUSSPOINTS) )
        allocate( GPShapeFunction_AllElements(nel_NURBS(IPatch_Temporary), ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        allocate( GPShapeFunctionDerivative_AllElements(nel_NURBS(IPatch_Temporary), ELEMENTGAUSSPOINTS, ELEMENTNODES, 1))! NDOFL) )
        
        !nen_NURBS
        
        !allocate( GPWeight_AllElements(nel_NURBS, ELEMENTGAUSSPOINTS) )
        !allocate( GPShapeFunction_AllElements(nel_NURBS, ELEMENTGAUSSPOINTS, ELEMENTNODES) )
        !allocate( GPShapeFunctionDerivative_AllElements(nel_NURBS, ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
        
        ! allocate shape function variables. If considering NURBS, this would be: 
        ! Xi direction. 
        !allocate( GPWeightBoundary(ELEMENTBOUNDARYGAUSSPOINTS) )
        !allocate( GPShapeFunctionBoundary(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES) )
        !allocate( GPShapeFunctionDerivativeBoundary(ELEMENTBOUNDARYGAUSSPOINTS, ELEMENTBOUNDARYNODES, NDOFL-1) )
    
        ! Xi direction 
        allocate( GPWeightBoundary(ELEMENTBOUNDARYGAUSSPOINTS_XI) )
        allocate( GPShapeFunctionBoundary(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI) )
        allocate( GPShapeFunctionDerivativeBoundary(ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI, 1))! NDOFL-1) )
                
    
        ! Eta direction 
        allocate( GPWeightBoundary_ETA(ELEMENTBOUNDARYGAUSSPOINTS_ETA) )
        allocate( GPShapeFunctionBoundary_ETA(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA) )
        allocate( GPShapeFunctionDerivativeBoundary_ETA(ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA, 1) )!NDOFL-1) )
        
        ! Zeta direction 
        allocate( GPWeightBoundary_ZETA(ELEMENTBOUNDARYGAUSSPOINTS_ZETA) )
        allocate( GPShapeFunctionBoundary_ZETA(ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA) )
        allocate( GPShapeFunctionDerivativeBoundary_ZETA(ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA, 1)) !NDOFL-1) )
        
        ! ------------------ all elements needed for traction implementation -------------------------------------------------
        
        !nel_NURBS = (nn_NURBS_NumberOfUnivariateXiKnots-NXiKnotOrder) * (mm_NURBS_NumberOfUnivariateEtaKnots-NEtaKnotOrder) 
        
        ! Xi direction - all elements 
        allocate( GPWeightBoundary_XI_AllElements((nn_NURBS_NumberOfUnivariateXiKnots(IPatch_Temporary)-NXiKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_XI) )
        allocate( GPShapeFunctionBoundary_XI_AllElements((nn_NURBS_NumberOfUnivariateXiKnots(IPatch_Temporary)-NXiKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI) )
        allocate( GPShapeFunctionDerivativeBoundary_XI_AllElements((nn_NURBS_NumberOfUnivariateXiKnots(IPatch_Temporary)-NXiKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_XI, ELEMENTBOUNDARYNODES_XI, 1)) !NDOFL-1) )
        
        
        ! Eta direction - all elements  
        allocate( GPWeightBoundary_ETA_AllElements((mm_NURBS_NumberOfUnivariateEtaKnots(IPatch_Temporary)-NEtaKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_ETA) )
        allocate( GPShapeFunctionBoundary_ETA_AllElements((mm_NURBS_NumberOfUnivariateEtaKnots(IPatch_Temporary)-NEtaKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA) )
        allocate( GPShapeFunctionDerivativeBoundary_ETA_AllElements((mm_NURBS_NumberOfUnivariateEtaKnots(IPatch_Temporary)-NEtaKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_ETA, ELEMENTBOUNDARYNODES_ETA, 1)) !NDOFL-1) )
        
        
        ! Zeta direction - all elements  
        allocate( GPWeightBoundary_ZETA_AllElements((oo_NURBS_NumberOfUnivariateZetaKnots(IPatch_Temporary)-NZetaKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_ZETA) )
        allocate( GPShapeFunctionBoundary_ZETA_AllElements((oo_NURBS_NumberOfUnivariateZetaKnots(IPatch_Temporary)-NZetaKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA) )
        allocate( GPShapeFunctionDerivativeBoundary_ZETA_AllElements((oo_NURBS_NumberOfUnivariateZetaKnots(IPatch_Temporary)-NZetaKnotOrder(IPatch_Temporary)), ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA, 1))!NDOFL-1) )
        
        !allocate( GPWeightBoundary_ZETA(ELEMENTBOUNDARYGAUSSPOINTS_ZETA) )
        !allocate( GPShapeFunctionBoundary_ZETA(ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA) )
        !allocate( GPShapeFunctionDerivativeBoundary_ZETA(ELEMENTBOUNDARYGAUSSPOINTS_ZETA, ELEMENTBOUNDARYNODES_ZETA, NDOFL-1) )
        
        ! ------------------ all elements needed for traction implementation -------------------------------------------------
        
        ! initialise shape function variables
        
        GPWeight = 0.0
        GPShapeFunction = 0.0 
        GPShapeFunctionDerivative = 0.0
        
        ! Xi direction 
        GPWeightBoundary = 0.0
        GPShapeFunctionBoundary = 0.0
        GPShapeFunctionDerivativeBoundary = 0.0
        
        ! Eta direction 
        GPWeightBoundary_ETA = 0.0
        GPShapeFunctionBoundary_ETA = 0.0
        GPShapeFunctionDerivativeBoundary_ETA = 0.0
        
        ! Zeta direction 
        GPWeightBoundary_ZETA = 0.0
        GPShapeFunctionBoundary_ZETA = 0.0
        GPShapeFunctionDerivativeBoundary_ZETA = 0.0
        
        ! note that all of these derivatives are evaluated at the center of the control points for a knot span...
        ! when we evaluate stresses we need to calculate them at the material points and not at the gauss points...
        ! this is a current bug in the evaulation of B when we calculate the internal forces...
        xi_tilde = 0.0
        eta_tilde = 0.0
        zeta_tilde = 0.0
        gauss_point_counter = 0 
        
        
        ! loop over elements here 
        GPShapeFunctionDerivative_AllElements = 0.0 
        GPShapeFunction_AllElements = 0.0 
        GPWeight_AllElements = 0.0
        
        GPShapeFunctionDerivativeBoundary_ETA_AllElements = 0
        GPShapeFunctionDerivativeBoundary_XI_AllElements = 0
        GPShapeFunctionDerivativeBoundary_ZETA_AllElements = 0
        
        !ElementConnectivities = IEN
        
        
        ! Xi direction 
        ! 2-noded 1D element -> not sure why we need this...
        
        
        
        !loop over the gauss points for each element 
        do gp_zeta = 1, NZetaGaussPoints 
            do gp_eta = 1, NEtaGaussPoints
                do gp_xi = 1, NXiGaussPoints
                
                    
                    xi_tilde(gp_xi) = XiGaussPoints(gp_xi)
                    eta_tilde(gp_eta) = EtaGaussPoints(gp_eta)
                    if (NDIM==3) then 
                    zeta_tilde(gp_zeta) = ZetaGaussPoints(gp_zeta)
                    end if
                
                    
                end do
            end do 
        end do 
        
            
        
        do IPatch = 1, NumberOfPatches ! loop over patches
        ! loop over NURBS elements 
        do ee = 1, nel_NURBS(IPatch)!(IPatch_Temporary) !<- why is this causing problems? for when ee = 1
        
        ! temporary fix 
        !ee = nel_NURBS
        
                   
        
            
        !do ii = nn_NURBS_NumberOfUnivariateXiKnots ! loop over the number of knot spans in xi direction   
        !    do jj = mm_NURBS_NumberOfUnivariateEtaKnots ! loop over the number of knot spans in eta direction 
           ni = INN(IEN(1,ee,IPatch),1,IPatch)    
           nj = INN(IEN(1,ee,IPatch),2,IPatch)
           !nk = INN(IEN(nen_NURBS,ee),3)
            
           if (NDIM == 3) then 
              nk = INN(IEN(1,ee,IPatch),3,IPatch)
          end if 
            
            Xi_ParametricDomain =  ( (XiKnotEntries(ni+1,IPatch) - XiKnotEntries(ni,IPatch) ) * xi_tilde &
                                    + (XiKnotEntries(ni+1,IPatch) + XiKnotEntries(ni,IPatch)) ) * 0.5
            
            Eta_ParametricDomain =  ( (EtaKnotEntries(nj+1,IPatch) - EtaKnotEntries(nj,IPatch) ) * eta_tilde &
                                    + (EtaKnotEntries(nj+1,IPatch) + EtaKnotEntries(nj,IPatch)) ) * 0.5 
            
            if (NDIM == 3) then 
                
                Zeta_ParametricDomain =  ( (ZetaKnotEntries(nk+1,IPatch) - ZetaKnotEntries(nk,IPatch) ) * zeta_tilde &
                                    + (ZetaKnotEntries(nk+1,IPatch) + ZetaKnotEntries(nk,IPatch)) ) * 0.5
                
            end if 
            
            
            
            !call InitialiseShapeFunctionsBoundaryPointer_NURBS(GPShapeFunctionBoundary, GPShapeFunctionDerivativeBoundary, GPWeightBoundary, &
            !                                                XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder)
            !
            !call InitialiseShapeFunctionsBoundaryPointer_NURBS(GPShapeFunctionBoundary_Eta, GPShapeFunctionDerivativeBoundary_Eta, GPWeightBoundary_Eta, &
            !                                                EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder)
            
            
            ! loop over gauss points now 
            !do gp_eta = 1, NEtaGaussPoints
            !    do gp_xi = 1, NXiGaussPoints
            !
            ! gauss_point_counter = gauss_point_counter + 1

            if (NDIM == 2) then 
            
            
            call InitialiseShapeFunctionsQUAD4_NURBS(GPShapeFunction, GPShapeFunctionDerivative, GPWeight, & !classic inout parameters
                                                    GPShapeFunctionBoundary, GPShapeFunctionDerivativeBoundary, GPWeightBoundary, &
                                                    GPShapeFunctionBoundary_ETA, GPShapeFunctionDerivativeBoundary_ETA, GPWeightBoundary_ETA, &
                                                    XiKnotEntries(:,IPatch), NXiKnotEntries(IPatch), Xi_ParametricDomain, NXiKnotOrder(IPatch), & !NURBS related inputs in the xi direction 
                                                    EtaKnotEntries(:,IPatch), NEtaKnotEntries(IPatch), Eta_ParametricDomain, NEtaKnotOrder(IPatch), &
                                                    ni, nj, &
                                                    IPatch) !NURBS related inputs in the eta direction 
            
            
            elseif (NDIM == 3) then 
                
                
                   call InitialiseShapeFunctionsHEXA_NURBS(GPShapeFunction, GPShapeFunctionDerivative, GPWeight, & !classic inout parameters
                                                    GPShapeFunctionBoundary, GPShapeFunctionDerivativeBoundary, GPWeightBoundary, &
                                                    GPShapeFunctionBoundary_ETA, GPShapeFunctionDerivativeBoundary_ETA, GPWeightBoundary_ETA, &
                                                    GPShapeFunctionBoundary_ZETA, GPShapeFunctionDerivativeBoundary_ZETA, GPWeightBoundary_ZETA, &
                                                    XiKnotEntries(:,IPatch), NXiKnotEntries(IPatch), Xi_ParametricDomain, NXiKnotOrder(IPatch), & !NURBS related inputs in the xi direction 
                                                    EtaKnotEntries(:,IPatch), NEtaKnotEntries(IPatch), Eta_ParametricDomain, NEtaKnotOrder(IPatch), &
                                                    ZetaKnotEntries(:,IPatch), NZetaKnotEntries(IPatch), Zeta_ParametricDomain, NZetaKnotOrder(IPatch), &
                                                    ni, nj, nk, &
                                                    ee, &
                                                    IPatch) !NURBS related inputs in the eta direction 
                
            end if 
            
            
            ! InitialiseShapeFunctionsQUAD4_NURBS
            ! Storing the information regarding the shape functions for all the elements 
            GPShapeFunctionDerivative_AllElements(ee, :, :, :) = GPShapeFunctionDerivative(:, :, :)  ! all of these are evaluated at the original gauss point of the knot span...
            
            ! allocate( GPShapeFunctionDerivative(nel_NURBS, ELEMENTGAUSSPOINTS, ELEMENTNODES, NDOFL) ) 
            ! (number of elements, number of gauss points, element nodes, whether xi(1) or eta(2) directions)
            ! after thought: we do not need to store this for every direction but we can store it at its tensor product evaluated value
            GPShapeFunction_AllElements(ee, :, :) = GPShapeFunction(:, :) ! this is the shape function derivative for all the material points 
            GPWeight_AllElements(ee, :) = GPWeight(:)
            
            ! Boundary elements for traction implementation 
            ! xi direction 
            if (ee <= (nn_NURBS_NumberOfUnivariateXiKnots(IPatch)-NXiKnotOrder(IPatch)) ) then 
            GPWeightBoundary_XI_AllElements(ee,:) = GPWeightBoundary
            GPShapeFunctionBoundary_XI_AllElements(ee,:,:) = GPShapeFunctionBoundary
            GPShapeFunctionDerivativeBoundary_XI_AllElements(ee,:,:,:) = GPShapeFunctionDerivativeBoundary
            end if
            
            ! eta direction ... eta direction is just not right
            if (ee <= (mm_NURBS_NumberOfUnivariateEtaKnots(IPatch)-NEtaKnotOrder(IPatch)) ) then 
            GPWeightBoundary_ETA_AllElements(ee,:) = GPWeightBoundary_ETA
            GPShapeFunctionBoundary_ETA_AllElements(ee,:,:) = GPShapeFunctionBoundary_ETA
            GPShapeFunctionDerivativeBoundary_ETA_AllElements(ee,:,:,:) = GPShapeFunctionDerivativeBoundary_ETA
            end if 
            
            ! zeta direction ...? 
            if (ee <= (oo_NURBS_NumberOfUnivariateZetaKnots(IPatch)-NZetaKnotOrder(IPatch)) ) then 
            GPWeightBoundary_ZETA_AllElements(ee,:) = GPWeightBoundary_ZETA
            GPShapeFunctionBoundary_ZETA_AllElements(ee,:,:) = GPShapeFunctionBoundary_ZETA
            GPShapeFunctionDerivativeBoundary_ZETA_AllElements(ee,:,:,:) = GPShapeFunctionDerivativeBoundary_ZETA
            end if 
             
            !GPShapeFunctionDerivative = 0.0
            !GPShapeFunction = 0.0
            !GPWeight = 0.0
            
            ! loop over the number of knot spans in xi direction 
        
        
        end do ! loop over the number of knot spans in eta direction 
        
        
        ! seperate these into two loops 
        ! one loop you loop over all the ni in all the elements
        ! Another loop is for you to loop over all the nj in all the element 
        ! or maybe just have them into one loop 
        
        
        
        
        
        !(HS, dHS, Wt, &               !classic inout parameters
        !                                                XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder) 
        !
        ! initialising the shape functions for one gauss point (2D/3D)
        !call InitialiseShapeFunctionsPointer_NURBS(GPShapeFunction, GPShapeFunctionDerivative, GPWeight)
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !
        !! Eta direction 
        !! 2-noded 1D element -> not sure why we need this...
        !call InitialiseShapeFunctionsBoundaryPointer_NURBS(GPShapeFunctionBoundary_ETA, GPShapeFunctionDerivativeBoundary_ETA, GPWeightBoundary_ETA)
        !! initialising the shape functions for one gauss point (2D/3D)
        !call InitialiseShapeFunctionsPointer_NURBS(GPShapeFunction_ETA, GPShapeFunctionDerivative_ETA, GPWeight_ETA)
        !
        !! Zeta direction 
        !! 2-noded 1D element -> not sure why we need this...
        !call InitialiseShapeFunctionsBoundaryPointer_NURBS(GPShapeFunctionBoundary_ZETA, GPShapeFunctionDerivativeBoundary_ZETA, GPWeightBoundary_ZETA)
        !! initialising the shape functions for one gauss point (2D/3D)
        !call InitialiseShapeFunctionsPointer_NURBS(GPShapeFunction_ZETA, GPShapeFunctionDerivative_ZETA, GPWeight_ZETA)
        !
        !!if (ELEMENTTYPE == QUAD4_NURBS) then 
        !!    
        !!    
        !!    
        !!    
        !!    ! initialise shape function variables
        !!    GPWeightBoundary_ETA = 0.0
        !!    GPShapeFunctionBoundary_ETA = 0.0
        !!    GPShapeFunctionDerivativeBoundary_ETA = 0.0
        !!    GPWeight_ETA = 0.0
        !!    GPShapeFunction_ETA = 0.0 
        !!    GPShapeFunctionDerivative_ETA = 0.0
        !!    
        !!    ! 2-noded 1D element -> not sure why we need this...
        !!    call InitialiseShapeFunctionsBoundaryPointer(GPShapeFunctionBoundary_ETA, GPShapeFunctionDerivativeBoundary_ETA, GPWeightBoundary_ETA)
        !!    ! initialising the shape functions for one gauss point (2D/3D)
        !!    call InitialiseShapeFunctionsPointer(GPShapeFunction_ETA, GPShapeFunctionDerivative_ETA, GPWeight_ETA)  
        !!    
        !!    
        !!    
        !!    
        !!    
        !!    
        !!    
        !!    
        !!    ! initialise shape function variables
        !!    GPWeightBoundary_ZETA = 0.0
        !!    GPShapeFunctionBoundary_ZETA = 0.0
        !!    GPShapeFunctionDerivativeBoundary_ZETA = 0.0
        !!    GPWeight_ZETA = 0.0
        !!    GPShapeFunction_ZETA = 0.0 
        !!    GPShapeFunctionDerivative_ZETA = 0.0
        !!    
        !!    ! 2-noded 1D element -> not sure why we need this...
        !!    call InitialiseShapeFunctionsBoundaryPointer(GPShapeFunctionBoundary_ZETA, GPShapeFunctionDerivativeBoundary_ZETA, GPWeightBoundary_ZETA)
        !!    ! initialising the shape functions for one gauss point (2D/3D)
        !!    call InitialiseShapeFunctionsPointer(GPShapeFunction_ZETA, GPShapeFunctionDerivative_ZETA, GPWeight_ZETA)  
        !!    
        !!    
        !!    
        !!    
        !!    
        !!    
        !!end if 
        !
        !    end do
        !end do
        
        end do ! loop over patches
        
        end subroutine !InitialiseShapeFunctions
              
      
      subroutine Normal_T3(Int, Co, IConL, IelTyp3, dHS, Vn, Vl) 
!***********************************************************************
!     Determine vector V normal to a plane based on 6-noded flat element
!     Three vertices (T3)
!     first determine A = (dx/dXi , dy/dXi , dz/dXi )
!                 and B = (dx/dEta, dy/dEta, dz/dEta)
!     V1 =  a2*b3 - a3*b2
!     V2 =  a3*b1 - a1*b3
!     V3 =  a1*b2 - a2*b1
!     VL is length of (V1,V2,V3)
!     Finally normalize V
!***********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: int, ieltyp3  
      real(REAL_TYPE), dimension(:, :, :), intent(in) :: dHS 
      real(REAL_TYPE), dimension(:, :), intent(in) :: Co 
      integer(INTEGER_TYPE), dimension(:), intent(in) :: IConL
      real(REAL_TYPE), dimension(:), intent(inout) :: Vn
      real(REAL_TYPE), intent(inout) :: Vl
      
      real(REAL_TYPE), dimension(3) :: A, B
      integer(INTEGER_TYPE) :: k, nn, j
      
      Vn = 0.0
      A = 0.0
      B = 0.0

      Do K=1,IelTyp3
        NN= IConL(K)
        Do J=1,3 
          A(J) = A(J) + dHS(Int,K,1)* Co(NN,J)  ! Sum d()/dXi
          B(J) = B(J) + dHS(Int,K,2)* Co(NN,J)  ! Sum d()/dEta
        End Do
      End Do

      ! Vn is cross product of A and B
      Vn = CrossProduct(A, B)

      ! Vl is vector-length of cross product Vn
      Vl = Length(Vn, 3)

      ! normalise vector Vn
      Vn = VectorNorm(Vn, 3)

      end subroutine Normal_T3 
      
      
      subroutine NormalOnLine(IntegPoint, Co, IConL, dHS, NormalVector, VectorLength, NumberOfTractionElements, NURBSTractionNodes, IElements)!, IConLElement)
      !***********************************************************************
      !     Determine size of the vector formed by a line based on 2-noded line element
      !     first determine A = (dx/dXi , dy/dXi , dz/dXi )
      !     VectorLength is length of A
      !***********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IntegPoint
      real(REAL_TYPE), dimension(:, :, :, :), intent(in) :: dHS ! input this for all the elements in the XI direction so that we can loop accross elements 
      real(REAL_TYPE), dimension(:, :), intent(in) :: Co
      integer(INTEGER_TYPE), dimension(:), intent(in) :: IConL
      real(REAL_TYPE), dimension(:), intent(inout) :: NormalVector
      real(REAL_TYPE), dimension(4), intent(inout) :: VectorLength !hardcoded... need to fix here
      
      ! Multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary =1
      
      
      ! local variables
      integer(INTEGER_TYPE) :: NodeNumber, DoF, K
      real(REAL_TYPE), dimension(NVECTOR) :: A
      real(REAL_TYPE) :: dxi_dtildexi
      
      integer(INTEGER_TYPE) :: ni, nj, ee, LeftNodeIndex, RightNodeIndex
      
      integer(INTEGER_TYPE), dimension(:) :: NURBSTractionNodes
      
      !integer(INTEGER_TYPE), dimension(NumberOfTractionElements, 2) :: IConLElement !, intent(inout)
      
      ! Abdelrahman Alsardi (1 August 2023): Hardcoded : I commented this line because IPatch needs to be an input argument and I 
      ! still have not done that yet. 
      !integer(INTEGER_TYPE), dimension(NumberOfTractionElements, NXiKnotOrder(IPatch_Temporary)+1) :: IConLElement !, intent(inout)

      
      integer(INTEGER_TYPE), intent(in) :: NumberOfTractionElements
      
      integer(INTEGER_TYPE), dimension(:), intent(in) :: IElements

      A = 0.0
      
      VectorLength(1) = 0.5*(Co(NURBSTractionNodes(1),1) - Co(NURBSTractionNodes(2),1))
      VectorLength(2) = 0.5*(Co(NURBSTractionNodes(1),1) - Co(NURBSTractionNodes(2),1)) + 0.5*(Co(NURBSTractionNodes(2),1) - Co(NURBSTractionNodes(3),1))
      VectorLength(3) = 0.5*(Co(NURBSTractionNodes(3),1) - Co(NURBSTractionNodes(4),1)) + 0.5*(Co(NURBSTractionNodes(2),1) - Co(NURBSTractionNodes(3),1))
      VectorLength(4) = 0.5*(Co(NURBSTractionNodes(3),1) - Co(NURBSTractionNodes(4),1))
      
      NormalVector(1) = 0
      NormalVector(2) = 1
      
      !dxi_dtildexi = ( XiKnotEntries(4) - XiKnotEntries(3) )/2 ! hardcoded --> need to know where you are in the element 
      
      !allocate(IConLElement(NumberOfTractionElements, 3))
      
      !do ee = 1, NumberOfTractionElements
      !LeftNodeIndex = ee
      !RightNodeIndex = ee+2
      !IConLElement(ee, 1:3) = NURBSTractionNodes(LeftNodeIndex:RightNodeIndex)
      !end do 
      !IConLElement(1,1) = 4
      !IConLElement(1,2) = 3
      !IConLElement(1,3) = 2
      
      !IConLElement(2,1) = 3
      !IConLElement(2,2) = 2
      !IConLElement(2,3) = 1
      
      !if (NXiKnotOrder == 3) then 
      !
      !    IConLElement(1,1) = 4
      !    IConLElement(1,2) = 3
      !    IConLElement(1,3) = 2
      !    IConLElement(1,4) = 1
      !
      !elseif (NXiKnotOrder == 2) then 
      !
      !    IConLElement(1,1) = 3
      !    IConLElement(1,2) = 2
      !    IConLElement(1,3) = 1
      !
      !elseif (NXiKnotOrder == 1) then
      !
      !    IConLElement(1,1) = 2
      !    IConLElement(1,2) = 1
      !
      !end if
      !
      !do ee = 1, NumberOfTractionElements
      !    A = 0.0
      !
      !    !LeftNodeIndex = ee
      !    !RightNodeIndex = ee+2
      !    !    
      !    !IConLElement(1:3) = NURBSTractionNodes(LeftNodeIndex:RightNodeIndex)
      !
      !    !IConElement = 
      !    do K = 1, NXiKnotOrder+1  ! 2D shape functions (hardcoded NURBS shapefunctions)
      !        
      !        
      !  
      !        NodeNumber = IConLElement(ee, K)
      !  
      !        !NodeNumber = ElementConnectivities
      !        
      !        
      !       !ni = INN(IEN(1,NodeNumber),1)
      !       !nj = INN(IEN(2,NodeNumber),2)
      !  
      !  ni = INN(IEN(1,IElements(ee)),1)
      !  nj = INN(IEN(2,IElements(ee)),2)
      !  
      !        
      !  ! applying this only in the xi direction ... need to generalize in eta direction too
      !  dxi_dtildexi = ( XiKnotEntries(ni+1) - XiKnotEntries(ni) )/2 ! hardcoded --> need to know where you are in the element 
      !
      !  !Need to include here the tributary length based on the node number 
      !  !i.e., 
      !  !between (1) and 2 only 
      !  !between 1,(2),3
      !  !between 2,(3),4
      !  !between (4), 1
      !  !note that in linear elements you only consider one half and this applies to get the forces on both nodes
      !  !but in NURBS, that is not the case... you don't know how much should be on either side. So it is not just simply 
      !  !a jacobian multiplication. 
      !  
      !  !VectorLength needs to be a vector where 
      !  !VectorLength(1) = (0.025/3 * 0.5) ! edge control point
      !  !VectorLength(2) = (0.025/3 * 0.5) + (0.025/3 * 0.5) ! middle control point
      !  !VectorLength(3) = (0.025/3 * 0.5) + (0.025/3 * 0.5) ! middle control point
      !  !VectorLength(4) = (0.025/3 * 0.5) ! edge control point
      !  
      !  ! dHS multiplied by the coordinates is going to get you a one. If you multiply this by dxi_dtildexi, you get dxi_dtildexi which is just half the knot span. 
      !  ! But it doesn't make sense when it comes 
      !  do DoF = 1, NVECTOR
      !      A(DoF) = A(DoF) + dHS(ee, IntegPoint, K, 1) * Co(NodeNumber, DoF) * dxi_dtildexi  ! Sum d()/dXi --> need to multiply by dxi/dtildexi
      !  end do
      !    end do
      !
      !
      !
      !NormalVector(1) = -A(2)
      !NormalVector(2) = A(1) !I added a minus here to make it point outwards... need to check 
      !
      !! length of normal vector
      !VectorLength = sqrt( NormalVector(1)*NormalVector(1) + NormalVector(2)*NormalVector(2) ) !* dxi_dtildexi !* dxi_dtildexi !* 0.25 ! hardcoded --> vector length 
      !
      !! normalising the normal vector
      !NormalVector = NormalVector / (VectorLength) ! hardcoded 2D implementation of this for quadratic element. Note that this is because of the node in the middle of each segment 
      !
      !end do

      end subroutine NormalOnLine
      
     subroutine DetermineSideNodes(ISide,SideNodes)
        !**********************************************************************
        !
        !    Function:  Returns the local nodeID of  ISide
        !               works for Tetrahedra and triangular element 
        !    ISide :      ID of the considered side 
        !    SideNodes : ID of the boundary node
        !
        !**********************************************************************

        implicit none

        integer(INTEGER_TYPE), intent(in) :: ISide
        integer(INTEGER_TYPE), dimension(:), intent(out):: SideNodes
        ! Local variables
        integer(INTEGER_TYPE):: J
      
        if ( NDIM == 3 ) then
          if ((ELEMENTTYPE == TETRAOLD).OR.(ELEMENTTYPE == TETRA4).OR.(ELEMENTTYPE == TETRA10)) then
            do J = 1,  ELEMENTBOUNDARYNODES
              SideNodes(J) = DetermineSideNodesTetrahedronHOE(ISide, J)
            end do
            else
               call GiveError('DetermineSideNodes not implemented for element type '//trim(String(ELEMENTTYPE))) 
            end if
        else if ( NDIM == 2 ) then
          if ((ELEMENTTYPE == TRI3).OR.(ELEMENTTYPE == TRI6)) then
            do J = 1,  ELEMENTBOUNDARYNODES  
              SideNodes(J) = DetermineSideNodesTRI6(ISide, J)
            end do            
          else
           call GiveError('determineSideNodes not implemented for element type '//trim(String(ELEMENTTYPE))) 
           end if
        end if
        
        end subroutine DetermineSideNodes
      
      end module ModElementEvaluation
      