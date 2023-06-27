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


      module ModElementEvaluationQUAD
      !**********************************************************************
      !
      !    DESCRIPTION:
      !    This module provides specific routines for evaluating quadrilateral elements.
      !
      !    $Revision: 7657 $
      !    $Date: 2018-08-15 09:54:30 +0200 (Wed, 15 Aug 2018) $
      !
      !**********************************************************************

      use ModGeometryMath
      use ModString
      use ModReadCalculationData
      use ModGlobalConstants
      use ModElementEvaluationTRI
      !use ModMeshInfo

      implicit none

      
      contains
    

        !**********************************************************************
        !
        !    SUBROUTINE: DetermineSideDataQUAD 
        !
        !    DESCRIPTION: 
        !>   Returns the normal vector of the side with SideID and a PointInitialLocalCoordinatesQUAD
        !>   on that plane. The normal vector is normalised and points inward.
        !    
        !>   @note : 2D element
        !
        !>   @param[in] ISide : Side ID of the element
        !
        !>   @param[out] PlaneNormal : Normal to the element side
        !>   @param[out] PlanePoint : Point on the element side
        !
        !**********************************************************************
        subroutine DetermineSideDataQUAD(ISide, PlaneNormal, PlanePoint)
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          
          integer(INTEGER_TYPE), intent(in) :: ISide
          real(REAL_TYPE), dimension(IDim), intent(out) :: PlaneNormal
          real(REAL_TYPE), dimension(IDim), intent(out) :: PlanePoint !-> I do not this being assigned at any point here
          
          PlaneNormal = 0.0 !initializing the vector with zeros 
          PlanePoint = 0.0 !initializing the vector with zeros
          
          ! 4 cases for each side
          select case(ISide)
            case(1)
              PlaneNormal(2) = 1.0 !horizontal side at bottom of square 
            case(2)
              PlaneNormal(1) = 1.0 !vertical side at left of square  
            case(3)
              PlaneNormal(2) = -1.0 !horizontal side at top of square
            case(4)
              PlaneNormal(1) = -1.0 !vertical side at right of square
            case default
              call GiveError("Undefined side number in [subroutine DetermineSideDataQUAD()].")
          end select
        
        end subroutine DetermineSideDataQUAD


        !**********************************************************************
        !
        !    FUNCTION: PointSideDistanceQUAD
        !
        !    DESCRIPTION:
        !>   Returns the minimum distance between side SideID and LocPos.
        !
        !>   @param[in] SideID : ID of the considered side (1..4)
        !>   @param[in] LocPos : Local coordinates of considered point
        !
        !>   @return PointSideDistanceQuadrilateral
        !
        !**********************************************************************
        real(REAL_TYPE) function PointSideDistanceQUAD(SideID, LocPos)
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element 
          
          integer(INTEGER_TYPE), intent(in) :: SideID
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          
          ! local variables
          real(REAL_TYPE), dimension(IDim) :: LineNormal, LinePoint

          call DetermineSideDataQUAD(SideID, LineNormal, LinePoint)
        
          PointSideDistanceQUAD = LinePointDistance(LineNormal, LinePoint, LocPos)
        
        end function PointSideDistanceQUAD


        !**********************************************************************
        !
        !    SUBROUTINE: InitialLocalMaterialPointCoordinatesQUAD
        !
        !    DESCRIPTION: 
        !>   Determines the local coordinates and integration weight assigned to material point with ID
        !>   IParticle which are returned through WeiGP and PosGP. Currently, all material points are placed at the same local positions.
        !
        !>   @param[in] IParticle : Number of the material point
        !>   @param[in] SolidPointsElement : Number of solid material points per element
        !>   @param[in] LiquidPointsElement : Number of liquid material points per element
        !>   @param[inout] WeiGP : Initial weight assigned to material point IParticle
        !>   @param[inout] PosGP : Initial local position of material point IParticle
        !
        !**********************************************************************
        subroutine InitialLocalMaterialPointCoordinatesQUAD(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP

          ! local variables
          integer(INTEGER_TYPE) :: ID

          ! first the solid material points are determined
          if ( (IParticle <= SolidPointsElement) .and. (SolidPointsElement > 0) ) then
            ID = IParticle
            select case(SolidPointsElement)
            case (1) ! 1 material point per element
              call InitialQUAD_MP1(ID, PosGP, WeiGP)
            case (4) ! 4 material points per element
              call InitialQUAD_MP4(ID, PosGP, WeiGP)
            case (9) ! 9 material points per element 
              call InitialQUAD_MP9(ID, PosGP, WeiGP) 
            case (16) ! 16 material points per element 
              call InitialQUAD_MP16(ID, PosGP, WeiGP)  
            case (25) ! 25 material points per element
                call InitialQUAD_MP25(ID, PosGP, WeiGP)
            case (36) ! 36 material points per element    
                call InitialQUAD_MP36(ID, PosGP, WeiGP)
            case (49) ! 49 material points per element    
                call InitialQUAD_MP49(ID, PosGP, WeiGP)
            case (64) ! 64 material points per element    
                call InitialQUAD_MP64(ID, PosGP, WeiGP)   
            case (81) ! 81 material points per element
                call InitialQUAD_MP81(ID, PosGP, WeiGP)
                
                ! unique cases of # rows and 3 columns 
            !case (3) ! 1 row x 3 columns = 3 material points per element  
            !    call InitialQUAD_MP_1ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (6)! 2 row x 3 columns = 6 material points per element      
            !    call InitialQUAD_MP_2ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (9) ! 3 row x 3 columns = 9 material points per element     
            !    call InitialQUAD_MP_3ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (12)! 4 row x 3 columns = 12 material points per element     
            !    call InitialQUAD_MP_4ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (15)! 5 row x 3 columns = 15 material points per element     
            !    call InitialQUAD_MP_5ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (63) ! 6 row x 3 columns = 18 material points per element     
            !    call InitialQUAD_MP_6ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (73) ! 7 row x 3 columns = 21 material points per element     
            !    call InitialQUAD_MP_7ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (83) ! 8 row x 3 columns = 24 material points per element     
            !    call InitialQUAD_MP_8ROW_3COLUMNS(ID, PosGP, WeiGP)
                
              
        
        
        
        
        
        
        
        
        
        
        
        
        
                
                
                
                
                
                
                
                
            !case (81) ! 81 material points per element    
            !    call InitialQUAD_MP81(ID, PosGP, WeiGP)
            !case (100) ! 100 material points per element    
            !    call InitialQUAD_MP100(ID, PosGP, WeiGP)
            !case (121) ! 121 material points per element    
            !    call InitialQUAD_MP121(ID, PosGP, WeiGP)
                
            !case (41) ! Nurbs 4 in xi direction 1 in eta direction (need to generalize this) 
            !case (42) ! Nurbs 4 in xi direction 2 in eta direction (need to generalize this)
            !case (43) ! Nurbs 4 in xi direction 3 in eta direction (need to generalize this)
            !case (44) ! Nurbs 4 in xi direction 4 in eta direction (need to generalize this)
            !case (45) ! Nurbs 4 in xi direction 5 in eta direction (need to generalize this)
            !case (46) ! Nurbs 4 in xi direction 6 in eta direction (need to generalize this)
            !case (47) ! Nurbs 4 in xi direction 7 in eta direction (need to generalize this)
            !case (48) ! Nurbs 4 in xi direction 8 in eta direction (need to generalize this)
                
                case default
              call GiveError("Number of solid material points, " // trim(String(SolidPointsElement)) // &
                              ", is not available for quadrilateral elements! Supported numbers are 1, 4, 9, and 16. Error in [subroutine InitialLocalMaterialPointCoordinatesQUAD()].")
            end select
          end if
          
          if ( (IParticle > SolidPointsElement) .and. (LiquidPointsElement > 0) ) then
            ID = IParticle - SolidPointsElement
            select case(LiquidPointsElement)
            case (1) ! 1 material point per element
              call InitialQUAD_MP1(ID, PosGP, WeiGP)
            case (4) ! 4 material points per element
              call InitialQUAD_MP4(ID, PosGP, WeiGP)
            case default
              call GiveError("Number of solid material points, " // trim(String(LiquidPointsElement)) // &
                              ", is not available for quadrilateral elements! Supported numbers are 1 and 4. Error in [subroutine InitialLocalMaterialPointCoordinatesQUAD()].")
            end select
          end if

        end subroutine InitialLocalMaterialPointCoordinatesQUAD

        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialLocalMaterialPointCoordinatesHEXA
        !
        !    DESCRIPTION: 
        !>   Determines the local coordinates and integration weight assigned to material point with ID
        !>   IParticle which are returned through WeiGP and PosGP. Currently, all material points are placed at the same local positions.
        !
        !>   @param[in] IParticle : Number of the material point
        !>   @param[in] SolidPointsElement : Number of solid material points per element
        !>   @param[in] LiquidPointsElement : Number of liquid material points per element
        !>   @param[inout] WeiGP : Initial weight assigned to material point IParticle
        !>   @param[inout] PosGP : Initial local position of material point IParticle
        !
        !**********************************************************************
        subroutine InitialLocalMaterialPointCoordinatesHEXA(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP

          ! local variables
          integer(INTEGER_TYPE) :: ID

          ! first the solid material points are determined
          if ( (IParticle <= SolidPointsElement) .and. (SolidPointsElement > 0) ) then
            ID = IParticle
            select case(SolidPointsElement)
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
            !case (4) ! 4 material points per element
            !  call InitialQUAD_MP4(ID, PosGP, WeiGP)
            !case (9) ! 9 material points per element 
            !  call InitialQUAD_MP9(ID, PosGP, WeiGP) 
            !case (16) ! 16 material points per element 
            !  call InitialQUAD_MP16(ID, PosGP, WeiGP)  
            !case (25) ! 25 material points per element
            !    call InitialQUAD_MP25(ID, PosGP, WeiGP)
            !case (36) ! 36 material points per element    
            !    call InitialQUAD_MP36(ID, PosGP, WeiGP)
            !case (49) ! 49 material points per element    
            !    call InitialQUAD_MP49(ID, PosGP, WeiGP)
            !case (64) ! 64 material points per element    
            !    call InitialQUAD_MP64(ID, PosGP, WeiGP)   
            !case (81) ! 81 material points per element
            !    call InitialQUAD_MP81(ID, PosGP, WeiGP)
                
                ! unique cases of # rows and 3 columns 
            !case (3) ! 1 row x 3 columns = 3 material points per element  
            !    call InitialQUAD_MP_1ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (6)! 2 row x 3 columns = 6 material points per element      
            !    call InitialQUAD_MP_2ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (9) ! 3 row x 3 columns = 9 material points per element     
            !    call InitialQUAD_MP_3ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (12)! 4 row x 3 columns = 12 material points per element     
            !    call InitialQUAD_MP_4ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (15)! 5 row x 3 columns = 15 material points per element     
            !    call InitialQUAD_MP_5ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (63) ! 6 row x 3 columns = 18 material points per element     
            !    call InitialQUAD_MP_6ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (73) ! 7 row x 3 columns = 21 material points per element     
            !    call InitialQUAD_MP_7ROW_3COLUMNS(ID, PosGP, WeiGP)
            !case (83) ! 8 row x 3 columns = 24 material points per element     
            !    call InitialQUAD_MP_8ROW_3COLUMNS(ID, PosGP, WeiGP)
                
              
        
        
        
        
        
        
        
        
        
        
        
        
        
                
                
                
                
                
                
                
                
            !case (81) ! 81 material points per element    
            !    call InitialQUAD_MP81(ID, PosGP, WeiGP)
            !case (100) ! 100 material points per element    
            !    call InitialQUAD_MP100(ID, PosGP, WeiGP)
            !case (121) ! 121 material points per element    
            !    call InitialQUAD_MP121(ID, PosGP, WeiGP)
                
            !case (41) ! Nurbs 4 in xi direction 1 in eta direction (need to generalize this) 
            !case (42) ! Nurbs 4 in xi direction 2 in eta direction (need to generalize this)
            !case (43) ! Nurbs 4 in xi direction 3 in eta direction (need to generalize this)
            !case (44) ! Nurbs 4 in xi direction 4 in eta direction (need to generalize this)
            !case (45) ! Nurbs 4 in xi direction 5 in eta direction (need to generalize this)
            !case (46) ! Nurbs 4 in xi direction 6 in eta direction (need to generalize this)
            !case (47) ! Nurbs 4 in xi direction 7 in eta direction (need to generalize this)
            !case (48) ! Nurbs 4 in xi direction 8 in eta direction (need to generalize this)
                
                case default
              call GiveError("Number of solid material points, " // trim(String(SolidPointsElement)) // &
                              ", is not available for Hexahedral elements! Supported numbers are 1, 8, 27, and 64. Error in [subroutine InitialLocalMaterialPointCoordinatesHEXA()].")
            end select
          end if
          
          if ( (IParticle > SolidPointsElement) .and. (LiquidPointsElement > 0) ) then
            ID = IParticle - SolidPointsElement
            select case(LiquidPointsElement)
            !case (1) ! 1 material point per element
            !  call InitialQUAD_MP1(ID, PosGP, WeiGP)
            !case (4) ! 4 material points per element
            !  call InitialQUAD_MP4(ID, PosGP, WeiGP)
            case default
              call GiveError("Number of solid material points, " // trim(String(LiquidPointsElement)) // &
                              ", is not available for quadrilateral elements! Supported numbers are 1 and 4. Error in [subroutine InitialLocalMaterialPointCoordinatesQUAD()].")
            end select
          end if

        end subroutine InitialLocalMaterialPointCoordinatesHEXA
        
          
                
                
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP_1ROW_3COLUMNS
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 row x 3 columns = 3 material points per element are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP_1ROW_3COLUMNS(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE) :: nMP
          
          
          
          ! local variables
          real(REAL_TYPE) :: a
          real(REAL_TYPE) :: b
          real(REAL_TYPE) :: row_1
  
          real(REAL_TYPE) :: weight
          nMP = 1*3
          a = 1/nMP
          b = 0.0
          row_1 = 0.0
          weight = 4/nMP
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = row_1
              WeiGP = weight
            case (2)
              PosGP(1) = b
              PosGP(2) = row_1
              WeiGP = weight
            case (3)
              PosGP(1) = +a
              PosGP(2) = row_1
              WeiGP = weight
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP_1ROW_3COLUMNS
                
                
                
                
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP_2ROW_3COLUMNS
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 row x 3 columns = 3 material points per element are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP_2ROW_3COLUMNS(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE) :: nMP_X
          real(REAL_TYPE) :: nMP_Y
          !real(REAL_TYPE) :: row_1
          !real(REAL_TYPE) :: row_2
          
          
          
          ! local variables
          real(REAL_TYPE) :: a
          real(REAL_TYPE) :: b
          real(REAL_TYPE) :: row_2
          real(REAL_TYPE) :: row_1
  
          real(REAL_TYPE) :: weight
          
          nMP_X = 3
          nMP_Y = 2
          a = 1/nMP_X
          b = 0.0
          row_2 = 1/3
          row_1 = -1/3
  
          weight = 4/(nMP_X*nMP_Y)
          
          
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = row_1
              WeiGP = weight
            case (2)
              PosGP(1) = b
              PosGP(2) = row_1
              WeiGP = weight
            case (3)
              PosGP(1) = +a
              PosGP(2) = row_1
              WeiGP = weight
              
              
            case (4)
              PosGP(1) = -a
              PosGP(2) = row_2
              WeiGP = weight
            case (5)
              PosGP(1) = b
              PosGP(2) = row_2
              WeiGP = weight
            case (6)
              PosGP(1) = a
              PosGP(2) = row_2
              WeiGP = weight
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP_2ROW_3COLUMNS
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP_1ROW_3COLUMNS
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 row x 3 columns = 3 material points per element are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP_3ROW_3COLUMNS(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE) :: nMP_X
          real(REAL_TYPE) :: nMP_Y
          
          
          
          ! local variables
          real(REAL_TYPE) :: a
          real(REAL_TYPE) :: b
          real(REAL_TYPE) :: row_3
          real(REAL_TYPE) :: row_2
          real(REAL_TYPE) :: row_1
          
          real(REAL_TYPE) :: weight
          nMP_X = 3
          nMP_Y = 3
          
          a = 1/nMP_X
          b = 0.0
          row_3 = 0.5
          row_2 = 0.0
          row_1 = -0.5
          
          weight = 4/(nMP_X*nMP_Y)
          
          
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = row_1
              WeiGP = weight
            case (2)
              PosGP(1) = b
              PosGP(2) = row_1
              WeiGP = weight
            case (3)
              PosGP(1) = +a
              PosGP(2) = row_1
              WeiGP = weight
              
              
            case (4)
              PosGP(1) = -a
              PosGP(2) = row_2
              WeiGP = weight
            case (5)
              PosGP(1) = b
              PosGP(2) = row_2
              WeiGP = weight
            case (6)
              PosGP(1) = +a
              PosGP(2) = row_2
              WeiGP = weight
              
              
            case (7)
              PosGP(1) = -a
              PosGP(2) = row_3
              WeiGP = weight
            case (8)
              PosGP(1) = b
              PosGP(2) = row_3
              WeiGP = weight
            case (9)
              PosGP(1) = +a
              PosGP(2) = row_3
              WeiGP = weight
              
              
            
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP_3ROW_3COLUMNS
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP_1ROW_3COLUMNS
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 row x 3 columns = 3 material points per element are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP_4ROW_3COLUMNS(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE) :: nMP_X
          real(REAL_TYPE) :: nMP_Y
          real(REAL_TYPE) :: a,b,row_1,row_2,row_3,row_4
          
          
          
          ! local variables
          !real(REAL_TYPE) :: a 
          !real(REAL_TYPE) :: b 
          !real(REAL_TYPE) :: row_4 
          !real(REAL_TYPE) :: row_3 
          !real(REAL_TYPE) :: row_2 
          !real(REAL_TYPE) :: row_1 
          
          real(REAL_TYPE) :: weight
          
          a = 1/nMP_X
          b = 0.0
          row_4 = 0.6
          row_3 = 0.2
          row_2 = -0.2
          row_1 = -0.6
          weight = 4/(nMP_X*nMP_Y)
          
          nMP_X = 3
          nMP_Y = 4
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = row_1
              WeiGP = weight
            case (2)
              PosGP(1) = b
              PosGP(2) = row_1
              WeiGP = weight
            case (3)
              PosGP(1) = +a
              PosGP(2) = row_1
              WeiGP = weight
              
              
            case (4)
              PosGP(1) = -a
              PosGP(2) = row_2
              WeiGP = weight
            case (5)
              PosGP(1) = b
              PosGP(2) = row_2
              WeiGP = weight
            case (6)
              PosGP(1) = +a
              PosGP(2) = row_2
              WeiGP = weight
              
              
            case (7)
              PosGP(1) = -a
              PosGP(2) = row_3
              WeiGP = weight
            case (8)
              PosGP(1) = b
              PosGP(2) = row_3
              WeiGP = weight
            case (9)
              PosGP(1) = +a
              PosGP(2) = row_3
              WeiGP = weight
              
              
              
              case (10)
              PosGP(1) = -a
              PosGP(2) = row_4
              WeiGP = weight
            case (11)
              PosGP(1) = b
              PosGP(2) = row_4
              WeiGP = weight
            case (12)
              PosGP(1) = +a
              PosGP(2) = row_4
              WeiGP = weight
            
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP_4ROW_3COLUMNS
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP_1ROW_3COLUMNS
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 row x 3 columns = 3 material points per element are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP_5ROW_3COLUMNS(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE) :: nMP_X
          real(REAL_TYPE) :: nMP_Y
          real(REAL_TYPE) :: a,b,row_1,row_2,row_3,row_4,row_5
          
          
          
          ! local variables
          !real(REAL_TYPE) :: a
          !real(REAL_TYPE) :: b
          !real(REAL_TYPE) :: row_5
          !real(REAL_TYPE) :: row_4
          !real(REAL_TYPE) :: row_3
          !real(REAL_TYPE) :: row_2
          !real(REAL_TYPE) :: row_1
          
          real(REAL_TYPE) :: weight
          
          
          nMP_X = 3
          nMP_Y = 5
          
          a = 1/nMP_X
          b = 0.0
          row_5 = 0.66666667
          row_4 = 0.33333333
          row_3 = 0.0
          row_2 = -0.33333333
          row_1 = -0.66666667
          
          weight = 4/(nMP_X*nMP_Y)
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = row_1
              WeiGP = weight
            case (2)
              PosGP(1) = b
              PosGP(2) = row_1
              WeiGP = weight
            case (3)
              PosGP(1) = +a
              PosGP(2) = row_1
              WeiGP = weight
              
              
            case (4)
              PosGP(1) = -a
              PosGP(2) = row_2
              WeiGP = weight
            case (5)
              PosGP(1) = b
              PosGP(2) = row_2
              WeiGP = weight
            case (6)
              PosGP(1) = +a
              PosGP(2) = row_2
              WeiGP = weight
              
              
            case (7)
              PosGP(1) = -a
              PosGP(2) = row_3
              WeiGP = weight
            case (8)
              PosGP(1) = b
              PosGP(2) = row_3
              WeiGP = weight
            case (9)
              PosGP(1) = +a
              PosGP(2) = row_3
              WeiGP = weight
              
              
              
              case (10)
              PosGP(1) = -a
              PosGP(2) = row_4
              WeiGP = weight
            case (11)
              PosGP(1) = b
              PosGP(2) = row_4
              WeiGP = weight
            case (12)
              PosGP(1) = +a
              PosGP(2) = row_4
              WeiGP = weight
            
              
              case (13)
              PosGP(1) = -a
              PosGP(2) = row_5
              WeiGP = weight
            case (14)
              PosGP(1) = b
              PosGP(2) = row_5
              WeiGP = weight
            case (15)
              PosGP(1) = +a
              PosGP(2) = row_5
              WeiGP = weight
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP_5ROW_3COLUMNS
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !3D
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP1
        !
        !    DESCRIPTION:  
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 material point is placed in each element, whose initial
        !>   location and weight is identical with those of the Gauss Point.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP1(IParticle, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.0
              PosGP(2) = 0.0 
              PosGP(3) = 0.0 !position right in the center of the local element
              WeiGP = 8.0 !why is the gauss point weight 4? This should be two
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP1()].")   
          end select
        
        end subroutine InitialHEXA_MP1
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP4
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   4 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP8(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.57735026918962576450914878050196 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 1.0
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !2D -> observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          !3D -> observe how the sum of all WeiGP is 8 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b 
            case (2)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b
            case (3)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b
            case (5)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = b 
            case (6)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = b
            case (7)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b
            case (8)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b  
            case default
              call GiveError("Undefined number of material points in [subroutine InitialHEXA_MP8()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialHEXA_MP8
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP9
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   9 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP27(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.77459666924 ! = sqrt(3/5)
          !real(REAL_TYPE) :: b = 0.0
          real(REAL_TYPE) :: c = 0.88888888888
          real(REAL_TYPE) :: d = 0.55555555555


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a          
              PosGP(3) = -a
              
              WeiGP = d*d *d
            case (2)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = -a
              
              WeiGP = d*c *d
            case (3)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              
              WeiGP = d*d *d
            case (4)
              PosGP(1) = -a
              PosGP(2) = 0.0
              PosGP(3) = -a

              WeiGP = c*d *d
            case (5)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              PosGP(3) = -a
              
              WeiGP = c*c *d
            case (6)
              PosGP(1) = a
              PosGP(2) = 0.0
              PosGP(3) = -a

              WeiGP = c*d *d
            case (7)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = d*d *d
            case (8)
              PosGP(1) = 0.0
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = d*c *d
            case (9)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = d*d *d
              
              
              
              
            case (10)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = 0.0

              WeiGP = d*d *c
            case (11)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = 0.0

              WeiGP = d*c *c
            case (12)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = 0.0

              WeiGP = d*d *c
            case (13)
              PosGP(1) = -a
              PosGP(2) = 0.0
              PosGP(3) = 0.0

              WeiGP = c*d *c
            case (14)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              PosGP(3) = 0.0

              WeiGP = c*c *c
            case (15)
              PosGP(1) = a
              PosGP(2) = 0.0
              PosGP(3) = 0.0

              WeiGP = c*d *c
            case (16)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = d*d *c
            case (17)
              PosGP(1) = 0.0
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = d*c *c
            case (18)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = d*d *c
              
              
              
              
              
              
              
              case (19)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = a

              WeiGP = d*d *d
            case (20)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = a

              WeiGP = d*c *d
            case (21)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = a

              WeiGP = d*d *d
            case (22)
              PosGP(1) = -a
              PosGP(2) = 0.0
              PosGP(3) = a 

              WeiGP = c*d *d
            case (23)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              PosGP(3) = a

              WeiGP = c*c *d
            case (24)
              PosGP(1) = a
              PosGP(2) = 0.0
              PosGP(3) = a

              WeiGP = c*d *d
            case (25)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = a 

              WeiGP = d*d *d
            case (26)
              PosGP(1) = 0.0
              PosGP(2) = a
              PosGP(3) = a 

              WeiGP = d*c *d
            case (27)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a 

              WeiGP = d*d *d
              
              
              
              
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP4()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialHEXA_MP27
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP16
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   16 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP64(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.86113631159
          real(REAL_TYPE) :: b = 0.33998104358 
          !real(REAL_TYPE) :: c = +0.33998104358
          !real(REAL_TYPE) :: d = +0.86113631159
          
          
          real(REAL_TYPE) :: e = 0.34785484513 ! --> a
          real(REAL_TYPE) :: f = 0.65214515486 ! --> b
          


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
              
            ! 3--> -a 
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -a
              
              WeiGP = e*e *e
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -a

              WeiGP = e*f *e
            case (3)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -a

              WeiGP = e*f *e
            case (4)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a

              WeiGP = e*e *e
            case (5)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = e*f *e
            case (6)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = f*f *e
            case (7)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = f*f *e
            case (8)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = e*f *e
            case (9)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = e*f *e              
            case (10)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = f*f *e
            case (11)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = f*f *e
            case (12)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = e*f *e
            case (13)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = e*e *e
            case (14)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = e*f *e
            case (15)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = e*f *e
            case (16) 
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = e*e *e
              
              
              
              
              
              ! 3--> -b 
              case (17)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -b
              
              WeiGP = e*e *f
            case (18)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -b

              WeiGP = e*f *f
            case (19)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -b

              WeiGP = e*f *f
            case (20)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -b

              WeiGP = e*e *f
            case (21)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = e*f *f
            case (22)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = f*f *f
            case (23)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = f*f *f
            case (24)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = e*f *f
            case (25)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = e*f *f              
            case (26)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = f*f *f
            case (27)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = f*f *f
            case (28)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = e*f *f
            case (29)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = e*e *f
            case (30)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = e*f *f
            case (31)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = e*f *f
            case (32) 
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = e*e *f
              
              
              
              ! 3--> b
              case (33)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = b
              
              WeiGP = e*e *f
            case (34)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = b

              WeiGP = e*f *f
            case (35)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = b

              WeiGP = e*f *f
            case (36)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = b

              WeiGP = e*e *f
            case (37)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = e*f *f
            case (38)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = f*f *f
            case (39)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = f*f *f
            case (40)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = e*f *f
            case (41)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = e*f *f              
            case (42)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = f*f *f
            case (43)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = f*f *f
            case (44)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = e*f *f
            case (45)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = e*e *f
            case (46)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = e*f *f
            case (47)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = e*f *f
            case (48) 
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = e*e *f
              
              
              
              
              
              
              ! 3--> a
            case (49)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = a
              
              WeiGP = e*e *e
            case (50)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = a

              WeiGP = e*f *e
            case (51)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = a

              WeiGP = e*f *e
            case (52)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = a

              WeiGP = e*e *e
            case (53)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = e*f *e
            case (54)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = f*f *e
            case (55)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = f*f *e
            case (56)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = e*f *e
            case (57)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = e*f *e              
            case (58)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = f*f *e
            case (59)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = f*f *e
            case (60)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = e*f *e
            case (61)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = e*e *e
            case (62)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = e*f *e
            case (63)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = e*f *e
            case (64) 
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = e*e *e
              
              
              
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP16()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialHEXA_MP64
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP25
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   25 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP125(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.90617984593
          real(REAL_TYPE) :: b = 0.5384693101
          
          real(REAL_TYPE) :: c = 0.23692688505 ! --> a
          real(REAL_TYPE) :: d = 0.56888888888 ! --> 0
          real(REAL_TYPE) :: e = 0.47862867049 ! --> b
          


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          
          !!!!!!!!!!!!!!!!!!!!!! --> -a
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c*c*c
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = e*c*c
            case (3)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = d*c*c
            case (4)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c*e*c
            case (5)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c*c*c
            case (6)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c*e*c
            case (7)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = e*e*c
            case (8)
              PosGP(1) = 0
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = d*e*c
            case (9)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = e*e*c          
            case (10)
              PosGP(1) = a 
              PosGP(2) = -b
              PosGP(3) = -a

              WeiGP = c*e*c 
            case (11)
              PosGP(1) = -a
              PosGP(2) = 0
              PosGP(3) = -a

              WeiGP = c*d*c
            case (12)
              PosGP(1) = -b
              PosGP(2) = 0
              PosGP(3) = -a

              WeiGP = d*e*c
            case (13)
              PosGP(1) = 0
              PosGP(2) = 0
              PosGP(3) = -a

              WeiGP = d*d*c
            case (14)
              PosGP(1) = b
              PosGP(2) = 0
              PosGP(3) = -a

              WeiGP = e*d*c
            case (15)
              PosGP(1) = a
              PosGP(2) = 0
              PosGP(3) = -a

              WeiGP = c*d*c
            case (16)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = c*e*c
              
            case (17)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = e*e*c
            case (18)
              PosGP(1) = 0.0
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = d*e*c
            case (19)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = e*e*c
            case (20)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -a

              WeiGP = c*e*c
            case (21)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = c*c*c
            case (22)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = e*c*c
            case (23)
              PosGP(1) = 0
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = d*c*c
            case (24)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = c*e*c
            case (25)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a

              WeiGP = c*c*c              
            
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --> -b
              
              case (26)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = c*c*e
            case (27)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = e*c*e
            case (28)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = d*c*e
            case (29)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = e*c*e
            case (30)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = c*c*e
            case (31)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = c*e*e
            case (32)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = e*e*e
            case (33)
              PosGP(1) = 0
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = d*e*e
            case (34)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = e*e*e          
            case (35)
              PosGP(1) = a 
              PosGP(2) = -b
              PosGP(3) = -b

              WeiGP = c*e*e 
            case (36)
              PosGP(1) = -a
              PosGP(2) = 0
              PosGP(3) = -b

              WeiGP = c*d*e
            case (37)
              PosGP(1) = -b
              PosGP(2) = 0
              PosGP(3) = -b

              WeiGP = d*e*e
            case (38)
              PosGP(1) = 0
              PosGP(2) = 0
              PosGP(3) = -b

              WeiGP = d*d*e
            case (39)
              PosGP(1) = b
              PosGP(2) = 0
              PosGP(3) = -b

              WeiGP = e*d*e
            case (40)
              PosGP(1) = a
              PosGP(2) = 0
              PosGP(3) = -b

              WeiGP = c*d*e
            case (41)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = c*e*e
              
            case (42)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = e*e*e
            case (43)
              PosGP(1) = 0.0
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = d*e*e
            case (44)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = e*e*e
            case (45)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -b

              WeiGP = c*e*e
            case (46)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = c*c*e
            case (47)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = e*c*e
            case (48)
              PosGP(1) = 0
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = d*c*e
            case (49)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = e*c*e
            case (50)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -b

              WeiGP = c*c*e              
            
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --> 0.0
              
              
              case (51)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = 0.0
              WeiGP = c*c*d
            case (52)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = 0.0
              WeiGP = e*c*d
            case (53)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = 0.0
              WeiGP = d*c*d
            case (54)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = 0.0
              WeiGP = c*e*d
            case (55)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = 0.0
              WeiGP = c*c*d
            case (56)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = 0.0
              WeiGP = c*e*d
            case (57)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = 0.0
              WeiGP = e*e*d
            case (58)
              PosGP(1) = 0
              PosGP(2) = -b
              PosGP(3) = 0.0

              WeiGP = d*e*d
            case (59)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = 0.0

              WeiGP = e*e*d          
            case (60)
              PosGP(1) = a 
              PosGP(2) = -b
              PosGP(3) = 0.0

              WeiGP = c*e*d 
            case (61)
              PosGP(1) = -a
              PosGP(2) = 0
              PosGP(3) = 0.0

              WeiGP = c*d*d
            case (62)
              PosGP(1) = -b
              PosGP(2) = 0
              PosGP(3) = 0.0

              WeiGP = d*e*d
            case (63)
              PosGP(1) = 0
              PosGP(2) = 0
              PosGP(3) = 0.0

              WeiGP = d*d*d
            case (64)
              PosGP(1) = b
              PosGP(2) = 0
              PosGP(3) = 0.0

              WeiGP = e*d*d
            case (65)
              PosGP(1) = a
              PosGP(2) = 0
              PosGP(3) = 0.0

              WeiGP = c*d*d
            case (66)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = 0.0

              WeiGP = c*e*d
              
            case (67)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = 0.0

              WeiGP = e*e*d
            case (68)
              PosGP(1) = 0.0
              PosGP(2) = b
              PosGP(3) = 0.0

              WeiGP = d*e*d
            case (69)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = 0.0

              WeiGP = e*e*d
            case (70)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = 0.0

              WeiGP = c*e*d
            case (71)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = c*c*d
            case (72)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = e*c*d
            case (73)
              PosGP(1) = 0
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = d*c*d
            case (74)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = c*e*d
            case (75)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = 0.0

              WeiGP = c*c*d              
            
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --> b
              
              case (76)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c*c*e
            case (77)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = e*c*e
            case (78)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = d*c*e
            case (79)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c*e*e
            case (80)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c*c*e
            case (81)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = c*e*e
            case (82)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = e*e*e
            case (83)
              PosGP(1) = 0
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = d*e*e
            case (84)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = e*e*e          
            case (85)
              PosGP(1) = a 
              PosGP(2) = -b
              PosGP(3) = b

              WeiGP = c*e*e 
            case (86)
              PosGP(1) = -a
              PosGP(2) = 0
              PosGP(3) = b

              WeiGP = c*d*e
            case (87)
              PosGP(1) = -b
              PosGP(2) = 0
              PosGP(3) = b

              WeiGP = d*e*e
            case (88)
              PosGP(1) = 0
              PosGP(2) = 0
              PosGP(3) = b

              WeiGP = d*d*e
            case (89)
              PosGP(1) = b
              PosGP(2) = 0
              PosGP(3) = b

              WeiGP = e*d*e
            case (90)
              PosGP(1) = a
              PosGP(2) = 0
              PosGP(3) = b

              WeiGP = c*d*e
            case (91)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = c*e*e
              
            case (92)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = e*e*e
            case (93)
              PosGP(1) = 0.0
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = d*e*e
            case (94)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = e*e*e
            case (95)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = b

              WeiGP = c*e*e
            case (96)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = c*c*e
            case (97)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = e*c*e
            case (98)
              PosGP(1) = 0
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = d*c*e
            case (99)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = c*e*e
            case (100)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = b

              WeiGP = c*c*e              
            
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --> a
              
              case (101)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = c*c*c
            case (102)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = e*c*c
            case (103)
              PosGP(1) = 0.0
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = d*c*c
            case (104)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = c*e*c
            case (105)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = c*c*c
            case (106)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = c*e*c
            case (107)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = e*e*c
            case (108)
              PosGP(1) = 0
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = d*e*c
            case (109)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = e*e*c          
            case (110)
              PosGP(1) = a 
              PosGP(2) = -b
              PosGP(3) = a

              WeiGP = c*e*c 
            case (111)
              PosGP(1) = -a
              PosGP(2) = 0
              PosGP(3) = a

              WeiGP = c*d*c
            case (112)
              PosGP(1) = -b
              PosGP(2) = 0
              PosGP(3) = a

              WeiGP = d*e*c
            case (113)
              PosGP(1) = 0
              PosGP(2) = 0
              PosGP(3) = a

              WeiGP = d*d*c
            case (114)
              PosGP(1) = b
              PosGP(2) = 0
              PosGP(3) = a

              WeiGP = e*d*c
            case (115)
              PosGP(1) = a
              PosGP(2) = 0
              PosGP(3) = a

              WeiGP = c*d*c
            case (116)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = c*e*c
              
            case (117)
              PosGP(1) = -b 
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = e*e*c
            case (118)
              PosGP(1) = 0.0
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = d*e*c
            case (119)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = e*e*c
            case (120)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = a

              WeiGP = c*e*c
            case (121)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = c*c*c
            case (122)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = e*c*c
            case (123)
              PosGP(1) = 0
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = d*c*c
            case (124)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = c*e*c
            case (125)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a

              WeiGP = c*c*c              
            
              
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
              
              
              
              
              
              
              
              
              
              
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialHEXA_MP125
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP36
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   36 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP216(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.9324695142031521
          real(REAL_TYPE) :: b = 0.6612093864662645
          real(REAL_TYPE) :: c = 0.2386191860831969
  
          real(REAL_TYPE) :: a_weight = 0.1713244923791704 
          real(REAL_TYPE) :: b_weight = 0.3607615730481386 
          real(REAL_TYPE) :: c_weight = 0.4679139345726910 
          


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight
            case (4)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight
            case (5)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight
            case (6)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
            case (7)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight
            case (8)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight
            case (9)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight    
            case (10)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight
            case (11)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight
            case (12)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight
            case (13)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (14)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (15)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (16)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (17)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (18)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (19)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (20)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (21)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (22)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (23)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (24)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (25)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight    
              case (26)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight  
              case (27)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight  
              case (28)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight  
              case (29)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight 
              case (30)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight  
              case (31)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight  
              case (32)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight  
              case (33)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight  
              case (34)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight  
              case (35)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight  
              case (36)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight  
              
              
              
              
              
              
              
              case (37)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight
            case (38)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight
            case (39)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight
            case (40)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight
            case (41)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight
            case (42)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight
            case (43)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight
            case (44)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight
            case (45)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight    
            case (46)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight
            case (47)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight
            case (48)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight
            case (49)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (50)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (51)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (52)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (53)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (54)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (55)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (56)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (57)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (58)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (59)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (60)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (61)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight    
              case (62)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight  
              case (63)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight  
              case (64)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight  
              case (65)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight 
              case (66)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight  
              case (67)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight  
              case (68)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight  
              case (69)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight  
              case (70)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight  
              case (71)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight  
              case (72)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight
              
              
              
              
              
              
              
              
              case (73)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight
            case (74)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight
            case (75)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight
            case (76)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight
            case (77)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight
            case (78)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight
            case (79)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight
            case (80)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight
            case (81)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight    
            case (82)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight
            case (83)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight
            case (84)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight
            case (85)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (86)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (87)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (88)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (89)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (90)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (91)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (92)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (93)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (94)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (95)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (96)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (97)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight    
              case (98)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*C_weight  
              case (99)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight  
              case (100)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight  
              case (101)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight 
              case (102)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight  
              case (103)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight  
              case (104)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight  
              case (105)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight  
              case (106)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight  
              case (107)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight  
              case (108)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight
              
              
              
              case (109)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight
            case (110)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight
            case (111)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight
            case (112)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight
            case (113)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight
            case (114)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight
            case (115)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight
            case (116)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight
            case (117)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight    
            case (118)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight
            case (119)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight
            case (120)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight
            case (121)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (122)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (123)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (124)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (125)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (126)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (127)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (128)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (129)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (130)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (131)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (132)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (133)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight    
              case (134)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight  
              case (135)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight  
              case (136)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight  
              case (137)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight 
              case (138)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight  
              case (139)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight  
              case (140)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight  
              case (141)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight  
              case (142)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight  
              case (143)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight  
              case (144)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight
              
              
              
              
              
              
              case (145)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight
            case (146)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight
            case (147)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight
            case (148)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight
            case (149)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight
            case (150)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight
            case (151)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight
            case (152)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight
            case (153)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight    
            case (154)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight
            case (155)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight
            case (156)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight
            case (157)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (158)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (159)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (160)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (161)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (162)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (163)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (164)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (165)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (166)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (167)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (168)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (169)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight    
              case (170)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight  
              case (171)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight  
              case (172)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight  
              case (173)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight 
              case (174)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight  
              case (175)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight  
              case (176)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight  
              case (177)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight  
              case (178)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight  
              case (179)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight  
              case (180)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight
              
              
              
              
              
              
              case (181)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = a_weight*a_weight*a_weight
            case (182)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = b_weight*a_weight*a_weight
            case (183)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = c_weight*a_weight*a_weight
            case (184)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = c_weight*a_weight*a_weight
            case (185)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = b_weight*a_weight*a_weight
            case (186)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = a_weight*a_weight*a_weight
            case (187)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = a_weight*b_weight*a_weight
            case (188)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = b_weight*b_weight*a_weight
            case (189)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = c_weight*b_weight*a_weight    
            case (190)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = c_weight*b_weight*a_weight
            case (191)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = b_weight*b_weight*a_weight
            case (192)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = a
              WeiGP = a_weight*b_weight*a_weight
            case (193)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = a
              WeiGP = a_weight*c_weight*a_weight
            case (194)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = a
              WeiGP = b_weight*c_weight*a_weight
            case (195)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = a
              WeiGP = c_weight*c_weight*a_weight
            case (196)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = a
              WeiGP = c_weight*c_weight*a_weight
            case (197)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = a
              WeiGP = b_weight*c_weight*a_weight
            case (198)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = a
              WeiGP = a_weight*c_weight*a_weight
            case (199)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = a
              WeiGP = a_weight*c_weight*a_weight
            case (200)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = a
              WeiGP = b_weight*c_weight*a_weight
            case (201)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = a
              WeiGP = c_weight*c_weight*a_weight
            case (202)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = a
              WeiGP = c_weight*c_weight*a_weight
            case (203)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = a
              WeiGP = b_weight*c_weight*a_weight
            case (204)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = a
              WeiGP = a_weight*c_weight*a_weight
            case (205)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = a
              WeiGP = a_weight*b_weight*a_weight    
              case (206)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = a
              WeiGP = b_weight*b_weight*a_weight  
              case (207)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = a
              WeiGP = c_weight*b_weight*a_weight  
              case (208)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = a
              WeiGP = c_weight*b_weight*a_weight  
              case (209)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = a
              WeiGP = b_weight*b_weight*a_weight 
              case (210)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = a
              WeiGP = a_weight*b_weight*a_weight  
              case (211)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = a_weight*a_weight*a_weight  
              case (212)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b_weight*a_weight*a_weight  
              case (213)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = c_weight*a_weight*a_weight  
              case (214)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = c_weight*a_weight*a_weight  
              case (215)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b_weight*a_weight*a_weight  
              case (216)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = a_weight*a_weight*a_weight
              
              
              
              
              
              
              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialHEXA_MP216
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP49
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   49 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialHEXA_MP343(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.9491079123427585
          real(REAL_TYPE) :: b = 0.7415311855993945
          real(REAL_TYPE) :: c = 0.4058451513773972
          real(REAL_TYPE) :: d = 0.0
  
          real(REAL_TYPE) :: a_weight = 0.1294849661688697
          real(REAL_TYPE) :: b_weight = 0.2797053914892766
          real(REAL_TYPE) :: c_weight = 0.3818300505051189
          real(REAL_TYPE) :: d_weight = 0.4179591836734694 


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight
            case (4)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = d_weight*a_weight*a_weight
            case (5)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight
            case (6)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight
            case (7)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
            case (8)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight
            case (9)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight
            case (10)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight
            case (11)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = d_weight*b_weight*a_weight      
            case (12)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight
            case (13)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight
            case (14)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight
            case (15)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (16)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (17)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (18)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = d_weight*c_weight*a_weight
            case (19)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (20)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (21)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (22)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = a_weight*d_weight*a_weight
            case (23)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = b_weight*d_weight*a_weight
            case (24)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = c_weight*d_weight*a_weight
            case (25)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = d_weight*d_weight*a_weight
            case (26)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = c_weight*d_weight*a_weight
            case (27)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = b_weight*d_weight*a_weight
            case (28)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = a_weight*d_weight*a_weight 
            case (29)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (30)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (31)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (32)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = d_weight*c_weight*a_weight
            case (33)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (34)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (35)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (36)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight    
              case (37)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight  
              case (38)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight 
              case (39)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = d_weight*b_weight*a_weight  
              case (40)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight  
              case (41)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight  
              case (42)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight  
              case (43)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight  
              case (44)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight  
              case (45)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight  
              case (46)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = d_weight*a_weight*a_weight  
              case (47)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight  
              case (48)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight  
              case (49)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight  
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              case (50)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight
            case (51)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight
            case (52)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight
            case (53)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = d_weight*a_weight*b_weight
            case (54)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight
            case (55)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight
            case (56)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight
            case (57)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight
            case (58)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight
            case (59)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight
            case (60)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = d_weight*b_weight*b_weight      
            case (61)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight
            case (62)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight
            case (63)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight
            case (64)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (65)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (66)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (67)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = d_weight*c_weight*b_weight
            case (68)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (69)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (70)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (71)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = a_weight*d_weight*b_weight
            case (72)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = b_weight*d_weight*b_weight
            case (73)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = c_weight*d_weight*b_weight
            case (74)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = d_weight*d_weight*b_weight
            case (75)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = c_weight*d_weight*b_weight
            case (76)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = b_weight*d_weight*b_weight
            case (77)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = -b
              WeiGP = a_weight*d_weight*b_weight 
            case (78)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (79)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (80)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (81)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = d_weight*c_weight*b_weight
            case (82)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = c_weight*c_weight*b_weight
            case (83)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = b_weight*c_weight*b_weight
            case (84)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -b
              WeiGP = a_weight*c_weight*b_weight
            case (85)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight    
              case (86)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight  
              case (87)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight 
              case (88)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = d_weight*b_weight*b_weight  
              case (89)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = c_weight*b_weight*b_weight  
              case (90)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = b_weight*b_weight*b_weight  
              case (91)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -b
              WeiGP = a_weight*b_weight*b_weight  
              case (92)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight  
              case (93)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight  
              case (94)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight  
              case (95)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = d_weight*a_weight*b_weight  
              case (96)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = c_weight*a_weight*b_weight  
              case (97)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = b_weight*a_weight*b_weight  
              case (98)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -b
              WeiGP = a_weight*a_weight*b_weight
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              case (99)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight
            case (100)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight
            case (101)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight
            case (102)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = d_weight*a_weight*c_weight
            case (103)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight
            case (104)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight
            case (105)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight
            case (106)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight
            case (107)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight
            case (108)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight
            case (109)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = d_weight*b_weight*c_weight      
            case (110)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight
            case (111)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight
            case (112)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight
            case (113)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (114)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (115)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (116)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = d_weight*c_weight*c_weight
            case (117)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (118)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (119)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (120)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = a_weight*d_weight*c_weight
            case (121)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = b_weight*d_weight*c_weight
            case (122)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = c_weight*d_weight*c_weight
            case (123)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = d_weight*d_weight*c_weight
            case (124)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = c_weight*d_weight*c_weight
            case (125)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = b_weight*d_weight*c_weight
            case (126)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = -c
              WeiGP = a_weight*d_weight*c_weight 
            case (127)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (128)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (129)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (130)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = d_weight*c_weight*c_weight
            case (131)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = c_weight*c_weight*c_weight
            case (132)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = b_weight*c_weight*c_weight
            case (133)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -c
              WeiGP = a_weight*c_weight*c_weight
            case (134)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight    
              case (135)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight  
              case (136)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight 
              case (137)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = d_weight*b_weight*c_weight  
              case (138)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = c_weight*b_weight*c_weight  
              case (139)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = b_weight*b_weight*c_weight  
              case (140)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -c
              WeiGP = a_weight*b_weight*c_weight  
              case (141)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight  
              case (142)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight  
              case (143)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight  
              case (144)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = d_weight*a_weight*c_weight  
              case (145)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = c_weight*a_weight*c_weight  
              case (146)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = b_weight*a_weight*c_weight  
              case (147)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -c
              WeiGP = a_weight*a_weight*c_weight
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              case (148)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = a_weight*a_weight*d_weight
            case (149)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = b_weight*a_weight*d_weight
            case (150)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = c_weight*a_weight*d_weight
            case (151)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = d_weight*a_weight*d_weight
            case (152)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = c_weight*a_weight*d_weight
            case (153)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = b_weight*a_weight*d_weight
            case (154)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = d
              WeiGP = a_weight*a_weight*d_weight
            case (155)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = a_weight*b_weight*d_weight
            case (156)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = b_weight*b_weight*d_weight
            case (157)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = c_weight*b_weight*d_weight
            case (158)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = d_weight*b_weight*d_weight      
            case (159)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = c_weight*b_weight*d_weight
            case (160)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = b_weight*b_weight*d_weight
            case (161)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = d
              WeiGP = a_weight*b_weight*d_weight
            case (162)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = a_weight*c_weight*d_weight
            case (163)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = b_weight*c_weight*d_weight
            case (164)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = c_weight*c_weight*d_weight
            case (165)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = d_weight*c_weight*d_weight
            case (166)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = c_weight*c_weight*d_weight
            case (167)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = b_weight*c_weight*d_weight
            case (168)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = d
              WeiGP = a_weight*c_weight*d_weight
            case (169)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = a_weight*d_weight*d_weight
            case (170)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = b_weight*d_weight*d_weight
            case (171)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = c_weight*d_weight*d_weight
            case (172)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = d_weight*d_weight*d_weight
            case (173)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = c_weight*d_weight*d_weight
            case (174)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = b_weight*d_weight*d_weight
            case (175)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = d
              WeiGP = a_weight*d_weight*d_weight 
            case (176)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = a_weight*c_weight*d_weight
            case (177)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = b_weight*c_weight*d_weight
            case (178)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = c_weight*c_weight*d_weight
            case (179)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = d_weight*c_weight*d_weight
            case (180)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = c_weight*c_weight*d_weight
            case (181)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = b_weight*c_weight*d_weight
            case (182)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = d
              WeiGP = a_weight*c_weight*d_weight
            case (183)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = a_weight*b_weight*d_weight    
              case (184)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = b_weight*b_weight*d_weight  
              case (185)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = c_weight*b_weight*d_weight 
              case (186)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = d_weight*b_weight*d_weight  
              case (187)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = c_weight*b_weight*d_weight  
              case (188)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = b_weight*b_weight*d_weight  
              case (189)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = d
              WeiGP = a_weight*b_weight*a_weight  
              case (190)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = a_weight*a_weight*d_weight  
              case (191)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = b_weight*a_weight*d_weight  
              case (192)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = c_weight*a_weight*d_weight  
              case (193)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = d_weight*a_weight*d_weight  
              case (194)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = c_weight*a_weight*d_weight  
              case (195)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = b_weight*a_weight*d_weight  
              case (196)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = d
              WeiGP = a_weight*a_weight*d_weight
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              case (197)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight
            case (198)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight
            case (199)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight
            case (200)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = d_weight*a_weight*c_weight
            case (201)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight
            case (202)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight
            case (203)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight
            case (204)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight
            case (205)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight
            case (206)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight
            case (207)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = d_weight*b_weight*c_weight      
            case (208)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight
            case (209)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight
            case (210)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight
            case (211)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (212)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (213)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (214)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = d_weight*c_weight*c_weight
            case (215)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (216)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (217)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (218)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = a_weight*d_weight*c_weight
            case (219)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = b_weight*d_weight*c_weight
            case (220)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = c_weight*d_weight*c_weight
            case (221)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = d_weight*d_weight*c_weight
            case (222)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = c_weight*d_weight*c_weight
            case (223)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = b_weight*d_weight*c_weight
            case (224)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = c
              WeiGP = a_weight*d_weight*c_weight 
            case (225)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (226)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (227)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (228)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = d_weight*c_weight*c_weight
            case (229)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = c_weight*c_weight*c_weight
            case (230)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = b_weight*c_weight*c_weight
            case (231)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = c
              WeiGP = a_weight*c_weight*c_weight
            case (232)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight    
              case (233)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight  
              case (234)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight 
              case (235)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = d_weight*b_weight*c_weight  
              case (236)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = c_weight*b_weight*c_weight  
              case (237)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = b_weight*b_weight*c_weight  
              case (238)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = c
              WeiGP = a_weight*b_weight*c_weight  
              case (239)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight  
              case (240)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight  
              case (241)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight  
              case (242)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = d_weight*a_weight*c_weight  
              case (243)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = c_weight*a_weight*c_weight  
              case (244)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = b_weight*a_weight*c_weight  
              case (245)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = c
              WeiGP = a_weight*a_weight*c_weight
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              case (246)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight
            case (247)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight
            case (248)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight
            case (249)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = d_weight*a_weight*b_weight
            case (250)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight
            case (251)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight
            case (252)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight
            case (253)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight
            case (254)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight
            case (255)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight
            case (256)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = d_weight*b_weight*b_weight      
            case (257)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight
            case (258)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight
            case (259)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight
            case (260)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (261)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (262)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (263)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = d_weight*c_weight*b_weight
            case (264)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (265)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (266)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (267)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = a_weight*d_weight*b_weight
            case (268)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = b_weight*d_weight*b_weight
            case (269)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = c_weight*d_weight*b_weight
            case (270)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = d_weight*d_weight*b_weight
            case (271)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = c_weight*d_weight*b_weight
            case (272)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = b_weight*d_weight*b_weight
            case (273)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = b
              WeiGP = a_weight*d_weight*b_weight 
            case (274)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (275)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (276)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (277)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = d_weight*c_weight*b_weight
            case (278)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = c_weight*c_weight*b_weight
            case (279)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = b_weight*c_weight*b_weight
            case (280)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = b
              WeiGP = a_weight*c_weight*b_weight
            case (281)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight    
              case (282)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight  
              case (283)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight 
              case (284)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = d_weight*b_weight*b_weight  
              case (285)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = c_weight*b_weight*b_weight  
              case (286)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = b_weight*b_weight*b_weight  
              case (287)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = b
              WeiGP = a_weight*b_weight*b_weight  
              case (288)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight  
              case (289)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight  
              case (290)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight  
              case (291)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = d_weight*a_weight*b_weight  
              case (292)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = c_weight*a_weight*b_weight  
              case (293)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = b_weight*a_weight*b_weight  
              case (294)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = b
              WeiGP = a_weight*a_weight*b_weight
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              case (295)
              PosGP(1) = -a 
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
            case (296)
              PosGP(1) = -b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight
            case (297)
              PosGP(1) = -c
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight
            case (298)
              PosGP(1) = d
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = d_weight*a_weight*a_weight
            case (299)
              PosGP(1) = c
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight
            case (300)
              PosGP(1) = b
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight
            case (301)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
            case (302)
              PosGP(1) = -a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight
            case (303)
              PosGP(1) = -b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight
            case (304)
              PosGP(1) = -c
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight
            case (305)
              PosGP(1) = d
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = d_weight*b_weight*a_weight      
            case (306)
              PosGP(1) = c
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight
            case (307)
              PosGP(1) = b
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight
            case (308)
              PosGP(1) = a
              PosGP(2) = -b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight
            case (309)
              PosGP(1) = -a
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (310)
              PosGP(1) = -b
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (311)
              PosGP(1) = -c
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (312)
              PosGP(1) = d
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = d_weight*c_weight*a_weight
            case (313)
              PosGP(1) = c
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (314)
              PosGP(1) = b 
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (315)
              PosGP(1) = a
              PosGP(2) = -c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (316)
              PosGP(1) = -a 
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = a_weight*d_weight*a_weight
            case (317)
              PosGP(1) = -b
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = b_weight*d_weight*a_weight
            case (318)
              PosGP(1) = -c
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = c_weight*d_weight*a_weight
            case (319)
              PosGP(1) = d
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = d_weight*d_weight*a_weight
            case (320)
              PosGP(1) = c
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = c_weight*d_weight*a_weight
            case (321)
              PosGP(1) = b
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = b_weight*d_weight*a_weight
            case (322)
              PosGP(1) = a
              PosGP(2) = d
              PosGP(3) = -a
              WeiGP = a_weight*d_weight*a_weight 
            case (323)
              PosGP(1) = -a
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (324)
              PosGP(1) = -b
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (325)
              PosGP(1) = -c
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (326)
              PosGP(1) = d
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = d_weight*c_weight*a_weight
            case (327)
              PosGP(1) = c
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = c_weight*c_weight*a_weight
            case (328)
              PosGP(1) = b
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = b_weight*c_weight*a_weight
            case (329)
              PosGP(1) = a
              PosGP(2) = c
              PosGP(3) = -a
              WeiGP = a_weight*c_weight*a_weight
            case (330)
              PosGP(1) = -a
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight    
              case (331)
              PosGP(1) = -b
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight  
              case (332)
              PosGP(1) = -c
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight 
              case (333)
              PosGP(1) = d
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = d_weight*b_weight*a_weight  
              case (334)
              PosGP(1) = c
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = c_weight*b_weight*a_weight  
              case (335)
              PosGP(1) = b
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = b_weight*b_weight*a_weight  
              case (336)
              PosGP(1) = a
              PosGP(2) = b
              PosGP(3) = -a
              WeiGP = a_weight*b_weight*a_weight  
              case (337)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight  
              case (338)
              PosGP(1) = -b
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight  
              case (339)
              PosGP(1) = -c
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight  
              case (340)
              PosGP(1) = d
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = d_weight*a_weight*a_weight  
              case (341)
              PosGP(1) = c
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = c_weight*a_weight*a_weight  
              case (342)
              PosGP(1) = b
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b_weight*a_weight*a_weight  
              case (343)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = a_weight*a_weight*a_weight
              
              
              
              
              
              
              
              
              
              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialHEXA_MP343
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP1
        !
        !    DESCRIPTION:  
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   1 material point is placed in each element, whose initial
        !>   location and weight is identical with those of the Gauss Point.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP1(IParticle, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.0
              PosGP(2) = 0.0 !position right in the center of the local element
              WeiGP = 4.0 !why is the gauss point weight 4? This should be two
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP1()].")   
          end select
        
        end subroutine InitialQUAD_MP1


        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP4
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   4 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP4(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.57735026918962576450914878050196 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 1.0
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = b 
            case (2)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = b
            case (3)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = b
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP4()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP4
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP9
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   9 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP9(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.77459666924 ! = sqrt(3/5)
          !real(REAL_TYPE) :: b = 0.0
          real(REAL_TYPE) :: c = 0.88888888888
          real(REAL_TYPE) :: d = 0.55555555555


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = d*d 
            case (2)
              PosGP(1) = 0.0
              PosGP(2) = -a
              WeiGP = d*c
            case (3)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = d*d
            case (4)
              PosGP(1) = -a
              PosGP(2) = 0.0
              WeiGP = c*d
            case (5)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = c*c
            case (6)
              PosGP(1) = a
              PosGP(2) = 0.0
              WeiGP = c*d
            case (7)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = d*d
            case (8)
              PosGP(1) = 0.0
              PosGP(2) = a
              WeiGP = d*c
            case (9)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = d*d
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP4()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP9
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP16
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   16 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP16(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.86113631159
          real(REAL_TYPE) :: b = 0.33998104358 
          !real(REAL_TYPE) :: c = +0.33998104358
          !real(REAL_TYPE) :: d = +0.86113631159
          
          
          real(REAL_TYPE) :: e = 0.34785484513 ! --> a
          real(REAL_TYPE) :: f = 0.65214515486 ! --> b
          


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = e*e
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = e*f
            case (3)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = e*f
            case (4)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = e*e
            case (5)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = e*f
            case (6)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = f*f
            case (7)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = f*f
            case (8)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = e*f
            case (9)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = e*f              
            case (10)
              PosGP(1) = -b 
              PosGP(2) = b
              WeiGP = f*f 
            case (11)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = f*f
            case (12)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = e*f
            case (13)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = e*e
            case (14)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = e*f
            case (15)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = e*f
            case (16)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = e*e

            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP16()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP16
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP81
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   81 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP81(IParticle, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.9681602395076261
          real(REAL_TYPE) :: b = 0.8360311073266358
          real(REAL_TYPE) :: c = 0.6133714327005904
          real(REAL_TYPE) :: d = 0.3242534234038089
        
          real(REAL_TYPE) :: a_weight = 0.0812743883615744
          real(REAL_TYPE) :: b_weight = 0.1806481606948574
          real(REAL_TYPE) :: c_weight = 0.2606106964029354
          real(REAL_TYPE) :: d_weight = 0.3123470770400029
        
        
          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            ! row 1
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (4)
              PosGP(1) = -d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
            case (5)
              PosGP(1) = 0.0
              PosGP(2) = -a
              WeiGP = 0.3302393550012598*a_weight
            case (6)
              PosGP(1) = d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
            case (7)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (8)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (9)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              
              ! row 2
            case (10)
              PosGP(1) = -a 
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
            case (11)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (12)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
            case (13)
              PosGP(1) = -d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight
            case (14)
              PosGP(1) = 0.0
              PosGP(2) = -b
              WeiGP = 0.3302393550012598*b_weight
            case (15)
              PosGP(1) = d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight
            case (16)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
            case (17)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (18)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              
              ! row 3
            case (19)
              PosGP(1) = -a 
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
            case (20)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (21)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (22)
              PosGP(1) = -d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
            case (23)
              PosGP(1) = 0.0
              PosGP(2) = -c
              WeiGP = 0.3302393550012598*c_weight
            case (24)
              PosGP(1) = d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
            case (25)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (26)
              PosGP(1) = b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (27)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              
              
              ! row 4
            case (28)
              PosGP(1) = -a 
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
            case (29)
              PosGP(1) = -b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
            case (30)
              PosGP(1) = -c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
            case (31)
              PosGP(1) = -d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
            case (32)
              PosGP(1) = 0.0
              PosGP(2) = -d
              WeiGP = 0.3302393550012598*d_weight
            case (33)
              PosGP(1) = d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
            case (34)
              PosGP(1) = c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
            case (35)
              PosGP(1) = b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
            case (36)
              PosGP(1) = a
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
            
              
              
              
              ! row 5
            case (37)
              PosGP(1) = -a 
              PosGP(2) = 0.0
              WeiGP = a_weight*0.3302393550012598
            case (38)
              PosGP(1) = -b
              PosGP(2) = 0.0
              WeiGP = b_weight*0.3302393550012598
            case (39)
              PosGP(1) = -c
              PosGP(2) = 0.0
              WeiGP = c_weight*0.3302393550012598
            case (40)
              PosGP(1) = -d
              PosGP(2) = 0.0
              WeiGP = d_weight*0.3302393550012598
            case (41)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = 0.3302393550012598*0.3302393550012598
            case (42)
              PosGP(1) = d
              PosGP(2) = 0.0
              WeiGP = d_weight*0.3302393550012598
            case (43)
              PosGP(1) = c
              PosGP(2) = 0.0
              WeiGP = c_weight*0.3302393550012598
            case (44)
              PosGP(1) = b
              PosGP(2) = 0.0
              WeiGP = b_weight*0.3302393550012598
            case (45)
              PosGP(1) = a
              PosGP(2) = 0.0
              WeiGP = a_weight*0.3302393550012598
              
              
              
              ! row 6
            case (46)
              PosGP(1) = -a 
              PosGP(2) = d
              WeiGP = a_weight*d_weight
            case (47)
              PosGP(1) = -b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
            case (48)
              PosGP(1) = -c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
            case (49)
              PosGP(1) = -d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
            case (50)
              PosGP(1) = 0.0
              PosGP(2) = d
              WeiGP = 0.3302393550012598*d_weight
            case (51)
              PosGP(1) = d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
            case (52)
              PosGP(1) = c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
            case (53)
              PosGP(1) = b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
            case (54)
              PosGP(1) = a
              PosGP(2) = d
              WeiGP = a_weight*d_weight
              
              
              
              
              ! row 7
            case (55)
              PosGP(1) = -a 
              PosGP(2) = c
              WeiGP = a_weight*c_weight
            case (56)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (57)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (58)
              PosGP(1) = -d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
            case (59)
              PosGP(1) = 0.0
              PosGP(2) = c
              WeiGP = 0.3302393550012598*c_weight
            case (60)
              PosGP(1) = d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
            case (61)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (62)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (63)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              
              ! row 8
            case (64)
              PosGP(1) = -a 
              PosGP(2) = b
              WeiGP = a_weight*b_weight
            case (65)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight
            case (66)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight
            case (67)
              PosGP(1) = -d
              PosGP(2) = b
              WeiGP = d_weight*b_weight
            case (68)
              PosGP(1) = 0.0
              PosGP(2) = b
              WeiGP = 0.3302393550012598*b_weight
            case (69)
              PosGP(1) = d
              PosGP(2) = b
              WeiGP = d_weight*b_weight
            case (70)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight
            case (71)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight
            case (72)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight
              
              
              
              ! row 9
            case (73)
              PosGP(1) = -a 
              PosGP(2) = a
              WeiGP = a_weight*a_weight
            case (74)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight
            case (75)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight
            case (76)
              PosGP(1) = -d
              PosGP(2) = a
              WeiGP = d_weight*a_weight
            case (77)
              PosGP(1) = 0.0
              PosGP(2) = a
              WeiGP = 0.3302393550012598*a_weight
            case (78)
              PosGP(1) = d
              PosGP(2) = a
              WeiGP = d_weight*a_weight
            case (79)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight
            case (80)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight
            case (81)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select
        
          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )
        
        end subroutine InitialQUAD_MP81
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP64
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   64 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP64(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.9602898564975363
          real(REAL_TYPE) :: b = 0.7966664774136267
          real(REAL_TYPE) :: c = 0.5255324099163290
          real(REAL_TYPE) :: d = 0.1834346424956498
  
          real(REAL_TYPE) :: a_weight = 0.1012285362903763
          real(REAL_TYPE) :: b_weight = 0.2223810344533745
          real(REAL_TYPE) :: c_weight = 0.3137066458778873
          real(REAL_TYPE) :: d_weight = 0.3626837833783620


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (4)
              PosGP(1) = -d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
            case (5)
              PosGP(1) = d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
            case (6)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (7)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (8)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              
              
              
            case (9)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
            case (10)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (11)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight      
            case (12)
              PosGP(1) = -d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight      
            case (13)
              PosGP(1) = d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight  
            case (14)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
            case (15)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (16)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              
            case (17)
              PosGP(1) = -a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
            case (18)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (19)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (20)
              PosGP(1) = -d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
            case (21)
              PosGP(1) = d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
            case (22)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (23)
              PosGP(1) = b 
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (24)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              
              
              case (25)
              PosGP(1) = -a 
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
            case (26)
              PosGP(1) = -b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
            case (27)
              PosGP(1) = -c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
            case (28)
              PosGP(1) = -d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
            case (29)
              PosGP(1) = d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
            case (30)
              PosGP(1) = c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
            case (31)
              PosGP(1) = b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
            case (32)
              PosGP(1) = a
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
              
              
            case (33)
              PosGP(1) = -a 
              PosGP(2) = d
              WeiGP = a_weight*d_weight
            case (34)
              PosGP(1) = -b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
            case (35)
              PosGP(1) = -c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
            case (36)
              PosGP(1) = -d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
            case (37)
              PosGP(1) = d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
            case (38)
              PosGP(1) = c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
            case (39)
              PosGP(1) = b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
            case (40)
              PosGP(1) = a
              PosGP(2) = d
              WeiGP = a_weight*d_weight
            
              
              
              
              
            case (41)
              PosGP(1) = -a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
            case (42)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (43)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (44)
              PosGP(1) = -d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
            case (45)
              PosGP(1) = d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
            case (46)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (47)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (48)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              
            case (49)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = a_weight*b_weight    
              case (50)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              case (51)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight 
              case (52)
              PosGP(1) = -d
              PosGP(2) = b
              WeiGP = d_weight*b_weight  
              case (53)
              PosGP(1) = d
              PosGP(2) = b
              WeiGP = d_weight*b_weight  
              case (54)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              case (55)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              case (56)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight  
              
              
              
              
              
              
              case (57)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              case (58)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              case (59)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              case (60)
              PosGP(1) = -d
              PosGP(2) = a
              WeiGP = d_weight*a_weight  
              case (61)
              PosGP(1) = d
              PosGP(2) = a
              WeiGP = d_weight*a_weight  
              case (62)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              case (63)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              case (64)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP64
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP49
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   49 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP49(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.9491079123427585
          real(REAL_TYPE) :: b = 0.7415311855993945
          real(REAL_TYPE) :: c = 0.4058451513773972
          real(REAL_TYPE) :: d = 0.0
  
          real(REAL_TYPE) :: a_weight = 0.1294849661688697
          real(REAL_TYPE) :: b_weight = 0.2797053914892766
          real(REAL_TYPE) :: c_weight = 0.3818300505051189
          real(REAL_TYPE) :: d_weight = 0.4179591836734694 


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (4)
              PosGP(1) = d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
            case (5)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (6)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (7)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              
              
              
            case (8)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
            case (9)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (10)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight      
            case (11)
              PosGP(1) = d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight      
            case (12)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
            case (13)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (14)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              
            case (15)
              PosGP(1) = -a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
            case (16)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (17)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (18)
              PosGP(1) = d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
            case (19)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (20)
              PosGP(1) = b 
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (21)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              
              
            case (22)
              PosGP(1) = -a 
              PosGP(2) = d
              WeiGP = a_weight*d_weight
            case (23)
              PosGP(1) = -b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
            case (24)
              PosGP(1) = -c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
            case (25)
              PosGP(1) = d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
            case (26)
              PosGP(1) = c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
            case (27)
              PosGP(1) = b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
            case (28)
              PosGP(1) = a
              PosGP(2) = d
              WeiGP = a_weight*d_weight
            
              
              
              
              
            case (29)
              PosGP(1) = -a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
            case (30)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (31)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (32)
              PosGP(1) = d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
            case (33)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (34)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (35)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              
            case (36)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = a_weight*b_weight    
              case (37)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              case (38)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight 
              case (39)
              PosGP(1) = d
              PosGP(2) = b
              WeiGP = d_weight*b_weight  
              case (40)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              case (41)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              case (42)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight  
              
              
              
              
              
              
              case (43)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              case (44)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              case (45)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              case (46)
              PosGP(1) = d
              PosGP(2) = a
              WeiGP = d_weight*a_weight  
              case (47)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              case (48)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              case (49)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP49
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP36
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   36 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP36(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.9324695142031521
          real(REAL_TYPE) :: b = 0.6612093864662645
          real(REAL_TYPE) :: c = 0.2386191860831969
  
          real(REAL_TYPE) :: a_weight = 0.1713244923791704 
          real(REAL_TYPE) :: b_weight = 0.3607615730481386 
          real(REAL_TYPE) :: c_weight = 0.4679139345726910 
          


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (4)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
            case (5)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            case (6)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
            case (7)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
            case (8)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (9)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight             
            case (10)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
            case (11)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
            case (12)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
            case (13)
              PosGP(1) = -a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
            case (14)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (15)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (16)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
            case (17)
              PosGP(1) = b 
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
            case (18)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
            case (19)
              PosGP(1) = -a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
            case (20)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (21)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (22)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
            case (23)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
            case (24)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
            case (25)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = a_weight*b_weight    
              case (26)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              case (27)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              case (28)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              case (29)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              case (30)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight  
              case (31)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              case (32)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              case (33)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              case (34)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              case (35)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              case (36)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP36
        
        
        
        
        
        !**********************************************************************
        !
        !    SUBROUTINE: InitialQUAD_MP25
        !
        !>   Returns the initial local coordinates and weight for IParticle.
        !>   25 material points are placed in each element, whose initial
        !>   locations and weights are identical with those of Gauss Points.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF, page 3
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !    
        !>   @param[in] IParticle : Local number of the considered particle inside an element
        !
        !>   @param[inout] PosGP : Returns the initial local coordinates of the material point with local number IParticle
        !>   @param[inout] WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        subroutine InitialQUAD_MP25(IParticle, PosGP, WeiGP)
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.90617984593
          real(REAL_TYPE) :: b = 0.5384693101
          
          real(REAL_TYPE) :: c = 0.23692688505 ! --> a
          real(REAL_TYPE) :: d = 0.56888888888 ! --> 0
          real(REAL_TYPE) :: e = 0.47862867049 ! --> b
          


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = c*c
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = e*c
            case (3)
              PosGP(1) = 0.0
              PosGP(2) = -a
              WeiGP = d*c
            case (4)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = c*e
            case (5)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = c*c
            case (6)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = c*e
            case (7)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = e*e
            case (8)
              PosGP(1) = 0
              PosGP(2) = -b
              WeiGP = d*e
            case (9)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = e*e              
            case (10)
              PosGP(1) = a 
              PosGP(2) = -b
              WeiGP = c*e 
            case (11)
              PosGP(1) = -a
              PosGP(2) = 0
              WeiGP = c*d
            case (12)
              PosGP(1) = -b
              PosGP(2) = 0
              WeiGP = d*e
            case (13)
              PosGP(1) = 0
              PosGP(2) = 0
              WeiGP = d*d
            case (14)
              PosGP(1) = b
              PosGP(2) = 0
              WeiGP = e*d
            case (15)
              PosGP(1) = a
              PosGP(2) = 0
              WeiGP = c*d
            case (16)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = c*e
              
            case (17)
              PosGP(1) = -b 
              PosGP(2) = b
              WeiGP = e*e
            case (18)
              PosGP(1) = 0.0
              PosGP(2) = b
              WeiGP = d*e
            case (19)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = e*e
            case (20)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = c*e
            case (21)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = c*c
            case (22)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = e*c
            case (23)
              PosGP(1) = 0
              PosGP(2) = a
              WeiGP = d*c
            case (24)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = c*e
            case (25)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = c*c              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP25
        
        
        
        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        !**********************************************************************
        !
        !    SUBROUTINE: GaussQUAD_Q1
        !
        !    DESCRIPTION:
        !>   Returns the local coordinates and weight for IGaussPoint.
        !>   1 Gauss Point is located in the quadrilaters element.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !
        !>   @param[in] IGaussPoint : Local number of the considered Gauss Point inside the element
        !
        !>   @param[inout] PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        !>   @param[inout] WeiGP : Returns the weight of IGaussPoint
        ! 
        !**********************************************************************
        subroutine GaussQUAD_Q1(IGaussPoint, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = 4.0 ! (2*2 because of each direction -> notice that the total is always 4 for 2D)
              
            !case(2)
            !    
            !    PosGP(1) = 
            !    
            !case(3)
            !    
            !case(4)
                
            case default
              call GiveError("Undefined number of Gauss points in [subroutine GaussQUAD_Q1()].")    
          end select

        end subroutine GaussQUAD_Q1


        !**********************************************************************
        !
        !    SUBROUTINE: GaussQUAD_Q4
        !
        !    DESCRIPTION:
        !>   Returns the local coordinates and weight for IGaussPoint.
        !>   4 Gauss Points are located in the quadrilateral element.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !
        !>   @param[in] IGaussPoint : Local number of the considered Gauss Point inside the element
        !
        !>   @param[inout] PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        !>   @param[inout] WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************
        subroutine GaussQUAD_Q4(IGaussPoint, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.57735026918962576450914878050196 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 1.0
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = -a
              PosGP(2) = -a
              WeiGP = b
            case (2)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = b
            case (3)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = b
            case default
              call GiveError("Undefined number of Gauss points in [subroutine GaussQUAD_Q4()].")  
          end select
         
        end subroutine GaussQUAD_Q4
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
!**********************************************************************
        !
        !    SUBROUTINE: GaussQUAD_Q4
        !
        !    DESCRIPTION:
        !>   Returns the local coordinates and weight for IGaussPoint.
        !>   4 Gauss Points are located in the quadrilateral element.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !
        !>   @param[in] IGaussPoint : Local number of the considered Gauss Point inside the element
        !
        !>   @param[inout] PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        !>   @param[inout] WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************
        subroutine GaussHEXA_Q8(IGaussPoint, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.57735026918962576450914878050196 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 1.0
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = -a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b
            case (2)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = -a
              WeiGP = b
            case (3)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b
            case (4)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = -a
              WeiGP = b  
            case (5)
              PosGP(1) = -a
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = b
            case (6)
              PosGP(1) = a
              PosGP(2) = -a
              PosGP(3) = a
              WeiGP = b
            case (7)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b
            case (8)
              PosGP(1) = -a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b
            case default
              call GiveError("Undefined number of Gauss points in [subroutine GaussQUAD_Q4()].")  
          end select
         
        end subroutine GaussHexa_Q8
        
        
        
        
        
        
        
        
        
        
        
                
        
        
!**********************************************************************
        !
        !    SUBROUTINE: GaussQUAD_Q1
        !
        !    DESCRIPTION:
        !>   Returns the local coordinates and weight for IGaussPoint.
        !>   4 Gauss Points are located in the quadrilateral element.
        !
        !>   @ note : 2D element
        !>   @ note : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !>   @ note : T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987
        !
        !>   @param[in] IGaussPoint : Local number of the considered Gauss Point inside the element
        !
        !>   @param[inout] PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        !>   @param[inout] WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************
        subroutine GaussHEXA_Q1(IGaussPoint, PosGP, WeiGP)
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          real(REAL_TYPE) :: a = 0.0 ! = 1 / sqrt(3)
          real(REAL_TYPE) :: b = 8.0
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = a
              PosGP(2) = a
              PosGP(3) = a
              WeiGP = b
            !case (2)
            !  PosGP(1) = a
            !  PosGP(2) = -a
            !  PosGP(3) = -a
            !  WeiGP = b
            !case (3)
            !  PosGP(1) = a
            !  PosGP(2) = a
            !  PosGP(3) = -a
            !  WeiGP = b
            !case (4)
            !  PosGP(1) = -a
            !  PosGP(2) = a
            !  PosGP(3) = -a
            !  WeiGP = b  
            !case (5)
            !  PosGP(1) = -a
            !  PosGP(2) = -a
            !  PosGP(3) = a
            !  WeiGP = b
            !case (6)
            !  PosGP(1) = a
            !  PosGP(2) = -a
            !  PosGP(3) = a
            !  WeiGP = b
            !case (7)
            !  PosGP(1) = a
            !  PosGP(2) = a
            !  PosGP(3) = a
            !  WeiGP = b
            !case (8)
            !  PosGP(1) = -a
            !  PosGP(2) = a
            !  PosGP(3) = a
            !  WeiGP = b
            case default
              call GiveError("Undefined number of Gauss points in [subroutine GaussQUAD_Q4()].")  
          end select
         
        end subroutine GaussHexa_Q1
        
        
        
        
        
        
        
        subroutine InitialiseShapeFunctionsQUAD4(HS, dHS, Wt)
        !**********************************************************************
        !
        !    SUBROUTINE: InitialiseShapeFunctionsQUAD4
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at  one Gaussian integration point for a 4-noded 2D quadrilateral element.
        !
        !>   @note : 2D element
        !>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983
        !
        !>   @param[in/out] HS(i,j) : Value of shape function j at integration point i
        !>   @param[in/out] dHS(i,j,k) : Value of derivative of shape function j at integration point i with respect to direction k
        !>   @param[in/out] Wt : Local weights for integration 
        !
        !             4) (-1,1)   ^ Eta    3) (1,1)
        !                 4       |
        !                +---------------+ 3
        !                |        |      |
        !                |        |      |
        !                |        |      |
        !                |        -------|---> Xi
        !                |               |
        !                |               |
        !                |1              | 2
        !                +---------------+-
        !             1) (-1,-1)           2) (-1,1)
        !**********************************************************************

        implicit none
      
          !real(REAL_TYPE), dimension(:), intent(inout) :: LocPos
          real(REAL_TYPE), dimension(:,:), intent(inout) :: HS
          real(REAL_TYPE), dimension(:,:,:), intent(inout) :: dHS
          real(REAL_TYPE), dimension(:), intent(inout) :: Wt

          ! local variables
          real(REAL_TYPE) :: Xi, Eta
          integer(INTEGER_TYPE) :: int, I1, Nint1
          
          ! Note this is two dimensional 
          
          Nint1=1 !number of gauss points
          Int = 0 !counter
          
          do I1 = 1, Nint1

              Xi = 0.0 !local position in Xi (local) direction 
              Eta = 0.0 !local position in Eta (local) direction
              
              Int = Int+1
              
              Wt(Int) = 2.0 !1d0 / Nint1 * 0.5 !This should be =2... double check!!!! 
           
              ! HS(i)
              HS(Int, 1) = (1.0 - Xi) * (1.0 - Eta) / 4.0 ! a=1
              HS(Int, 2) = (1.0 + Xi) * (1.0 - Eta) / 4.0 ! a=2
              HS(Int, 3) = (1.0 + Xi) * (1.0 + Eta) / 4.0 ! a=3
              HS(Int, 4) = (1.0 - Xi) * (1.0 + Eta) / 4.0 ! a=4
         
              ! dHS(i,1) = dHS / dXi
              dHS(Int,1,1) =  - (1.0 - Eta) / 4.0 ! a=1
              dHS(Int,2,1) =    (1.0 - Eta) / 4.0 ! a=2
              dHS(Int,3,1) =    (1.0 + Eta) / 4.0 ! a=3
              dHS(Int,4,1) =  - (1.0 + Eta) / 4.0 ! a=4
        
              ! dHS(i,2) = dHS / dEta
              dHS(Int,1,2) =  - (1.0 - Xi) / 4.0 ! a=1
              dHS(Int,2,2) =  - (1.0 + Xi) / 4.0 ! a=2
              dHS(Int,3,2) =    (1.0 + Xi) / 4.0 ! a=3
              dHS(Int,4,2) =    (1.0 - Xi) / 4.0 ! a=4
          
          end do
          

        end subroutine InitialiseShapeFunctionsQUAD4


        !**********************************************************************
        !
        !    SUBROUTINE: ShapeLocPosQUAD4
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at LocPos for a 4-noded 2D quadrilateral element.
        !
        !>   @note : 2D element
        !>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983
        !
        !>   @param[in] LocPos : Local coordinates of a point inside an element
        !
        !>   @param[out] HS(i) : Value of shape function i at LocPos
        !>   @param[out] dHS(i,j) : Value of derivative of shape function i at LocPos with respect to direction j
        !
        !             4) (-1,1)   ^ Eta    3) (1,1)
        !                 4       |
        !                +---------------+ 3
        !                |        |      |
        !                |        |      |
        !                |        |      |
        !                |        -------|---> Xi
        !                |               |
        !                |               |
        !                |1              | 2
        !                +---------------+-
        !             1) (-1,-1)           2) (-1,1)
        !**********************************************************************
        subroutine ShapeLocPosQUAD4(LocPos, HS, dHS)

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          real(REAL_TYPE), dimension(:), intent(out) :: HS
          real(REAL_TYPE), dimension(:, :), intent(out) :: dHS

          ! local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi = LocPos(1) !local position in Xi (local) direction 
          Eta = LocPos(2) !local position in Eta (local) direction

          ! HS(i)
          HS(1) = (1.0 - Xi) * (1.0 - Eta) / 4.0 ! a=1
          HS(2) = (1.0 + Xi) * (1.0 - Eta) / 4.0 ! a=2
          HS(3) = (1.0 + Xi) * (1.0 + Eta) / 4.0 ! a=3
          HS(4) = (1.0 - Xi) * (1.0 + Eta) / 4.0 ! a=4

          ! dHS(i,1) = dHS / dXi
          dHS(1,1) =  - (1.0 - Eta) / 4.0 ! a=1
          dHS(2,1) =    (1.0 - Eta) / 4.0 ! a=2
          dHS(3,1) =    (1.0 + Eta) / 4.0 ! a=3
          dHS(4,1) =  - (1.0 + Eta) / 4.0 ! a=4

          ! dHS(i,2) = dHS / dEta
          dHS(1,2) =  - (1.0 - Xi) / 4.0 ! a=1
          dHS(2,2) =  - (1.0 + Xi) / 4.0 ! a=2
          dHS(3,2) =    (1.0 + Xi) / 4.0 ! a=3
          dHS(4,2) =    (1.0 - Xi) / 4.0 ! a=4

        end subroutine ShapeLocPosQUAD4


        !**********************************************************************
        !
        !    SUBROUTINE: ShapeLocPosQUAD8
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at LocPos for a 8-noded 2D quadrilateral element.
        !
        !>   @note : 2D element
        !>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983, p.81
        !
        !>   @param[in] LocPos : Local coordinates of a point inside an element
        !
        !>   @param[out] HS(i) : Value of shape function i at LocPos
        !>   @param[out] dHS(i,j) : Value of derivative of shape function i at LocPos with respect to direction j
        !
        !                         ^ Eta
        !                 4       | 7
        !                +--------+-------+ 3
        !                |        |       |
        !                |        |       |
        !                | 8      |       | 6
        !                +        --------|---> Xi
        !                |                |
        !                |                |
        !                |1       5       | 2
        !                +--------+-------+
        !
        !**********************************************************************
        subroutine ShapeLocPosQUAD8(LocPos, HS, dHS)

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          real(REAL_TYPE), dimension(:), intent(out) :: HS
          real(REAL_TYPE), dimension(:, :), intent(out) :: dHS

          ! local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi = LocPos(1)
          Eta = LocPos(2)

          ! HS(i)
          ! 1..4 corner nodes, 5..8 middle nodes
          HS(1) = - (1.0 - Xi) * (1.0 - Eta) * (1.0 + Xi + Eta) / 4.0   
          HS(2) = - (1.0 + Xi) * (1.0 - Eta) * (1.0 - Xi + Eta) / 4.0 
          HS(3) = - (1.0 + Xi) * (1.0 + Eta) * (1.0 - Xi - Eta) / 4.0
          HS(4) = - (1.0 - Xi) * (1.0 + Eta) * (1.0 + Xi - Eta) / 4.0
          HS(5) =   (1.0 - Xi) * (1.0 + Xi) * (1.0 - Eta) / 2.0
          HS(6) =   (1.0 + Xi) * (1.0 + Eta) * (1.0 - Eta) / 2.0
          HS(7) =   (1.0 - Xi) * (1.0 + Xi) * (1.0 + Eta) / 2.0
          HS(8) =   (1.0 - Xi) * (1.0 + Eta) * (1.0 - Eta) / 2.0

          ! dHS(i,1) = dHS / dXi
          ! 1..4 corner nodess, 5..8 middle nodes
          dHS(1, 1) = - (-1.0) * (1.0 - Eta) * (1.0 + Xi + Eta) / 4.0  - (1.0 - Xi) * (1.0 - Eta) * (1.0) / 4.0
          dHS(2, 1) = - ( 1.0) * (1.0 - Eta) * (1.0 - Xi + Eta) / 4.0  - (1.0 + Xi) * (1.0 - Eta) * (-1.0) / 4.0
          dHS(3, 1) = - ( 1.0) * (1.0 + Eta) * (1.0 - Xi - Eta) / 4.0  - (1.0 + Xi) * (1.0 + Eta) * (-1.0) / 4.0
          dHS(4, 1) = - (-1.0) * (1.0 + Eta) * (1.0 + Xi - Eta) / 4.0  - (1.0 - Xi) * (1.0 + Eta) * (1.0) / 4.0
          dHS(5, 1) = (-1.0) * (1.0 + Xi) * (1.0 - Eta) / 2.0 + (1.0 - Xi) * (1.0) * (1.0 - Eta) / 2.0
          dHS(6, 1) = (1.0) * (1.0 + Eta) * (1.0 - Eta) / 2.0 
          dHS(7, 1) = (-1.0) * (1.0 + Xi) * (1.0 + Eta) / 2.0 + (1.0 - Xi) * (1.0) * (1.0 + Eta) / 2.0
          dHS(8, 1) = (-1.0) * (1.0 + Eta) * (1.0 - Eta) / 2.0 

          ! dHS(i,2) = dHS / dEta
          ! 1..4 corner nodes, 5..8 middle nodes
          dHS(1, 2) = - (1.0 - Xi) * (-1.0) * (1.0 + Xi + Eta) / 4.0 - (1.0 - Xi) * (1.0 - Eta) * (1.0) / 4.0 
          dHS(2, 2) = - (1.0 + Xi) * (-1.0) * (1.0 - Xi + Eta) / 4.0 - (1.0 + Xi) * (1.0 - Eta) * (1.0) / 4.0 
          dHS(3, 2) = - (1.0 + Xi) * (1.0) * (1.0 - Xi - Eta) / 4.0 - (1.0 + Xi) * (1.0 + Eta) * (-1.0) / 4.0
          dHS(4, 2) = - (1.0 - Xi) * (1.0) * (1.0 + Xi - Eta) / 4.0 - (1.0 - Xi) * (1.0 + Eta) * (-1.0) / 4.0
          dHS(5, 2) = (1.0 - Xi) * (1.0 + Xi) * (-1.0) / 2.0
          dHS(6, 2) = (1.0 + Xi) * (1.0) * (1.0 - Eta) / 2.0 + (1.0 + Xi) * (1.0 + Eta) * (- 1.0) / 2.0
          dHS(7, 2) = (1.0 - Xi) * (1.0 + Xi) * (1.0) / 2.0
          dHS(8, 2) = (1.0 - Xi) * (1.0) * (1.0 - Eta) / 2.0 + (1.0 - Xi) * (1.0 + Eta) * (-1.0) / 2.0

        end subroutine ShapeLocPosQUAD8


        !**********************************************************************
        !
        !    SUBROUTINE: GetLocalCoordinatesQUAD4
        !
        !    DESCRIPTION:
        !>   Determination of local coordinates LocPos from global coordinates GlobPos, 
        !>   assuming that the point lies inside the quadrilateral element.
        !>   OutsideElement returns .false. if the local position is in the element.
        !>   CrossedSide returns the number of the side that has been crossed.
        !
        !>   @param[in] GlobPos : Global coordinates of a point inside the element
        !>   @param[in] MInv : Element matrix
        !>   @param[in] MIX1 : Element vector of first element node
        !
        !>   @param[out] LocPos : Local coordinates of the considered point inside IElement
        !>   @param[out] OutsideElement : True, if the local coordinate lie outside IElement
        !>   @param[out] CrossedSide : ID of crossed side, if local coordinate lie outside IElement
        !
        !**********************************************************************
        subroutine GetLocalCoordinatesQUAD4(GlobPos, LocPos, OutsideElement, MInv, MIX1, CrossedSide)

        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(IDim), intent(out):: LocPos
          logical, intent(out) :: OutsideElement
          real(REAL_TYPE), dimension(IDim, IDim), intent(in) :: MInv
          real(REAL_TYPE), dimension(IDim), intent(in) :: MIX1
          integer(INTEGER_TYPE), intent(out) :: CrossedSide
          
          ! local variables
          integer(INTEGER_TYPE) :: J
          
          OutsideElement = .true.

          LocPos = 0.0
          do J = 1, 2
            LocPos(1) = LocPos(1) + MInv(1, J) *  GlobPos(J)
            LocPos(2) = LocPos(2) + MInv(2, J) *  GlobPos(J)
          end do
         
          LocPos(1) = LocPos(1) - MIX1(1)
          LocPos(2) = LocPos(2) - MIX1(2)

          CrossedSide = -1
          if (LocPos(1) > 1.0) then
            CrossedSide = 2
          elseif (LocPos(1) < -1.0) then
            CrossedSide = 4
          elseif (LocPos(2) > 1.0) then
            CrossedSide = 3
          elseif (LocPos(2) < -1.0) then
            CrossedSide = 1
          else
            OutsideElement = .false.
          end if
        
        end subroutine GetLocalCoordinatesQUAD4

        
        
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 3D Hexes 
        !**********************************************************************
        !
        !    FUNCTION: IsInsideElementLocPosQUAD
        !  
        !    DESCRIPTION:
        !>   Returns .true. if LocPos (local coordinates) lies inside the 
        !>   area of the quadrilateral element.
        !
        !>   @param[in] LocPos : Local coordinates of the considered point inside the quadrilateral element
        !
        !>   @return IsInsideElementLocPosQUAD : True, if the point lies inside the quadrilateral element
        !
        !**********************************************************************
        logical function IsInsideElementLocPosHEXA(LocPos)
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        
          if ( (LocPos(1) > 1.0) .or. (LocPos(1) < -1.0) .or. (LocPos(2) > 1.0) .or. (LocPos(2) < -1.0) .or. (LocPos(3) > 1.0) .or. (LocPos(3) < -1.0) ) then
            IsInsideElementLocPosHEXA = .false.
          else
            IsInsideElementLocPosHEXA = .true.
          end if
        
        end function IsInsideElementLocPosHEXA

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        !**********************************************************************
        !
        !    FUNCTION: IsInsideElementLocPosQUAD
        !  
        !    DESCRIPTION:
        !>   Returns .true. if LocPos (local coordinates) lies inside the 
        !>   area of the quadrilateral element.
        !
        !>   @param[in] LocPos : Local coordinates of the considered point inside the quadrilateral element
        !
        !>   @return IsInsideElementLocPosQUAD : True, if the point lies inside the quadrilateral element
        !
        !**********************************************************************
        logical function IsInsideElementLocPosQUAD(LocPos)
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        
          if ( (LocPos(1) > 1.0) .or. (LocPos(1) < -1.0) .or. (LocPos(2) > 1.0) .or. (LocPos(2) < -1.0) ) then
            IsInsideElementLocPosQUAD = .false.
          else
            IsInsideElementLocPosQUAD = .true.
          end if
        
        end function IsInsideElementLocPosQUAD


        !**********************************************************************
        !
        !    FUNCTION: IsInsideElementGlobPosQUAD
        !
        !    DESCRIPTION:
        !>   Returns .true. if GlobPos (global coordinates) lies inside the 
        !>   area of the quadrilateral element.
        !
        !>   @param[in] GlobPos : Global coordinates of the considered point inside an element
        !>   @param[in] ElementID : ID of the considered element
        !>   @param[in] NodTot : Total number of nodes
        !>   @param[in] IElTyp : Number of nodes per element
        !>   @param[in] NEl : Total number of elements
        !>   @param[in] NodeCoord : Nodal coordinates
        !>   @param[in] ICon : Element connectivities
        !
        !>   @return IsInsideElementGlobPosQUAD : True, if the point lies inside the element
        !
        !**********************************************************************
        logical function IsInsideElementGlobPosQUAD(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          integer(INTEGER_TYPE), intent(in) :: NodTot, IElTyp, NEl
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          
          ! local variables
          real(REAL_TYPE), dimension(IDim) :: A, B, C, D ! vertices of quadrilateral
          
          A = NodeCoord(ICon(1, ElementID), 1:2)
          B = NodeCoord(ICon(2, ElementID), 1:2)
          C = NodeCoord(ICon(3, ElementID), 1:2)
          D = NodeCoord(ICon(4, ElementID), 1:2)

          IsInsideElementGlobPosQUAD = CheckInsideQuadrilateral(A, B, C, D, GlobPos) 

        end function IsInsideElementGlobPosQUAD


        !**********************************************************************
        !
        !    SUBROUTINE: DetermineCheckEdgeNodesQUAD8
        !
        !    DESCRIPTION:
        !>   Determines which edge node of each side to check in order
        !>   to detect an adjacent element for each side.
        !
        !>   @param[inout] CheckEdgeNodes : Array containing for each side of the element, the local number of the edge node
        !
        !**********************************************************************
        subroutine DetermineCheckEdgeNodesQUAD8(CheckEdgeNodes)
                
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: CheckEdgeNodes ! size(nCheckNodes, NumberOfElementSides)

          CheckEdgeNodes(1,1) = 5 ! Edge node to check for side 1 spanned by nodes 1-2
          CheckEdgeNodes(1,2) = 6 ! Edge node to check for side 2 spanned by nodes 2-3
          CheckEdgeNodes(1,3) = 7 ! Edge node to check for side 3 spanned by nodes 3-4
          CheckEdgeNodes(1,4) = 8 ! Edge node to check for side 3 spanned by nodes 4-1

        end subroutine DetermineCheckEdgeNodesQUAD8
        
      
        !**********************************************************************
        !
        !    FUNCTION: DetermineSideNodesQUAD8
        !
        !    DESCRIPTION:
        !>   Returns the element connectivity of LocalNodeID of side SideID.
        !>     Side 1 is spanned by nodes 1, 2, 5
        !>     Side 2 is spanned by nodes 2, 3, 6
        !>     Side 3 is spanned by nodes 3, 4, 7
        !>     Side 4 is spanned by nodes 4, 1, 8
        !
        !>   @param[in] SideID : ID of the considered side (1 .. 4)
        !>   @param[in] LocalNodeID : ID of the side node  (1 .. 8)
        !
        !>   @return DetermineSideNodesQUAD8 : Local ID of the considered node (1 .. 8)
        !
        !**********************************************************************
        integer function DetermineSideNodesQUAD8(SideID, LocalNodeID)
        
        implicit none

        integer(INTEGER_TYPE), intent(in) :: SideID, LocalNodeID
        
        ! local variables
        integer(INTEGER_TYPE), dimension(3, 4) :: SideConnectivities

        SideConnectivities = reshape( (/  1, 2, 5, &
                                          2, 3, 6, &
                                          3, 4, 7, &
                                          4, 1, 8/), &
                                       (/ 3, 4 /) )

        DetermineSideNodesQUAD8 = SideConnectivities(LocalNodeID, SideID)

      end function DetermineSideNodesQUAD8

      
      subroutine CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !          CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodTot, IElTyp, NEl, NodeCoord, ICon, CrossedSide, IsInside)
      !**********************************************************************
      !
      !    SUBROUTINE: CheckQUADForGlobPos
      !
      !    DESCRIPTION:
      !>   Determines whether GlobPos lies inside the element with ElementID
      !>   (result written to IsInside) and, which side of the quadrilateral is
      !>   crossed by the line between the centrepoint of ElementID and GlobPos
      !>   if GlobPos lies in another element (CrossedSide).
      !>     Side 1 (nodes 1-2 & 5) at xi line
      !>     Side 2 (nodes 2-3 & 6) at eta line
      !>     Side 3 (nodes 3-4 & 7) at xi line
      !>     Side 4 (nodes 4-1 & 8) at eta line
      !
      !>   @param[in] GlobPos : Global coordinates of a point inside the mesh
      !>   @param[in] ElementID : Considered element
      !>   @param[in] CentrePoint : Centrepoint of ElementID
      !>   @param[in] NodTot : Total number of nodes
      !>   @param[in] IElTyp : Number of node connectivities of IElement
      !>   @param[in] NEl : Number of elements
      !>   @param[in] NodeCoord : Global nodal coordinates
      !>   @param[in] ICon : Element connectivities ICon(I, J): global node number of local node I in element J
      !
      !>   @param[out] CrossedSide : Side which contains the intersection point of the above mentioned line
      !>   @param[out] IsInside : True, if GlobPos lies inside ElementID
      !
      !**********************************************************************

      ! DummyCheckForGlobPosPointer(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      
      implicit none
      
        integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as only 2D element
        real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
        integer(INTEGER_TYPE), intent(in) :: ElementID !NodTot, IElTyp, NEl
        real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
        real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
        integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        integer(INTEGER_TYPE), intent(out) :: CrossedSide
        logical, intent(out) :: IsInside
        
        ! local variables
        integer(INTEGER_TYPE), dimension(2, 4) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
        real(REAL_TYPE), dimension(IDim) :: A, B, C, D ! coordinates of nodes of one side of a quadrilateral
        integer(INTEGER_TYPE) :: Success, I, II
        real(REAL_TYPE) :: distanceP_CentrePoint

        ! variables for NURBS
        integer(INTEGER_TYPE), dimension(ELEMENTSIDES) :: KnotConn
        real(REAL_TYPE), dimension(ELEMENTSIDES, NDIM) :: KnotBezierMesh_Elements
        
        KnotConn = 0
        KnotBezierMesh_Elements = 0
        
        !VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
        !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
        !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
        !                         4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
        !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
        
        VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 2, 4,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 4, 3,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 3, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

        
        
        ! generating the knot connectivity for the knot we are currently in 
        do I = 1, ELEMENTSIDES
           KnotConn(I) = KnotBezierMeshConnectivity(ElementID,I) ! 4 by 1
        end do 
        
        ! generating the knot span 'elements'
        do I = 1, ELEMENTSIDES
            do II = 1, NDIM 
                KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
            end do
        end do 
        
        
        !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotBezierMeshConnectivity(ElementID,I) , II) ! 4 by 2
        
        ! KnotBezierMeshConnectivity = it tells you the knot connectivity in a structured mesh
        ! KnotBezierMesh = it tells you the knot coordinates in the physical domain    
        ! KnotBezierMesh_Elements = 
        ! KnotBezierMesh
        
        ! 
        
            
            
            
        
        
        
        do I = 1, ELEMENTSIDES ! loop over sides of quadrilateral
            
            
            !do II = 1, NDIM 
             
            !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
            !end do
            
            
            !VerticesID(1, I)
            
            A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
            
            B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
            
            
            
          !  
          !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
          !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)
          !C = NodeCoord(ICon(VerticesID(1, 3), ElementID), 1:NDIM)
          !D = NodeCoord(ICon(VerticesID(1, 4), ElementID), 1:NDIM)

          ! compute distance P-CentrePoint
          distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
          
          if(distanceP_CentrePoint > SMALL) then
              Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)!CheckInsideSubQuadrilateral(A, B, C, D, GlobPos)
          else
              Success = 2 ! point is exactly on top of the center node 
          endif
          if (Success == 1) then
            CrossedSide = I
            EXIT
          else if (Success == 2) then
            IsInside = .true.
            EXIT
          end if

        end do

      end subroutine CheckQUADForGlobPos

      
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Hexahedron 
      !subroutine CheckHexahedronForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !!          CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodTot, IElTyp, NEl, NodeCoord, ICon, CrossedSide, IsInside)
      !!**********************************************************************
      !!
      !!    SUBROUTINE: CheckQUADForGlobPos
      !!
      !!    DESCRIPTION:
      !!>   Determines whether GlobPos lies inside the element with ElementID
      !!>   (result written to IsInside) and, which side of the quadrilateral is
      !!>   crossed by the line between the centrepoint of ElementID and GlobPos
      !!>   if GlobPos lies in another element (CrossedSide).
      !!>     Side 1 (nodes 1-2 & 5) at xi line
      !!>     Side 2 (nodes 2-3 & 6) at eta line
      !!>     Side 3 (nodes 3-4 & 7) at xi line
      !!>     Side 4 (nodes 4-1 & 8) at eta line
      !!
      !!>   @param[in] GlobPos : Global coordinates of a point inside the mesh
      !!>   @param[in] ElementID : Considered element
      !!>   @param[in] CentrePoint : Centrepoint of ElementID
      !!>   @param[in] NodTot : Total number of nodes
      !!>   @param[in] IElTyp : Number of node connectivities of IElement
      !!>   @param[in] NEl : Number of elements
      !!>   @param[in] NodeCoord : Global nodal coordinates
      !!>   @param[in] ICon : Element connectivities ICon(I, J): global node number of local node I in element J
      !!
      !!>   @param[out] CrossedSide : Side which contains the intersection point of the above mentioned line
      !!>   @param[out] IsInside : True, if GlobPos lies inside ElementID
      !!
      !!**********************************************************************
      !
      !! DummyCheckForGlobPosPointer(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !
      !implicit none
      !
      !  integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 2D element
      !  real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
      !  integer(INTEGER_TYPE), intent(in) :: ElementID !NodTot, IElTyp, NEl
      !  real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
      !  real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
      !  integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
      !  integer(INTEGER_TYPE), intent(out) :: CrossedSide
      !  logical, intent(out) :: IsInside
      !  
      !  ! local variables
      !  integer(INTEGER_TYPE), dimension(3, ELEMENTSIDES) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
      !  !real(REAL_TYPE) :: A, B, C, D, E, F, G, H ! coordinates of nodes of one side of a quadrilateral
      !  integer(INTEGER_TYPE) :: A, B, C, D
      !  real(REAL_TYPE), dimension(NDIM) :: CenterPoint_Surface
      !  integer(INTEGER_TYPE) :: Success, I, II
      !  real(REAL_TYPE) :: distanceP_CentrePoint
      !
      !  ! variables for NURBS
      !  integer(INTEGER_TYPE), dimension(8) :: KnotConn
      !  real(REAL_TYPE), dimension(8, NDIM) :: KnotBezierMesh_Elements
      !  integer(INTEGER_TYPE), dimension(4, 6) :: SideConnectivities
      !  KnotConn = 0
      !  KnotBezierMesh_Elements = 0
      !  
      !  !VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
      !  !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
      !  !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
      !  !                         4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
      !  !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
      !  
      !  SideConnectivities = &
      !        reshape( (/ 1, 2, 4, 3, &! 6, 7, & ! side 1
      !                   2, 6, 8, 4, &! 9, 5, & ! side 2 
      !                   1, 5, 7, 3, &!10, 8, & ! side 3
      !                   5, 6, 8, 7, & ! side 4
      !                   4, 3, 7, 8, & ! side 5
      !                   1, 2, 6, 5/), & ! side 6
      !                  (/ 4, 6 /) )
      !  
      !  
      !  
      !  
      !  
      !  
      !  !! I need to modify this for Hex elements
      !  !VerticesID = reshape( (/ 1, 2, 3, &
      !  !                         1, 3, 5, &
      !  !                         1, 5, 2, &
      !  !                         5, 2, 1, &
      !  !                         !&
      !  !    
      !  !                        5, 6, 8, &
      !  !                        2, 8, 6, &
      !  !                        2, 6, 5, &
      !  !                        5, 2, 6, &
      !  !    
      !  !                        !
      !  !                        !2, 6, 8, &
      !  !                        !6, 8, 5, &
      !  !                        !2, 5, 6, &
      !  !                        !2, 8, 6, &
      !  !                        !!&
      !  !                       3, 4, 2, &
      !  !                       3, 8, 5, &
      !  !                       2, 8, 4, &
      !  !                       2, 3, 4, &
      !  !                      !&
      !  !                      3, 7, 8, &
      !  !                      7, 3, 5, &
      !  !                      7, 8, 5, &
      !  !                      8, 7, 5 /), & 
      !  !                    (/ 3, 12/) ) 
      !  
      !  
      !  
      !  
      !  
      !  
      !  
      !  !reshape( (/ 1, 2, & 
      !  !                2, 6, & 
      !  !                6, 5, & 
      !  !                5, 1, &
      !  !                1, 3, &
      !  !                2, 4, &
      !  !                6, 8, &
      !  !                5, 7, &
      !  !                3, 4, &
      !  !                4, 8, &
      !  !                8, 7, &
      !  !                7, 3 /), &
      !  !                (/ 2, 12 /) ) ! 
      !  
      !  
      !  
      !  
      !  
      !  !reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
      !  !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
      !  !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
      !  !                         4, 5,   &
      !  !                         5, 6,   &
      !  !                         6, 7,   &
      !  !                         7, 8,   &
      !  !                         8, 1/),   &
      !  !                        (/ 2, 8 /))
      !  !                         !6, 5,   &
      !  !                         !5, 7,   &
      !  !                         !7, 8,   &
      !  !                         !8, 6,   &
      !  !                         !6, 2,   &
      !  !                         !2, 4,   &
      !  !                         !4, 8,   &
      !  !                         !8, 6/), &
      !                      
      !  
      !  
      !  !/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
      !  !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
      !
      !  
      !  
      !  ! generating the knot connectivity for the knot we are currently in 
      !  ! 12 sides for hex elements 
      !  do I = 1, 8!ELEMENTSIDES
      !     KnotConn(I) = KnotBezierMeshConnectivity(ElementID,I) ! 4 by 1
      !  end do 
      !  
      !  ! generating the knot span 'elements'
      !  ! 12 sides for Hex elements
      !  do I = 1, 8!ELEMENTSIDES --> vertices of the Bezier mesh 
      !      do II = 1, NDIM 
      !          KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
      !      end do
      !  end do 
      !  
      !  ! Loop over the surfaces of the hexahedron 1 to 6 
      !  
      !  
      !  Success = 0
      !  
      !  do I = 1, 6 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
      !      
      !      
      !
      !      A = KnotConn(SideConnectivities(1, I)) ! coordinate of one side 
      !      B = KnotConn(SideConnectivities(2, I)) ! coordinate of another side 
      !      C = KnotConn(SideConnectivities(3, I)) ! coordinate of another side 
      !      D = KnotConn(SideConnectivities(4, I)) ! coordinate of another side 
      !      !
      !      !
      !      !!A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
      !      !!B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
      !      !!C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
      !      !
      !      !!Success = CheckInsideSubTetrahedron(A, B, C, CentrePoint, GlobPos)
      !      !!
      !      !!Succcess = CheckInsideHexahedronSide (A, B, C, D, GlobPos) 
      !      !
      !      !
      !      ! get center point of surface 
      !      !CenterPoint_Surface(1) = &
      !      !    (KnotBezierMesh(A, 1) + KnotBezierMesh(B, 1) + KnotBezierMesh(C, 1) + KnotBezierMesh(D, 1))/4 
      !      !CenterPoint_Surface(2) = &
      !      !    (KnotBezierMesh(A, 2) + KnotBezierMesh(B, 2) + KnotBezierMesh(C, 2) + KnotBezierMesh(D, 2))/4
      !      !CenterPoint_Surface(3) = &
      !      !    (KnotBezierMesh(A, 3) + KnotBezierMesh(B, 3) + KnotBezierMesh(C, 3) + KnotBezierMesh(D, 3))/4
      !      
      !      
      !      
      !      ! side 1 
      !      
      !      if ((GlobPos(3)<KnotBezierMesh(A, 3)) .and. (I ==1)) then 
      !      CrossedSide = I 
      !      Success=1
      !  else if ((GlobPos(1)>KnotBezierMesh(A, 1)) .and. (I ==2)) then 
      !      CrossedSide = I 
      !      Success=1
      !  else if ((GlobPos(1)<KnotBezierMesh(A, 1)) .and. (I ==3)) then 
      !      CrossedSide = I 
      !      Success=1
      !  else if ((GlobPos(3)>KnotBezierMesh(A, 3)) .and. (I ==4)) then 
      !      CrossedSide = I 
      !      Success=1
      !  else if ((GlobPos(2)>KnotBezierMesh(A, 2)) .and. (I ==5)) then 
      !      CrossedSide = I 
      !      Success=1
      !  else if ((GlobPos(2)<KnotBezierMesh(A, 2)) .and. (I ==6)) then 
      !      CrossedSide = I 
      !      Success=1
      !  end if 
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
      !      if (Success==1) then
      !        CrossedSide = I
      !        EXIT
      !      else if (Success==2) then
      !        IsInside = .true.
      !        EXIT
      !      else 
      !          !do nothing
      !      end if
      !      
      !    end do
      !  
      !  
      !  
      !  
      !  
      !  
      !  
      !  !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotBezierMeshConnectivity(ElementID,I) , II) ! 4 by 2
      !  
      !  ! KnotBezierMeshConnectivity = it tells you the knot connectivity in a structured mesh
      !  ! KnotBezierMesh = it tells you the knot coordinates in the physical domain    
      !  ! KnotBezierMesh_Elements = 
      !  ! KnotBezierMesh
      !  
      !  ! 
      !  
      !      
      !      
      !      
      !  
      !  
      !  
      !  !do I = 1, ELEMENTSIDES ! loop over sides of quadrilateral
      !  !    
      !  !    
      !  !    !do II = 1, NDIM 
      !  !     
      !  !    !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
      !  !    !end do
      !  !    
      !  !    
      !  !    !VerticesID(1, I)
      !  !    
      !  !    A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
      !  !    
      !  !    B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
      !  !    
      !  !    
      !  !    
      !  !  !  
      !  !  !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
      !  !  !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)
      !  !  !C = NodeCoord(ICon(VerticesID(1, 3), ElementID), 1:NDIM)
      !  !  !D = NodeCoord(ICon(VerticesID(1, 4), ElementID), 1:NDIM)
      !  !
      !  !  ! compute distance P-CentrePoint
      !  !  distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
      !  !  
      !  !  if(distanceP_CentrePoint > SMALL) then
      !  !      Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)!CheckInsideSubQuadrilateral(A, B, C, D, GlobPos)
      !  !  else
      !  !      Success = 2 ! point is exactly on top of the center node 
      !  !  endif
      !  !  if (Success == 1) then
      !  !    CrossedSide = I
      !  !    EXIT
      !  !  else if (Success == 2) then
      !  !    IsInside = .true.
      !  !    EXIT
      !  !  end if
      !  !
      !  !end do
      !
      !end subroutine CheckHexahedronForGlobPos
      
      
      
      
      
      subroutine CheckHexahedronForGlobPos_Structured(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !          CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodTot, IElTyp, NEl, NodeCoord, ICon, CrossedSide, IsInside)
      !**********************************************************************
      !
      !    SUBROUTINE: CheckQUADForGlobPos
      !
      !    DESCRIPTION:
      !>   Determines whether GlobPos lies inside the element with ElementID
      !>   (result written to IsInside) and, which side of the quadrilateral is
      !>   crossed by the line between the centrepoint of ElementID and GlobPos
      !>   if GlobPos lies in another element (CrossedSide).
      !>     Side 1 (nodes 1-2 & 5) at xi line
      !>     Side 2 (nodes 2-3 & 6) at eta line
      !>     Side 3 (nodes 3-4 & 7) at xi line
      !>     Side 4 (nodes 4-1 & 8) at eta line
      !
      !>   @param[in] GlobPos : Global coordinates of a point inside the mesh
      !>   @param[in] ElementID : Considered element
      !>   @param[in] CentrePoint : Centrepoint of ElementID
      !>   @param[in] NodTot : Total number of nodes
      !>   @param[in] IElTyp : Number of node connectivities of IElement
      !>   @param[in] NEl : Number of elements
      !>   @param[in] NodeCoord : Global nodal coordinates
      !>   @param[in] ICon : Element connectivities ICon(I, J): global node number of local node I in element J
      !
      !>   @param[out] CrossedSide : Side which contains the intersection point of the above mentioned line
      !>   @param[out] IsInside : True, if GlobPos lies inside ElementID
      !
      !**********************************************************************

      ! DummyCheckForGlobPosPointer(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      
      implicit none
      
        integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 2D element
        real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
        integer(INTEGER_TYPE), intent(in) :: ElementID !NodTot, IElTyp, NEl
        real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
        real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
        integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        integer(INTEGER_TYPE), intent(out) :: CrossedSide
        logical, intent(out) :: IsInside
        
        ! local variables
        integer(INTEGER_TYPE), dimension(3, ELEMENTSIDES) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
        real(REAL_TYPE), dimension(IDim) :: Vec_A, Vec_B, Vec_C, Vec_D, Vec_E, Vec_F, Vec_G, Vec_H ! coordinates of nodes of one side of a quadrilateral
        
        integer(INTEGER_TYPE), dimension(4, 6) :: SideConnectivities
        
        integer(INTEGER_TYPE) :: A, B, C, D
        integer(INTEGER_TYPE) :: Success, I, II
        real(REAL_TYPE) :: distanceP_CentrePoint

        ! variables for NURBS
        integer(INTEGER_TYPE), dimension(8) :: KnotConn
        real(REAL_TYPE), dimension(8, NDIM) :: KnotBezierMesh_Elements
        
        KnotConn = 0
        KnotBezierMesh_Elements = 0
        
        !VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
        !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
        !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
        !                         4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
        !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
        
        ! I need to modify this for Hex elements
        
        SideConnectivities = &
              reshape( (/ 1, 2, 4, 3, &! 6, 7, & ! side 1
                         2, 6, 8, 4, &! 9, 5, & ! side 2 
                         1, 5, 7, 3, &!10, 8, & ! side 3
                         5, 6, 8, 7, & ! side 4
                         4, 3, 7, 8, & ! side 5
                         1, 2, 6, 5/), & ! side 6
                        (/ 4, 6 /) )
        
        
        
        
        
        
        
        
        
        
        VerticesID = reshape( (/ 1, 2, 3, &
                                 1, 3, 5, &
                                 1, 5, 2, &
                                 5, 2, 1, &
                                 !&
            
                                5, 6, 8, &
                                2, 8, 6, &
                                2, 6, 5, &
                                5, 2, 6, &
            
                                !
                                !2, 6, 8, &
                                !6, 8, 5, &
                                !2, 5, 6, &
                                !2, 8, 6, &
                                !!&
                               3, 4, 2, &
                               3, 8, 5, &
                               2, 8, 4, &
                               2, 3, 4, &
                              !&
                              3, 7, 8, &
                              7, 3, 5, &
                              7, 8, 5, &
                              8, 7, 5 /), & 
                            (/ 3, 12/) ) 
        
        
        
        
        
        
        
        !reshape( (/ 1, 2, & 
        !                2, 6, & 
        !                6, 5, & 
        !                5, 1, &
        !                1, 3, &
        !                2, 4, &
        !                6, 8, &
        !                5, 7, &
        !                3, 4, &
        !                4, 8, &
        !                8, 7, &
        !                7, 3 /), &
        !                (/ 2, 12 /) ) ! 
        
        
        
        
        
        !reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
        !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
        !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
        !                         4, 5,   &
        !                         5, 6,   &
        !                         6, 7,   &
        !                         7, 8,   &
        !                         8, 1/),   &
        !                        (/ 2, 8 /))
        !                         !6, 5,   &
        !                         !5, 7,   &
        !                         !7, 8,   &
        !                         !8, 6,   &
        !                         !6, 2,   &
        !                         !2, 4,   &
        !                         !4, 8,   &
        !                         !8, 6/), &
                            
        
        
        !/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
        !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

        
        
        ! generating the knot connectivity for the knot we are currently in 
        ! 12 sides for hex elements 
        do I = 1, 8!ELEMENTSIDES
           KnotConn(I) = KnotBezierMeshConnectivity(ElementID,I) ! 4 by 1
        end do 
        
        ! generating the knot span 'elements'
        ! 12 sides for Hex elements
        do I = 1, 8!ELEMENTSIDES --> vertices of the Bezier mesh 
            do II = 1, NDIM 
                KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
            end do
        end do 
        
        ! final check if it is inside or not the cell 
        
        
        
        
        
        do I = 1, 12 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
            
            

            Vec_A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
            Vec_B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
            Vec_C = KnotBezierMesh_Elements(VerticesID(3, I), 1:NDIM) ! coordinate of another side 


            
            !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
            !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
            !C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
            
            Success = CheckInsideSubTetrahedron(Vec_A, Vec_B, Vec_C, CentrePoint, GlobPos)
            
            if (Success==1) then
              CrossedSide = I
              EXIT
            else if (Success==2) then
              IsInside = .true.
              EXIT
            end if
            
        end do
        
        
        !Check if inside
        do I = 1, 6 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
            
            
      
            A = KnotConn(SideConnectivities(1, I)) ! coordinate of one side 
            B = KnotConn(SideConnectivities(2, I)) ! coordinate of another side 
            C = KnotConn(SideConnectivities(3, I)) ! coordinate of another side 
            D = KnotConn(SideConnectivities(4, I)) ! coordinate of another side 
            !
            !
            !!A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
            !!B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
            !!C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
            !
            !!Success = CheckInsideSubTetrahedron(A, B, C, CentrePoint, GlobPos)
            !!
            !!Succcess = CheckInsideHexahedronSide (A, B, C, D, GlobPos) 
            !
            !
            ! get center point of surface 
            !CenterPoint_Surface(1) = &
            !    (KnotBezierMesh(A, 1) + KnotBezierMesh(B, 1) + KnotBezierMesh(C, 1) + KnotBezierMesh(D, 1))/4 
            !CenterPoint_Surface(2) = &
            !    (KnotBezierMesh(A, 2) + KnotBezierMesh(B, 2) + KnotBezierMesh(C, 2) + KnotBezierMesh(D, 2))/4
            !CenterPoint_Surface(3) = &
            !    (KnotBezierMesh(A, 3) + KnotBezierMesh(B, 3) + KnotBezierMesh(C, 3) + KnotBezierMesh(D, 3))/4
            
            
            
            ! side 1 
            
            if ((GlobPos(3)>KnotBezierMesh(A, 3)) .and. (I ==1)) then 
            !CrossedSide = I 
            Success=2
            IsInside = .true.
        else if ((GlobPos(1)<KnotBezierMesh(A, 1)) .and. (I ==2)) then 
            !CrossedSide = I 
            Success=2
            IsInside = .true.
        else if ((GlobPos(1)>KnotBezierMesh(A, 1)) .and. (I ==3)) then 
            !CrossedSide = I 
            Success=2
            IsInside = .true.
        else if ((GlobPos(3)<KnotBezierMesh(A, 3)) .and. (I ==4)) then 
            !CrossedSide = I 
            Success=2
            IsInside = .true.
        else if ((GlobPos(2)<KnotBezierMesh(A, 2)) .and. (I ==5)) then 
            !CrossedSide = I 
            Success=2
            IsInside = .true.
        else if ((GlobPos(2)>KnotBezierMesh(A, 2)) .and. (I ==6)) then 
            !CrossedSide = I 
            Success=2
            IsInside = .true.
        end if 
        
        end do 
        
        
        
        
        
        
        
        
        !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotBezierMeshConnectivity(ElementID,I) , II) ! 4 by 2
        
        ! KnotBezierMeshConnectivity = it tells you the knot connectivity in a structured mesh
        ! KnotBezierMesh = it tells you the knot coordinates in the physical domain    
        ! KnotBezierMesh_Elements = 
        ! KnotBezierMesh
        
        ! 
        
            
            
            
        
        
        
        !do I = 1, ELEMENTSIDES ! loop over sides of quadrilateral
        !    
        !    
        !    !do II = 1, NDIM 
        !     
        !    !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
        !    !end do
        !    
        !    
        !    !VerticesID(1, I)
        !    
        !    A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
        !    
        !    B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
        !    
        !    
        !    
        !  !  
        !  !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
        !  !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)
        !  !C = NodeCoord(ICon(VerticesID(1, 3), ElementID), 1:NDIM)
        !  !D = NodeCoord(ICon(VerticesID(1, 4), ElementID), 1:NDIM)
        !
        !  ! compute distance P-CentrePoint
        !  distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
        !  
        !  if(distanceP_CentrePoint > SMALL) then
        !      Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)!CheckInsideSubQuadrilateral(A, B, C, D, GlobPos)
        !  else
        !      Success = 2 ! point is exactly on top of the center node 
        !  endif
        !  if (Success == 1) then
        !    CrossedSide = I
        !    EXIT
        !  else if (Success == 2) then
        !    IsInside = .true.
        !    EXIT
        !  end if
        !
        !end do
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

      end subroutine CheckHexahedronForGlobPos_Structured
                                         
      
      
      subroutine CheckHexahedronForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !          CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodTot, IElTyp, NEl, NodeCoord, ICon, CrossedSide, IsInside)
      !**********************************************************************
      !
      !    SUBROUTINE: CheckQUADForGlobPos
      !
      !    DESCRIPTION:
      !>   Determines whether GlobPos lies inside the element with ElementID
      !>   (result written to IsInside) and, which side of the quadrilateral is
      !>   crossed by the line between the centrepoint of ElementID and GlobPos
      !>   if GlobPos lies in another element (CrossedSide).
      !>     Side 1 (nodes 1-2 & 5) at xi line
      !>     Side 2 (nodes 2-3 & 6) at eta line
      !>     Side 3 (nodes 3-4 & 7) at xi line
      !>     Side 4 (nodes 4-1 & 8) at eta line
      !
      !>   @param[in] GlobPos : Global coordinates of a point inside the mesh
      !>   @param[in] ElementID : Considered element
      !>   @param[in] CentrePoint : Centrepoint of ElementID
      !>   @param[in] NodTot : Total number of nodes
      !>   @param[in] IElTyp : Number of node connectivities of IElement
      !>   @param[in] NEl : Number of elements
      !>   @param[in] NodeCoord : Global nodal coordinates
      !>   @param[in] ICon : Element connectivities ICon(I, J): global node number of local node I in element J
      !
      !>   @param[out] CrossedSide : Side which contains the intersection point of the above mentioned line
      !>   @param[out] IsInside : True, if GlobPos lies inside ElementID
      !
      !**********************************************************************

      ! DummyCheckForGlobPosPointer(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      
      implicit none
      
        integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 2D element
        real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
        integer(INTEGER_TYPE), intent(in) :: ElementID !NodTot, IElTyp, NEl
        real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
        real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
        integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        integer(INTEGER_TYPE), intent(out) :: CrossedSide
        logical, intent(out) :: IsInside
        
        ! local variables
        integer(INTEGER_TYPE), dimension(3, ELEMENTSIDES) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
        real(REAL_TYPE), dimension(IDim) :: Vec_A, Vec_B, Vec_C, Vec_D, Vec_E, Vec_F, Vec_G, Vec_H ! coordinates of nodes of one side of a quadrilateral
        
        integer(INTEGER_TYPE), dimension(4, 6) :: SideConnectivities
        
        integer(INTEGER_TYPE) :: A, B, C, D
        integer(INTEGER_TYPE) :: Success, I, II
        real(REAL_TYPE) :: distanceP_CentrePoint

        ! variables for NURBS
        integer(INTEGER_TYPE), dimension(8) :: KnotConn
        real(REAL_TYPE), dimension(8, NDIM) :: KnotBezierMesh_Elements
        
        KnotConn = 0
        KnotBezierMesh_Elements = 0
        
        !VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
        !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
        !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
        !                         4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
        !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
        
        ! I need to modify this for Hex elements
        
        SideConnectivities = &
              reshape( (/ 1, 2, 4, 3, &! 6, 7, & ! side 1
                         2, 6, 8, 4, &! 9, 5, & ! side 2 
                         1, 5, 7, 3, &!10, 8, & ! side 3
                         5, 6, 8, 7, & ! side 4
                         4, 3, 7, 8, & ! side 5
                         1, 2, 6, 5/), & ! side 6
                        (/ 4, 6 /) )
        
        
        
        
        
        
        
        
        
        
        VerticesID = reshape( (/ 1, 2, 3, &
                                 1, 3, 5, &
                                 1, 5, 2, &
                                 5, 2, 1, &
                                 !&
            
                                5, 6, 8, &
                                2, 8, 6, &
                                2, 6, 5, &
                                5, 2, 6, &
            
                                !
                                !2, 6, 8, &
                                !6, 8, 5, &
                                !2, 5, 6, &
                                !2, 8, 6, &
                                !!&
                               3, 4, 2, &
                               3, 8, 5, &
                               2, 8, 4, &
                               2, 3, 4, &
                              !&
                              3, 7, 8, &
                              7, 3, 5, &
                              7, 8, 5, &
                              8, 7, 5 /), & 
                            (/ 3, 12/) ) 
        
        
        
        
        
        
        
        !reshape( (/ 1, 2, & 
        !                2, 6, & 
        !                6, 5, & 
        !                5, 1, &
        !                1, 3, &
        !                2, 4, &
        !                6, 8, &
        !                5, 7, &
        !                3, 4, &
        !                4, 8, &
        !                8, 7, &
        !                7, 3 /), &
        !                (/ 2, 12 /) ) ! 
        
        
        
        
        
        !reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
        !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
        !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
        !                         4, 5,   &
        !                         5, 6,   &
        !                         6, 7,   &
        !                         7, 8,   &
        !                         8, 1/),   &
        !                        (/ 2, 8 /))
        !                         !6, 5,   &
        !                         !5, 7,   &
        !                         !7, 8,   &
        !                         !8, 6,   &
        !                         !6, 2,   &
        !                         !2, 4,   &
        !                         !4, 8,   &
        !                         !8, 6/), &
                            
        
        
        !/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
        !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

        
        
        ! generating the knot connectivity for the knot we are currently in 
        ! 12 sides for hex elements 
        do I = 1, 8!ELEMENTSIDES
           KnotConn(I) = KnotBezierMeshConnectivity(ElementID,I) ! 4 by 1
        end do 
        
        ! generating the knot span 'elements'
        ! 12 sides for Hex elements
        do I = 1, 8!ELEMENTSIDES --> vertices of the Bezier mesh 
            do II = 1, NDIM 
                KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
            end do
        end do 
        
        ! final check if it is inside or not the cell 
        
        
        
        
        
        do I = 1, 12 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
            
            

            Vec_A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
            Vec_B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
            Vec_C = KnotBezierMesh_Elements(VerticesID(3, I), 1:NDIM) ! coordinate of another side 


            
            !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
            !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
            !C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
            
            Success = CheckInsideSubTetrahedron(Vec_A, Vec_B, Vec_C, CentrePoint, GlobPos)
            
            if (Success==1) then
              CrossedSide = I
              EXIT
            else if (Success==2) then
              IsInside = .true.
              EXIT
            end if
            
        end do
        
        
        !Check if inside
        do I = 1, 6 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
            
            
      
            A = KnotConn(SideConnectivities(1, I)) ! coordinate of one side 
            B = KnotConn(SideConnectivities(2, I)) ! coordinate of another side 
            C = KnotConn(SideConnectivities(3, I)) ! coordinate of another side 
            D = KnotConn(SideConnectivities(4, I)) ! coordinate of another side 
            !
            !
            !!A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
            !!B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
            !!C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
            !
            !!Success = CheckInsideSubTetrahedron(A, B, C, CentrePoint, GlobPos)
            !!
            !!Succcess = CheckInsideHexahedronSide (A, B, C, D, GlobPos) 
            !
            !
            ! get center point of surface 
            !CenterPoint_Surface(1) = &
            !    (KnotBezierMesh(A, 1) + KnotBezierMesh(B, 1) + KnotBezierMesh(C, 1) + KnotBezierMesh(D, 1))/4 
            !CenterPoint_Surface(2) = &
            !    (KnotBezierMesh(A, 2) + KnotBezierMesh(B, 2) + KnotBezierMesh(C, 2) + KnotBezierMesh(D, 2))/4
            !CenterPoint_Surface(3) = &
            !    (KnotBezierMesh(A, 3) + KnotBezierMesh(B, 3) + KnotBezierMesh(C, 3) + KnotBezierMesh(D, 3))/4
            
            
            
            ! side 1 
            
        !    if ((GlobPos(3)>KnotBezierMesh(A, 3)) .and. (I ==1)) then 
        !    !CrossedSide = I 
        !    Success=2
        !    IsInside = .true.
        !else if ((GlobPos(1)<KnotBezierMesh(A, 1)) .and. (I ==2)) then 
        !    !CrossedSide = I 
        !    Success=2
        !    IsInside = .true.
        !else if ((GlobPos(1)>KnotBezierMesh(A, 1)) .and. (I ==3)) then 
        !    !CrossedSide = I 
        !    Success=2
        !    IsInside = .true.
        !else if ((GlobPos(3)<KnotBezierMesh(A, 3)) .and. (I ==4)) then 
        !    !CrossedSide = I 
        !    Success=2
        !    IsInside = .true.
        !else if ((GlobPos(2)<KnotBezierMesh(A, 2)) .and. (I ==5)) then 
        !    !CrossedSide = I 
        !    Success=2
        !    IsInside = .true.
        !else if ((GlobPos(2)>KnotBezierMesh(A, 2)) .and. (I ==6)) then 
        !    !CrossedSide = I 
        !    Success=2
        !    IsInside = .true.
        !end if 
            
        !    if (NewLocPos(1) <-1 ) then 
        !    CrossedSide = 4!side 3
        !elseif (NewLocPos(1) >1 ) then
        !    CrossedSide = 2!side 2
        !elseif (NewLocPos(2) <-1 ) then 
        !    CrossedSide = 1!side 6
        !elseif (NewLocPos(2) >1 ) then
        !    CrossedSide = 9!side 5
        !elseif (NewLocPos(3) <-1 ) then
        !    CrossedSide = 1!side 1
        ! elseif (NewLocPos(3) >1 ) then
        !    CrossedSide = 3!side 4
        ! end if 
        
        end do 
        
        
        
        
        
        
        
        
        !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotBezierMeshConnectivity(ElementID,I) , II) ! 4 by 2
        
        ! KnotBezierMeshConnectivity = it tells you the knot connectivity in a structured mesh
        ! KnotBezierMesh = it tells you the knot coordinates in the physical domain    
        ! KnotBezierMesh_Elements = 
        ! KnotBezierMesh
        
        ! 
        
            
            
            
        
        
        
        !do I = 1, ELEMENTSIDES ! loop over sides of quadrilateral
        !    
        !    
        !    !do II = 1, NDIM 
        !     
        !    !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
        !    !end do
        !    
        !    
        !    !VerticesID(1, I)
        !    
        !    A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
        !    
        !    B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
        !    
        !    
        !    
        !  !  
        !  !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
        !  !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)
        !  !C = NodeCoord(ICon(VerticesID(1, 3), ElementID), 1:NDIM)
        !  !D = NodeCoord(ICon(VerticesID(1, 4), ElementID), 1:NDIM)
        !
        !  ! compute distance P-CentrePoint
        !  distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
        !  
        !  if(distanceP_CentrePoint > SMALL) then
        !      Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)!CheckInsideSubQuadrilateral(A, B, C, D, GlobPos)
        !  else
        !      Success = 2 ! point is exactly on top of the center node 
        !  endif
        !  if (Success == 1) then
        !    CrossedSide = I
        !    EXIT
        !  else if (Success == 2) then
        !    IsInside = .true.
        !    EXIT
        !  end if
        !
        !end do
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

      end subroutine CheckHexahedronForGlobPos
      
      !subroutine CheckHexahedronForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside) !CentrePoint, !IsInside,
      !!          CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodTot, IElTyp, NEl, NodeCoord, ICon, CrossedSide, IsInside)
      !!**********************************************************************
      !!
      !!    SUBROUTINE: CheckQUADForGlobPos
      !!
      !!    DESCRIPTION:
      !!>   Determines whether GlobPos lies inside the element with ElementID
      !!>   (result written to IsInside) and, which side of the quadrilateral is
      !!>   crossed by the line between the centrepoint of ElementID and GlobPos
      !!>   if GlobPos lies in another element (CrossedSide).
      !!>     Side 1 (nodes 1-2 & 5) at xi line
      !!>     Side 2 (nodes 2-3 & 6) at eta line
      !!>     Side 3 (nodes 3-4 & 7) at xi line
      !!>     Side 4 (nodes 4-1 & 8) at eta line
      !!
      !!>   @param[in] GlobPos : Global coordinates of a point inside the mesh
      !!>   @param[in] ElementID : Considered element
      !!>   @param[in] CentrePoint : Centrepoint of ElementID
      !!>   @param[in] NodTot : Total number of nodes
      !!>   @param[in] IElTyp : Number of node connectivities of IElement
      !!>   @param[in] NEl : Number of elements
      !!>   @param[in] NodeCoord : Global nodal coordinates
      !!>   @param[in] ICon : Element connectivities ICon(I, J): global node number of local node I in element J
      !!
      !!>   @param[out] CrossedSide : Side which contains the intersection point of the above mentioned line
      !!>   @param[out] IsInside : True, if GlobPos lies inside ElementID
      !!
      !!**********************************************************************
      !
      !! DummyCheckForGlobPosPointer(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
      !
      !implicit none
      !
      !  integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 2D element
      !  real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
      !  integer(INTEGER_TYPE), intent(in) :: ElementID !NodTot, IElTyp, NEl
      !  !real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
      !  real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
      !  integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
      !  integer(INTEGER_TYPE), intent(out) :: CrossedSide
      !  !logical, intent(out) :: IsInside
      !  
      !   real(REAL_TYPE), dimension(:), intent(in) :: NewLocPos
      !  
      !  ! local variables
      !  integer(INTEGER_TYPE), dimension(3, ELEMENTSIDES) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
      !  real(REAL_TYPE), dimension(IDim) :: Vec_A, Vec_B, Vec_C, Vec_D, Vec_E, Vec_F, Vec_G, Vec_H ! coordinates of nodes of one side of a quadrilateral
      !  
      !  integer(INTEGER_TYPE), dimension(4, 6) :: SideConnectivities
      !  
      !  integer(INTEGER_TYPE) :: A, B, C, D
      !  integer(INTEGER_TYPE) :: Success, I, II
      !  real(REAL_TYPE) :: distanceP_CentrePoint
      !
      !  ! variables for NURBS
      !  integer(INTEGER_TYPE), dimension(8) :: KnotConn
      !  real(REAL_TYPE), dimension(8, NDIM) :: KnotBezierMesh_Elements
      !  
      !  !KnotConn = 0
      !  !KnotBezierMesh_Elements = 0
      !  
      !  !VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
      !  !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
      !  !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
      !  !                         4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
      !  !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
      !  
      !  ! I need to modify this for Hex elements
      !  
      !  SideConnectivities = &
      !        reshape( (/ 1, 2, 4, 3, &! 6, 7, & ! side 1
      !                   2, 6, 8, 4, &! 9, 5, & ! side 2 
      !                   1, 5, 7, 3, &!10, 8, & ! side 3
      !                   5, 6, 8, 7, & ! side 4
      !                   4, 3, 7, 8, & ! side 5
      !                   1, 2, 6, 5/), & ! side 6
      !                  (/ 4, 6 /) )
      !  
      !  
      !  
      !  !IsInside = .false.
      !  ! note that crossed side here is crossed EDGE (one of 8) and not side (one of 6)
      !  if (NewLocPos(1) <-1 ) then 
      !      CrossedSide = 4!side 3
      !  elseif (NewLocPos(1) >1 ) then
      !      CrossedSide = 2!side 2
      !  elseif (NewLocPos(2) <-1 ) then 
      !      CrossedSide = 1!side 6
      !  elseif (NewLocPos(2) >1 ) then
      !      CrossedSide = 9!side 5
      !  elseif (NewLocPos(3) <-1 ) then
      !      CrossedSide = 1!side 1
      !   elseif (NewLocPos(3) >1 ) then
      !      CrossedSide = 3!side 4
      !   end if 
      !      
      !      
      !  
      !  
      !  
      !  
      !  
      !  
      !  !VerticesID = reshape( (/ 1, 2, 3, &
      !  !                         1, 3, 5, &
      !  !                         1, 5, 2, &
      !  !                         5, 2, 1, &
      !  !                         !&
      !  !    
      !  !                        5, 6, 8, &
      !  !                        2, 8, 6, &
      !  !                        2, 6, 5, &
      !  !                        5, 2, 6, &
      !  !    
      !  !                        !
      !  !                        !2, 6, 8, &
      !  !                        !6, 8, 5, &
      !  !                        !2, 5, 6, &
      !  !                        !2, 8, 6, &
      !  !                        !!&
      !  !                       3, 4, 2, &
      !  !                       3, 8, 5, &
      !  !                       2, 8, 4, &
      !  !                       2, 3, 4, &
      !  !                      !&
      !  !                      3, 7, 8, &
      !  !                      7, 3, 5, &
      !  !                      7, 8, 5, &
      !  !                      8, 7, 5 /), & 
      !  !                    (/ 3, 12/) ) 
      !  
      !  
      !  
      !  
      !  
      !  
      !  
      !  !reshape( (/ 1, 2, & 
      !  !                2, 6, & 
      !  !                6, 5, & 
      !  !                5, 1, &
      !  !                1, 3, &
      !  !                2, 4, &
      !  !                6, 8, &
      !  !                5, 7, &
      !  !                3, 4, &
      !  !                4, 8, &
      !  !                8, 7, &
      !  !                7, 3 /), &
      !  !                (/ 2, 12 /) ) ! 
      !  
      !  
      !  
      !  
      !  
      !  !reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
      !  !                         2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
      !  !                         3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
      !  !                         4, 5,   &
      !  !                         5, 6,   &
      !  !                         6, 7,   &
      !  !                         7, 8,   &
      !  !                         8, 1/),   &
      !  !                        (/ 2, 8 /))
      !  !                         !6, 5,   &
      !  !                         !5, 7,   &
      !  !                         !7, 8,   &
      !  !                         !8, 6,   &
      !  !                         !6, 2,   &
      !  !                         !2, 4,   &
      !  !                         !4, 8,   &
      !  !                         !8, 6/), &
      !                      
      !  
      !  
      !  !/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
      !  !                      (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.
      !
      !  
      !  
      !  ! generating the knot connectivity for the knot we are currently in 
      !  ! 12 sides for hex elements 
      !  !do I = 1, 8!ELEMENTSIDES
      !  !   KnotConn(I) = KnotBezierMeshConnectivity(ElementID,I) ! 4 by 1
      !  !end do 
      !  !
      !  !! generating the knot span 'elements'
      !  !! 12 sides for Hex elements
      !  !do I = 1, 8!ELEMENTSIDES --> vertices of the Bezier mesh 
      !  !    do II = 1, NDIM 
      !  !        KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
      !  !    end do
      !  !end do 
      !  
      !  ! final check if it is inside or not the cell 
      !  
      !  
      !  
      !  
      !  
      !  !do I = 1, 12 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
      !  !    
      !  !    
      !  !
      !  !    Vec_A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
      !  !    Vec_B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
      !  !    Vec_C = KnotBezierMesh_Elements(VerticesID(3, I), 1:NDIM) ! coordinate of another side 
      !  !
      !  !
      !  !    
      !  !    !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
      !  !    !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
      !  !    !C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
      !  !    
      !  !    Success = CheckInsideSubTetrahedron(Vec_A, Vec_B, Vec_C, CentrePoint, GlobPos)
      !  !    
      !  !    if (Success==1) then
      !  !      CrossedSide = I
      !  !      EXIT
      !  !    else if (Success==2) then
      !  !      IsInside = .true.
      !  !      EXIT
      !  !    end if
      !  !    
      !  !end do
      !  
      !  
      !  !Check if inside
      !  !do I = 1, 6 ! Loop over sides of tetrahedron --> hexahedron has 12 sides 
      !  !    
      !  !    
      !  !
      !  !    A = KnotConn(SideConnectivities(1, I)) ! coordinate of one side 
      !  !    B = KnotConn(SideConnectivities(2, I)) ! coordinate of another side 
      !  !    C = KnotConn(SideConnectivities(3, I)) ! coordinate of another side 
      !  !    D = KnotConn(SideConnectivities(4, I)) ! coordinate of another side 
      !      !
      !      !
      !      !!A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
      !      !!B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
      !      !!C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
      !      !
      !      !!Success = CheckInsideSubTetrahedron(A, B, C, CentrePoint, GlobPos)
      !      !!
      !      !!Succcess = CheckInsideHexahedronSide (A, B, C, D, GlobPos) 
      !      !
      !      !
      !      ! get center point of surface 
      !      !CenterPoint_Surface(1) = &
      !      !    (KnotBezierMesh(A, 1) + KnotBezierMesh(B, 1) + KnotBezierMesh(C, 1) + KnotBezierMesh(D, 1))/4 
      !      !CenterPoint_Surface(2) = &
      !      !    (KnotBezierMesh(A, 2) + KnotBezierMesh(B, 2) + KnotBezierMesh(C, 2) + KnotBezierMesh(D, 2))/4
      !      !CenterPoint_Surface(3) = &
      !      !    (KnotBezierMesh(A, 3) + KnotBezierMesh(B, 3) + KnotBezierMesh(C, 3) + KnotBezierMesh(D, 3))/4
      !      
      !      
      !      
      !      ! side 1 
      !      
      !  !    if ((GlobPos(3)>KnotBezierMesh(A, 3)) .and. (I ==1)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(1)<KnotBezierMesh(A, 1)) .and. (I ==2)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(1)>KnotBezierMesh(A, 1)) .and. (I ==3)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(3)<KnotBezierMesh(A, 3)) .and. (I ==4)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(2)<KnotBezierMesh(A, 2)) .and. (I ==5)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(2)>KnotBezierMesh(A, 2)) .and. (I ==6)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !end if 
      !  
      !  !end do 
      !  
      !  
      !  
      !  
      !  
      !  
      !  
      !  
      !  !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotBezierMeshConnectivity(ElementID,I) , II) ! 4 by 2
      !  
      !  ! KnotBezierMeshConnectivity = it tells you the knot connectivity in a structured mesh
      !  ! KnotBezierMesh = it tells you the knot coordinates in the physical domain    
      !  ! KnotBezierMesh_Elements = 
      !  ! KnotBezierMesh
      !  
      !  ! 
      !  
      !      
      !      
      !      
      !  
      !  
      !  
      !  !do I = 1, ELEMENTSIDES ! loop over sides of quadrilateral
      !  !    
      !  !    
      !  !    !do II = 1, NDIM 
      !  !     
      !  !    !KnotBezierMesh_Elements(I,II) = KnotBezierMesh(KnotConn(I) , II) ! 4 by 2
      !  !    !end do
      !  !    
      !  !    
      !  !    !VerticesID(1, I)
      !  !    
      !  !    A = KnotBezierMesh_Elements(VerticesID(1, I), 1:NDIM) ! coordinate of one side 
      !  !    
      !  !    B = KnotBezierMesh_Elements(VerticesID(2, I), 1:NDIM) ! coordinate of another side 
      !  !    
      !  !    
      !  !    
      !  !  !  
      !  !  !A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
      !  !  !B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)
      !  !  !C = NodeCoord(ICon(VerticesID(1, 3), ElementID), 1:NDIM)
      !  !  !D = NodeCoord(ICon(VerticesID(1, 4), ElementID), 1:NDIM)
      !  !
      !  !  ! compute distance P-CentrePoint
      !  !  distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
      !  !  
      !  !  if(distanceP_CentrePoint > SMALL) then
      !  !      Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)!CheckInsideSubQuadrilateral(A, B, C, D, GlobPos)
      !  !  else
      !  !      Success = 2 ! point is exactly on top of the center node 
      !  !  endif
      !  !  if (Success == 1) then
      !  !    CrossedSide = I
      !  !    EXIT
      !  !  else if (Success == 2) then
      !  !    IsInside = .true.
      !  !    EXIT
      !  !  end if
      !  !
      !  !end do
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
      !  !if ((GlobPos(3)>KnotBezierMesh(A, 3)) .and. (I ==1)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(1)<KnotBezierMesh(A, 1)) .and. (I ==2)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(1)>KnotBezierMesh(A, 1)) .and. (I ==3)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(3)<KnotBezierMesh(A, 3)) .and. (I ==4)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(2)<KnotBezierMesh(A, 2)) .and. (I ==5)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !else if ((GlobPos(2)>KnotBezierMesh(A, 2)) .and. (I ==6)) then 
      !  !    !CrossedSide = I 
      !  !    Success=2
      !  !    IsInside = .true.
      !  !end if 
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
      !  
      !  
      !  
      !  
      !  
      !  
      !  
      !
      !end subroutine CheckHexahedronForGlobPos
      
      
      
      
      end module ModElementEvaluationQUAD