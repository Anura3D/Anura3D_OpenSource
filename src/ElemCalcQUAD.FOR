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
          real(REAL_TYPE), dimension(IDim), intent(out) :: PlanePoint
          
          PlaneNormal = 0.0
          PlanePoint = 0.0
          
          select case(ISide)
            case(1)
              PlaneNormal(2) = 1.0
            case(2)
              PlaneNormal(1) = 1.0
            case(3)
              PlaneNormal(2) = -1.0
            case(4)
              PlaneNormal(1) = -1.0
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
        !>   @param[in] SideID : ID of the considered side (1..3)
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
            case default
              call GiveError("Number of solid material points, " // trim(String(SolidPointsElement)) // &
                              ", is not available for quadrilateral elements! Supported numbers are 1 and 4. Error in [subroutine InitialLocalMaterialPointCoordinatesQUAD()].")
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
              PosGP(2) = 0.0
              WeiGP = 4.0
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
              WeiGP = 4.0
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
        !                         ^ Eta
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
        !
        !**********************************************************************
        subroutine ShapeLocPosQUAD4(LocPos, HS, dHS)

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
          real(REAL_TYPE), dimension(:), intent(out) :: HS
          real(REAL_TYPE), dimension(:, :), intent(out) :: dHS

          ! local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi = LocPos(1)
          Eta = LocPos(2)

          ! HS(i)
          HS(1) = (1.0 - Xi) * (1.0 - Eta) / 4.0
          HS(2) = (1.0 + Xi) * (1.0 - Eta) / 4.0
          HS(3) = (1.0 + Xi) * (1.0 + Eta) / 4.0
          HS(4) = (1.0 - Xi) * (1.0 + Eta) / 4.0

          ! dHS(i,1) = dHS / dXi
          dHS(1,1) =  - (1.0 - Eta) / 4.0
          dHS(2,1) =    (1.0 - Eta) / 4.0
          dHS(3,1) =    (1.0 + Eta) / 4.0
          dHS(4,1) =  - (1.0 + Eta) / 4.0

          ! dHS(i,2) = dHS / dEta
          dHS(1,2) =  - (1.0 - Xi) / 4.0
          dHS(2,2) =  - (1.0 + Xi) / 4.0
          dHS(3,2) =    (1.0 + Xi) / 4.0
          dHS(4,2) =    (1.0 - Xi) / 4.0

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

      
      subroutine CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodTot, IElTyp, NEl, NodeCoord, ICon, CrossedSide, IsInside)
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

      implicit none
      
        integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as only 2D element
        
        real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
        integer(INTEGER_TYPE), intent(in) :: NodTot, IElTyp, NEl, ElementID
        real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
        real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
        integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
        integer(INTEGER_TYPE), intent(out) :: CrossedSide
        logical, intent(out) :: IsInside
        
        ! local variables
        integer(INTEGER_TYPE), dimension(2, 4) :: VerticesID ! size (Number of node on each side, Number of sides of the quadrilateral)
        real(REAL_TYPE), dimension(IDim) :: A, B ! coordinates of nodes of one side of a quadrilateral
        integer(INTEGER_TYPE) :: Success, I
        real(REAL_TYPE) :: distanceP_CentrePoint

        VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                 2, 3,   & ! Corner nodes to check for side 2 spanned by nodes 2-3
                                 3, 4,   & ! Corner nodes to check for side 3 spanned by nodes 3-4
                                 4, 1/), & ! Corner nodes to check for side 4 spanned by nodes 4-1
                              (/ 2, 4 /) ) ! 8 elements which should be reshaped to a 2x4 matrix. N_SIDE is always 4 for a quadrilateral.

        do I = 1, 4 ! loop over sides of quadrilateral
            
          A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
          B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)

          ! compute distance P-CentrePoint
          distanceP_CentrePoint = Distance(GlobPos, CentrePoint, NDIM)
          
          if(distanceP_CentrePoint > SMALL) then
              Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)
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

                                         
      end module ModElementEvaluationQUAD