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

	  module ModElementEvaluationTETRA
      !**********************************************************************
      !
      !    Function:  This module provides specific routines for evaluating
      !               tetrahedral elements.
      !
      !    Implemented in the frame of the MPM project.
      !
      !    $Revision: 8842 $
      !    $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************

      use ModGeometryMath
      use ModGlobalConstants

      implicit none

      contains 

        subroutine DetermineSidePlaneDataTETRA(ISide, PlaneNormal, PlanePoint)
        !**********************************************************************
        !
        !    Function: Returns a normal vector of the side with SideID and a point
        !              on that plane. The normal vector is normalised and points inward.
        !> @note : 3D element
        !
        ! I   ISide : Side ID of the element
        !
        ! O   PlaneNormal : Normal to the element side
        ! O   PlanePoint : Point on the element side
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as 3D element
          integer(INTEGER_TYPE), intent(in) :: ISide
          real(REAL_TYPE), dimension(IDim), intent(out) :: PlaneNormal, PlanePoint

          PlaneNormal = 0.0
          PlanePoint = 0.0

          select case(ISide)
            case(1)
              PlaneNormal(2) = 1.0
            case(2)
              PlaneNormal(1) = 1.0
            case(3)
              PlaneNormal(3) = 1.0
            case(4)
              PlaneNormal(1) = -1.0
              PlaneNormal(2) = -1.0
              PlaneNormal(3) = -1.0
              PlaneNormal = VectorNorm(PlaneNormal, IDim)
              PlanePoint(1) = 1.0
          end select
        
        end subroutine DetermineSidePlaneDataTETRA

        
        real(REAL_TYPE) function PointSideDistanceTETRA(SideID, LocPos)
        !**********************************************************************
        !
        !    Function:  Returns the minimum distance between side SideID and LocPos.
        !> @note : 3D element
        !
        !  I   SideID : ID of the considered side (1..4)
        !  I   LocPos : Local coordinates of considered point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as 3D element
          integer(INTEGER_TYPE), intent(in) :: SideID
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          ! Local variables
          real(REAL_TYPE), dimension(IDim) :: PlaneNormal, PlanePoint

          call DetermineSidePlaneDataTETRA(SideID, PlaneNormal, PlanePoint)
        
          PointSideDistanceTETRA = PlanePointDistance(PlaneNormal, PlanePoint, LocPos)
        
        end function PointSideDistanceTETRA

     
        subroutine InitialLocalMaterialPointCoordinatesTETRA(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)
        !**********************************************************************
        !
        !  Function : Determines the local coordinates and integration weight assigned to material point with ID
        !             IParticle which are returned through WeiGP and PosGP.
        !             Currently, all material points are placed at the same local positions.
        !> @note : 3D element
        !
        !  I  IParticle : Number of the material point
        !  I  SolidPointsElement : Number of solid material points per element
        !  I  LiquidPointsElement : Number of liquid material points per element
        !  O  WeiGP : Initial weight assigned to material point IParticle
        !  O  PosGP : Initial local position of material point IParticle
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          ! local variables
          integer(INTEGER_TYPE) :: ID

          ! first the solid material points are determined
          if ((IParticle<=SolidPointsElement).and.(SolidPointsElement>0)) then
            ID = IParticle
            select case(SolidPointsElement)
            case(1) ! 1 material points per element
              call InitialTETRA_MP1(PosGP, WeiGP)
            case(4) ! 4 material points per element
              call InitialTETRA_MP4(ID, PosGP, WeiGP)
            case(7) ! 7 material points per element
              call InitialTETRA_MP7(ID, PosGP, WeiGP)
            case(8) ! 8 material points per element
              call InitialTETRA_MP8(ID, PosGP, WeiGP)  
            case(10) ! 10 material points per element
              call InitialTETRA_MP10(ID, PosGP, WeiGP)
            case(13) ! 13 material points per element
              call InitialTETRA_MP13(ID, PosGP, WeiGP)
            case(20) ! 20 material points per element
              call InitialTETRA_MP20(ID, PosGP, WeiGP)
            end select
          end if
          
          if ((IParticle>SolidPointsElement).and.(LiquidPointsElement>0)) then
            ID = IParticle - SolidPointsElement
            select case(LiquidPointsElement)
            case(1) ! 1 material points per element
              call InitialTETRA_MP1(PosGP, WeiGP)
            case(4) ! 4 material points per element
              call InitialTETRA_MP4(ID, PosGP, WeiGP)
            case(7) ! 7 material points per element
              call InitialTETRA_MP7(ID, PosGP, WeiGP)
            case(8) ! 8 material points per element
              call InitialTETRA_MP8(ID, PosGP, WeiGP)  
            case(10) ! 10 material points per element
              call InitialTETRA_MP10(ID, PosGP, WeiGP)
            case(13) ! 13 material points per element
              call InitialTETRA_MP13(ID, PosGP, WeiGP)
            case(20) ! 20 material points per element
              call InitialTETRA_MP20(ID, PosGP, WeiGP)
            end select
          end if
          
        end subroutine InitialLocalMaterialPointCoordinatesTETRA
      
     
        subroutine InitialTETRA_MP1(PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               1 material point is placed in each element
        !> @note : 3D element
        !
        ! I  IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as 3D element
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: A, B, C

          A = (5 -     sqrt(5d0) ) / 20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) / 20 ! ~ 0.585 410
          C = ((3*A)+B)/4.0
          
          do I = 1, IDim 
            PosGP(I) = C
          end do
          WeiGP = 1.0 / 6.0
                   
        end subroutine InitialTETRA_MP1

        
        subroutine InitialTETRA_MP4(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               4 material points are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !> @note: 3D element
        !
        ! I  IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as 3D element
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A, B

          A = (5 -     sqrt(5d0) ) / 20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) / 20 ! ~ 0.585 410

          select case (IParticle)
            case (1)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 24.0
            case (2)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = B
              WeiGP = 1.0 / 24.0
            case (3)
              PosGP(1) = B
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 24.0
            case (4)
              PosGP(1) = A
              PosGP(2) = B
              PosGP(3) = A
              WeiGP = 1.0 / 24.0
          end select
        
        end subroutine InitialTETRA_MP4
        
        
        subroutine InitialTETRA_MP7(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               7 material points are placed in each element.
        !    Note : 3D element
        !
        ! I IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A, B, C

          A = (5 -     sqrt(5d0) ) / 20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) / 20 ! ~ 0.585 410
          C = (B + A) / 2

          select case (IParticle)
            case (1)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 42.0
            case (2)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = B
              WeiGP = 1.0 / 42.0
            case (3)
              PosGP(1) = B
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 42.0
            case (4)
              PosGP(1) = A
              PosGP(2) = B
              PosGP(3) = A
              WeiGP = 1.0 / 42.0
            case (5)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = C
              WeiGP = 1.0 / 42.0
            case (6)
              PosGP(1) = C
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 42.0
            case (7)
              PosGP(1) = A
              PosGP(2) = C
              PosGP(3) = A
              WeiGP = 1.0 / 42.0
          end select
        
        end subroutine InitialTETRA_MP7
        
        
        subroutine InitialTETRA_MP8(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               8 material points are placed in each element
        !> @note : 3D element
        !
        ! I  IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as 3D element
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP


          select case (IParticle)
            case (1)
              PosGP(1) = 0.125d0
              PosGP(2) = 0.125d0
              PosGP(3) = 0.125d0
              WeiGP = 1.0 / 48.0
            case (2)
              PosGP(1) = 0.125d0
              PosGP(2) = 0.125d0
              PosGP(3) = 0.625d0
              WeiGP = 1.0 / 48.0
            case (3)
              PosGP(1) = 0.625d0
              PosGP(2) = 0.125d0
              PosGP(3) = 0.125d0
              WeiGP = 1.0 / 48.0
            case (4)
              PosGP(1) = 0.125d0
              PosGP(2) = 0.625d0
              PosGP(3) = 0.125d0
              WeiGP = 1.0 / 48.0
            case (5)
              PosGP(1) = 0.25d0
              PosGP(2) = 0.125d0
              PosGP(3) = 0.25d0
              WeiGP = 1.0 / 48.0
            case (6)
              PosGP(1) = 0.25d0
              PosGP(2) = 0.375d0
              PosGP(3) = 0.25d0
              WeiGP = 1.0 / 48.0
            case (7)
              PosGP(1) = 0.125d0
              PosGP(2) = 0.25d0
              PosGP(3) = 0.375d0
              WeiGP = 1.0 / 48.0
            case (8)
              PosGP(1) = 0.375d0
              PosGP(2) = 0.25d0
              PosGP(3) = 0.125d0
              WeiGP = 1.0 / 48.0
          end select
        
        end subroutine InitialTETRA_MP8
        
        
        subroutine InitialTETRA_MP10(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               10 material points are placed in each element.
        !> @note : 3D element
        !
        ! I IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as 3D element
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A, B, C

          A = (5 -     sqrt(5d0) ) / 20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) / 20 ! ~ 0.585 410
          C = (B + A) / 2

          select case (IParticle)
            case (1)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 60.0
            case (2)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = B
              WeiGP = 1.0 / 60.0
            case (3)
              PosGP(1) = B
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 60.0
            case (4)
              PosGP(1) = A
              PosGP(2) = B
              PosGP(3) = A
              WeiGP = 1.0 / 60.0
            case (5)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = C
              WeiGP = 1.0 / 60.0
            case (6)
              PosGP(1) = C
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 60.0
            case (7)
              PosGP(1) = A
              PosGP(2) = C
              PosGP(3) = A
              WeiGP = 1.0 / 60.0
            case (8)
              PosGP(1) = C
              PosGP(2) = A
              PosGP(3) = C
              WeiGP = 1.0 / 60.0
            case (9)
              PosGP(1) = A
              PosGP(2) = C
              PosGP(3) = C
              WeiGP = 1.0 / 60.0
            case (10)
              PosGP(1) = C
              PosGP(2) = C
              PosGP(3) = A
              WeiGP = 1.0 / 60.0
          end select
        
        end subroutine InitialTETRA_MP10
        
        
         subroutine InitialTETRA_MP20(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               10 material points are placed in each element.
        !> @note : 3D element
        !
        ! I  IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A, B, C

          A = (5 -     sqrt(5d0) ) / 20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) / 20 ! ~ 0.585 410
          C = (B + A) / 2

          select case (IParticle)
            case (1)
              PosGP(1) = 0.05
              PosGP(2) = 0.075
              PosGP(3) = 0.075
              WeiGP = 1.0 / 120.0
            case (2)
              PosGP(1) = 0.3
              PosGP(2) = 0.1
              PosGP(3) = 0.05
              WeiGP = 1.0 / 120.0
            case (3)
              PosGP(1) = 0.55
              PosGP(2) = 0.1
              PosGP(3) = 0.075
              WeiGP = 1.0 / 120.0
            case (4)
              PosGP(1) = 0.775
              PosGP(2) = 0.1
              PosGP(3) = 0.05
              WeiGP = 1.0 / 120.0
            case (5)
              PosGP(1) = 0.075
              PosGP(2) = 0.35
              PosGP(3) = 0.05
              WeiGP = 1.0 / 120.0
            case (6)
              PosGP(1) = 0.35
              PosGP(2) = 0.35
              PosGP(3) = 0.075
              WeiGP = 1.0 / 120.0
            case (7)
              PosGP(1) = 0.525
              PosGP(2) = 0.35
              PosGP(3) = 0.05
              WeiGP = 1.0 / 120.0
            case (8)
              PosGP(1) = 0.05
              PosGP(2) = 0.6
              PosGP(3) = 0.075
              WeiGP = 1.0 / 120.0
            case (9)
              PosGP(1) = 0.25
              PosGP(2) = 0.6
              PosGP(3) = 0.1
              WeiGP = 1.0 / 120.0
            case (10)
              PosGP(1) = 0.075
              PosGP(2) = 0.8
              PosGP(3) = 0.1
              WeiGP = 1.0 / 120.0
             case (11)
              PosGP(1) = 0.075
              PosGP(2) = 0.1
              PosGP(3) = 0.325
              WeiGP = 1.0 / 120.0
            case (12)
              PosGP(1) = 0.25
              PosGP(2) = 0.1
              PosGP(3) = 0.35
              WeiGP = 1.0 / 120.0
            case (13)
              PosGP(1) = 0.4
              PosGP(2) = 0.1
              PosGP(3) = 0.325
              WeiGP = 1.0 / 120.0
            case (14)
              PosGP(1) = 0.1
              PosGP(2) = 0.3
              PosGP(3) = 0.35
              WeiGP = 1.0 / 120.0
            case (15)
              PosGP(1) = 0.2
              PosGP(2) = 0.3
              PosGP(3) = 0.35
              WeiGP = 1.0 / 120.0
            case (16)
              PosGP(1) = 0.075
              PosGP(2) = 0.45
              PosGP(3) = 0.325
              WeiGP = 1.0 / 120.0
            case (17)
              PosGP(1) = 0.075
              PosGP(2) = 0.075
              PosGP(3) = 0.525
              WeiGP =1.0 / 120.0
            case (18)
              PosGP(1) = 0.175
              PosGP(2) = 0.075
              PosGP(3) = 0.55
              WeiGP = 1.0 / 120.0
            case (19)
              PosGP(1) = 0.075
              PosGP(2) = 0.175
              PosGP(3) = 0.55
              WeiGP = 1.0 / 120.0
            case (20)
              PosGP(1) = 0.075
              PosGP(2) = 0.075
              PosGP(3) = 0.75
              WeiGP = 1.0 / 120.0
          end select
        
        end subroutine InitialTETRA_MP20
        
        
        subroutine InitialTETRA_MP13(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               13 material points are placed in each element.
        !> @note : 3D element
        !
        ! I  IParticle : Local number of the considered material point inside an element
        !
        ! O PosGP : Returns the initial local coordinates of the material point with local number IParticle
        ! O WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! fixed dimension as only 3D functionality
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(IDim), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A, B, C, D, E

          A = (5 -     sqrt(5d0) ) / 20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) / 20 ! ~ 0.585 410
          C = (B + A) / 2
          D = (3 * A + B) / 4
          E = (A + 3 * B) / 4

          select case (IParticle)
            case (1)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (2)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = B
              WeiGP = 1.0 / 78.0
            case (3)
              PosGP(1) = B
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (4)
              PosGP(1) = A
              PosGP(2) = B
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (5)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = C
              WeiGP = 1.0 / 78.0
            case (6)
              PosGP(1) = C
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (7)
              PosGP(1) = A
              PosGP(2) = C
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (8)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = E
              WeiGP = 1.0 / 78.0
            case (9)
              PosGP(1) = E
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (10)
              PosGP(1) = A
              PosGP(2) = E
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (11)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = D
              WeiGP = 1.0 / 78.0
            case (12)
              PosGP(1) = D
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
            case (13)
              PosGP(1) = A
              PosGP(2) = D
              PosGP(3) = A
              WeiGP = 1.0 / 78.0
          end select
        
        end subroutine InitialTETRA_MP13
        
        
        subroutine GaussTETRA_Q4(IGaussPoint, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the local coordinates and weight for IGaussPoint.
        !               4 Gauss Points are located in each tetrahedral element.
        !> @note : 3D element
        !
        ! I  IGaussPoint : Local number of the considered Gauss Point inside an element
        !
        ! O PosGP : Returns the xi, eta, zeta local coordinates of the Gauss Point with local number IGaussPoint
        ! O WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(IDim), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A, B

          A = (5 -     sqrt(5d0) ) /20 ! ~ 0.138 1966
          B = (5 + 3 * sqrt(5d0) ) /20 ! ~ 0.585 410

          select case (IGaussPoint)
            case (1)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 24.0
            case (2)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = B
              WeiGP = 1.0 / 24.0
            case (3)
              PosGP(1) = B
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 24.0
            case (4)
              PosGP(1) = A
              PosGP(2) = B
              PosGP(3) = A
              WeiGP = 1.0 / 24.0
          end select
        
        end subroutine GaussTETRA_Q4

        
        subroutine GaussTETRA_Q1(IGaussPoint, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the local coordinates and weight for IGaussPoint.
        !               1 Gauss Point is located in each 4-noded tetrahedral element.
        !> @note : 3D element
        !
        ! I  IGaussPoint : Local number of the considered Gauss Point inside an element
        !
        ! O PosGP : Returns the xi, eta, zeta local coordinates of the Gauss Point with local number IGaussPoint
        ! O WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          ! Local variables
          real(REAL_TYPE) :: A
          
          A =  0.25d0

          select case (IGaussPoint)
            case (1)
              PosGP(1) = A
              PosGP(2) = A
              PosGP(3) = A
              WeiGP = 1.0 / 6.0
          end select
        
        end subroutine GaussTETRA_Q1

        
        subroutine InitialiseShapeFunctionsTETRA4(HS, dHS, Wt)
    !**********************************************************************
    !
    !    FUNCTION:  InitialiseShapeFunctionsTETRA4
    ! 
    !    DESCRIPTION:
    !>   To calculate the values of shape functions and their
    !>   derivatives at 1 Gaussian integration points in a
    !>   4-noded 3D tetrahedral elements
    !
    !>   @note: 3D element
    !
    !>    @param[inout] Wt : Local weights for integration
    !>    @param[inout] HS(i,j): Shape function j at integration point i
    !>    @param[inout] dHS(i,j,k) : Derivative of shape function j at integration point i wrt direction k
    !
    !                ^ Eta
    !                |
    !                + 4
    !               /|  \
    !                |     \
    !             /  |
    !                |        \
    !            /   |          \
    !                |             \
    !           /    |1              3
    !                +---------------+----> xi
    !           /   /             /
    !              /          /
    !             /       /
    !          / /     /
    !           /  /
    !          +2
    !         / Zeta
    !**********************************************************************
      implicit none

      real(REAL_TYPE), dimension(:, :), intent(inout) :: HS
      real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS
      real(REAL_TYPE), dimension(:), intent(inout) :: Wt
      ! local variables
      real(REAL_TYPE) :: a, b, xi, eta, zeta
      integer(INTEGER_TYPE) :: int, j

          a = 0.25d0 ! ( 5d0) ) /20 
          b = 0.25d0 ! ( 5d0) ) /20 
          
          do int = 1, 1
              
            Wt(Int) = 1d0 / 6
            
            do j = 1, 4
              hs(int,j) = 0.123456
            end do
            
            xi   = a
            eta  = a
            zeta = a
            
            If (int == 2) Zeta = b
            If (int == 3) Xi   = b
            If (int == 4) Eta  = b

            ! 1..4 corners
            HS(Int,1) = (1-Xi-Eta-Zeta)
            HS(Int,2) = Zeta
            HS(Int,3) = Xi  
            HS(Int,4) = Eta 

            ! dHS(,,1) = d HS / dXi   [1..4 : Corners]
            dHS(Int,1,1) = -1
            dHS(Int,2,1) = 0
            dHS(Int,3,1) = 1  
            dHS(Int,4,1) = 0

            ! dHS(,,2) = d HS / dEta  [1..4 : Corners]
            dHS(Int,1,2) = -1
            dHS(Int,2,2) =  0
            dHS(Int,3,2) =  0
            dHS(Int,4,2) =  1 

            ! dHS(,,3) = d HS / dZeta  [1..4 : Corners]
            dHS(Int,1,3) = -1
            dHS(Int,2,3) = 1 
            dHS(Int,3,3) = 0
            dHS(Int,4,3) = 0
            
          end do
          
        end subroutine InitialiseShapeFunctionsTETRA4
        
        
        subroutine ShapeLocPosTETRA10(LocPos, HS, dHS)
        !**********************************************************************
        !
        !    Function:  To calculate the values of shape functions and their
        !               derivatives at LocPos for 10-noded 3D tetrahedral elements.
        !> @note : 3D element
        !
        ! I    LocPos : Local coordinates of a point inside an element
        !
        ! O   HS(i,j) : Shape function j at integration point i
        ! O   dHS(i,j,k) : Derivative of shape function j at integration point i
        !                  wrt direction k
        !
        !                ^ Eta
        !                |
        !                + 4
        !               /|  \
        !                |     \
        !             /  |
        !                +8      +10
        !            /   |          \
        !                |             \
        !          9+    |1      7        3
        !                +-------+-------+----> xi
        !           /   /             /
        !              /          /
        !             +5      +6
        !          / /     /
        !           /  /
        !          +2
        !         / Zeta
        !
        !**********************************************************************

        implicit none
    
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos ! Local coordinates
          real(REAL_TYPE), dimension(:), intent(out) :: HS ! Value of shape functions at LocPos
          real(REAL_TYPE), dimension(:,:), intent(out) :: dHS ! derivative of shape function i with respect to direction j
          ! Local variables
          real(REAL_TYPE) :: Xi, Eta, Zeta

          Xi = LocPos(1)
          Eta = LocPos(2)
          Zeta = LocPos(3)

          ! 1..4 corners
          HS(1) = (1 -Xi -Eta -Zeta) * (0.5 -Xi -Eta -Zeta) / 0.5
          HS(2) = Zeta * (2 * Zeta -1)
          HS(3) = Xi * (2 * Xi -1)
          HS(4) = Eta * (2 * Eta -1)
          ! 5..10 mids
          HS(5) = 4 * Zeta * (1 -Xi -Eta -Zeta)
          HS(6) = 4 * Xi * Zeta
          HS(7) = 4 * Xi * (1 -Xi -Eta -Zeta)
          HS(8) = 4 * Eta * (1 -Xi -Eta -Zeta)
          HS(9) = 4 * Eta * Zeta
          HS(10) = 4 * Xi * Eta

          ! dHS(,,1) = d HS / dXi
          ! 1..4 corners
          dHS(1, 1) = (-1) * (0.5 -Xi -Eta -Zeta) / 0.5 + (1 -Xi -Eta -Zeta) * (-1) / 0.5
          dHS(2, 1) = 0
          dHS(3, 1) = 1 * (2 * Xi -1) + Xi * (2)
          dHS(4, 1) = 0
          ! 5..10 mids
          dHS(5, 1) = 4 * Zeta * (-1)
          dHS(6, 1) = 4 * Zeta
          dHS(7, 1) = 4 * (1 -Xi -Eta -Zeta) + 4 * Xi * (-1)
          dHS(8, 1) = 4 * Eta * (-1)
          dHS(9, 1) = 0
          dHS(10, 1) = 4 * Eta

          ! dHS(,,2) = d HS / dEta
          ! 1..4 corners
          dHS(1, 2) = (-1) * (0.5 -Xi -Eta -Zeta) / 0.5 + (1 -Xi -Eta -Zeta) * (-1) / 0.5
          dHS(2, 2) = 0
          dHS(3, 2) = 0
          dHS(4, 2) = (2 * Eta -1) + Eta * (2)
          ! 5..10 mids
          dHS(5, 2) = 4 * Zeta * (-1)
          dHS(6, 2) = 0
          dHS(7, 2) = 4 * Xi * (-1)
          dHS(8, 2) = 4 * (1 -Xi -Eta -Zeta) + 4 * Eta * (-1)
          dHS(9, 2) = 4 * Zeta
          dHS(10, 2)= 4 * Xi

          ! dHS(,,3) = d HS / dZeta
          ! 1..4 corners
          dHS(1, 3) = (-1) * (0.5 -Xi -Eta -Zeta) / 0.5 + (1 -Xi -Eta -Zeta) * (-1) / 0.5
          dHS(2, 3) = (2 * Zeta -1) + Zeta * (2)
          dHS(3, 3) = 0
          dHS(4, 3) = 0
          ! 5..10 mids
          dHS(5, 3) = 4 * (1 -Xi -Eta -Zeta) + 4 * Zeta * (-1)
          dHS(6, 3) = 4 * Xi
          dHS(7, 3) = 4 * Xi * (-1)
          dHS(8, 3) = 4 * Eta * (-1)
          dHS(9, 3) = 4 * Eta
          dHS(10, 3) = 0

        end subroutine ShapeLocPosTETRA10
        
        
        subroutine ShapeLocPosTETRA4(LocPos, HS, dHS) 
        !**********************************************************************
        !
        !    Function:  To calculate the values of shape functions and their
        !               derivatives at LocPos for 4-noded 3D tetrahedral elements.
        !> @note : 3D element
        !
        ! I   LocPos : Local coordinates of a point inside an element
        !
        ! O   HS(i,j) : Shape function j at integration point i
        ! O   dHS(i,j,k) : Derivative of shape function j at integration point i
        !                  wrt direction k
        !
        !                ^ Eta
        !                |
        !                + 4
        !               /|  \
        !                |     \
        !             /  |
        !                |        \
        !            /   |          \
        !                |             \
        !           /    |1              3
        !                +---------------+----> xi
        !          /    /             /
        !              /          /
        !         /   /        /
        !         /  /     /
        !           /  /
        !          +2
        !         / Zeta
        !
        !**********************************************************************

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos ! Local coordinates
          real(REAL_TYPE), dimension(:), intent(out) :: HS ! Value of shape functions at LocPos
          ! Value of shape function derivatives at LocPos
          ! dHS(i, j) : derivative of shape function i with respect to direction j
          real(REAL_TYPE), dimension(:, :), intent(out) :: dHS

          ! Local variables
          real(REAL_TYPE) :: Xi, Eta, Zeta

          Xi = LocPos(1)
          Eta = LocPos(2)
          Zeta = LocPos(3)

          HS(1) = (1-Xi-Eta-Zeta)
          HS(2) = Zeta
          HS(3) = Xi  
          HS(4) = Eta 

          !dHS(,,1) = d HS / dXi
          dHS(1,1) = -1
          dHS(2,1) =  0
          dHS(3,1) =  1
          dHS(4,1) =  0

          !dHS(,,2) = d HS / dEta
          dHS(1,2) = -1
          dHS(2,2) =  0
          dHS(3,2) =  0
          dHS(4,2) =  1

          !dHS(,,3) = d HS / dZeta
          dHS(1,3) = -1
          dHS(2,3) =  1
          dHS(3,3) =  0
          dHS(4,3) =  0

        end subroutine ShapeLocPosTETRA4
        

        logical function IsInsideElementLocPosTETRA(LocPos)
        !**********************************************************************
        !
        !    Function:  Returns .true. if LocPos (local coordinates) lies inside the 
        !               volume of the tetrahedral element.
        !> @note : 3D element
        !
        ! I   LocPos : Local coordinates of the considered point inside an element
        !
        ! O   IsInsideElementLocPosTETRA : True, if the point lies inside the element.
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        
          if ( (LocPos(1) < 0.0) .or. (LocPos(2 )< 0.0) .or. (LocPos(3) < 0.0) .or. (LocPos(1) + LocPos(2) + LocPos(3) > 1.00) ) then
            IsInsideElementLocPosTETRA = .false.
          else
            IsInsideElementLocPosTETRA = .true.
          end if
        
        end function IsInsideElementLocPosTETRA

        
        logical function IsInsideElementGlobPosTETRA(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon) 
        !**********************************************************************
        !
        !    Function:  Returns .true. if GlobPos (global coordinates) lies inside the 
        !               volume of the tetrahedral element.
        !    Note : 3D element
        !
        !     GlobPos : Global coordinates of the considered point inside an element
        !     ElementID : ID of the considered element
        !     NodTot : Total number of nodes
        !     IElTyp : Number of nodes per element
        !     NEl : Total number of elements
        !     NodeCoord : Nodal coordinates
        !     ICon : Element connectivities
        !
        ! O   IsInsideElementGlobPosTETRA : True, if the point lies inside the element
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 3 
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          integer(INTEGER_TYPE), intent(in) :: NodTot, IElTyp, NEl
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          ! Local variables
          real(REAL_TYPE), dimension(IDim) :: A, B, C, D ! Vertices of tetrahedron
          
          A = NodeCoord(ICon(1, ElementID), 1:3)
          B = NodeCoord(ICon(2, ElementID), 1:3)
          C = NodeCoord(ICon(3, ElementID), 1:3)
          D = NodeCoord(ICon(4, ElementID), 1:3)
        
          IsInsideElementGlobPosTETRA = CheckInsideTetrahedron(A, B, C, D, GlobPos)
        
        end function IsInsideElementGlobPosTETRA
        

        subroutine DetermineAdjacentParticlesTETRA_MP1(NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles lie next to ISide
        !               (tetrahedral element with initially 1 material point).
        !> @note : 3D element
        !               
        !> @note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage 
        !
        ! I    NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

          ParticleStatus(1) = .true.

        end subroutine DetermineAdjacentParticlesTETRA_MP1
        

        subroutine DetermineAdjacentParticlesTETRA_MP4(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which material points lie next to ISide
        !               (tetrahedral element with initially 4 material points).
        !               Side 1 (nodes 1-2-3 & 5-6-7) in xi-zeta plane
        !               Side 2 (nodes 1-4-2 & 8-9-5) in eta-zeta plane
        !               Side 3 (nodes 1-3-4 & 7-10-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3 & 9-10-6) in 'inclined' plane
        !               Particle 1: (0.138,0.138,0.138) near node 1
        !               Particle 2: (0.138,0.138,0.585) near node 2
        !               Particle 3: (0.585,0.138,0.138) near node 3
        !               Particle 4: (0.138,0.585,0.138) near node 4
        !> @note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage
        !> @note : 3D element
        !
        ! I    ISide : Local number of considered element side
        ! I    NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
     
          select case (ISide)
            case (1) ! Particles 1, 2, 3
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
            case (2) ! Particles 1, 2, 4
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
            case (3) ! Particles 1, 3, 4
              ParticleStatus(1) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
            case (4) ! Particles 2, 3, 4
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
          end select
                                 
        end subroutine DetermineAdjacentParticlesTETRA_MP4
        
        
        subroutine DetermineAdjacentParticlesTETRA_MP7(ISide, NElementParticles, ParticleStatus) 
        !**********************************************************************
        !
        !    Function:  Determines which material points lie next to ISide
        !               (tetrahedral element with initially 7 material points).
        !               Side 1 (nodes 1-2-3 & 5-6-7) in xi-zeta plane
        !               Side 2 (nodes 1-4-2 & 8-9-5) in eta-zeta plane
        !               Side 3 (nodes 1-3-4 & 7-10-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3 & 9-10-6) in 'inclined' plane
        !> @note : ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage 
        !> @note : 3D element
        !
        !     ISide : Local number of considered element side
        !     NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
          
         select case (ISide)
            case (1) ! Particles 1, 2, 3, 5, 6
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(6) = .true.
            case (2) ! Particles 1, 2, 4, 5, 7
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(7) = .true.
            case (3) ! Particles 1, 3, 4, 6, 7
              ParticleStatus(1) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(7) = .true.
            case (4) ! Particles 2, 3, 4
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
          end select

        end subroutine DetermineAdjacentParticlesTETRA_MP7
        
        
        subroutine DetermineAdjacentParticlesTETRA_MP8(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which material points lie next to ISide
        !               (tetrahedral element with initially 8 material points).
        !               Side 1 (nodes 1-2-3-5) in xi-zeta plane
        !               Side 2 (nodes 1-4-2-7) in eta-zeta plane
        !               Side 3 (nodes 1-3-4-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3-6) in 'inclined' plane
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage
        !> @note : 3D element
        !
        ! I    ISide : Local number of considered element side
        ! I    NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
          
         select case (ISide)
            case (1) ! Particles 1, 2, 3, 5
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(5) = .true.
            case (2) ! Particles 1, 2, 4, 7
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(7) = .true.
            case (3) ! Particles 1, 3, 4, 8
              ParticleStatus(1) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(8) = .true.
            case (4) ! Particles 2, 3, 4, 6
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
          end select

        end subroutine DetermineAdjacentParticlesTETRA_MP8
        
        
        subroutine DetermineAdjacentParticlesTETRA_MP10(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which material points lie next to ISide
        !               (tetrahedral element with initially 10 material points).
        !               Side 1 (nodes 1-2-3 & 5-6-7) in xi-zeta plane
        !               Side 2 (nodes 1-4-2 & 8-9-5) in eta-zeta plane
        !               Side 3 (nodes 1-3-4 & 7-10-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3 & 9-10-6) in 'inclined' plane
        !> @note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage
        !> @note : 3D element
        !
        ! I    ISide : Local number of considered element side
        ! I    NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
          
         select case (ISide)
            case (1) ! Particles 1, 2, 3, 5, 6, 8
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(8) = .true.
            case (2) ! Particles 1, 2, 4, 5, 7, 9
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(9) = .true.
            case (3) ! Particles 1, 3, 4, 6, 7, 10
              ParticleStatus(1) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(10) = .true.
            case (4) ! Particles 2, 3, 4, 8, 9, 10
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(10) = .true.
          end select

        end subroutine DetermineAdjacentParticlesTETRA_MP10
        
        
        subroutine DetermineAdjacentParticlesTETRA_MP13(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which material points lie next to ISide
        !               (tetrahedral element with initially 13 material points).
        !               Side 1 (nodes 1-2-3 & 5-6-7) in xi-zeta plane
        !               Side 2 (nodes 1-4-2 & 8-9-5) in eta-zeta plane
        !               Side 3 (nodes 1-3-4 & 7-10-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3 & 9-10-6) in 'inclined' plane
        !> @note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage
        !> @note : 3D element
        !
        ! I    ISide : Local number of considered element side
        ! I    NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
          
         select case (ISide)
            case (1) ! Particles 1, 2, 3, 5, 6, 8, 9, 11, 12
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(11) = .true.
              ParticleStatus(12) = .true.
            case (2) ! Particles 1, 2, 4, 5, 7, 8, 10, 11, 13
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(11) = .true.
              ParticleStatus(13) = .true.
            case (3) ! Particles 1, 3, 4, 6, 7, 9, 10, 12, 13
              ParticleStatus(1) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(12) = .true.
              ParticleStatus(13) = .true.
            case (4) ! Particles 2, 3, 4
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
          end select

        end subroutine DetermineAdjacentParticlesTETRA_MP13

        
        subroutine DetermineAdjacentParticlesTETRA_MP20(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which material points lie next to ISide
        !               (tetrahedral element with initially 13 material points).
        !               Side 1 (nodes 1-2-3 & 5-6-7) in xi-zeta plane
        !               Side 2 (nodes 1-4-2 & 8-9-5) in eta-zeta plane
        !               Side 3 (nodes 1-3-4 & 7-10-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3 & 9-10-6) in 'inclined' plane
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage
        !> @note : 3D element
        !
        ! I    ISide : Local number of considered element side
        ! I    NElementParticles : Initial number of material points per elment
        !
        ! O   ParticleStatus : Set to .true., if the material point lies next to ISide
        !
        !**********************************************************************
     
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
          
         select case (ISide)
            case (1) ! Particles 1, 2, 3, 4, 11, 12, 13, 17, 18, 20
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(11) = .true.
              ParticleStatus(12) = .true.
              ParticleStatus(13) = .true.
              ParticleStatus(17) = .true.
              ParticleStatus(18) = .true.
              ParticleStatus(20) = .true.
            case (2) ! Particles 1, 5, 8, 10, 11, 16, 17, 19, 20
              ParticleStatus(1) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(11) = .true.
              ParticleStatus(16) = .true.
              ParticleStatus(17) = .true.
              ParticleStatus(19) = .true.
              ParticleStatus(20) = .true.
            case (3) ! Particles 1, 2, 3, 4, 5, 6, 7, 8
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(8) = .true.
            case (4) ! Particles 4, 7, 9, 10
              ParticleStatus(4) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(10) = .true.
          end select

        end subroutine DetermineAdjacentParticlesTETRA_MP20
        

        subroutine DetermineCheckEdgeNodesTETRA10(CheckEdgeNodes) 
        !**********************************************************************
        !
        !    Function:  Determines which edge nodes of each side to check in order
        !               to detect an adjacent element for each side.
        !     Note : 3D function
        !
        ! O   CheckEdgeNodes : Array containing for each side of the element, 2 local number of edge nodes
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: CheckEdgeNodes ! size(2, 4), NumberOfEdgeNodes=2, NumberOfSides=4
        
          CheckEdgeNodes = &
            reshape( (/ 5, 6, & 
                        8, 9, & 
                        7, 8, & 
                        9, 6 /), & 
                     (/ 2, 4 /) ) ! Only two edge nodes need to be checked on each of the four sides
        
        end subroutine DetermineCheckEdgeNodesTETRA10
        

        integer(INTEGER_TYPE) function DetermineSideNodesTetrahedronHOE(SideID, LocalNodeID) 
        !**********************************************************************
        !
        !    Function:  Returns the element connectivity of LocalNodeID of side SideID.
        !               Side 1 is spanned by nodes 1, 2, 3, 5, 6, 7
        !               Side 2 is spanned by nodes 1, 4, 2, 8, 9, 5
        !               Side 3 is spanned by nodes 1, 3, 4, 7, 10, 8
        !               Side 4 is spanned by nodes 2, 4, 3, 9, 10, 6
        !     Note : 3D function
        !
        !     SideID : ID of the considered side (1 .. 4)
        !     LocalNodeID : ID of the side node (1 .. 6)
        !
        ! O   DetermineSideNodesTetrahedronHOE : Local ID of the considered node (1 .. 10)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: SideID, LocalNodeID
          ! Local variables
          integer(INTEGER_TYPE), dimension(6, 4) :: SideConnectivities
        
          SideConnectivities = &
            reshape( (/ 1, 2, 3, 5, 6, 7, &
                         1, 4, 2, 8, 9, 5, &
                         1, 3, 4, 7, 10, 8, &
                         2, 4, 3, 9, 10, 6 /), &
                     (/ 6, 4 /) )
        
          DetermineSideNodesTetrahedronHOE = SideConnectivities(LocalNodeID, SideID)
        
        end function DetermineSideNodesTetrahedronHOE
  
        
        subroutine DetermineSideNodeLocPosTetrahedronHOE(INode, Xi, Eta) 
        !**********************************************************************
        !
        !    Function:  Returns the local coordinates of node INode in local coordinate
        !               system of the element side.
        !               Side 1 is spanned by nodes 1, 2, 3, 5, 6, 7
        !               Side 2 is spanned by nodes 1, 4, 2, 8, 9, 5
        !               Side 3 is spanned by nodes 1, 3, 4, 7, 10, 8
        !               Side 4 is spanned by nodes 2, 4, 3, 9, 10, 6
        !     Note : 3D function
        !
        !     INode : Node number in counter-clockwise side node numbering
        !
        ! O   Xi, Eta : Local coordinates of nodes
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: INode
          real(REAL_TYPE), intent(out) :: Xi, Eta
        
          select case(INode)
            case(1)
              Xi = 0.0
              Eta = 0.0
            case(2)
              Xi = 0.5
              Eta = 0.0
            case(3)
              Xi = 1.0
              Eta = 0.0
            case(4)
              Xi = 0.5
              Eta = 0.5
            case(5)
              Xi = 0.0
              Eta = 1.0
            case(6)
              Xi = 0.0
              Eta = 0.5
          end select
        
        end subroutine DetermineSideNodeLocPosTetrahedronHOE

        
        subroutine ShapeXiEtaTetrahedronHOE(Xi, Eta, HS, DHS) 
        !**********************************************************************
        !
        !    Function:  Determines the nodal shape function values HS and derivatives DHS for
        !               (Xi, Eta).
        !     Note : 3D function
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

          real(REAL_TYPE), intent(in) :: Xi, Eta
          real(REAL_TYPE), dimension(6), intent(out) :: HS
          real(REAL_TYPE), dimension(6, 2), intent(out) :: DHS

          ! Shape functions HS
          HS(1) = (1 - Xi - Eta) * (1 - 2 * Xi - 2 * Eta)
          HS(2) = -Xi * (1 - 2 * Xi)
          HS(3) = -Eta * (1 - 2 * Eta)
          HS(4) = 4 * Xi * (1 - Xi - Eta)
          HS(5) = 4 * Xi * Eta
          HS(6) = 4 * Eta * (1 - Xi - Eta)

          ! Shape function derivatives DHS/DXi
          DHS(1, 1) = 4 * (Xi + Eta) - 3
          DHS(2, 1) = 4 * Xi -1
          DHS(3, 1) = 0
          DHS(4, 1) = 4 * (1 - 2 * Xi - Eta)
          DHS(5, 1) = 4 * Eta
          DHS(6, 1) = -4 * Eta

          ! Shape function derivatives DHS/DEta
          DHS(1, 2) = 4 * (Xi + Eta) - 3
          DHS(2, 2) = 0
          DHS(3, 2) = 4 * Eta - 1
          DHS(4, 2) = -4 * Xi
          DHS(5, 2) = 4 * Xi
          DHS(6, 2) = 4 * (1 - Xi - 2*Eta)

        end subroutine ShapeXiEtaTetrahedronHOE

        
        subroutine ShapeXiEtaTetrahedronLOE(Xi, Eta, HS, DHS)
        !**********************************************************************
        !
        !    Function:  Determines the nodal shape function values HS and derivatives DHS for
        !               (Xi, Eta).
        !     Note : 3D function
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

          real(REAL_TYPE), intent(in) :: Xi, Eta
          real(REAL_TYPE), dimension(3), intent(out) :: HS
          real(REAL_TYPE), dimension(3, 2), intent(out) :: DHS

          ! Shape functions HS
          HS(1) = (1 - Xi - Eta)
          HS(2) = Xi
          HS(3) = Eta

          ! Shape function derivatives DHS/DXi
          DHS(1, 1) = -1
          DHS(2, 1) = +1
          DHS(3, 1) = 0

          ! Shape function derivatives DHS/DEta
          DHS(1, 2) = -1
          DHS(2, 2) = 0
          DHS(3, 2) = +1

        end subroutine ShapeXiEtaTetrahedronLOE

        
        subroutine CheckTetrahedronForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside) 
        !**********************************************************************
        !
        !    Function:  Determines whether GlobPos lies inside the element with ElementID
        !               (result written to IsInside) and, which side of the tetrahedron is
        !               crossed by the line between the centrepoint of ElementID and GlobPos
        !               if GlobPos lies in another element (CrossedSide).
        !
        !               Side 1 (nodes 1-2-3 & 5-6-7) in xi-zeta plane
        !               Side 2 (nodes 1-4-2 & 8-9-5) in eta-zeta plane
        !               Side 3 (nodes 1-3-4 & 7-10-8) in xi-eta plane
        !               Side 4 (nodes 2-4-3 & 9-10-6) in 'inclined' plane
        !
        !     GlobPos : Global coordinates of a point inside the mesh
        !     ElementID : Considered element
        !     CentrePoint : Centrepoint of ElementID
        !     NodeCoord : Global nodal coordinates
        !     ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        !
        ! O   CrossedSide : Side which contains the intersection point of the above mentioned line
        ! O   IsInside : True, if GlobPos lies inside ElementID
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 3 
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          integer(INTEGER_TYPE), intent(out) :: CrossedSide
          logical, intent(out) :: IsInside
          ! Local variables
          integer(INTEGER_TYPE), dimension(IDim, 4) :: VerticesID
          real(REAL_TYPE), dimension(IDim) :: A, B, C ! nodes of one side of a tetrahedron
          integer(INTEGER_TYPE) :: Success, I
          
          VerticesID = &
            reshape( (/ 1, 2, 3, & ! Corner nodes to check for side 1 spanned by nodes 1-2-3
                        1, 4, 2, & ! Corner nodes to check for side 2 spanned by nodes 1-4-2
                        1, 3, 4, & ! Corner nodes to check for side 3 spanned by nodes 1-3-4
                        2, 4, 3 /), & ! Corner nodes to check for side 4 spanned by nodes 2-4-3
                     (/ 3, 4 /) )
          
          do I = 1, 4 ! Loop over sides of tetrahedron
            
            A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:3)
            B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:3)
            C = NodeCoord(ICon(VerticesID(3, I), ElementID), 1:3)
            
            Success = CheckInsideSubTetrahedron(A, B, C, CentrePoint, GlobPos)
            
            if (Success==1) then
              CrossedSide = I
              EXIT
            else if (Success==2) then
              IsInside = .true.
              EXIT
            end if
            
          end do

        end subroutine CheckTetrahedronForGlobPos
        

        subroutine DataInterpolationTetrahedronHOE(Xi, Eta, Zeta, NGP, InterpolationFunction) 
        !**********************************************************************
        !
        !    Function: Returns the values of the linear interpolation function
        !              evaluated at (Xi, Eta, Zeta) for 10-noded
        !              tetrahedral element with 4 Gauss points for interpolation
        !              of data from Gauss points to (Xi, Eta, Zeta).
        !     Note : 3D function
        !
        !     Xi, Eta, Zeta : Local coordinates of a point inside the element
        !     NGP : Number of Gauss points per element ( != 4)
        !
        ! O   InterpolationFunction : Values of interpolation function at provided local coordinates
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          real(REAL_TYPE), intent(in) :: Xi, Eta, Zeta
          integer(INTEGER_TYPE), intent(in) :: NGP
          real(REAL_TYPE), dimension(NGP), intent(out) :: InterpolationFunction
          ! Local variables
          real(REAL_TYPE) :: A, B
        
          A = (1 -     sqrt(5d0) ) / 4 
          B = (1 + 3 * sqrt(5d0) ) / 4 
        
          InterpolationFunction(1) = B - sqrt(5d0) * Xi - sqrt(5d0) * Eta - sqrt(5d0) * Zeta
          InterpolationFunction(2) = A + sqrt(5d0) * Zeta
          InterpolationFunction(3) = A + sqrt(5d0) * Xi
          InterpolationFunction(4) = A + sqrt(5d0) * Eta
        
        end subroutine DataInterpolationTetrahedronHOE
        
        
        subroutine DataInterpolationTetrahedronLOE(NGP, InterpolationFunction) 
        !**********************************************************************
        !
        !    Function: Returns the averaging value 1 / NGP as interpolation function.
        !     Note : 3D function
        !
        !     NGP : Number of Gauss points per element
        !
        ! O   InterpolationFunction : Averaging values
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: NGP
          real(REAL_TYPE), dimension(NGP), intent(out) :: InterpolationFunction
          ! Local variables
          integer(INTEGER_TYPE) :: I
          
          do I = 1, NGP
            InterpolationFunction(I)= 1.0 / NGP
          end do
         
        end subroutine DataInterpolationTetrahedronLOE
        

        subroutine InterpolateDataFromGP(IElTyp, NInt1, NComponents, Xi, Eta, Zeta, GPData, MPData) 
        !**********************************************************************
        !
        !    Function: Returns data at point (Xi, Eta, Zeta) inside a
        !              4-noded or 10-noded tetrahedral element interpolated
        !              from NInt1 Gauss points. In case of the 4-noded
        !              element, data is average. In case of the 10-noded
        !              element, linear interpolation is performed.
        !     Note : 3D function
        !
        !     IElTyp : Number of nodes per element
        !     NInt1 : Number of Gauss points per element
        !     NComponents : Number of components of GPData (f.e. 6 in case of stresses)
        !     Xi : Local coordinate of point inside element
        !     Eta : Local coordinate of point inside element
        !     Zeta : Local coordinate of point inside element
        !     GPData : Data at Gauss points
        !
        ! O   MPData : Interpolated data at (Xi, Eta, Zeta)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElTyp, NInt1, NComponents
          real(REAL_TYPE), intent(in) :: Xi, Eta, Zeta
          real(REAL_TYPE), dimension(NInt1, NComponents), intent(in) :: GPData
          real(REAL_TYPE), dimension(NComponents), intent(out) :: MPData
          ! Local variables
          integer(INTEGER_TYPE) :: I, J
          real(REAL_TYPE), dimension(NInt1) :: InterpolationFunction

          ! Determine stress interpolation values for particle
          select case(IElTyp)
            case(10) ! 10-noded tetrahedral element
              call DataInterpolationTetrahedronHOE(Xi, Eta, Zeta, NInt1, InterpolationFunction)   
            case(4) ! 4-noded tetrahedral element
              call DataInterpolationTetrahedronLOE(NInt1, InterpolationFunction)  
          end select 
                   
          ! Interpolate particle data from Gauss point data
          MPData = 0.0
          do I = 1, NInt1 ! Loop over the Gauss points per element
            do J = 1, NComponents ! Loop over the components of each Gauss point
                
              MPData(J) = MPData(J) + GPData(I, J) * InterpolationFunction(I)
            end do
          end do
        
        end subroutine InterpolateDataFromGP
        
        
      end module ModElementEvaluationTETRA