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
         PosGP(1) = -a
         PosGP(2) = a
         WeiGP = b
       case (4)
         PosGP(1) = a
         PosGP(2) = a
         WeiGP = b
       case default
         call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP4()].")
      end select

      PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

   end subroutine InitialQUAD_MP4
   
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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
          ! local variables
          real(REAL_TYPE) :: a = 0.77459666924 ! = sqrt(3/5) ! --> d
          !real(REAL_TYPE) :: b = 0.0
          real(REAL_TYPE) :: c = 0.88888888888 ! --> c
          real(REAL_TYPE) :: d = 0.55555555555


          
          ! note that these gauss points are the same as 1D
          ! refer to https://en.wikipedia.org/wiki/Gaussian_quadrature
          
          !observe how the sum of all WeiGP is 4 with each GP having a weight of 1
          
          select case (IParticle)
            case (1)
              PosGP(1) = -a 
              PosGP(2) = -a
              WeiGP = d*d 
              
              WeiGP_Xi = d
              WeiGP_Eta = d
            case (2)
              PosGP(1) = 0.0
              PosGP(2) = -a
              WeiGP = d*c
              
              WeiGP_Xi = c
              WeiGP_Eta = d
            case (3)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = d*d
              
              WeiGP_Xi = d
              WeiGP_Eta = d
            case (4)
              PosGP(1) = -a
              PosGP(2) = 0.0
              WeiGP = c*d
              
              WeiGP_Xi = d
              WeiGP_Eta = c
            case (5)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = c*c
              
              WeiGP_Xi = c
              WeiGP_Eta = c
            case (6)
              PosGP(1) = a
              PosGP(2) = 0.0
              WeiGP = c*d
              
              WeiGP_Xi = d
              WeiGP_Eta = c
            case (7)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = d*d
              
              WeiGP_Xi = d
              WeiGP_Eta = d
            case (8)
              PosGP(1) = 0.0
              PosGP(2) = a
              WeiGP = d*c
              
              WeiGP_Xi = c
              WeiGP_Eta = d
            case (9)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = d*d
              
              WeiGP_Xi = d
              WeiGP_Eta = d
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP4()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP9
        
        
        
        
        
        
        
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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
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
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = e*f
              
              WeiGP_Xi = f
              WeiGP_Eta = e
            case (3)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = e*f
              
              WeiGP_Xi = f
              WeiGP_Eta = e
            case (4)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = e*e
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (5)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = e*f
              
              WeiGP_Xi = e
              WeiGP_Eta = f
            case (6)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = f*f
              
              WeiGP_Xi = f
              WeiGP_Eta = f
            case (7)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = f*f
              
              WeiGP_Xi = f
              WeiGP_Eta = f
            case (8)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = e*f
              
              WeiGP_Xi = e
              WeiGP_Eta = f
            case (9)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = e*f        
              
              WeiGP_Xi = e
              WeiGP_Eta = f
            case (10)
              PosGP(1) = -b 
              PosGP(2) = b
              WeiGP = f*f 
              
              WeiGP_Xi = f
              WeiGP_Eta = f
            case (11)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = f*f
              
              WeiGP_Xi = f
              WeiGP_Eta = f
            case (12)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = e*f
              
              WeiGP_Xi = e
              WeiGP_Eta = f
            case (13)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = e*e
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (14)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = e*f
              
              WeiGP_Xi = f
              WeiGP_Eta = e
            case (15)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = e*f
              
              WeiGP_Xi = f
              WeiGP_Eta = e
            case (16)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = e*e

              WeiGP_Xi = e
              WeiGP_Eta = e
              
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP16()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP16
        
        
        
        
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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
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
              
              WeiGP_Xi = c
              WeiGP_Eta = c
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = e*c
              
              WeiGP_Xi = e
              WeiGP_Eta = c
            case (3)
              PosGP(1) = 0.0
              PosGP(2) = -a
              WeiGP = d*c
              
              WeiGP_Xi = d
              WeiGP_Eta = c
            case (4)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = c*e
              
              WeiGP_Xi = e
              WeiGP_Eta = c
            case (5)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = c*c
              
              WeiGP_Xi = c
              WeiGP_Eta = c
            case (6)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = c*e
              
              WeiGP_Xi = c
              WeiGP_Eta = e
            case (7)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = e*e
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (8)
              PosGP(1) = 0
              PosGP(2) = -b
              WeiGP = d*e
              
              WeiGP_Xi = d
              WeiGP_Eta = e
            case (9)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = e*e    
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (10)
              PosGP(1) = a 
              PosGP(2) = -b
              WeiGP = c*e 
              
              WeiGP_Xi = c
              WeiGP_Eta = e
            case (11)
              PosGP(1) = -a
              PosGP(2) = 0
              WeiGP = c*d
              
              WeiGP_Xi = c
              WeiGP_Eta = d
            case (12)
              PosGP(1) = -b
              PosGP(2) = 0
              WeiGP = d*e
              
              WeiGP_Xi = e
              WeiGP_Eta = d
            case (13)
              PosGP(1) = 0
              PosGP(2) = 0
              WeiGP = d*d
              
              WeiGP_Xi = d
              WeiGP_Eta = d
            case (14)
              PosGP(1) = b
              PosGP(2) = 0
              WeiGP = e*d
              
              WeiGP_Xi = e
              WeiGP_Eta = d
            case (15)
              PosGP(1) = a
              PosGP(2) = 0
              WeiGP = c*d
              
              WeiGP_Xi = c
              WeiGP_Eta = d
            case (16)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = c*e
              
              WeiGP_Xi = c
              WeiGP_Eta = e
              
            case (17)
              PosGP(1) = -b 
              PosGP(2) = b
              WeiGP = e*e
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (18)
              PosGP(1) = 0.0
              PosGP(2) = b
              WeiGP = d*e
              
              WeiGP_Xi = d
              WeiGP_Eta = e
            case (19)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = e*e
              
              WeiGP_Xi = e
              WeiGP_Eta = e
            case (20)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = c*e
              
              WeiGP_Xi = c
              WeiGP_Eta = e
            case (21)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = c*c
              
              WeiGP_Xi = c
              WeiGP_Eta = c
            case (22)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = e*c
              
              WeiGP_Xi = e
              WeiGP_Eta = c
            case (23)
              PosGP(1) = 0
              PosGP(2) = a
              WeiGP = d*c
              
              WeiGP_Xi = d
              WeiGP_Eta = c
            case (24)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = c*e
              
              WeiGP_Xi = e
              WeiGP_Eta = c
            case (25)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = c*c              
              
              WeiGP_Xi = c
              WeiGP_Eta = c
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP25
        
        
        
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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
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
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
                          
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (4)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (5)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (6)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case (7)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
            case (8)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (9)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight        
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (10)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (11)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (12)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
            case (13)
              PosGP(1) = -a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (14)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (15)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (16)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (17)
              PosGP(1) = b 
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (18)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (19)
              PosGP(1) = -a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (20)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (21)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (22)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (23)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (24)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (25)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = a_weight*b_weight    
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              case (26)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
              case (27)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
              case (28)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
              case (29)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
              case (30)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              case (31)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
              case (32)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              case (33)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
              case (34)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
              case (35)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              case (36)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
              WeiGP_Xi = a_weight  
              WeiGP_Eta = a_weight  
              
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP36
        
        
        
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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
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
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight 
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (4)
              PosGP(1) = d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
            case (5)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (6)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (7)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
              
              
              
              
            case (8)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
            case (9)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (10)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight    
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (11)
              PosGP(1) = d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight 
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (12)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (13)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (14)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              
              
            case (15)
              PosGP(1) = -a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (16)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (17)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (18)
              PosGP(1) = d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (19)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (20)
              PosGP(1) = b 
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (21)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
              
              
              
            case (22)
              PosGP(1) = -a 
              PosGP(2) = d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            case (23)
              PosGP(1) = -b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (24)
              PosGP(1) = -c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (25)
              PosGP(1) = d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (26)
              PosGP(1) = c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (27)
              PosGP(1) = b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (28)
              PosGP(1) = a
              PosGP(2) = d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            
              
              
              
              
            case (29)
              PosGP(1) = -a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (30)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (31)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (32)
              PosGP(1) = d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (33)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (34)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (35)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
              
              
            case (36)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = a_weight*b_weight    
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              case (37)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
              case (38)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight 
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
              case (39)
              PosGP(1) = d
              PosGP(2) = b
              WeiGP = d_weight*b_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
              case (40)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
              case (41)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
              case (42)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              
              
              
              
              
              
              case (43)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
              case (44)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              case (45)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
              case (46)
              PosGP(1) = d
              PosGP(2) = a
              WeiGP = d_weight*a_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
              case (47)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
              
              
              case (48)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              case (49)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP49
        
        
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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
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
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (4)
              PosGP(1) = -d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight 
            case (5)
              PosGP(1) = d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
            case (6)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (7)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (8)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
              
              
            case (9)
              PosGP(1) = -a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
            case (10)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (11)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight      
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (12)
              PosGP(1) = -d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight      
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (13)
              PosGP(1) = d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (14)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (15)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (16)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              
            case (17)
              PosGP(1) = -a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (18)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (19)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (20)
              PosGP(1) = -d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (21)
              PosGP(1) = d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (22)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (23)
              PosGP(1) = b 
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (24)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
              
              
              
              case (25)
              PosGP(1) = -a 
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            case (26)
              PosGP(1) = -b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (27)
              PosGP(1) = -c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (28)
              PosGP(1) = -d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (29)
              PosGP(1) = d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (30)
              PosGP(1) = c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (31)
              PosGP(1) = b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (32)
              PosGP(1) = a
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
              
              
            case (33)
              PosGP(1) = -a 
              PosGP(2) = d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            case (34)
              PosGP(1) = -b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (35)
              PosGP(1) = -c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (36)
              PosGP(1) = -d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (37)
              PosGP(1) = d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (38)
              PosGP(1) = c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (39)
              PosGP(1) = b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight 
            case (40)
              PosGP(1) = a
              PosGP(2) = d
              WeiGP = a_weight*d_weight
            
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
              
              
              
            case (41)
              PosGP(1) = -a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight 
            case (42)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (43)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (44)
              PosGP(1) = -d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (45)
              PosGP(1) = d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (46)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (47)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (48)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
              
              
            case (49)
              PosGP(1) = -a
              PosGP(2) = b
              WeiGP = a_weight*b_weight    
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              case (50)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
              case (51)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight 
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
              case (52)
              PosGP(1) = -d
              PosGP(2) = b
              WeiGP = d_weight*b_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
              case (53)
              PosGP(1) = d
              PosGP(2) = b
              WeiGP = d_weight*b_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
              case (54)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
              case (55)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
              case (56)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              
              
              
              
              case (57)
              PosGP(1) = -a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
              case (58)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              case (59)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
              case (60)
              PosGP(1) = -d
              PosGP(2) = a
              WeiGP = d_weight*a_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
              case (61)
              PosGP(1) = d
              PosGP(2) = a
              WeiGP = d_weight*a_weight  
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
              case (62)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight  
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
              case (63)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight  
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              case (64)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight  
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select

          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )

        end subroutine InitialQUAD_MP64
        

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
          
          real(REAL_TYPE) :: WeiGP_Xi
          real(REAL_TYPE) :: WeiGP_Eta
          
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
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case (2)
              PosGP(1) = -b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (3)
              PosGP(1) = -c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (4)
              PosGP(1) = -d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
            case (5)
              PosGP(1) = 0.0
              PosGP(2) = -a
              WeiGP = 0.3302393550012598*a_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = a_weight
            case (6)
              PosGP(1) = d
              PosGP(2) = -a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
            case (7)
              PosGP(1) = c
              PosGP(2) = -a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (8)
              PosGP(1) = b
              PosGP(2) = -a
              WeiGP = b_weight*a_weight
            
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
              
            case (9)
              PosGP(1) = a
              PosGP(2) = -a
              WeiGP = a_weight*a_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
              
              
              ! row 2
            case (10)
              PosGP(1) = -a 
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
            case (11)
              PosGP(1) = -b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (12)
              PosGP(1) = -c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (13)
              PosGP(1) = -d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (14)
              PosGP(1) = 0.0
              PosGP(2) = -b
              WeiGP = 0.3302393550012598*b_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = b_weight
            case (15)
              PosGP(1) = d
              PosGP(2) = -b
              WeiGP = d_weight*b_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (16)
              PosGP(1) = c
              PosGP(2) = -b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (17)
              PosGP(1) = b
              PosGP(2) = -b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (18)
              PosGP(1) = a
              PosGP(2) = -b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              
              
              ! row 3
            case (19)
              PosGP(1) = -a 
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (20)
              PosGP(1) = -b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (21)
              PosGP(1) = -c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (22)
              PosGP(1) = -d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (23)
              PosGP(1) = 0.0
              PosGP(2) = -c
              WeiGP = 0.3302393550012598*c_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = c_weight
            case (24)
              PosGP(1) = d
              PosGP(2) = -c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (25)
              PosGP(1) = c
              PosGP(2) = -c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (26)
              PosGP(1) = b
              PosGP(2) = -c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (27)
              PosGP(1) = a
              PosGP(2) = -c
              WeiGP = a_weight*c_weight
              
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
              
              
              ! row 4
            case (28)
              PosGP(1) = -a 
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            case (29)
              PosGP(1) = -b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (30)
              PosGP(1) = -c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (31)
              PosGP(1) = -d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (32)
              PosGP(1) = 0.0
              PosGP(2) = -d
              WeiGP = 0.3302393550012598*d_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = d_weight
            case (33)
              PosGP(1) = d
              PosGP(2) = -d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (34)
              PosGP(1) = c
              PosGP(2) = -d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (35)
              PosGP(1) = b
              PosGP(2) = -d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (36)
              PosGP(1) = a
              PosGP(2) = -d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            
              
              
              
              ! row 5
            case (37)
              PosGP(1) = -a 
              PosGP(2) = 0.0
              WeiGP = a_weight*0.3302393550012598
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = 0.3302393550012598
            case (38)
              PosGP(1) = -b
              PosGP(2) = 0.0
              WeiGP = b_weight*0.3302393550012598
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = 0.3302393550012598
            case (39)
              PosGP(1) = -c
              PosGP(2) = 0.0
              WeiGP = c_weight*0.3302393550012598
              
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = 0.3302393550012598
            case (40)
              PosGP(1) = -d
              PosGP(2) = 0.0
              WeiGP = d_weight*0.3302393550012598
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = 0.3302393550012598
            case (41)
              PosGP(1) = 0.0
              PosGP(2) = 0.0
              WeiGP = 0.3302393550012598*0.3302393550012598
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = 0.3302393550012598
            case (42)
              PosGP(1) = d
              PosGP(2) = 0.0
              WeiGP = d_weight*0.3302393550012598
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = 0.3302393550012598
            case (43)
              PosGP(1) = c
              PosGP(2) = 0.0
              WeiGP = c_weight*0.3302393550012598
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = 0.3302393550012598
            case (44)
              PosGP(1) = b
              PosGP(2) = 0.0
              WeiGP = b_weight*0.3302393550012598
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = 0.3302393550012598
            case (45)
              PosGP(1) = a
              PosGP(2) = 0.0
              WeiGP = a_weight*0.3302393550012598
              WeiGP_Xi = a_weight
              WeiGP_Eta = 0.3302393550012598
              
              
              ! row 6
            case (46)
              PosGP(1) = -a 
              PosGP(2) = d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
            case (47)
              PosGP(1) = -b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (48)
              PosGP(1) = -c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (49)
              PosGP(1) = -d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (50)
              PosGP(1) = 0.0
              PosGP(2) = d
              WeiGP = 0.3302393550012598*d_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = d_weight
            case (51)
              PosGP(1) = d
              PosGP(2) = d
              WeiGP = d_weight*d_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = d_weight
            case (52)
              PosGP(1) = c
              PosGP(2) = d
              WeiGP = c_weight*d_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = d_weight
            case (53)
              PosGP(1) = b
              PosGP(2) = d
              WeiGP = b_weight*d_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = d_weight
            case (54)
              PosGP(1) = a
              PosGP(2) = d
              WeiGP = a_weight*d_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = d_weight
              
              
              
              
              ! row 7
            case (55)
              PosGP(1) = -a 
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
            case (56)
              PosGP(1) = -b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (57)
              PosGP(1) = -c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (58)
              PosGP(1) = -d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (59)
              PosGP(1) = 0.0
              PosGP(2) = c
              WeiGP = 0.3302393550012598*c_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = c_weight
            case (60)
              PosGP(1) = d
              PosGP(2) = c
              WeiGP = d_weight*c_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = c_weight
            case (61)
              PosGP(1) = c
              PosGP(2) = c
              WeiGP = c_weight*c_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = c_weight
            case (62)
              PosGP(1) = b
              PosGP(2) = c
              WeiGP = b_weight*c_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = c_weight
            case (63)
              PosGP(1) = a
              PosGP(2) = c
              WeiGP = a_weight*c_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = c_weight
              
              
              ! row 8
            case (64)
              PosGP(1) = -a 
              PosGP(2) = b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
            case (65)
              PosGP(1) = -b
              PosGP(2) = b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (66)
              PosGP(1) = -c
              PosGP(2) = b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (67)
              PosGP(1) = -d
              PosGP(2) = b
              WeiGP = d_weight*b_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (68)
              PosGP(1) = 0.0
              PosGP(2) = b
              WeiGP = 0.3302393550012598*b_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = b_weight
            case (69)
              PosGP(1) = d
              PosGP(2) = b
              WeiGP = d_weight*b_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = b_weight
            case (70)
              PosGP(1) = c
              PosGP(2) = b
              WeiGP = c_weight*b_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = b_weight
            case (71)
              PosGP(1) = b
              PosGP(2) = b
              WeiGP = b_weight*b_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = b_weight
            case (72)
              PosGP(1) = a
              PosGP(2) = b
              WeiGP = a_weight*b_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = b_weight
              
              
              
              ! row 9
            case (73)
              PosGP(1) = -a 
              PosGP(2) = a
              WeiGP = a_weight*a_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case (74)
              PosGP(1) = -b
              PosGP(2) = a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (75)
              PosGP(1) = -c
              PosGP(2) = a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (76)
              PosGP(1) = -d
              PosGP(2) = a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
            case (77)
              PosGP(1) = 0.0
              PosGP(2) = a
              WeiGP = 0.3302393550012598*a_weight
              
              WeiGP_Xi = 0.3302393550012598
              WeiGP_Eta = a_weight
            case (78)
              PosGP(1) = d
              PosGP(2) = a
              WeiGP = d_weight*a_weight
              
              WeiGP_Xi = d_weight
              WeiGP_Eta = a_weight
            case (79)
              PosGP(1) = c
              PosGP(2) = a
              WeiGP = c_weight*a_weight
              
              WeiGP_Xi = c_weight
              WeiGP_Eta = a_weight
            case (80)
              PosGP(1) = b
              PosGP(2) = a
              WeiGP = b_weight*a_weight
              
              WeiGP_Xi = b_weight
              WeiGP_Eta = a_weight
            case (81)
              PosGP(1) = a
              PosGP(2) = a
              WeiGP = a_weight*a_weight
              
              WeiGP_Xi = a_weight
              WeiGP_Eta = a_weight
            case default
              call GiveError("Undefined number of material points in [subroutine InitialQUAD_MP25()].")      
          end select
        
          PosGP = PosGP * ( 1.0 - CalParams%shrinkageMateriaPointPositionFactor )
        
        end subroutine InitialQUAD_MP81

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
         PosGP(1) = -a
         PosGP(2) = a
         WeiGP = b
       case (4)
         PosGP(1) = a
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
      dHS(1,1) =- (1.0 - Eta) / 4.0
      dHS(2,1) =  (1.0 - Eta) / 4.0
      dHS(3,1) =  (1.0 + Eta) / 4.0
      dHS(4,1) =- (1.0 + Eta) / 4.0

      ! dHS(i,2) = dHS / dEta
      dHS(1,2) =- (1.0 - Xi) / 4.0
      dHS(2,2) =- (1.0 + Xi) / 4.0
      dHS(3,2) =  (1.0 + Xi) / 4.0
      dHS(4,2) =  (1.0 - Xi) / 4.0

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
   subroutine CheckQUADForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
   

      implicit none

      integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as only 2D element

      real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
      integer(INTEGER_TYPE), intent(in) :: ElementID
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
   subroutine InitialiseShapeFunctionsQUAD4(HS, dHS, Wt)

      implicit none

      !real(REAL_TYPE), dimension(:), intent(inout) :: LocPos
      real(REAL_TYPE), dimension(:,:), intent(inout) :: HS
      real(REAL_TYPE), dimension(:,:,:), intent(inout) :: dHS
      real(REAL_TYPE), dimension(:), intent(inout) :: Wt

      ! local variables
      real(REAL_TYPE) :: Xi, Eta
      integer(INTEGER_TYPE) :: int, I1, Nint1
      real(REAL_TYPE), dimension(NDIM) :: LocPosGP
      
      !allocate(LocPosGP(NDIM))

      ! Note this is two dimensional

      Nint1=ELEMENTGAUSSPOINTS !number of gauss points
      Int = 0 !counter

      do I1 = 1, Nint1

          
          ! need to get the local position of the gauss points
          call InitialLocalMaterialPointCoordinatesPointer(I1, Nint1, 0, Wt(I1), LocPosGP)
          
         Xi = LocPosGP(1)!0.0 !local position in Xi (local) direction
         Eta = LocPosGP(2)!0.0 !local position in Eta (local) direction

         Int = Int+1

         !Wt(Int) = 2.0 !1d0 / Nint1 * 0.5 !This should be =2... double check!!!!

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
   !    FUNCTION:  GetMinAltitudeQUAD
   !
   !    DESCRIPTION:
   !>   Determines the minimum altitude of a quadilateral element...
   !    I would not call it an altitude as it is practically the smallest side of the rectangle.
   !    WARNING: Only works for 4 noded elements
   !>   @param[not defined] NodeNr : node number - See fixme below for info on why intent is not defined
   !    @param[in] NodeCoord       : Global nodal coordinates
   !    @param[out] Lmin           : minimum altitude of the element
   !
   !**********************************************************************
   subroutine GetMinAltitudeQUAD(NodeNr, NodeCoord, Lmin)
      ! idea is to input all nodes of a rectangle in the order within which they are supposed to be
      ! e.g.
      ! 1) node 1, node 2, node 3, node 4 in an anticlockwise fashion
      ! 2) form lines between node 1 and node 2
      ! 3) find the magnitude of these lines
      ! 4) repeat this process with nodes 2-3, nodes 3-4
      ! 5) select the minimum out of all of these values.
      implicit none
      real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
      real(REAL_TYPE), intent (inout):: Lmin

      !FIXME: For an unknown reason the pointer GetminAltitudePoint which links to DummyGetMinAltitudePointer doesn't set an intent
      !for NodeNr. The pointer and the referenced function have to have matching inputs so NodeNr can't be set with an intent in this
      !subroutine. This isn't correct, you should be able to set an intent for NodeNr. This will probably require changing the other
      !element types
      integer(INTEGER_TYPE), dimension(:) :: NodeNr


      ! local variables
      integer(INTEGER_TYPE):: i, second_node, num_nodes
      real(REAL_TYPE), allocatable :: d(:)

      select case(ELEMENTTYPE)

       case(QUAD4) ! quadrilateral_4-noded
         num_nodes = 4
         allocate(d(num_nodes))

         ! Loop over the nodes
         do i= 1, num_nodes
            ! Use the mod to cycle around the array- needed for the last distance calc of node 4 - node 1
            second_node = mod(i, num_nodes) +1

            ! Calc the distance between each node
            ! The distance formula for a vector a and b is equivalent to the 2-norm of (a-b) ie. norm2(a-b)
            d(i) = norm2(NodeCoord(NodeNr(i), :) - NodeCoord(NodeNr(second_node), :))
         end do

       case default
         call GiveError("In GetMinAltitudeQUAD, calculating the length of quadrilateral elements is only implemented for 4-noded elements")
      end select
      ! Select and return the minimum element width
      Lmin = minval(d)

   end subroutine GetMinAltitudeQUAD


   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP1(NElementParticles, ParticleStatus)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.

   end subroutine DetermineAdjacentParticlesQUAD4_MP1

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP4(NElementParticles, ParticleStatus)
      
      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.

   end subroutine DetermineAdjacentParticlesQUAD4_MP4

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP9(NElementParticles, ParticleStatus)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.

   end subroutine DetermineAdjacentParticlesQUAD4_MP9

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP16(NElementParticles, ParticleStatus)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.
      ParticleStatus(10) = .true.
      ParticleStatus(11) = .true.
      ParticleStatus(12) = .true.
      ParticleStatus(13) = .true.
      ParticleStatus(14) = .true.
      ParticleStatus(15) = .true.
      ParticleStatus(16) = .true.


   end subroutine DetermineAdjacentParticlesQUAD4_MP16

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP25(NElementParticles, ParticleStatus)
      
      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.
      ParticleStatus(10) = .true.
      ParticleStatus(11) = .true.
      ParticleStatus(12) = .true.
      ParticleStatus(13) = .true.
      ParticleStatus(14) = .true.
      ParticleStatus(15) = .true.
      ParticleStatus(16) = .true.
      ParticleStatus(17) = .true.
      ParticleStatus(18) = .true.
      ParticleStatus(19) = .true.
      ParticleStatus(20) = .true.
      ParticleStatus(21) = .true.
      ParticleStatus(22) = .true.
      ParticleStatus(23) = .true.
      ParticleStatus(24) = .true.
      ParticleStatus(25) = .true.


   end subroutine DetermineAdjacentParticlesQUAD4_MP25

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP36(NElementParticles, ParticleStatus)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.
      ParticleStatus(10) = .true.
      ParticleStatus(11) = .true.
      ParticleStatus(12) = .true.
      ParticleStatus(13) = .true.
      ParticleStatus(14) = .true.
      ParticleStatus(15) = .true.
      ParticleStatus(16) = .true.
      ParticleStatus(17) = .true.
      ParticleStatus(18) = .true.
      ParticleStatus(19) = .true.
      ParticleStatus(20) = .true.
      ParticleStatus(21) = .true.
      ParticleStatus(22) = .true.
      ParticleStatus(23) = .true.
      ParticleStatus(24) = .true.
      ParticleStatus(25) = .true.
      ParticleStatus(26) = .true.
      ParticleStatus(27) = .true.
      ParticleStatus(28) = .true.
      ParticleStatus(29) = .true.
      ParticleStatus(30) = .true.
      ParticleStatus(31) = .true.
      ParticleStatus(32) = .true.
      ParticleStatus(33) = .true.
      ParticleStatus(34) = .true.
      ParticleStatus(35) = .true.
      ParticleStatus(36) = .true.


   end subroutine DetermineAdjacentParticlesQUAD4_MP36

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP49(NElementParticles, ParticleStatus)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.
      ParticleStatus(10) = .true.
      ParticleStatus(11) = .true.
      ParticleStatus(12) = .true.
      ParticleStatus(13) = .true.
      ParticleStatus(14) = .true.
      ParticleStatus(15) = .true.
      ParticleStatus(16) = .true.
      ParticleStatus(17) = .true.
      ParticleStatus(18) = .true.
      ParticleStatus(19) = .true.
      ParticleStatus(20) = .true.
      ParticleStatus(21) = .true.
      ParticleStatus(22) = .true.
      ParticleStatus(23) = .true.
      ParticleStatus(24) = .true.
      ParticleStatus(25) = .true.
      ParticleStatus(26) = .true.
      ParticleStatus(27) = .true.
      ParticleStatus(28) = .true.
      ParticleStatus(29) = .true.
      ParticleStatus(30) = .true.
      ParticleStatus(31) = .true.
      ParticleStatus(32) = .true.
      ParticleStatus(33) = .true.
      ParticleStatus(34) = .true.
      ParticleStatus(35) = .true.
      ParticleStatus(36) = .true.
      ParticleStatus(37) = .true.
      ParticleStatus(38) = .true.
      ParticleStatus(39) = .true.
      ParticleStatus(40) = .true.
      ParticleStatus(41) = .true.
      ParticleStatus(42) = .true.
      ParticleStatus(43) = .true.
      ParticleStatus(44) = .true.
      ParticleStatus(45) = .true.
      ParticleStatus(46) = .true.
      ParticleStatus(47) = .true.
      ParticleStatus(48) = .true.
      ParticleStatus(49) = .true.


   end subroutine DetermineAdjacentParticlesQUAD4_MP49

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP64(NElementParticles, ParticleStatus)
      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.
      ParticleStatus(10) = .true.
      ParticleStatus(11) = .true.
      ParticleStatus(12) = .true.
      ParticleStatus(13) = .true.
      ParticleStatus(14) = .true.
      ParticleStatus(15) = .true.
      ParticleStatus(16) = .true.
      ParticleStatus(17) = .true.
      ParticleStatus(18) = .true.
      ParticleStatus(19) = .true.
      ParticleStatus(20) = .true.
      ParticleStatus(21) = .true.
      ParticleStatus(22) = .true.
      ParticleStatus(23) = .true.
      ParticleStatus(24) = .true.
      ParticleStatus(25) = .true.
      ParticleStatus(26) = .true.
      ParticleStatus(27) = .true.
      ParticleStatus(28) = .true.
      ParticleStatus(29) = .true.
      ParticleStatus(30) = .true.
      ParticleStatus(31) = .true.
      ParticleStatus(32) = .true.
      ParticleStatus(33) = .true.
      ParticleStatus(34) = .true.
      ParticleStatus(35) = .true.
      ParticleStatus(36) = .true.
      ParticleStatus(37) = .true.
      ParticleStatus(38) = .true.
      ParticleStatus(39) = .true.
      ParticleStatus(40) = .true.
      ParticleStatus(41) = .true.
      ParticleStatus(42) = .true.
      ParticleStatus(43) = .true.
      ParticleStatus(44) = .true.
      ParticleStatus(45) = .true.
      ParticleStatus(46) = .true.
      ParticleStatus(47) = .true.
      ParticleStatus(48) = .true.
      ParticleStatus(49) = .true.
      ParticleStatus(50) = .true.
      ParticleStatus(51) = .true.
      ParticleStatus(52) = .true.
      ParticleStatus(53) = .true.
      ParticleStatus(54) = .true.
      ParticleStatus(55) = .true.
      ParticleStatus(56) = .true.
      ParticleStatus(57) = .true.
      ParticleStatus(58) = .true.
      ParticleStatus(59) = .true.
      ParticleStatus(60) = .true.
      ParticleStatus(61) = .true.
      ParticleStatus(62) = .true.
      ParticleStatus(63) = .true.
      ParticleStatus(64) = .true.



   end subroutine DetermineAdjacentParticlesQUAD4_MP64

   !**********************************************************************
   !
   !    Function:  Determines which particles of an element lie next to side ISide
   !               (linear quadrilateral element with initially 1 material point).
   !               Note: ParticleStatus is not initialised to .false. in order
   !                     to allow for a more flexible usage!
   !
   ! I  NElementParticles : Initial number of particles per element
   ! O  ParticleStatus : Set to .true.
   !
   !**********************************************************************
   subroutine DetermineAdjacentParticlesQUAD4_MP81(NElementParticles, ParticleStatus)
      
      implicit none

      integer(INTEGER_TYPE), intent(in) :: NElementParticles
      logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

      ParticleStatus(1) = .true.
      ParticleStatus(2) = .true.
      ParticleStatus(3) = .true.
      ParticleStatus(4) = .true.
      ParticleStatus(5) = .true.
      ParticleStatus(6) = .true.
      ParticleStatus(7) = .true.
      ParticleStatus(8) = .true.
      ParticleStatus(9) = .true.
      ParticleStatus(10) = .true.
      ParticleStatus(11) = .true.
      ParticleStatus(12) = .true.
      ParticleStatus(13) = .true.
      ParticleStatus(14) = .true.
      ParticleStatus(15) = .true.
      ParticleStatus(16) = .true.
      ParticleStatus(17) = .true.
      ParticleStatus(18) = .true.
      ParticleStatus(19) = .true.
      ParticleStatus(20) = .true.
      ParticleStatus(21) = .true.
      ParticleStatus(22) = .true.
      ParticleStatus(23) = .true.
      ParticleStatus(24) = .true.
      ParticleStatus(25) = .true.
      ParticleStatus(26) = .true.
      ParticleStatus(27) = .true.
      ParticleStatus(28) = .true.
      ParticleStatus(29) = .true.
      ParticleStatus(30) = .true.
      ParticleStatus(31) = .true.
      ParticleStatus(32) = .true.
      ParticleStatus(33) = .true.
      ParticleStatus(34) = .true.
      ParticleStatus(35) = .true.
      ParticleStatus(36) = .true.
      ParticleStatus(37) = .true.
      ParticleStatus(38) = .true.
      ParticleStatus(39) = .true.
      ParticleStatus(40) = .true.
      ParticleStatus(41) = .true.
      ParticleStatus(42) = .true.
      ParticleStatus(43) = .true.
      ParticleStatus(44) = .true.
      ParticleStatus(45) = .true.
      ParticleStatus(46) = .true.
      ParticleStatus(47) = .true.
      ParticleStatus(48) = .true.
      ParticleStatus(49) = .true.
      ParticleStatus(50) = .true.
      ParticleStatus(51) = .true.
      ParticleStatus(52) = .true.
      ParticleStatus(53) = .true.
      ParticleStatus(54) = .true.
      ParticleStatus(55) = .true.
      ParticleStatus(56) = .true.
      ParticleStatus(57) = .true.
      ParticleStatus(58) = .true.
      ParticleStatus(59) = .true.
      ParticleStatus(60) = .true.
      ParticleStatus(61) = .true.
      ParticleStatus(62) = .true.
      ParticleStatus(63) = .true.
      ParticleStatus(64) = .true.
      ParticleStatus(65) = .true.
      ParticleStatus(66) = .true.
      ParticleStatus(67) = .true.
      ParticleStatus(68) = .true.
      ParticleStatus(69) = .true.
      ParticleStatus(70) = .true.
      ParticleStatus(71) = .true.
      ParticleStatus(72) = .true.
      ParticleStatus(73) = .true.
      ParticleStatus(74) = .true.
      ParticleStatus(75) = .true.
      ParticleStatus(76) = .true.
      ParticleStatus(77) = .true.
      ParticleStatus(78) = .true.
      ParticleStatus(79) = .true.
      ParticleStatus(80) = .true.
      ParticleStatus(81) = .true.

   end subroutine DetermineAdjacentParticlesQUAD4_MP81

end module ModElementEvaluationQUAD

