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
	   
	   
	  module ModElementEvaluationTRI
      !**********************************************************************
      !
      !    Function:  This module provides specific routines for evaluating triangular elements.
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************

      use ModGeometryMath
      use ModString
      use ModReadCalculationData
      use ModGlobalConstants

      implicit none

      contains ! routines of this module

        subroutine DetermineSideDataTRI(ISide, PlaneNormal, PlanePoint)
        !**********************************************************************
        !
        !    Function: Returns a normal vector of the side with SideID and a 
        !              pointInitialLocalCoordinatesTriangle on that plane. 
        !              The normal vector is normalised and points inward.
        !    @note : 2D element
        !
        ! I  ISide : Side ID of the element
        ! O  PlaneNormal : Normal to the element side
        ! O  PlanePoint : Point on the element side
        !
        !**********************************************************************
        
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
              PlaneNormal(1) = -1.0
              PlaneNormal(2) = -1.0
              PlaneNormal = VectorNorm(PlaneNormal, IDim)
              PlanePoint(1) = 1.0
            case default
              call GiveError("DetermineSideDataTriangle: Side number " // trim(String(ISide)) // " does not exist with triangular elements")
          end select
        
        end subroutine DetermineSideDataTRI


        real(REAL_TYPE) function PointSideDistanceTRI(SideID, LocPos)
        !**********************************************************************
        !
        !    Function:  Returns the minimum distance between side SideID and LocPos.
        !    @note : 2D element
        !
        ! I  SideID : ID of the considered side (1..3)
        ! I  LocPos : Local coordinates of considered point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element       
          integer(INTEGER_TYPE), intent(in) :: SideID
          real(REAL_TYPE), dimension(IDim), intent(in) :: LocPos
          
          ! Local variables
          real(REAL_TYPE), dimension(IDim) :: LineNormal, LinePoint

          call DetermineSideDataTRI(SideID, LineNormal, LinePoint)
        
          PointSideDistanceTRI = LinePointDistance(LineNormal, LinePoint, LocPos)
        
        end function PointSideDistanceTRI


        subroutine InitialLocalMaterialPointCoordinatesTRI(IParticle, SolidPointsElement, LiquidPointsElement, WeiGP, PosGP)
        !**********************************************************************
        !
        !    Function : Determines the local coordinates and integration weight assigned to material point with ID
        !               IParticle which are returned through WeiGP and PosGP.
        !               Currently, all material points are placed at the same local positions.
        !    @note : 2D element
        !
        ! I  IParticle : Number of the material point
        ! I  SolidPointsElement : Number of solid material points per element
        ! I  LiquidPointsElement : Number of liquid material points per element
        ! O  WeiGP : Initial weight assigned to material point IParticle
        ! O  PosGP : Initial local position of material point IParticle
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle, SolidPointsElement, LiquidPointsElement
          real(REAL_TYPE), intent(inout) :: WeiGP
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP

          ! local variables
          integer :: ID

          ! first the solid material points are determined
          if ( (IParticle <= SolidPointsElement) .and. (SolidPointsElement > 0) ) then
            ID = IParticle
            select case(SolidPointsElement)
            case (1)! 1 material points per element
              call InitialTRI_MP1(ID, PosGP, WeiGP)
            case (3)! 3 material points per element
              call InitialTRI_MP3(ID, PosGP, WeiGP)
            case (6)! 6 material points per element
              call InitialTRI_MP6(ID, PosGP, WeiGP)
            case (12)! 12 material points per element
              call InitialTRI_MP12(ID, PosGP, WeiGP)
            case (16)! 16 material points per element
              call InitialTRI_MP16(ID, PosGP, WeiGP)
            case (25)! 25 material points per element
              call InitialTRI_MP25(ID, PosGP, WeiGP)
            case(46)! 46 material points per element
              call InitialTRI_MP46(ID, PosGP, WeiGP)
            case(88)! 88 material points per element
              call InitialTRI_MP88(ID, PosGP, WeiGP)
            case default
              call GiveError("InitialLocalCoordinatesTriangle: Number of solid material points, " // &
                              trim(String(SolidPointsElement)) // ", is not allowed for triangular elements!" // &
                              "Supported numbers: 1, 3, 6, 12, 16, 25, 46 & 88")
            end select
          end if
          
          if ( (IParticle > SolidPointsElement) .and. (LiquidPointsElement > 0) ) then
            ID = IParticle - SolidPointsElement
            select case(LiquidPointsElement)
            case (1)! 1 material points per element
              call InitialTRI_MP1(ID, PosGP, WeiGP)
            case (3)! 3 material points per element
              call InitialTRI_MP3(ID, PosGP, WeiGP)
            case (6)! 6 material points per element
              call InitialTRI_MP6(ID, PosGP, WeiGP)
            case (12)! 12 material points per element
              call InitialTRI_MP12(ID, PosGP, WeiGP)
            case (16)! 16 material points per element
              call InitialTRI_MP16(ID, PosGP, WeiGP)
            case (25)! 25 material points per element
              call InitialTRI_MP25(ID, PosGP, WeiGP)
            case(46)! 46 material points per element
              call InitialTRI_MP46(ID, PosGP, WeiGP)
            case(88)! 88 material points per element
              call InitialTRI_MP88(ID, PosGP, WeiGP)
            case default
              call GiveError("InitialLocalCoordinatesTriangle: Number of liquid material points, " // &
                              trim(String(LiquidPointsElement)) // ", is not allowed for triangular elements!" // &
                              "Supported numbers: 1, 3, 6, 12, 16, 25, 46 & 88")
            end select
            
          end if

        end subroutine InitialLocalMaterialPointCoordinatesTRI


        subroutine InitialTRI_MP1(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               1 particle is placed in each element.
        !    @note : 2D element
        !    @note : ref: http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !          
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.33333333333333
              PosGP(2) = 0.33333333333333
              WeiGP    = 0.5
          end select
        
        end subroutine InitialTRI_MP1


        subroutine InitialTRI_MP3(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               3 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !               [T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite
        !               Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987]
        !    @note : 2D element
        !    @note : ref: http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF   page 18
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        !                ^ Eta
        !                |
        !              3 + 
        !                | \
        !                |   \
        !                | +3  \
        !                |       \
        !                |         \
        !                |           \
        !                | +1     +2   \
        !                +---------------+----> xi
        !                1               2
        ! 
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.16666666666667
              PosGP(2) = 0.16666666666667
              WeiGP    = 0.33333333333333 * 0.5
            case (2)
              PosGP(1) = 0.66666666666667
              PosGP(2) = 0.16666666666667
              WeiGP    = 0.33333333333333 * 0.5
            case (3)
              PosGP(1) = 0.16666666666667
              PosGP(2) = 0.66666666666667
              WeiGP    = 0.33333333333333 * 0.5
          end select

          PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP3

        
        subroutine InitialTRI_MP6(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               6 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !               [T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite
        !               Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987]
        !    @note : 2D element
        !    @note : ref. http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP

          select case (IParticle)
            case (1)
              PosGP(1) = 0.44594849091597
              PosGP(2) = 0.44594849091597
              WeiGP    = 0.22338158967801 * 0.5
            case (2)
              PosGP(1) = 0.44594849091597
              PosGP(2) = 0.10810301816807
              WeiGP    = 0.22338158967801 * 0.5
            case (3)
              PosGP(1) = 0.10810301816807
              PosGP(2) = 0.44594849091597
              WeiGP    = 0.22338158967801 * 0.5
            case (4)
              PosGP(1) = 0.09157621350977
              PosGP(2) = 0.09157621350977
              WeiGP    = 0.10995174365532 * 0.5
            case (5)
              PosGP(1) = 0.09157621350977
              PosGP(2) = 0.81684757298046
              WeiGP    = 0.10995174365532 * 0.5
            case (6)
              PosGP(1) = 0.81684757298046
              PosGP(2) = 0.09157621350977
              WeiGP   =  0.10995174365532 * 0.5
          end select

          if (CalParams%UseUniformMaterialPointWeight) then
            WeiGP = (1.0/6.0) * 0.5
          end if

          PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP6
        

        subroutine InitialTRI_MP12(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               12 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !               [T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite
        !               Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987]
        !    @note : 2D element
        !    @note : ref. http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.24928674517091
              PosGP(2) = 0.24928674517091
              WeiGP    = 0.11678627572638 * 0.5
            case (2)
              PosGP(1) = 0.24928674517091
              PosGP(2) = 0.50142650965818
              WeiGP    = 0.11678627572638 * 0.5
            case (3)
              PosGP(1) = 0.50142650965818
              PosGP(2) = 0.24928674517091
              WeiGP    = 0.11678627572638 * 0.5
            case (4)
              PosGP(1) = 0.06308901449150
              PosGP(2) = 0.06308901449150
              WeiGP    = 0.05084490637021 * 0.5
            case (5)
              PosGP(1) = 0.06308901449150
              PosGP(2) = 0.87382197101700
              WeiGP    = 0.05084490637021 * 0.5
            case (6)
              PosGP(1) = 0.87382197101700
              PosGP(2) = 0.06308901449150
              WeiGP   =  0.05084490637021 * 0.5
            case (7)
              PosGP(1) = 0.31035245103378
              PosGP(2) = 0.63650249912140
              WeiGP   =  0.08285107561837 * 0.5
            case (8)
              PosGP(1) = 0.63650249912140
              PosGP(2) = 0.05314504984482
              WeiGP   =  0.08285107561837 * 0.5
            case (9)
              PosGP(1) = 0.05314504984482
              PosGP(2) = 0.31035245103378
              WeiGP   =  0.08285107561837 * 0.5
            case (10)
              PosGP(1) = 0.63650249912140
              PosGP(2) = 0.31035245103378
              WeiGP   =  0.08285107561837 * 0.5
            case (11)
              PosGP(1) = 0.31035245103378
              PosGP(2) = 0.05314504984482
              WeiGP   =  0.08285107561837 * 0.5
            case (12)
              PosGP(1) = 0.05314504984482
              PosGP(2) = 0.63650249912140
              WeiGP   =  0.08285107561837 * 0.5
          end select

          if (CalParams%UseUniformMaterialPointWeight) then
            WeiGP = (1.0/12.0) * 0.5
          end if

          PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP12
        

        subroutine InitialTRI_MP16(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               16 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !               [T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite
        !               Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987]
        !    @note : 2D element
        !    @note : ref. http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.33333333333333
              PosGP(2) = 0.33333333333333
              WeiGP    = 0.14431560767779 * 0.5
            case (2)
              PosGP(1) = 0.45929258829272
              PosGP(2) = 0.45929258829272
              WeiGP    = 0.09509163426728 * 0.5
            case (3)
              PosGP(1) = 0.45929258829272
              PosGP(2) = 0.08141482341455
              WeiGP    = 0.09509163426728 * 0.5
            case (4)
              PosGP(1) = 0.08141482341455
              PosGP(2) = 0.45929258829272
              WeiGP    = 0.09509163426728 * 0.5
            case (5)
              PosGP(1) = 0.17056930775176
              PosGP(2) = 0.17056930775176
              WeiGP    = 0.10321737053472 * 0.5
            case (6)
              PosGP(1) = 0.17056930775176
              PosGP(2) = 0.65886138449648
              WeiGP    = 0.10321737053472 * 0.5
            case (7)
              PosGP(1) = 0.65886138449648
              PosGP(2) = 0.17056930775176
              WeiGP    = 0.10321737053472 * 0.5
            case (8)
              PosGP(1) = 0.05054722831703
              PosGP(2) = 0.05054722831703
              WeiGP    = 0.03245849762320 * 0.5
            case (9)
              PosGP(1) = 0.05054722831703
              PosGP(2) = 0.89890554336594
              WeiGP    = 0.03245849762320 * 0.5
            case (10)
              PosGP(1) = 0.89890554336594
              PosGP(2) = 0.05054722831703
              WeiGP    = 0.03245849762320 * 0.5
            case (11)
              PosGP(1) = 0.26311282963464
              PosGP(2) = 0.72849239295540
              WeiGP    = 0.02723031417443 * 0.5
            case (12)
              PosGP(1) = 0.72849239295540
              PosGP(2) = 0.00839477740996
              WeiGP    = 0.02723031417443 * 0.5
            case (13)
              PosGP(1) = 0.00839477740996
              PosGP(2) = 0.26311282963464
              WeiGP    = 0.02723031417443 * 0.5
            case (14)
              PosGP(1) = 0.72849239295540
              PosGP(2) = 0.26311282963464
              WeiGP    = 0.02723031417443 * 0.5
            case (15)
              PosGP(1) = 0.26311282963464
              PosGP(2) = 0.00839477740996
              WeiGP    = 0.02723031417443 * 0.5
            case (16)
              PosGP(1) = 0.00839477740996
              PosGP(2) = 0.72849239295540
              WeiGP    = 0.02723031417443 * 0.5
          end select

          if (CalParams%UseUniformMaterialPointWeight) then
            WeiGP = (1.0/16.0) * 0.5
          end if

          PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP16
        

        subroutine InitialTRI_MP25(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               25 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !               [T.J.R. Hughes. The Finite Element Method: Linear Static and Dynamic Finite
        !               Element Analysis. Prentice-Hall, Inc., Englewood Cliffs, New Jersey, 1987]
        !    @note : 2D element
        !    @note : ref. https://www.clear.rice.edu/mech517/Books/FEEA/Chap_10.pdf
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IParticle)
            case (1)
              PosGP(1) = 0.333333333333333
              PosGP(2) = 0.333333333333333
              WeiGP    = 0.090817990382754 * 0.5
            case (2)
              PosGP(1) = 0.485577633383657
              PosGP(2) = 0.485577633383657
              WeiGP    = 0.036725957756467 * 0.5
            case (3)
              PosGP(1) = 0.485577633383657
              PosGP(2) = 0.028844733232685
              WeiGP    = 0.036725957756467 * 0.5
            case (4)
              PosGP(1) = 0.028844733232685
              PosGP(2) = 0.485577633383657
              WeiGP    = 0.036725957756467 * 0.5
            case (5)
              PosGP(1) = 0.109481575485037
              PosGP(2) = 0.109481575485037
              WeiGP    = 0.045321059435528 * 0.5
            case (6)
              PosGP(1) = 0.109481575485037
              PosGP(2) = 0.781036849029926
              WeiGP    = 0.045321059435528 * 0.5
            case (7)
              PosGP(1) = 0.781036849029926
              PosGP(2) = 0.109481575485037
              WeiGP    = 0.045321059435528 * 0.5
            case (8)
              PosGP(1) = 0.307939838764121
              PosGP(2) = 0.550352941820999
              WeiGP    = 0.072757916845420 * 0.5
            case (9)
              PosGP(1) = 0.550352941820999
              PosGP(2) = 0.141707219414880
              WeiGP    = 0.072757916845420 * 0.5
            case (10)
              PosGP(1) = 0.141707219414880
              PosGP(2) = 0.307939838764121
              WeiGP    = 0.072757916845420 * 0.5
            case (11)
              PosGP(1) = 0.550352941820999
              PosGP(2) = 0.307939838764121
              WeiGP    = 0.072757916845420 * 0.5
            case (12)
              PosGP(1) = 0.141707219414880
              PosGP(2) = 0.550352941820999
              WeiGP    = 0.072757916845420 * 0.5
            case (13)
              PosGP(1) = 0.307939838764121
              PosGP(2) = 0.141707219414880
              WeiGP    = 0.072757916845420 * 0.5
            case (14)
              PosGP(1) = 0.246672560639903
              PosGP(2) = 0.728323904597411
              WeiGP    = 0.028327242531057 * 0.5
            case (15)
              PosGP(1) = 0.728323904597411
              PosGP(2) = 0.025003534762686
              WeiGP    = 0.028327242531057 * 0.5
            case (16)
              PosGP(1) = 0.025003534762686
              PosGP(2) = 0.246672560639903
              WeiGP    = 0.028327242531057 * 0.5
            case (17)
              PosGP(1) = 0.728323904597411
              PosGP(2) = 0.246672560639903
              WeiGP    = 0.028327242531057 * 0.5
            case (18)
              PosGP(1) = 0.025003534762686
              PosGP(2) = 0.728323904597411
              WeiGP    = 0.028327242531057 * 0.5
            case (19)
              PosGP(1) = 0.246672560639903
              PosGP(2) = 0.025003534762686
              WeiGP    = 0.028327242531057 * 0.5
            case (20)
              PosGP(1) = 0.066803251012200
              PosGP(2) = 0.923655933587500
              WeiGP    = 0.009421666963733 * 0.5
            case (21)
              PosGP(1) = 0.923655933587500
              PosGP(2) = 0.009540815400299
              WeiGP    = 0.009421666963733 * 0.5
            case (22)
              PosGP(1) = 0.009540815400299
              PosGP(2) = 0.066803251012200
              WeiGP    = 0.009421666963733 * 0.5
            case (23)
              PosGP(1) = 0.923655933587500
              PosGP(2) = 0.066803251012200
              WeiGP    = 0.009421666963733 * 0.5
            case (24)
              PosGP(1) = 0.009540815400299
              PosGP(2) = 0.923655933587500
              WeiGP    = 0.009421666963733 * 0.5
            case (25)
              PosGP(1) = 0.066803251012200
              PosGP(2) = 0.009540815400299
              WeiGP    = 0.009421666963733 * 0.5
          end select

          if (CalParams%UseUniformMaterialPointWeight) then
            WeiGP = (1.0/25.0) * 0.5
          endif

          PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP25
        

        subroutine InitialTRI_MP46(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               46 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !    @note : 2D element
        !    @note : ref. http://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
        
          ! local variables
          integer(INTEGER_TYPE), parameter:: N_SIZE_21  = 7
          integer(INTEGER_TYPE), parameter:: N_SIZE_111 = 4
          integer(INTEGER_TYPE), parameter:: N_COMBINATION_21  = 3
          integer(INTEGER_TYPE), parameter:: N_COMBINATION_111  = 6
        
          real(REAL_TYPE), dimension(N_SIZE_21), parameter:: S21 = (/ 0.0099797608064584324152935295820524, &
                                                                      0.4799778935211883898105528650883899, &
                                                                      0.1538119591769669000000000000000000, &
                                                                      0.0740234771169878100000000000000000, &
                                                                      0.1303546825033300000000000000000000, &
                                                                      0.2306172260266531342996053700983831, &
                                                                      0.4223320834191478241144087137913939  /)

          real(REAL_TYPE), dimension(N_SIZE_21), parameter:: W21 = (/ 0.0017351512297252675680618638808094, &
                                                                      0.0261637825586145217778288591819783, &
                                                                      0.0039197292424018290965208275701454, &
                                                                      0.0122473597569408660972869899262505, &
                                                                      0.0281996285032579601073663071515657, &
                                                                      0.0508870871859594852960348275454540, &
                                                                      0.0504534399016035991910208971341189  /)

          real(REAL_TYPE), dimension(N_SIZE_111), parameter:: SA111 = (/ 0.7862373859346610033296221140330900, &
                                                                         0.6305521436606074416224090755688129, &
                                                                         0.6265773298563063142335123137534265, &
                                                                         0.9142099849296254122399670993850469  /)

          real(REAL_TYPE), dimension(N_SIZE_111), parameter:: SB111 = (/ 0.1906163600319009042461432828653034, &
                                                                         0.3623231377435471446183267343597729, &
                                                                         0.2907712058836674150248168174816732, &
                                                                         0.0711657108777507625475924502924336  /)

          real(REAL_TYPE), dimension(N_SIZE_111), parameter:: W111 = (/ 0.0170636442122334512900253993849472, &
                                                                        0.0096834664255066004075209630934194, &
                                                                        0.0363857559284850056220113277642717, &
                                                                        0.0069646633735184124253997225042413  /)

          integer(INTEGER_TYPE) :: iOffSet, iPack, index
          real(REAL_TYPE) :: a, b

          select case (IParticle)
            case (1)
              PosGP(1) = 0.333333333333333
              PosGP(2) = 0.333333333333333
              WeiGP    = 0.0585962852260285941278938063477560 * 0.5
            case (2:22)
              iOffSet = 1
              iPack = (((IParticle -1) - iOffSet) / N_COMBINATION_21) + 1
              index = (IParticle - iOffSet) - ((iPack - 1) * N_COMBINATION_21)
              select case (index)
              case (1)
                PosGP(1) = S21(iPack)
                PosGP(2) = S21(iPack)
                WeiGP    = W21(iPack) * 0.5
              case (2)
                PosGP(1) = S21(iPack)
                PosGP(2) = 1.0 - 2.0 * S21(iPack)
                WeiGP    = W21(iPack) * 0.5
              case (3)
                PosGP(1) = 1.0 - 2.0 * S21(iPack)
                PosGP(2) = S21(iPack)
                WeiGP    = W21(iPack) * 0.5
              case default
                call giveError('Error in initialisation of particle positions and weight. Undefined index number S12.')
              end select
            case (23:46)
              iOffSet = N_SIZE_21 * N_COMBINATION_21 + 1
              iPack = (((IParticle -1) - iOffSet) / N_COMBINATION_111) + 1
              index = (IParticle - iOffSet) - ((iPack - 1) * N_COMBINATION_111)
              a = SA111(iPack)
              b = SB111(iPack)
              WeiGP = W111(iPack) * 0.5
              select case (index)
              case (1)
                PosGP(1) = a
                PosGP(2) = b
              case (2)
                PosGP(1) = b
                PosGP(2) = a
              case (3)
                PosGP(1) = a
                PosGP(2) = 1.0 - a - b
              case (4)
                PosGP(1) = 1.0 - a - b
                PosGP(2) = a
              case (5)
                PosGP(1) = b
                PosGP(2) = 1.0 - a - b
              case (6)
                PosGP(1) = 1.0 - a - b
                PosGP(2) = b
              case default
                call giveError('Error in initialisation of particle positions and weight. Undefined index number S111.')
              end select
            case default
              call giveError('Error in initialisation of particle positions and weight. Undefined particle number.')
          end select

        if (CalParams%UseUniformMaterialPointWeight) then
          WeiGP = (1.0/46.0) * 0.5
        endif

        PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP46

        
        subroutine InitialTRI_MP88(IParticle, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the initial local coordinates and weight for IParticle.
        !               46 particles are placed in each element, whose initial
        !               locations and weights are identical with those of Gauss Points.
        !    @note : 2D element
        !    @note : ref. http://lsec.cc.ac.cn/~tcui/myinfo/paper/quad.pdf
        !
        ! I  IParticle : Local number of the considered particle inside an element
        ! O  PosGP : Returns the initial local coordinates of the particle with local number IParticle
        ! O  WeiGP : Returns the initial weight of IParticle
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          real(REAL_TYPE), dimension(:), intent(out) :: PosGP
          real(REAL_TYPE), intent(out) :: WeiGP
          
          ! Local variables
          integer(INTEGER_TYPE), parameter:: N_SIZE_21  = 5
          integer(INTEGER_TYPE), parameter:: N_SIZE_111 = 12
          integer(INTEGER_TYPE), parameter:: N_COMBINATION_21  = 3
          integer(INTEGER_TYPE), parameter:: N_COMBINATION_111  = 6
        
          real(REAL_TYPE), dimension(N_SIZE_21), parameter:: S21 = (/ 0.2158743059329919731902545438401828, &
                                                                      0.0753767665297472780972854309459163, &
                                                                      0.0103008281372217921136862160096969, &
                                                                      0.4936022112987001655119208321450536, &
                                                                      0.4615509381069252967410487102915180  /)
          real(REAL_TYPE), dimension(N_SIZE_21), parameter:: W21 = (/ 0.0274718698764242137484535496073598, &
                                                                      0.0097652722770514230413646914294237, &
                                                                      0.0013984195353918235239233631597867, &
                                                                      0.0092921026251851826304282034030330, &
                                                                      0.0165778760323669253260236250351840  /)
          real(REAL_TYPE), dimension(N_SIZE_111), parameter:: SA111 = (/ 0.3286214064242369933034974609509133, &
                                                                         0.2604803617865687564195930170811535, &
                                                                         0.1370742358464553000000000000000000, &
                                                                         0.1467269458722997843041609884874530, &
                                                                         0.0269989777425532900000000000000000, &
                                                                         0.0618717859336170268417124700122339, &
                                                                         0.0477243674276219962083526801042934, &
                                                                         0.1206005151863643799672337870400794, &
                                                                         0.0026971477967097876716489145012827, &
                                                                         0.0030156332779423626572762598234710, &
                                                                         0.0299053757884570188069287738643386, &
                                                                         0.0067566542224609885399458175192278  /)
          real(REAL_TYPE), dimension(N_SIZE_111), parameter:: SB111 = (/ 0.4293405702582103752139588004663984, &
                                                                         0.1015775342809694461687550061961797, &
                                                                         0.7100659730011301599879040745464079, &
                                                                         0.4985454776784148493896226967076119, &
                                                                         0.0491867226725820016197037125775872, &
                                                                         0.7796601465405693953603506190768108, &
                                                                         0.3704915391495476369201496202567388, &
                                                                         0.8633469487547526484979879960925217, &
                                                                         0.0561949381877455029878923019865887, &
                                                                         0.2086750067484213509575944630613577, &
                                                                         0.7211512409120340910281041502050941, &
                                                                         0.6400554419405418899040536682721647  /)
          real(REAL_TYPE), dimension(N_SIZE_111), parameter:: W111 = (/ 0.0206677623486650769614219700129729, &
                                                                        0.0208222355211545073068785561993297, &
                                                                        0.0095686384198490606888758450458320, &
                                                                        0.0244527709689724638856439207024089, &
                                                                        0.0031557306306305340038264003207296, &
                                                                        0.0121367963653212969370133090807574, &
                                                                        0.0149664801438864490365249118515707, &
                                                                        0.0063275933217777395693240327504398, &
                                                                        0.0013425603120636958849798512981433, &
                                                                        0.0027760769163475540677293561558015, &
                                                                        0.0107398444741849415551734474479517, &
                                                                        0.0053678057381874532052474100212697  /)
          integer :: iOffSet, iPack, index
          real(REAL_TYPE) :: a, b

          select case (IParticle)
            case (1)
              PosGP(1) = 0.333333333333333
              PosGP(2) = 0.333333333333333
              WeiGP    = 0.0125376079944966565735856367723948 * 0.5
            case (2:16)
              iOffSet = 1
              iPack = (((IParticle -1) - iOffSet) / N_COMBINATION_21) + 1
              index = (IParticle - iOffSet) - ((iPack - 1) * N_COMBINATION_21)
              select case (index)
              case (1)
                PosGP(1) = S21(iPack)
                PosGP(2) = S21(iPack)
                WeiGP    = W21(iPack) * 0.5
              case (2)
                PosGP(1) = S21(iPack)
                PosGP(2) = 1.0 - 2.0 * S21(iPack)
                WeiGP    = W21(iPack) * 0.5
              case (3)
                PosGP(1) = 1.0 - 2.0 * S21(iPack)
                PosGP(2) = S21(iPack)
                WeiGP    = W21(iPack) * 0.5
              case default
                call giveError('Error in initialisation of particle positions and weight. Undefined index number S12.')
              end select
            case (17:88)
              iOffSet = N_SIZE_21 * N_COMBINATION_21 + 1
              iPack = (((IParticle -1) - iOffSet) / N_COMBINATION_111) + 1
              index = (IParticle - iOffSet) - ((iPack - 1) * N_COMBINATION_111)
              a = SA111(iPack)
              b = SB111(iPack)
              WeiGP = W111(iPack) * 0.5
              select case (index)
              case (1)
                PosGP(1) = a
                PosGP(2) = b
              case (2)
                PosGP(1) = b
                PosGP(2) = a
              case (3)
                PosGP(1) = a
                PosGP(2) = 1.0 - a - b
              case (4)
                PosGP(1) = 1.0 - a - b
                PosGP(2) = a
              case (5)
                PosGP(1) = b
                PosGP(2) = 1.0 - a - b
              case (6)
                PosGP(1) = 1.0 - a - b
                PosGP(2) = b
              case default
                call giveError('Error in initialisation of particle positions and weight. Undefined index number S111.')
              end select
            case default
              call giveError('Error in initialisation of particle positions and weight. Undefined particle number.')
          end select

          if (CalParams%UseUniformMaterialPointWeight) then
            WeiGP = (1.0/88.0) * 0.5
          endif

          PosGP = (PosGP - 0.33333333333333d0) * (1.0 - CalParams%ShrinkageMateriaPointPositionFactor) +  0.33333333333333d0

        end subroutine InitialTRI_MP88
        

        subroutine GaussTRI_Q1(IGaussPoint, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the local coordinates and weight for IGaussPoint.
        !               1 Gauss Point is located in each triangular element.
        !    @note : 2D element
        !    @note : ref: http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !
        ! I  IGaussPoint : Local number of the considered Gauss Point inside an element
        ! O  PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        ! O  WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IGaussPoint)
              
            case (1)
              PosGP(1) = 0.33333333333333
              PosGP(2) = 0.33333333333333
              WeiGP    = 0.5
              
            case default
              call GiveError("Number of Gauss points not defined in [subroutine GaussTRI_Q1()].")
              
          end select

        end subroutine GaussTRI_Q1


        subroutine GaussTRI_Q3(IGaussPoint, PosGP, WeiGP)
        !**********************************************************************
        !
        !    Function:  Returns the local coordinates and weight for IGaussPoint.
        !               3 Gauss Points are located in each triangular element.
        !> @note : 2D element
        !> @note : ref: http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
        !
        ! I  IGaussPoint : Local number of the considered Gauss Point inside an element
        ! O  PosGP : Returns the xi, eta local coordinates of the Gauss Point with local number IGaussPoint
        ! O  WeiGP : Return the weight of IGaussPoint
        ! 
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IGaussPoint
          real(REAL_TYPE), dimension(:), intent(inout) :: PosGP
          real(REAL_TYPE), intent(inout) :: WeiGP
          
          select case (IGaussPoint)
            case (1)
              PosGP(1) = 0.16666666666667
              PosGP(2) = 0.16666666666667
              WeiGP    = 0.33333333333333 * 0.5
            case (2)
              PosGP(1) = 0.16666666666667
              PosGP(2) = 0.66666666666667
              WeiGP =    0.33333333333333 * 0.5
            case (3)
              PosGP(1) = 0.66666666666667
              PosGP(2) = 0.16666666666667
              WeiGP    = 0.33333333333333 * 0.5

            case default
              call GiveError("Number of Gauss points not defined in [subroutine GaussTRI_Q3()].")
              
          end select
         
        end subroutine GaussTRI_Q3

     
        subroutine InitialiseShapeFunctionsTRI3(HS, dHS, Wt)
        !**********************************************************************
        !
        !    Function:  Calculates the values of shape functions and their
        !                derivatives at one Gaussian integration point for
        !                3-noded 2D triangular elements.
        !    @note: 2D element
        !
        ! I/O  HS(i,j) : Shape function j at integration point i
        ! I/O  dHS(i,j,k) : Derivative of shape function j at integration point i 
        !                   with respect to direction k
        ! I/O  Wt : Local weights for integration
        !
        !
        !                ^ Eta
        !                |
        !              3 + 
        !                | \
        !                |   \
        !                |     \
        !                |       \
        !                |         \
        !                |           \
        !                |             \
        !                +---------------+----> xi
        !                1               2
        !**********************************************************************      
        implicit none

        real(REAL_TYPE), dimension(:, :), intent(inout) :: HS
        real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS
        real(REAL_TYPE), dimension(:), intent(inout) :: Wt
      
        ! local variables
        real(REAL_TYPE) :: xi, eta
        integer(INTEGER_TYPE) :: int, I1, Nint1
      
        Nint1=1
        Int = 0
      
        do I1 = 1, Nint1

          Xi = 1 / 3d0
          Eta = 1 / 3d0

          Int = Int+1

          Wt(Int) = 1d0 / Nint1 * 0.5
               
          Hs(Int,1) = (1-Xi-Eta) 
          Hs(Int,2) = Xi
          Hs(Int,3) = Eta

          ! dHs(Int,i,1) = dHS(i)/dXi
          dHs(Int,1,1) = -1
          dHs(Int,2,1) = 1
          dHs(Int,3,1) = 0

          ! dHs(i,2) = dHS(i)/dEta
          dHs(Int,1,2) = -1
          dHs(Int,2,2) = 0
          dHs(Int,3,2) = 1
        
        end do

      end subroutine InitialiseShapeFunctionsTRI3
       
       
        subroutine ShapeLocPosTRI3(LocPos, HS, dHS)
        !**********************************************************************
        !
        !    Function:  To calculate the values of shape functions and their
        !               derivatives at LocPos for 3-noded 2D triangular elements.
        !    @note : 2D element
        !
        ! I  LocPos : Local coordinates of a point inside an element
        ! O  HS(i) : Shape function i
        ! O  dHS(i,j) : Derivative of shape function i with respect to direction j
        !
        !                ^ Eta
        !                |
        !                + 3
        !                |  \
        !                |     \
        !                |
        !                |        \
        !                |          \
        !                |             \
        !                |1              2
        !                +---------------+----> xi
        !**********************************************************************

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos 
          real(REAL_TYPE), dimension(:), intent(out) :: HS 
          real(REAL_TYPE), dimension(:,:), intent(out) :: dHS 

          ! Local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi  = LocPos(1)
          Eta = LocPos(2)
          
          HS(1) = 1.0 - Xi - Eta
          HS(2) = Xi
          HS(3) = Eta
          
          dHS(1,1) = -1.0
          dHS(2,1) = 1.0
          dHS(3,1) = 0.0
          
          dHS(1,2) = -1.0
          dHS(2,2) = 0.0
          dHS(3,2) = 1.0

        end subroutine ShapeLocPosTRI3


        subroutine ShapeLocPosTRI6(LocPos, HS, dHS)
        !**********************************************************************
        !
        !    Function:  To calculate the values of shape functions and their
        !               derivatives at LocPos for 6-noded 2D triangular elements.
        !               [R. K. Livesley, Finite Elements: An Introduction for 
        !                Engineers, CUP Archive 1983, p. 81]
        !    @note : 2D element
        !    @note : ref http://www.ippt.pan.pl/Repository/o2210.pdf
        !
        ! I  LocPos : Local coordinates of a point inside an element
        ! O  HS(i) : Shape function i
        ! O  dHS(i,j) : Derivative of shape function i with respect to direction j
        !
        !                ^ Eta
        !                |
        !                + 3
        !                |  \
        !                |     \
        !                |
        !                +6      +5
        !                |          \
        !                |             \
        !                |1      4        2
        !                +-------+-------+----> xi
        !**********************************************************************

        implicit none
      
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos 
          real(REAL_TYPE), dimension(:), intent(out) :: HS 
          real(REAL_TYPE), dimension(:,:), intent(out) :: dHS 

          ! Local variables
          real(REAL_TYPE) :: Xi, Eta

          Xi = LocPos(1)
          Eta = LocPos(2)

          ! 1..3 corners
          HS(1) = Xi  * (2.0 * Xi - 1.0)                          ! 2 -> 1
          HS(2) = Eta * (2.0 * Eta - 1.0)                         ! 3 -> 2
          HS(3) = (1.0 - Xi - Eta) * (1.0 - 2.0 * Xi - 2.0 * Eta) ! 1 -> 3
          ! 4..6 mids
          HS(4) = 4.0 * Xi * Eta                                  ! 5 -> 4
          HS(5) = 4.0 * Eta * (1.0 - Xi - Eta)                    ! 6 -> 5
          HS(6) = 4.0 * Xi * (1.0 - Xi - Eta)                     ! 4 -> 6

          ! dHS(,,1) = d HS / dXi
          ! 1..3 corners
          dHS(1, 1) = (1.0) * (2.0 * Xi - 1.0) + Xi * (2.0)
          dHS(2, 1) = 0.0
          dHS(3, 1) = (-1.0) * (1.0 - 2.0 * Xi - 2.0 * Eta) + (1.0 - Xi - Eta) * (-2.0)
          ! 5..10 mids
          dHS(4, 1) = 4.0 * Eta
          dHS(5, 1) = 4.0 * Eta * (-1.0)
          dHS(6, 1) = (4.0) * (1.0 - Xi - Eta) + 4.0 * Xi * (-1.0)

          ! dHS(,,2) = d HS / dEta
          ! 1..3 corners
          dHS(1, 2) = 0.0
          dHS(2, 2) = (1.0) * (2.0 * Eta - 1.0) + Eta * (2.0)
          dHS(3, 2) = (-1.0) * (1.0 - 2.0 * Xi - 2.0 * Eta) + (1.0 - Xi - Eta) * (-2.0)
          ! 5..10 mids
          dHS(4, 2) = 4.0 * Xi
          dHS(5, 2) = (4.0) * (1.0 - Xi - Eta) + 4.0 * Eta * (-1.0)
          dHS(6, 2) = 4.0 * Xi * (-1.0)

        end subroutine ShapeLocPosTRI6


        subroutine GetLocalCoordinatesTRI3(GlobPos, LocPos, OutsideElement, MInv, MIX1, CrossedSide)
        !**********************************************************************
        !
        !    Function:  Determination of local coordinates from global coordinates,
        !               assuming that the point with global coordinates lies inside 
        !               the triangular element IElement.
        !               OutsideElement returns .false. if the local position is in the element.
        !               CrossedSide returns the number of the side that has been crossed.
        !
        ! I  GlobPos : Global coordinates of a point inside IElement
        ! I  MInv : Element matrix
        ! I  MIX1 : Element vector of first element node
        ! O  LocPos : Local coordinates of the considered point inside IElement
        ! O  OutsideElement : True, if the local coordinate lie outside IElement
        ! O  CrossedSide : number of crossed side, if local coordinate lie outside IElement
        !                ^ Eta
        !                |
        !                + 3
        !                |  \
        !                |    \
        !                |      \ side 3
        !         side 2 |        \
        !                |          \
        !                |            \
        !                |1              2
        !                +---------------+----> xi
        !                    side 1
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          real(REAL_TYPE), dimension(IDim), intent(out):: LocPos
          logical, intent(out) :: OutsideElement
          real(REAL_TYPE), dimension(IDim, IDim), intent(in) :: MInv
          real(REAL_TYPE), dimension(IDim), intent(in) :: MIX1
          integer(INTEGER_TYPE), intent(out) :: CrossedSide
          
          ! Local variables
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
          if     (LocPos(1) < 0) then
            CrossedSide = 2
          elseif (LocPos(2) < 0) then
            CrossedSide = 1
          elseif ((LocPos(1) + LocPos(2)) > 1) then
            CrossedSide = 3
          else
            OutsideElement = .false.
          end if
        
        end subroutine GetLocalCoordinatesTRI3


        logical function IsInsideElementLocPosTRI(LocPos)
        !**********************************************************************
        !
        !    Function:  Returns .true. if LocPos (local coordinates) lies inside the 
        !               area of the triangular element.
        !
        ! I  LocPos : Local coordinates of the considered point inside an element
        ! O  IsInsideElementLocPosTriangle : True, if the point lies inside the element
        !
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        
          if ( (LocPos(1) < 0.0) .or. (LocPos(2) < 0.0) .or. (LocPos(1) + LocPos(2) > 1.0) ) then
            IsInsideElementLocPosTRI = .false.
          else
            IsInsideElementLocPosTRI = .true.
          end if
        
        end function IsInsideElementLocPosTRI


        logical function IsInsideElementGlobPosTRI(GlobPos, ElementID, NodTot, IElTyp, NEl, NodeCoord, ICon)
        !**********************************************************************
        !
        !    Function:  Returns .true. if GlobPos (global coordinates) lies inside the 
        !               area of the triangular element.
        !
        ! I  GlobPos : Global coordinates of the considered point inside an element
        ! I  ElementID : ID of the considered element
        ! I  NodTot : Total number of nodes
        ! I  IElTyp : Number of nodes per element
        ! I  NEl : Total number of elements
        ! I  NodeCoord : Nodal coordinates
        ! I  ICon : Element connectivities
        ! O  IsInsideElementGlobPosTriangle : True, if the point lies inside the element
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as 2D element
          real(REAL_TYPE), dimension(IDim), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          integer(INTEGER_TYPE), intent(in) :: NodTot, IElTyp, NEl
          real(REAL_TYPE), dimension(NodTot, IDim), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(IElTyp, NEl), intent(in) :: ICon
          
          ! Local variables
          real(REAL_TYPE), dimension(IDim) :: A, B, C ! Vertices of triangle
          
          A = NodeCoord(ICon(1, ElementID), 1:2)
          B = NodeCoord(ICon(2, ElementID), 1:2)
          C = NodeCoord(ICon(3, ElementID), 1:2)

          IsInsideElementGlobPosTRI = CheckInsideTriangle(A, B, C, GlobPos)

        end function IsInsideElementGlobPosTRI


        subroutine DetermineAdjacentParticlesTRI_MP1(NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 1 material point).
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true.
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

          ParticleStatus(1) = .true.

        end subroutine DetermineAdjacentParticlesTRI_MP1


        subroutine DetermineAdjacentParticlesTRI_MP3(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 3 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               Particle 1: (0.167, 0.167) near node 1
        !               Particle 2: (0.667, 0.167) near node 2
        !               Particle 3: (0.167, 0.667) near node 3
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 1, 2
              ParticleStatus(1) = .true.
              ParticleStatus(2) = .true.
            case (2) ! Particles 1, 3
              ParticleStatus(1) = .true.
              ParticleStatus(3) = .true.
            case (3) ! Particles 2, 3
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP3
        
        
        subroutine DetermineAdjacentParticlesTRI_MP6(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 6 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 2, 4, 6
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
            case (2) ! Particles 3, 4, 5
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
            case (3) ! Particles 1, 5, 6
              ParticleStatus(1) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(6) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP6
        
        
        subroutine DetermineAdjacentParticlesTRI_MP12(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 12 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 4, 6, 8, 11            
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(11) = .true.
            case (2) ! Particles 4, 5, 9, 12
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(12) = .true.
            case (3) ! Particles 5, 6, 7, 10
              ParticleStatus(5) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(10) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP12
        
        
        subroutine DetermineAdjacentParticlesTRI_MP16(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 16 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 3, 8, 10, 12, 15
              ParticleStatus(3) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(12) = .true.
              ParticleStatus(15) = .true.
            case (2) ! Particles 4, 8, 9, 13, 16
              ParticleStatus(4) = .true.
              ParticleStatus(8) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(13) = .true.
              ParticleStatus(16) = .true.
            case (3) ! Particles 2, 9, 10, 11, 14
              ParticleStatus(2) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(11) = .true.
              ParticleStatus(14) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP16
        
        
        subroutine DetermineAdjacentParticlesTRI_MP25(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 25 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 3, 15, 19, 21, 25
              ParticleStatus(3) = .true.
              ParticleStatus(15) = .true.
              ParticleStatus(19) = .true.
              ParticleStatus(21) = .true.
              ParticleStatus(25) = .true.
            case (2) ! Particles 4, 16, 18, 22, 24
              ParticleStatus(4) = .true.
              ParticleStatus(16) = .true.
              ParticleStatus(18) = .true.
              ParticleStatus(22) = .true.
              ParticleStatus(24) = .true.
            case (3) ! Particles 2, 14, 17, 20, 23
              ParticleStatus(2) = .true.
              ParticleStatus(14) = .true.
              ParticleStatus(17) = .true.
              ParticleStatus(20) = .true.
              ParticleStatus(23) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP25
        
        
        subroutine DetermineAdjacentParticlesTRI_MP46(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 46 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 2, 4, 6, 25, 27, 31, 33, 43, 45
              ParticleStatus(2) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(6) = .true.
              ParticleStatus(25) = .true.
              ParticleStatus(27) = .true.
              ParticleStatus(31) = .true.
              ParticleStatus(33) = .true.
              ParticleStatus(43) = .true.
              ParticleStatus(45) = .true.
            case (2) ! Particles 2, 3, 7, 26, 28, 32, 34, 44, 46
              ParticleStatus(2) = .true.
              ParticleStatus(3) = .true.
              ParticleStatus(7) = .true.
              ParticleStatus(26) = .true.
              ParticleStatus(28) = .true.
              ParticleStatus(32) = .true.
              ParticleStatus(34) = .true.
              ParticleStatus(44) = .true.
              ParticleStatus(46) = .true.
            case (3) ! Particles 3, 4, 5, 23, 24, 29, 30, 41, 42
              ParticleStatus(3) = .true.
              ParticleStatus(4) = .true.
              ParticleStatus(5) = .true.
              ParticleStatus(23) = .true.
              ParticleStatus(24) = .true.
              ParticleStatus(29) = .true.
              ParticleStatus(30) = .true.
              ParticleStatus(41) = .true.
              ParticleStatus(42) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP46
               
        
        subroutine DetermineAdjacentParticlesTRI_MP88(ISide, NElementParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles of an element lie next to side ISide
        !               (triangular element with initially 88 material points).
        !               Side 1 (nodes 1-2 & 4) at xi line
        !               Side 2 (nodes 1-3 & 6) at eta line
        !               Side 3 (nodes 2-3 & 5) at 'inclined' line
        !               
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage!
        !
        ! I  ISide : Local number of considered element side
        ! I  NElementParticles : Initial number of particles per element
        ! O  ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: ISide, NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus
        
          select case (ISide)
            case (1) ! Particles 8, 10, 12, 61, 63, 66, 68, 72, 74, 84, 86
              ParticleStatus(8) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(12) = .true.
              ParticleStatus(61) = .true.
              ParticleStatus(63) = .true.
              ParticleStatus(66) = .true.
              ParticleStatus(68) = .true.
              ParticleStatus(72) = .true.
              ParticleStatus(74) = .true.
              ParticleStatus(84) = .true.
              ParticleStatus(86) = .true.
            case (2) ! Particles 8, 9, 13, 62, 64, 65, 67, 71, 73, 83, 85
              ParticleStatus(8) = .true.
              ParticleStatus(9) = .true.
              ParticleStatus(13) = .true.
              ParticleStatus(62) = .true.
              ParticleStatus(64) = .true.
              ParticleStatus(65) = .true.
              ParticleStatus(67) = .true.
              ParticleStatus(71) = .true.
              ParticleStatus(73) = .true.
              ParticleStatus(83) = .true.
              ParticleStatus(85) = .true.
            case (3) ! Particles 9, 10, 11, 59, 60, 69, 75, 76, 81, 82, 87
              ParticleStatus(9) = .true.
              ParticleStatus(10) = .true.
              ParticleStatus(11) = .true.
              ParticleStatus(59) = .true.
              ParticleStatus(60) = .true.
              ParticleStatus(69) = .true.
              ParticleStatus(75) = .true.
              ParticleStatus(76) = .true.
              ParticleStatus(81) = .true.
              ParticleStatus(82) = .true.
              ParticleStatus(87) = .true.
          end select
        
        end subroutine DetermineAdjacentParticlesTRI_MP88


        subroutine DetermineCheckEdgeNodesTRI6(CheckEdgeNodes)
        !**********************************************************************
        !
        !    Function:  Determines which edge node of each side to check in order
        !               to detect an adjacent element for each side.
        !
        ! O  CheckEdgeNodes : Array containing for each side of the element, the local number of the edge node
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), dimension(:, :), intent(inout) :: CheckEdgeNodes ! size(nCheckNodes, NumberOfElementSides)

          CheckEdgeNodes(1,1) = 4 ! Edge node to check for side 1 spanned by nodes 1-2
          CheckEdgeNodes(1,2) = 6 ! Edge node to check for side 2 spanned by nodes 1-3
          CheckEdgeNodes(1,3) = 5 ! Edge node to check for side 3 spanned by nodes 2-3

        end subroutine DetermineCheckEdgeNodesTRI6
      
        
        integer function DetermineSideNodesTRI6(SideID, LocalNodeID)
        !**********************************************************************
        !
        !    Function:  Returns the element connectivity of LocalNodeID of side SideID.
        !               Side 1 is spanned by nodes 1, 2, 4
        !               Side 2 is spanned by nodes 1, 3, 6
        !               Side 3 is spanned by nodes 2, 3, 5
        !
        ! I  SideID :      ID of the considered side (1 .. 3)
        ! I  LocalNodeID : ID of the side node  (1 .. 6)
        ! O  DetermineSideNodesTRI6 : Local ID of the considered node (1 .. 6)
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: SideID, LocalNodeID
        
          ! Local variables
          integer(INTEGER_TYPE), dimension(3, 3) :: SideConnectivities

          SideConnectivities = reshape( (/  1, 2, 4, &
                                            3, 1, 6, &
                                          2, 3, 5/), &
                                        (/ 3, 3 /) )

          DetermineSideNodesTRI6 = SideConnectivities(LocalNodeID, SideID)

        end function DetermineSideNodesTRI6


        subroutine CheckTRIForGlobPos(GlobPos, ElementID, CentrePoint, NodeCoord, ICon, CrossedSide, IsInside)
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
        !    @note : 2D element
        !
        ! I  GlobPos : Global coordinates of a point inside the mesh
        ! I  ElementID : Considered element
        ! I  CentrePoint : Centrepoint of ElementID
        ! I  NodeCoord : Global nodal coordinates
        ! I  ICon : Element connectivities ICon(I, J): global node number of local node I in element J
        ! O  CrossedSide : Side which contains the intersection point of the above mentioned line
        ! O  IsInside : True, if GlobPos lies inside ElementID
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), parameter :: IDim = 2 ! fixed dimension as only 2D functionality
          real(REAL_TYPE), dimension(:), intent(in) :: GlobPos
          integer(INTEGER_TYPE), intent(in) :: ElementID
          real(REAL_TYPE), dimension(:), intent(in) :: CentrePoint
          real(REAL_TYPE), dimension(:, :), intent(in) :: NodeCoord
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          integer(INTEGER_TYPE), intent(out) :: CrossedSide
          logical, intent(out) :: IsInside
        
          ! Local variables
          integer(INTEGER_TYPE), dimension(2, 3) :: VerticesID ! size(N_BOUND_TRIANGLE_ELEMENT_NODES_LOE, N_SIDE_TRIANGLE)
          real(REAL_TYPE), dimension(IDim) :: A, B ! coordinates of nodes of one side of a triangle
          integer :: Success, I

          VerticesID = reshape( (/ 1, 2,   & ! Corner nodes to check for side 1 spanned by nodes 1-2
                                   3, 1,   & ! Corner nodes to check for side 2 spanned by nodes 3-1
                                   2, 3/), & ! Corner nodes to check for side 3 spanned by nodes 2-3
                                (/ 2, 3 /) ) ! 6 elements which should be reshaped to a 2x3 matrix. N_SIDE is always 3 for a triangle. 

          do I = 1, 3 ! Loop over sides of triangle
            A = NodeCoord(ICon(VerticesID(1, I), ElementID), 1:NDIM)
            B = NodeCoord(ICon(VerticesID(2, I), ElementID), 1:NDIM)

            Success = CheckInsideSubTriangle(A, B, CentrePoint, GlobPos)

            if (Success == 1) then
              CrossedSide = I
              EXIT
            else if (Success == 2) then
              IsInside = .true.
              EXIT
            end if
          end do

        end subroutine CheckTRIForGlobPos


        subroutine RearrangeConnectivitiesBoundaryLine3(IConGlobal, ValuesGlobal, IConLocal,  ValuesLocal, NSurfaceNodes)
        !**********************************************************************
        !
        !    Function:  Rearranges the values of the two arrays 
        !               IConGlobal and ValuesGlobal and returns them
        !               in the arrays IConLocal and ValuesLocal
        !               for a 3-noded line element.
        !               The connectivities are rearranged such that corner nodes
        !               are located in the array first, afterwards the mid-nodes.
        !
        ! I  IConGlobal : Node connectivities
        ! I  ValuesGlobal : Values belonging to the connectivities of IConGlobal
        ! O  IConLocal : Rearranged node connectivities
        ! O  ValuesLocal : Rearranged node values
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NSurfaceNodes
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(in) :: IConGlobal
          real(REAL_TYPE), dimension(NSurfaceNodes, 2), intent(in) :: ValuesGlobal
          integer(INTEGER_TYPE), dimension(NSurfaceNodes), intent(out) :: IConLocal
          real(REAL_TYPE), dimension(NSurfaceNodes, 2), intent(out) :: ValuesLocal

          IConLocal(1) = IConGlobal(1)
          IConLocal(2) = IConGlobal(3)
          IConLocal(3) = IConGlobal(2)

          ValuesLocal(1, :) = ValuesGlobal(1, :)
          ValuesLocal(2, :) = ValuesGlobal(3, :)
          ValuesLocal(3, :) = ValuesGlobal(2, :)

        end subroutine RearrangeConnectivitiesBoundaryLine3
      
      
        subroutine InitialiseShapeFunctionsLINE2(HS, dHS, Wt)
        !**********************************************************************
        !
        !    Function:  To calculate the values of shape functions and their
        !               derivatives at one Gaussian integration point for
        !               2-noded 1D element.
        !
        ! O  HS(i,j)    : Shape function j at integration point i
        ! O  dHS(i,j,k) : Derivative of shape function j at integration point i
        !                 with respect to direction k
        ! O  Wt         : Local weights for integration
        !
        !                +---------------+----> xi
        !                1               2
        !
        !**********************************************************************
        implicit none
      
          real(REAL_TYPE), dimension(:, :), intent(inout) :: HS
          real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS
          real(REAL_TYPE), dimension(:), intent(inout) :: Wt
        
          ! local variables 
          integer :: I, NInt, Int
          real(REAL_TYPE) :: xi

          NInt = 1
          Int = 0

          do I = 1, NInt

            Xi = 0.0

            Int = Int + 1

            Wt(Int) = 2.0

            HS(Int, 1) = ( 1 - Xi ) / 2.0
            HS(Int, 2) = ( 1 + Xi ) / 2.0

            dHS(Int, 1, 1) = -0.5
            dHS(Int, 2, 1) = 0.5

          end do
      
        end subroutine InitialiseShapeFunctionsLINE2
      
      
      end module ModElementEvaluationTRI