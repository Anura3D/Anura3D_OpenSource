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


	! Module BuildDElastic
	!**********************************************************************
	!    Should this be a module 
	!
	!     $Revision: 8842 $
	!     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
	!
	!**********************************************************************       

    !**********************************************************************
    !
    !    FUNCTION:  BuildDElastic
    ! 
    !    DESCRIPTION:
    !>   To form the elastic material stiffness matrix (Hooke)
    !
    !>   @param[in] G : Shear modulus
    !    @param[in] xNu : Poisson ratio
    !    @param[out] D : Resulting matrix
    !
    !
    !                             | D1  D2  D2 o  o  o |          | D1 D2 o  o |
    !     Structure of      in 3D | D2  D1  D2 o  o  o |   in 2D  | D2 D1 o  o |
    !     elastic D matrix        | D2  D2  D1 o  o  o |          | o  o  D1 o |
    !                             | o   o   o  G  o  o |          | o  o  o  G | 
    !                             | o   o   o  o  G  o |
    !                             | o   o   o  o  o  G |
    !**********************************************************************  
    Subroutine BuildDElastic(G,xNu,D) 

      use ModGlobalConstants
      
      implicit none
      
      real(REAL_TYPE), intent(in) :: G, xNu
      real(REAL_TYPE), dimension(NTENSOR, NTENSOR) :: D
      ! local variables
      real(REAL_TYPE) :: FAC, D1, D2
      
      FAC= 2*G / (1D0 - 2*xNu)
      D1 = FAC * (1D0 - xNu)
      D2 = FAC * xNu

      D = 0.0 
      
      select case(NTENSOR)
      
        case(4)
          D(1,1)=D1
          D(2,2)=D1
          D(3,3)=D1
          D(1,2)=D2
          D(2,1)=D2
          D(4,4)=G
      
        case(6)
          D(1:3,1:3)=D2
          D(1,1)=D1
          D(2,2)=D1
          D(3,3)=D1
          D(4,4)=G
          D(5,5)=G
          D(6,6)=G
      
        end select
 
      end subroutine BuildDElastic
      
    !**********************************************************************
    !
    !    FUNCTION:  BuildDElasticInverse
    ! 
    !    DESCRIPTION:
    !>   To form the inverse of the elastic material stiffness matrix (Hooke)
    !
    !>   @param[in] G : Shear modulus
    !    @param[in] xNu : Poisson ratio
    !    @param[out] DI : Resulting matrix
    !
    !**********************************************************************  
    Subroutine BuildDElasticInverse(G,xNu,DI)
      use ModGlobalConstants
      
      implicit none
      
      real(REAL_TYPE), intent(in) :: G, XNU
      real(REAL_TYPE), dimension(NTENSOR, NTENSOR) :: DI
      ! local variables
      real(REAL_TYPE) :: E, D1, D2, D3
      
      E = 2d0*G*(1D0+xNu)
      D1 = 1d0/E
      D2 = -xNu/E
      D3 = 1d0/G
      
      DI = 0.0
      select case(NTENSOR)
      
        case(4)
          DI(1,1)=D1
          DI(2,2)=D1
          DI(3,3)=D1
          DI(1,2)=D2
          DI(2,1)=D2
          DI(4,4)=D3
          
        case(6)
          DI(1:3,1:3)=D2
          DI(1,1)=D1
          DI(2,2)=D1
          DI(3,3)=D1
          DI(4,4)=D3
          DI(5,5)=D3
          DI(6,6)=D3
      
      end select

      end subroutine BuildDElasticInverse