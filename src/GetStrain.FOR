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


! Module GetStrain
!**********************************************************************
!    Should this be a module 
!
!     $Revision: 8842 $
!     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
!
!**********************************************************************

  
  

    subroutine Get_Strain(IEl, IPoint, ICon, B, Disp, NDofEx, Eps)
    !*************************************************************************************   
    !    subroutine:     Get_Strain
    ! 
    !    DESCRIPTION:        
    !>   Calculates the strains in an integration point
    !
    !>   @param[in] IEl : Element number
    !    @param[in] IMP : Material point number
    !    @param[in] Icon : Element conectivities
    !    @param[in] B : B matrix
    !    @param[in] disp : Displacements of global dof
    !    @param[in] NDofEx : Reduced dof
    !
    !>   @param[out] Eps : Strain
    !
    !*************************************************************************************
      use ModCounters
      use ModMeshInfo
      use ModMPMData
      
      implicit none
      
      integer(INTEGER_TYPE), intent(in) :: IEl, IPoint
      integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%Nel), intent(in) :: ICon
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: B
      real(REAL_TYPE), dimension(Counters%N), intent(in) :: Disp
      integer(INTEGER_TYPE), dimension(Counters%NodTot+1), intent(in) :: NDofEx
      real(REAL_TYPE), dimension(NTENSOR), intent(out) :: Eps
      
      ! local variables
      integer(INTEGER_TYPE) :: J, NN, ParticleIndex
      real(REAL_TYPE) :: Position
      real(REAL_TYPE), dimension(NVECTOR) :: U
      real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
      
      Eps = 0.0

      select case(NDIM)
      
        case(2)
          
          if ( ISAXISYMMETRIC ) then
            if ( IsParticleIntegration(IEl) ) then ! MP-integration
              ParticleIndex = GetParticleIndex(IPoint, IEl) ! MP global ID
              Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
              ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
            else ! GP-integration
              Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
              ShapeValues(:) = GPShapeFunction(IPoint, :)
            end if
          end if  

          do J = 1, ELEMENTNODES 
            NN = Icon(J, Iel)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            Eps(1) = Eps(1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
            Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
            if ( ISAXISYMMETRIC ) then
              Eps(3) = Eps(3) + U(1) * ShapeValues(J) / Position ! Eps_tt = u_r*N/r
            else               
              Eps(3) = 0.0                                     ! E_zz = 0 in 2D
            end if
            Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
          end do   
            
        case(3)
            
          do J = 1, ELEMENTNODES 
            NN = Icon(J, Iel)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            U(3) = Disp(NDofEx(NN)+3)
            Eps(1) = Eps(1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
            Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
            Eps(3) = Eps(3) + B(3,J) * U(3)                  ! Eps_ZZ = dUz/dZ
            Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
            Eps(5) = Eps(5) + B(3,J) * U(2) + B(2,J) * U(3)  ! Gam_YZ = dUy/dZ+dUz/dY
            Eps(6) = Eps(6) + B(3,J) * U(1) + B(1,J) * U(3)  ! Gam_ZX = dUz/dX+dUx/dZ
          end do  
          
      end select
          
      end subroutine Get_Strain
