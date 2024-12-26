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
    !	Copyright (C) 2024  Members of the Anura3D MPM Research Community 
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

  
  

    subroutine Get_Strain(IEl, IPoint, ICon, B, Disp, NDofEx, Eps, IPatch)
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
      !integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%Nel), intent(in) :: ICon    
      !integer(INTEGER_TYPE), dimension(ELEMENTNODES, nel_NURBS(IPatch)), intent(in) :: ICon
      integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%NEl(IPatch)), intent(in) :: ICon


      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: B
      real(REAL_TYPE), dimension(Counters%N), intent(in) :: Disp
      integer(INTEGER_TYPE), dimension(NumberOfGlobalControlPointsUniqueMultipatch+1), intent(in) :: NDofEx !Counters%NodTot
      real(REAL_TYPE), dimension(NTENSOR), intent(out) :: Eps
      
      ! local variables
      integer(INTEGER_TYPE) :: J, NN, ParticleIndex
      real(REAL_TYPE) :: Position
      real(REAL_TYPE), dimension(NVECTOR) :: U
      real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
      
      ! Multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary
      integer(INTEGER_TYPE), intent(in) :: IPatch
      
      Eps = 0.0

      select case(NDIM)
      
      case(2)
          
          ! Axisymmetry is a 2D problem
          if ( ISAXISYMMETRIC ) then
            if ( IsParticleIntegration(IEl, IPatch) ) then ! MP-integration
              ParticleIndex = GetParticleIndex(IPoint, IEl, IPatch) ! MP global ID
              Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
              ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
            else ! GP-integration
              Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
              ShapeValues(:) = GPShapeFunction(IPoint, :) ! where is this GP shape function declared 
            end if
          end if  

          do J = 1, ELEMENTNODES 
            NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
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
          
            
          do J = 1, ELEMENTNODES!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES 
            NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            U(3) = Disp(NDofEx(NN)+3)
            Eps(1) = Eps(1) + ( B(1,J) * U(1) )                      ! Eps_XX = dUx/dX
            Eps(2) = Eps(2) + ( B(2,J) * U(2) )                      ! Eps_YY = dUy/dY
            Eps(3) = Eps(3) + ( B(3,J) * U(3) )                      ! Eps_ZZ = dUz/dZ
            Eps(4) = Eps(4) + ( B(2,J) * U(1) ) + ( B(1,J) * U(2) )  ! Gam_XY = dUx/dY+dUy/dX
            Eps(5) = Eps(5) + ( B(3,J) * U(2) ) + ( B(2,J) * U(3) )  ! Gam_YZ = dUy/dZ+dUz/dY
            Eps(6) = Eps(6) + ( B(3,J) * U(1) ) + ( B(1,J) * U(3) )  ! Gam_ZX = dUz/dX+dUx/dZ
          end do  
          
      end select
          
    end subroutine Get_Strain

    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    subroutine Get_Strain_BBar(IEl, IPoint, ICon, B, BBar, Disp, NDofEx, Eps, IPatch)
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
      !integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%Nel), intent(in) :: ICon    
      !integer(INTEGER_TYPE), dimension(ELEMENTNODES, nel_NURBS(IPatch)), intent(in) :: ICon
      integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%NEl(IPatch)), intent(in) :: ICon


      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: B
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: BBar
      real(REAL_TYPE), dimension(Counters%N), intent(in) :: Disp
      integer(INTEGER_TYPE), dimension(NumberOfGlobalControlPointsUniqueMultipatch+1), intent(in) :: NDofEx !Counters%NodTot
      real(REAL_TYPE), dimension(NTENSOR), intent(out) :: Eps
      
      ! local variables
      integer(INTEGER_TYPE) :: J, NN, ParticleIndex
      real(REAL_TYPE) :: Position
      real(REAL_TYPE), dimension(NVECTOR) :: U
      real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
      REAL(REAL_TYPE) :: NDIM_Factor
      
      ! Multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary
      integer(INTEGER_TYPE), intent(in) :: IPatch
      
      Eps = 0.0

      
      NDIM_Factor = REAL(1)/REAL(NDIM)
      
      
      select case(NDIM)
      
      case(2)
          
          ! Axisymmetry is a 2D problem
          if ( ISAXISYMMETRIC ) then
            if ( IsParticleIntegration(IEl, IPatch) ) then ! MP-integration
              ParticleIndex = GetParticleIndex(IPoint, IEl, IPatch) ! MP global ID
              Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
              ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
            else ! GP-integration
              Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
              ShapeValues(:) = GPShapeFunction(IPoint, :) ! where is this GP shape function declared 
            end if
          end if  

          do J = 1, ELEMENTNODES 
            NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            
            Eps(1) = Eps(1) + ((B(1,J) + (NDIM_Factor)*(BBar(1,J)-B(1,J))) * U(1))    &
                            + ((         (NDIM_Factor)*(BBar(2,J)-B(2,J))) * U(2))    ! Eps_XX = dUx/dX
            
            Eps(2) = Eps(2) + ((         (NDIM_Factor)*(BBar(1,J)-B(1,J))) * U(1))    &
                            + ((B(2,J) + (NDIM_Factor)*(BBar(2,J)-B(2,J))) * U(2))
            
            !Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
            if ( ISAXISYMMETRIC ) then
              Eps(3) = Eps(3) + U(1) * ShapeValues(J) / Position ! Eps_tt = u_r*N/r
            else               
              Eps(3) = 0.0                                     ! E_zz = 0 in 2D
            end if
            Eps(4) = Eps(4) + (B(2,J) * U(1)) + (B(1,J) * U(2))  ! Gam_XY = dUx/dY+dUy/dX
          end do   
            
      case(3)
          
            
          do J = 1, ELEMENTNODES!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES 
            NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            U(3) = Disp(NDofEx(NN)+3)
            Eps(1) = Eps(1) + ( (B(1,J) + (NDIM_Factor)*(BBar(1,J)-B(1,J) )) * U(1) )  &
                            + ( (         (NDIM_Factor)*(BBar(2,J)-B(2,J) )) * U(2) )  &
                            + ( (         (NDIM_Factor)*(BBar(3,J)-B(3,J) )) * U(3) )  ! Eps_XX = dUx/dX
            
            
            Eps(2) = Eps(2) + ( (         (NDIM_Factor)*(BBar(1,J)-B(1,J) )) * U(1) )  &
                            + ( (B(2,J) + (NDIM_Factor)*(BBar(2,J)-B(2,J) )) * U(2) )  &
                            + ( (         (NDIM_Factor)*(BBar(3,J)-B(3,J) )) * U(3) )  
            !Eps(2) = Eps(2) + ( B(2,J) * U(2) )                      ! Eps_YY = dUy/dY
            
            
            Eps(3) = Eps(3) + ( (         (NDIM_Factor)*(BBar(1,J)-B(1,J) )) * U(1) )  &
                            + ( (         (NDIM_Factor)*(BBar(2,J)-B(2,J) )) * U(2) )  &
                            + ( (B(3,J) + (NDIM_Factor)*(BBar(3,J)-B(3,J) )) * U(3) )  
            
            
            !Eps(3) = Eps(3) + ( B(3,J) * U(3) )                      ! Eps_ZZ = dUz/dZ
            Eps(4) = Eps(4) + ( B(2,J) * U(1) ) + ( B(1,J) * U(2) )  ! Gam_XY = dUx/dY+dUy/dX
            Eps(5) = Eps(5) + ( B(3,J) * U(2) ) + ( B(2,J) * U(3) )  ! Gam_YZ = dUy/dZ+dUz/dY
            Eps(6) = Eps(6) + ( B(3,J) * U(1) ) + ( B(1,J) * U(3) )  ! Gam_ZX = dUz/dX+dUx/dZ
          end do  
          
      end select
          
    end subroutine Get_Strain_BBar

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    subroutine Get_RelativeDeformationGradient(IEl, IPoint, ICon, B, Disp, NDofEx, RelativeDeformationGradient, IPatch)
    !-----------------------------------------------------------------
    ! calculate relative deformation gradient from the velocity field. 
    !-----------------------------------------------------------------
      use ModCounters
      use ModMeshInfo
      use ModMPMData
      
      implicit none
      
      integer(INTEGER_TYPE), intent(in) :: IEl, IPoint
      integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%NEl(IPatch)), intent(in) :: ICon

      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: B
      real(REAL_TYPE), dimension(Counters%N), intent(in) :: Disp
      integer(INTEGER_TYPE), dimension(Counters%NodTot(IPatch)+1), intent(in) :: NDofEx
      real(REAL_TYPE), dimension(NDIM,NDIM), intent(out) :: RelativeDeformationGradient
      
      ! local variables
      integer(INTEGER_TYPE) :: J, NN, ParticleIndex
      real(REAL_TYPE) :: Position
      real(REAL_TYPE), dimension(NVECTOR) :: U
      real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
      
      ! Multipatch variables 
      integer(INTEGER_TYPE) :: IPatch_Temporary
      integer(INTEGER_TYPE), intent(in) :: IPatch 
      integer(INTEGER_TYPE) :: IDim
    
      RelativeDeformationGradient = 0.0
      
      do IDim = 1, NDIM 
          RelativeDeformationGradient(IDim,IDim) = 1.0
      end do 
      

      select case(NDIM)
      
      case(2)
          
          ! Axisymmetry is a 2D problem
          if ( ISAXISYMMETRIC ) then
              print *, "Error: F-bar is not implemented for axisymmetric problems." 
            !if ( IsParticleIntegration(IEl, IPatch) ) then ! MP-integration
            !  ParticleIndex = GetParticleIndex(IPoint, IEl, IPatch) ! MP global ID
            !  Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
            !  ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
            !else ! GP-integration
            !  Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
            !  ShapeValues(:) = GPShapeFunction(IPoint, :) ! where is this GP shape function declared 
            !end if
          end if  

          do J = 1, ELEMENTNODES 

            NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            
            RelativeDeformationGradient(1,1) = RelativeDeformationGradient(1,1) + B(1,J) * U(1)                  ! F_XX = dUx/dX
            RelativeDeformationGradient(1,2) = RelativeDeformationGradient(1,2) + B(2,J) * U(1)                  ! F_XY = dUx/dY
            
            RelativeDeformationGradient(2,1) = RelativeDeformationGradient(2,1) + B(1,J) * U(2)                  ! F_YX = dUy/dX
            RelativeDeformationGradient(2,2) = RelativeDeformationGradient(2,2) + B(2,J) * U(2)                  ! F_YY = dUy/dY
          
          end do   
            
      case(3)
          
            
          do J = 1, ELEMENTNODES
            NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
            U(1) = Disp(NDofEx(NN)+1)
            U(2) = Disp(NDofEx(NN)+2)
            U(3) = Disp(NDofEx(NN)+3)
            
            !Eps(1) = Eps(1) + ( B(1,J) * U(1) )                      ! Eps_XX = dUx/dX
            !Eps(2) = Eps(2) + ( B(2,J) * U(2) )                      ! Eps_YY = dUy/dY
            !Eps(3) = Eps(3) + ( B(3,J) * U(3) )                      ! Eps_ZZ = dUz/dZ
            !Eps(4) = Eps(4) + ( B(2,J) * U(1) ) + ( B(1,J) * U(2) )  ! Gam_XY = dUx/dY+dUy/dX
            !Eps(5) = Eps(5) + ( B(3,J) * U(2) ) + ( B(2,J) * U(3) )  ! Gam_YZ = dUy/dZ+dUz/dY
            !Eps(6) = Eps(6) + ( B(3,J) * U(1) ) + ( B(1,J) * U(3) )  ! Gam_ZX = dUz/dX+dUx/dZ
            
            RelativeDeformationGradient(1,1) = RelativeDeformationGradient(1,1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
            RelativeDeformationGradient(1,2) = RelativeDeformationGradient(1,2) + B(2,J) * U(1)                  ! Eps_YY = dUy/dY
            RelativeDeformationGradient(1,3) = RelativeDeformationGradient(1,3) + B(3,J) * U(1)                  ! Eps_YY = dUy/dY
            
            RelativeDeformationGradient(2,1) = RelativeDeformationGradient(2,1) + B(1,J) * U(2)                  ! Eps_YY = dUy/dY
            RelativeDeformationGradient(2,2) = RelativeDeformationGradient(2,2) + B(2,J) * U(2)  
            RelativeDeformationGradient(2,3) = RelativeDeformationGradient(2,3) + B(3,J) * U(2)  
            
            RelativeDeformationGradient(3,1) = RelativeDeformationGradient(3,1) + B(1,J) * U(3)                  ! Eps_YY = dUy/dY
            RelativeDeformationGradient(3,2) = RelativeDeformationGradient(3,2) + B(2,J) * U(3)  
            RelativeDeformationGradient(3,3) = RelativeDeformationGradient(3,3) + B(3,J) * U(3)  
          
          end do  
          
      end select
    
    
    
    end subroutine Get_RelativeDeformationGradient
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !subroutine Get_Strain_VolLockSmooth(IEl, IPoint, ICon, B, Disp, NDofEx, Eps, IPatch)
    !!*************************************************************************************   
    !!    subroutine:     Get_Strain
    !! 
    !!    DESCRIPTION:        
    !!>   Calculates the strains in an integration point
    !!
    !!>   @param[in] IEl : Element number
    !!    @param[in] IMP : Material point number
    !!    @param[in] Icon : Element conectivities
    !!    @param[in] B : B matrix
    !!    @param[in] disp : Displacements of global dof
    !!    @param[in] NDofEx : Reduced dof
    !!
    !!>   @param[out] Eps : Strain
    !!
    !!*************************************************************************************
    !  use ModCounters
    !  use ModMeshInfo
    !  use ModMPMData
    !  
    !  implicit none
    !  
    !  integer(INTEGER_TYPE), intent(in) :: IEl, IPoint
    !  !integer(INTEGER_TYPE), dimension(ELEMENTNODES, Counters%Nel), intent(in) :: ICon    
    !  !integer(INTEGER_TYPE), dimension(ELEMENTNODES, nel_NURBS(IPatch)), intent(in) :: ICon
    !  integer(INTEGER_TYPE), dimension(ELEMENTNODES_VolLockSmooth, nel_NURBS(IPatch)), intent(in) :: ICon
    !
    !  real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES_VolLockSmooth), intent(in) :: B
    !  real(REAL_TYPE), dimension(Counters%N_VolLockSmooth), intent(in) :: Disp
    !  !integer(INTEGER_TYPE), dimension(Counters%NodTot_VolLockSmooth+1), intent(in) :: NDofEx
    !  integer(INTEGER_TYPE), dimension(NControlPoints_VolLockSmooth(1)+1), intent(in) :: NDofEx ! hardcoded for patch 1
    !
    !  real(REAL_TYPE), dimension(NTENSOR), intent(out) :: Eps
    !  
    !  ! local variables
    !  integer(INTEGER_TYPE) :: J, NN, ParticleIndex
    !  real(REAL_TYPE) :: Position
    !  real(REAL_TYPE), dimension(NVECTOR) :: U
    !  real(REAL_TYPE), dimension(ELEMENTNODES_VolLockSmooth) :: ShapeValues
    !  
    !  ! Multipatch variables 
    !  integer(INTEGER_TYPE) :: IPatch_Temporary
    !  integer(INTEGER_TYPE), intent(in) :: IPatch
    !  
    !  Eps = 0.0
    !
    !  select case(NDIM)
    !  
    !  case(2)
    !      
    !      ! Axisymmetry is a 2D problem
    !      if ( ISAXISYMMETRIC ) then
    !        if ( IsParticleIntegration(IEl, IPatch) ) then ! MP-integration
    !          ParticleIndex = GetParticleIndex(IPoint, IEl, IPatch) ! MP global ID
    !          Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
    !          ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
    !        else ! GP-integration
    !          Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
    !          ShapeValues(:) = GPShapeFunction(IPoint, :) ! where is this GP shape function declared 
    !        end if
    !      end if  
    !
    !      do J = 1, ELEMENTNODES_VolLockSmooth 
    !        !NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
    !        NN = Icon(J, Iel)
    !
    !        U(1) = Disp(NDofEx(NN)+1)
    !        U(2) = Disp(NDofEx(NN)+2)
    !        Eps(1) = Eps(1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
    !        Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
    !        if ( ISAXISYMMETRIC ) then
    !          Eps(3) = Eps(3) + U(1) * ShapeValues(J) / Position ! Eps_tt = u_r*N/r
    !        else               
    !          Eps(3) = 0.0                                     ! E_zz = 0 in 2D
    !        end if
    !        Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
    !      end do   
    !        
    !  case(3)
    !      
    !        
    !      do J = 1, ELEMENTNODES_VolLockSmooth!(NXiKnotOrder(IPatch)+1)*(NEtaKnotOrder(IPatch)+1)*(NZetaKnotOrder(IPatch)+1)!ELEMENTNODES 
    !        !NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
    !        NN = Icon(J, Iel)
    !
    !        U(1) = Disp(NDofEx(NN)+1)
    !        U(2) = Disp(NDofEx(NN)+2)
    !        U(3) = Disp(NDofEx(NN)+3)
    !        Eps(1) = Eps(1) + ( B(1,J) * U(1) )                      ! Eps_XX = dUx/dX
    !        Eps(2) = Eps(2) + ( B(2,J) * U(2) )                      ! Eps_YY = dUy/dY
    !        Eps(3) = Eps(3) + ( B(3,J) * U(3) )                      ! Eps_ZZ = dUz/dZ
    !        Eps(4) = Eps(4) + ( B(2,J) * U(1) ) + ( B(1,J) * U(2) )  ! Gam_XY = dUx/dY+dUy/dX
    !        Eps(5) = Eps(5) + ( B(3,J) * U(2) ) + ( B(2,J) * U(3) )  ! Gam_YZ = dUy/dZ+dUz/dY
    !        Eps(6) = Eps(6) + ( B(3,J) * U(1) ) + ( B(1,J) * U(3) )  ! Gam_ZX = dUz/dX+dUx/dZ
    !      end do  
    !      
    !  end select
    !      
    !end subroutine Get_Strain_VolLockSmooth
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !NodalStressArray
    
    
    
    
    !subroutine Get_StressGradient(IEl, IPoint, ICon, B,  NDofEx, IPatch, MPStressGradient) !Eps, !Disp,
    !!*************************************************************************************   
    !!    subroutine:     Get_StressGradient
    !! 
    !!    DESCRIPTION:        
    !!>   Calculates the stress gradient for an integration point
    !!
    !!>   @param[in] IEl : Element number
    !!    @param[in] IMP : Material point number
    !!    @param[in] Icon : Element conectivities
    !!    @param[in] B : B matrix
    !!    @param[in] disp : Displacements of global dof
    !!    @param[in] NDofEx : Reduced dof
    !!
    !!>   @param[out] stress gradient --> actually no output for stress gradient subroutine as 
    !!                       we pour directly into the nodal array. 
    !!
    !!*************************************************************************************
    !  use ModCounters
    !  use ModMeshInfo
    !  use ModMPMData
    !  
    !  implicit none
    !  
    !  integer(INTEGER_TYPE), intent(in) :: IEl, IPoint
    !  integer(INTEGER_TYPE), dimension(ELEMENTNODES, nel_NURBS(IPatch)), intent(in) :: ICon
    !
    !
    !  real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(in) :: B
    !  integer(INTEGER_TYPE), dimension(Counters%NodTot+1), intent(in) :: NDofEx
    !  
    !  ! local variables
    !  integer(INTEGER_TYPE) :: J, NN, ParticleIndex
    !  real(REAL_TYPE) :: Position
    !  real(REAL_TYPE), dimension(NVECTOR), intent(out) :: MPStressGradient
    !  real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues
    !  
    !  ! Multipatch variables 
    !  integer(INTEGER_TYPE) :: IPatch_Temporary
    !  integer(INTEGER_TYPE) :: aa
    !  integer(INTEGER_TYPE), intent(in) :: IPatch
    !  
    !  ! volumetric locking variables 
    !  real(REAL_TYPE), dimension(NVECTOR) :: Stress
    !  integer(INTEGER_TYPE) :: iEntity_Temporary
    !  
    !  MPStressGradient = 0.0
    !  iEntity_Temporary = 1
    !
    !  select case(NDIM)
    !  
    !  case(2)
    !      
    !      ! TODO: Axisymmetry is a 2D problem
    !      !if ( ISAXISYMMETRIC ) then
    !      !  if ( IsParticleIntegration(IEl, IPatch) ) then ! MP-integration
    !      !    ParticleIndex = GetParticleIndex(IPoint, IEl, IPatch) ! MP global ID
    !      !    Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
    !      !    ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
    !      !  else ! GP-integration
    !      !    Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
    !      !    ShapeValues(:) = GPShapeFunction(IPoint, :) ! where is this GP shape function declared 
    !      !  end if
    !      !end if  
    !
    !
    !      do J = 1, ELEMENTNODES 
    !        
    !          NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
    !        
    !          
    !          ! we only need the volumetric strains but not the shear strains for the 
    !          ! calculation of the divergence of the stress
    !          !U(1) = Disp(NDofEx(NN)+1)
    !          !U(2) = Disp(NDofEx(NN)+2)
    !          
    !          ! SigmaEffArray needs to be replaced by the nodal stress
    !          Stress(1) = NodalStress_VolLockSmooth(NDofEx(NN)+1, iEntity_Temporary) ! Sigma_XX
    !          Stress(2) = NodalStress_VolLockSmooth(NDofEx(NN)+2, iEntity_Temporary) ! Sigma_YY
    !        
    !          if (abs(Stress(1))>0.0 .or. abs(Stress(1))>0.0 ) then 
    !              aa = 1
    !          end if 
    !      
    !          
    !          MPStressGradient(1) = MPStressGradient(1) + ( B(1,J) * Stress(1) )                 ! dSigma_XX/dX
    !          MPStressGradient(2) = MPStressGradient(2) + ( B(2,J) * Stress(2) )                 ! dSigma_YY/dY
    !              
    !        if ( ISAXISYMMETRIC ) then
    !            ! TODO: needs completion for axisymmetric 
    !            !Eps(3) = Eps(3) + U(1) * ShapeValues(J) / Position ! Eps_tt = u_r*N/r
    !        else               
    !            !MPStressGradient(3) = 0.0                                     ! dSigma_ZZ/dZ = 0 in 2D
    !        end if
    !        !Eps(4) = Eps(4) + B(2,J) * U(1) + B(1,J) * U(2)  ! Gam_XY = dUx/dY+dUy/dX
    !      
    !      end do   
    !      
    !      ! Pour the nodal stress gradient in the global array
    !      !MPStressGradient(NDofEx(NN)+1, IPatch) = NodalStressGradientArray(NDofEx(NN)+1, IPatch) + StressGradient(1)
    !      !NodalStressGradientArray(NDofEx(NN)+2, IPatch) = NodalStressGradientArray(NDofEx(NN)+2, IPatch) + StressGradient(2)
    !        
    !  case(3)
    !      
    !        
    !      do J = 1, ELEMENTNODES
    !        NN = Multipatch_Connecting_Local_To_Global_ControlPoints(Icon(J, Iel), IPatch)
    !        
    !        ! we only need the volumetric strains but not the shear strains for the 
    !        ! calculation of the divergence of the stress
    !        Stress(1) = NodalStress_VolLockSmooth(NDofEx(NN)+1, iEntity_Temporary) ! Sigma_XX
    !        Stress(2) = NodalStress_VolLockSmooth(NDofEx(NN)+2, iEntity_Temporary) ! Sigma_YY
    !        Stress(3) = NodalStress_VolLockSmooth(NDofEx(NN)+3, iEntity_Temporary) ! Sigma_YY
    !        
    !        MPStressGradient(1) = MPStressGradient(1) + ( B(1,J) * Stress(1) )                      ! dSigma_XX/dX
    !        MPStressGradient(2) = MPStressGradient(2) + ( B(2,J) * Stress(2) )                      ! dSigma_YY/dY
    !        MPStressGradient(3) = MPStressGradient(3) + ( B(3,J) * Stress(3) )                      ! dSigma_ZZ/dZ
    !        
    !
    !        ! Pour the nodal stress gradient in the global array
    !        !NodalStressGradientArray(NDofEx(NN)+1, IPatch) = NodalStressGradientArray(NDofEx(NN)+1, IPatch) + StressGradient(1)
    !        !NodalStressGradientArray(NDofEx(NN)+2, IPatch) = NodalStressGradientArray(NDofEx(NN)+2, IPatch) + StressGradient(2)
    !        !NodalStressGradientArray(NDofEx(NN)+3, IPatch) = NodalStressGradientArray(NDofEx(NN)+3, IPatch) + StressGradient(3)
    !      end do  
    !      
    !  end select
    !      
    !  end subroutine Get_StressGradient