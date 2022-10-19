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


      ! Module BuildBJacDet
      !**********************************************************************
      !
      !     $Revision: 8842 $
      !     $Date: 2020-07-30 07:58:40 -0400 (Thu, 30 Jul 2020) $
      !
      !**********************************************************************
      
      subroutine FormB3(Int, IEl, ICon, Co, B, Det, WtN, DShapeValuesArray, maxParticle, maxEl, GetParticleIndex)
!********************************************************************
!
!    Function:  
!
!    Int : Local integration point number 
!      --> note that this is only one in the code and we never actually evaluate this for more than one gauss point 
!    IEl : Element number
!    ICon : Connectivities     
!    Co : Coordinates
!    B : B-matrix
!    Det : Determinant of Jacobi matrix
!    WtN : Weight = Local weight*Det    
!
!********************************************************************
      use ModCounters
      use ModElementEvaluation
      use ModString
      use ModFeedback
      use ModGlobalConstants
      !use ModMeshInfo
      !use ModMPMData
      !use ModAdjustParticleDiscretisation
            
      implicit none
      
      integer(INTEGER_TYPE), intent(in) :: Int
      integer(INTEGER_TYPE) :: IEl, ICon(ELEMENTNODES, Counters%NEL) 
      real(REAL_TYPE) :: Co(Counters%NodTot, NDIM), B(NDIM, ELEMENTNODES), Det, WtN !-> Det is not an input, it is calculated in this routine 
      ! local variables
      integer(INTEGER_TYPE) :: I, J, K, NN
      real(REAL_TYPE) :: RJAC(NDIM, NDIM), RJAC1(NDIM, NDIM), DET1, H
      real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
      
      integer(INTEGER_TYPE) :: ParticleIndex, IPoint
      integer(INTEGER_TYPE) :: NumberOfParticlesForInputElement
      
      integer(INTEGER_TYPE), intent(in) :: maxParticle, maxEl
      
      integer(INTEGER_TYPE), dimension(maxParticle, maxEl), intent(in) :: GetParticleIndex !hardcoded
      
      ! allocation of variables 
      real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
      real(REAL_TYPE), dimension(NDIM,NDIM) :: dx_dtildexi

      integer(INTEGER_TYPE) :: ni, nj  
      
      integer(INTEGER_TYPE) :: counter 
      
      counter = 1
      
       !integer(INTEGER_TYPE), intent(in)  :: NElemPart 
      !(INTEGER_TYPE)
      !maxParticle, maxEl
      ! There is a bug in this part of the code. 
      ! Note that we need to use the derivative evaluated at the material point and not at the original gauss point.  
      ! 
      
      ! evaluate shape function derivative here 
      !call DetJacob(Particles(ParticleIndex)%LocPos, Counters%NEl, Counters%NodTot, NVECTOR, IElement, ElementConnectivities, NodalCoordinates, RJac, InvRJac, DetJacUndeformed)
      
      !use the shape function derivative in the below subroutine 

      !do ii = 1, NElemPart ! loop over particles in the element 
          
      ParticleIndex = GetParticleIndex(Int, IEl) !Get Particle Index... Material Point number 1 locally 
      !NumberOfParticlesForInputElement = NumberOfIntegrationPoints(IEl)
      RJac = 0.0 
      do K = 1, ELEMENTNODES
        NN = ICon(K, IEl)
        do I = 1, NDIM
          do J = 1, NDIM
            !RJac(I,J) = RJac(I,J) + GPShapeFunctionDerivative(Int,K,I) * Co(NN,J)
          
              
              !RJac(I,J) = RJac(I,J) + GPShapeFunctionDerivative_AllElements(IEl,Int,K,I) * Co(NN,J)
              
              !if ( (CalParams%TimeStep ==1) .and. (CalParams%IStep ==1) ) then 
              !    
              !    RJac(I,J) = RJac(I,J) + GPShapeFunctionDerivative_AllElements(IEl,Int,K,I) * Co(NN,J)
              !    
              !
              !else 
              
              ! which particle for which element are we considering?
              ! 
              
                  RJac(I,J) = RJac(I,J) + ( DShapeValuesArray(ParticleIndex, K, I) * Co(NN,J) ) !update Jacob -> dx/dxi

            !RJac(I,J) = RJac(I,J) !+ ( DShapeValuesArray(1, K, I) * Co(NN,J) )
          
          !GPShapeFunctionDerivative_AllElements(IEl,Int,K,I)
          
              !end if 
              
          end do
        end do
      end do
      
      !if (RJac(1,1) == 0.0) then 
      !    return 
      !end if 
      
      
      

      ! Invert Jacobi-matrix
      call RJacInv(NDIM, RJac, RJac1, Det, Det1) !update Jacob -> dxi/dx
      counter = counter +1
      if (Det < 0) Then
        call WriteInLogFile('iEl : ' // trim(String(iEl)))
        call WriteInLogFile('int : ' // trim(String(int)))
        call WriteInLogFile('Error DET<0')
        do K = 1, ELEMENTNODES
          NN = iCon(K, IEl)
          do J = 1, NDIM
            call WriteInLogFile(trim(String(nn)) //' '// trim(String(co(nn,j))))
          end do
        end Do
        call GiveError('Determinant less than zero. [subroutine FormB3()].')
      end If

      do K = 1, ELEMENTNODES
        do I = 1, NDIM
          H = 0.0
          do J = 1, NDIM
              
              
              !if ( (CalParams%TimeStep ==1) .and. (CalParams%IStep ==1) ) then 
              !
              !
              !    H = H + RJac1(I,J) * GPShapeFunctionDerivative_AllElements(IEl,Int,K,J)
              !
              !else 
                  
                  H = H + RJac1(I,J) * DShapeValuesArray(ParticleIndex, K, J) !GPShapeFunctionDerivative_AllElements(IEl,Int,K,J) !update Jacob --> no need to change anything -> dR/dx

            !H = H + RJac1(I,J) * DShapeValuesArray(1, K, I) !GPShapeFunctionDerivative_AllElements(IEl,Int,K,J)
          
              !end if
              
              
              
          end do
          B(I,K) = H
        end do
      end do

      !!!!!!!!!!!!!!!!!!! PARAMETRIC TO PARENT DOMAIN !!!!!!!!!!!!!!!!!!!


      ! based on element number, IEl --> find ni, nj, nk  
      ! need to calculate dxi_dtildexi here... but we need to input ni, nj, nk 

      ni = INN(IEN(1,IEl),1)
      nj = INN(IEN(2,IEl),2)

      call Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries, XiKnotEntries, ni, &
                                                                         NEtaKnotEntries, EtaKnotEntries, nj, &   
                                                                         dxi_dtildexi)
          
      call Get_dx_dtildexi_MappingFromParentElementToParameterSpace(RJac, dxi_dtildexi, dx_dtildexi) ! -> RJac = dx/dxi * dxi/dtildexi
          
      ! reassign the jacobian here with respect to the parent coordinates 
      RJac = dx_dtildexi 
      
               
      !if (RJac(1,1) == 0.0) then 
      !    
      !
      !    return 
      !        
      !end if 
          
      ! recalculate jacobian determinant -> recalculate the deteminant and its inverse 
      call RJacInv(NDIM, RJac, RJac1, Det, Det1) ! -> we will not need the inverse from here 

      !call RJacInv(NDIM, RJac, RJac1, Det, Det1) !update Jacob
      !!!!!!!!!!!!!!!!!!! PARAMETRIC TO PARENT DOMAIN !!!!!!!!!!!!!!!!!!!
      
      
      
      WtN = GPWeight(1) * Det !hardcoded weight of Material Point !update Jacob
      
      !WtN = GPWeight(Int) * Det
      
      !end do 
      
      ! hardcoded... need a better generalization 
      ! need to multiply the first row with dxi/dxi_tilde 
      ! need to multiply the second row with deta/deta_tilde 
      !
      !B(1,:) = B(1,:)*dxi_dtildexi(1,1)
      !B(2,:) = B(2,:)*dxi_dtildexi(2,2)
      

    end subroutine FormB3
    
    
    !integer(INTEGER_TYPE) function NumberOfIntegrationPoints(IElement)
    !    !**********************************************************************
    !    !
    !    !    Function: Returns the number of particles in element IElement if the element is 
    !    !              partially filled and the number of Gauss points per element if the 
    !    !              element is fully filled.
    !    !
    !    !     IElement : ID of the considered element
    !    !
    !    ! O   NumberOfIntegrationPoints : Number of Gauss points or particles
    !    !
    !    ! Implemented in the frame of the MPM project.
    !    !
    !    !**********************************************************************
    ! 
    !    implicit none
    !    
    !      integer(INTEGER_TYPE) :: IElement
    !    
    !      if (IsParticleIntegration(IElement) ) then ! Partially filled element
    !        NumberOfIntegrationPoints = NPartEle(IElement) ! Number of particles in IElement
    !      else ! Fully filled element
    !        NumberOfIntegrationPoints = ELEMENTGAUSSPOINTS ! Number of Gauss points in element
    !      end if
    !    
    !    end function NumberOfIntegrationPoints
 

      subroutine RJacInv(IDimJ, RJac, RJac1, Det, Det1) 
!********************************************************************
!
!    Function:  Find inverse of Jacobi matrix, either 2*2 (2D) or 3*3 (3D) matrix
!
! I   IDimJ      : Dimension of Jacobian matrix (2 = 2D, 3 = 3D)
! I   RJac(i,j)  : Jacobi matrix = dXGi/dXLj
!                                   XG : Global coord. system
!                                   XL : Local  coord. system
! O   RJac1(i,j) : Inverse of Jacobi matrix (= dXLi/dXGj)
! O   Det        : Determinant of Jacobi matrix
! O   Det1       : Determinant of inverse of Jacobi matrix
!
!********************************************************************
      use ModString
      use ModFeedback
      use ModGlobalConstants
      
      implicit none
      
      integer(INTEGER_TYPE), intent(in) :: IDimJ
      real(REAL_TYPE), dimension(IDimJ, IDimJ), intent(in) :: RJac
      real(REAL_TYPE), dimension(IDimJ, IDimJ), intent(out) :: RJac1
      real(REAL_TYPE), intent(out) :: Det, Det1

      Det = 0.0 
        
      select case (IDimJ)
        
        case(2)  
          Det = RJac(1,1) * RJac(2,2) - RJac(1,2) * RJac(2,1) !determinant of the jacobian matrix for 2D case
          
          if (det < SMALL) Then
            call WriteInLogFile('det:' // trim(String(det)))
            call WriteInLogFile('Error DET<0')
          end If
          Det1 = 1/Det
          !calculating jacobian inverse for 2D case
          RJac1(1,1) =   RJac(2,2) * Det1
          RJac1(1,2) = - RJac(1,2) * Det1
          
          RJac1(2,1) = - RJac(2,1) * Det1
          RJac1(2,2) =   RJac(1,1) * Det1
          
        case(3)
          Det =       RJac(1,1) * (RJac(2,2) * RJac(3,3) - RJac(3,2) * RJac(2,3))
          Det = Det - RJac(1,2) * (RJac(2,1) * RJac(3,3) - RJac(3,1) * RJac(2,3))
          Det = Det + RJac(1,3) * (RJac(2,1) * RJac(3,2) - RJac(3,1) * RJac(2,2))
       
          if (det < SMALL) Then
            call WriteInLogFile('det:' // trim(String(det)))
            call WriteInLogFile('Error DET<0')
          end If
          Det1 = 1/Det

          RJac1(1,1) =   (RJac(2,2) * RJac(3,3) - RJac(3,2) * RJac(2,3)) * Det1
          RJac1(2,1) = - (RJac(2,1) * RJac(3,3) - RJac(3,1) * RJac(2,3)) * Det1
          RJac1(3,1) =   (RJac(2,1) * RJac(3,2) - RJac(3,1) * RJac(2,2)) * Det1
          
          RJac1(1,2) = - (RJac(1,2) * RJac(3,3) - RJac(3,2) * RJac(1,3)) * Det1
          RJac1(2,2) =   (RJac(1,1) * RJac(3,3) - RJac(3,1) * RJac(1,3)) * Det1
          RJac1(3,2) = - (RJac(1,1) * RJac(3,2) - RJac(3,1) * RJac(1,2)) * Det1
          
          RJac1(1,3) =   (RJac(1,2) * RJac(2,3) - RJac(2,2) * RJac(1,3)) * Det1
          RJac1(2,3) = - (RJac(1,1) * RJac(2,3) - RJac(2,1) * RJac(1,3)) * Det1
          RJac1(3,3) =   (RJac(1,1) * RJac(2,2) - RJac(2,1) * RJac(1,2)) * Det1
          
        case default
          call GiveError('Dimension not defined in subroutine RJacInv().')
          
        end select  
      
      end subroutine RJacInv

    
      real(REAL_TYPE) function GetElementDeterminant(IEl, ICon, Co, DShapeValuesArray) result(Det)
      !**********************************************************************
      !
      !    Function:  returns determinant of a specific element (IEl)
      !     
      !  I  IEl : Element number
      !  I  ICon : Connectivities
      !  I  Co : Coordinates
      !  O  Det
      !
      !**********************************************************************
      use ModCounters
      use ModElementEvaluation
      use ModGlobalConstants

      implicit none

      integer(INTEGER_TYPE), intent(in) :: IEl
      integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
      real(REAL_TYPE), dimension(:, :), intent(in) :: co
      ! local variables
      integer(INTEGER_TYPE) :: I, J, K, NN
      real(REAL_TYPE), dimension(NDIM, NDIM) :: RJac
      integer(INTEGER_TYPE), parameter :: int = 1
      
      real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray

      ! allocation of variables
      real(REAL_TYPE), dimension(NDIM,NDIM) :: dxi_dtildexi
      real(REAL_TYPE), dimension(NDIM,NDIM) :: dx_dtildexi

      integer(INTEGER_TYPE) :: ni, nj
      
      

      RJac = 0.0
      do K = 1, ELEMENTNODES
        NN = ICon(k, IEl)
        do I = 1, NDIM
          do J = 1, NDIM
            !RJac(i,j) = RJac(i,j) + GPShapeFunctionDerivative(Int,K,I) * Co(NN,J)
              
              !real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(in) :: DShapeValuesArray
              !
              !(1, K, I)
            
              RJac(i,j) = RJac(i,j) + DShapeValuesArray(IEl, K, I) * Co(NN,J) !update Jacob -> dx/dxi = x^e * dR/dxi
            

            !RJac(i,j) = RJac(i,j) + GPShapeFunctionDerivative_AllElements(IEl,Int,K,I) * Co(NN,J)

            ! this needs to use DShapeValuesArray instead of GPShapeFunctionDerivative_AllElements(IEl, Int, K, I)
          end do
        end do
      end do

      Det = 0.0
      
      
      !!!!!!!!!!!!!!!!!!! PARAMETRIC TO PARENT DOMAIN !!!!!!!!!!!!!!!!!!!


      ! based on element number, IEl --> find ni, nj, nk
      ! need to calculate dxi_dtildexi here... but we need to input ni, nj, nk

      ni = INN(IEN(1,IEl),1)
      nj = INN(IEN(2,IEl),2)

      call Get_dxi_dtildexi_MappingFromParentElementToParameterSpace(NXiKnotEntries, XiKnotEntries, ni, &
                                                                         NEtaKnotEntries, EtaKnotEntries, nj, &
                                                                         dxi_dtildexi)

      call Get_dx_dtildexi_MappingFromParentElementToParameterSpace(RJac, dxi_dtildexi, dx_dtildexi) ! -> RJac = dxi/dtildexi * dx/dxi

      ! reassign the jacobian here with respect to the parent coordinates
      RJac = dx_dtildexi ! -> reassign as Jacobain 

      !! recalculate jacobian determinant -> recalculate the deteminant and its inverse
      !call RJacInv(NDIM, RJac, RJac1, Det, Det1)

      !call RJacInv(NDIM, RJac, RJac1, Det, Det1) !update Jacob
      !!!!!!!!!!!!!!!!!!! PARAMETRIC TO PARENT DOMAIN !!!!!!!!!!!!!!!!!!!
      
      
      select case(NDIM)
        
        case(2)  
          Det = RJac(1,1) * RJac(2,2) - RJac(1,2) * RJac(2,1)
          
        case(3)
          Det =       RJac(1,1) * (RJac(2,2)*RJac(3,3) - RJac(3,2)*RJac(2,3))
          Det = Det - RJac(1,2) * (RJac(2,1)*RJac(3,3) - RJac(3,1)*RJac(2,3))
          Det = Det + RJac(1,3) * (RJac(2,1)*RJac(3,2) - RJac(3,1)*RJac(2,2))
          
        case default
          call GiveError('Dimension not defined. [function GetElementDeterminant()].')
          
      end select

      end function GetElementDeterminant
