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
	!
    !   This file incorporates work covered by the following copyright and 
    !   permission notice:
    !  
    !   Copyright © 2017-2020 Bentley Systems, Incorporated. All rights reserved.
    !
    !   Permission is hereby granted, free of charge, to any person obtaining a copy 
    !   of this software and associated documentation files (the "Software"), to deal 
    !   in the Software without restriction, including without limitation the rights 
    !   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
    !   copies of the Software, and to permit persons to whom the Software is furnished 
    !   to do so, subject to the following conditions:
    !
    !   The above copyright notice and this permission notice shall be included in all 
    !   copies or substantial portions of the Software.
    !
    !   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    !   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
    !   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
    !   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
    !   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
    !   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    !
    !*****************************************************************************   
	  
	  
	  module ModMatrixMath
      !**********************************************************************
      !
      !    Function:  Contains routines related to matrix calculations:
      !               - Determining the determinante of a quadratic matrix
      !               - Solving a linear system of equations by means of a Gauss solver
      !
      !    Implemented in the frame of the MPM project.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      
      implicit none

      contains ! Routines of this module

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getEpsV(strain) result (volumetricstrain)
      !********************************************************************
      ! Function: calculate volumetric strain
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: strain(:)

        volumetricstrain = strain(1) + strain(2) + strain(3)
        ! scalar, valid for 1D/2D/3D

      end function getEpsV

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getEpsq(strain) result (deviatoricstrain)
      !********************************************************************
      ! Function: calculate deviatoric strain
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: strain(:)
      integer(INTEGER_TYPE) :: IDim
      
      IDim = size(strain)
      select case(IDim)
        case(6)  
          deviatoricstrain = (sqrt(2.0d0)/3.0d0)*sqrt( (strain(1)-strain(2))**2 + (strain(1)-strain(3))**2 + (strain(2)-strain(3))**2 + &
                                                        6*(strain(4))**2 + 6*(strain(5))**2 + 6*(strain(6))**2 )
        case(4)  
          deviatoricstrain = (sqrt(2.0d0)/3.0d0)*sqrt( (strain(1)-strain(2))**2 + (strain(1)-strain(3))**2 + (strain(2)-strain(3))**2 + &
                                                       6*(strain(4))**2 )
        case default
          call GiveError('Dimension not defined in function getEpsq().')
          
        end select

      end function getEpsq

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getI1(principalvector) result (firstinvariant)
      !********************************************************************
      ! Function: calculate first invariant I1
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: principalvector(:)

        firstinvariant = (principalvector(1) + principalvector(2) + principalvector(3))
        ! scalar, valid for 1D/2D/3D
 
      end function getI1

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getI2(principalvector) result (secondinvariant)
      !********************************************************************
      ! Function: calculate second invariant I2
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: principalvector(:)

        secondinvariant = principalvector(1) * principalvector(2) &
                        + principalvector(2) * principalvector(3) &
                        + principalvector(3) * principalvector(1)
        ! scalar, valid for 1D/2D/3D

      end function getI2

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getI3(principalvector) result (thirdinvariant)
       !********************************************************************
      ! Function: calculate third invariant I3
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: principalvector(:)

        thirdinvariant = (principalvector(1) * principalvector(2) * principalvector(3))
        ! scalar, valid for 1D/2D/3D

      end function getI3

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getp(principalstress) result (meanstress)
       !********************************************************************
      ! Function: calculate mean stress
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in):: principalstress(:)

        meanstress = (principalstress(1) + principalstress(2) + principalstress(3))/3.0d0
        ! scalar, valid for 1D/2D/3D

      end function getp

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getSigma3(principalstress) result (minimumstress)
      !********************************************************************
      ! Function: calculate minimum principal stress
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: principalstress(:)

        minimumstress = minVal(principalstress)
        ! scalar, valid for 1D/2D/3D

      end function getSigma3

      !-----------------------------------------------------------------
      real(REAL_TYPE) function getq(principalstress) result (deviatoricstress)
       !********************************************************************
      ! Function: calculate deviatoric stess
      !********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: principalstress(:)

        deviatoricstress = sqrt( 0.5*( (principalstress(1)-principalstress(2))**2 &
                                     + (principalstress(2)-principalstress(3))**2 &
                                     + (principalstress(1)-principalstress(3))**2) )
        ! scalar, valid for 1D/2D/3D

      end function getq
      
      
      
        subroutine CoFactor(CheckMatrix, MatrixDimension, CoFactorMatrix)
        !****************************************************************************
        !
        !     Function:  Calculate CoFactorMatrix of the quadratric matrix CheckMatrix
        !                 whose dimensions are provided by MatrixDimension.
        !
        !****************************************************************************
        implicit none
      
          integer(INTEGER_TYPE), intent (in) :: MatrixDimension
          real(REAL_TYPE), dimension(MatrixDimension, MatrixDimension), intent (in) :: CheckMatrix
          real(REAL_TYPE), dimension(MatrixDimension, MatrixDimension), intent (out) :: CoFactorMatrix
          
          ! Local variables          
          real(REAL_TYPE), dimension(3, 3) :: MatrixTemp
          real(REAL_TYPE), dimension(3, 3) :: CoFactorMatrixTemp
          integer(INTEGER_TYPE) :: I, J
      
          call Assert(MatrixDimension <= 3, "CoFactor only available for dimension<=3")
          
          ! initialise identity matrix
          MatrixTemp = 0.0
          MatrixTemp(1,1) = 1.0
          MatrixTemp(2,2) = 1.0
          MatrixTemp(3,3) = 1.0
      
          ! overwrite (part of) identity matrix with CheckMatrix
          do i = 1, MatrixDimension
            do j = 1, MatrixDimension
              MatrixTemp(i,j) = CheckMatrix(i,j)
            end do
          end do
          
          CoFactorMatrixTemp(1,1) =   MatrixTemp(2,2) * MatrixTemp(3,3) - MatrixTemp(2,3) * MatrixTemp(3,2)
          CoFactorMatrixTemp(1,2) = - MatrixTemp(2,1) * MatrixTemp(3,3) + MatrixTemp(2,3) * MatrixTemp(3,1)
          CoFactorMatrixTemp(1,3) =   MatrixTemp(2,1) * MatrixTemp(3,2) - MatrixTemp(2,2) * MatrixTemp(3,1)

          CoFactorMatrixTemp(2,1) = - MatrixTemp(1,2) * MatrixTemp(3,3) + MatrixTemp(1,3) * MatrixTemp(3,2)
          CoFactorMatrixTemp(2,2) =   MatrixTemp(1,1) * MatrixTemp(3,3) - MatrixTemp(1,3) * MatrixTemp(3,1)
          CoFactorMatrixTemp(2,3) = - MatrixTemp(1,1) * MatrixTemp(3,2) + MatrixTemp(1,2) * MatrixTemp(3,1)

          CoFactorMatrixTemp(3,1) =   MatrixTemp(1,2) * MatrixTemp(2,3) - MatrixTemp(1,3) * MatrixTemp(2,2)
          CoFactorMatrixTemp(3,2) = - MatrixTemp(1,1) * MatrixTemp(2,3) + MatrixTemp(1,3) * MatrixTemp(2,1)
          CoFactorMatrixTemp(3,3) =   MatrixTemp(1,1) * MatrixTemp(2,2) - MatrixTemp(1,2) * MatrixTemp(2,1)
      
          CoFactorMatrix = 0.0
          do I = 1, MatrixDimension
            do J = 1, MatrixDimension
              CoFactorMatrix(I,J) = CoFactorMatrixTemp(I,J)
            end do
          end do

        end subroutine CoFactor
        
        
        real(REAL_TYPE) function CalculateDeterminant(CheckMatrix, MatrixDimension)
        !**********************************************************************
        !
        !    Function: Returns the determinante of the quadratic matrix CheckMatrix
        !              whose dimensions are provided by MatrixDimension.
        !     
        !     CheckMatrix : Matrix whose determinante is calculated
        !     MatrixDimension : Dimensions of the quadratic matrix CheckMatrix
        !
        ! O   CalculateDeterminant : Determinant of CheckMatrix
        !
        !**********************************************************************
        
        implicit none
        
          integer, intent(in) :: MatrixDimension
          real(REAL_TYPE), dimension (MatrixDimension, MatrixDimension), intent(in) :: CheckMatrix
          
          ! Local variables
          real(REAL_TYPE), dimension(3,3) :: MatrixTemp
          real(REAL_TYPE) :: a, b, c
          integer :: i, j
          
          call Assert(MatrixDimension <= 3, "CalculateDeterminant only available for dimension<=3")
          
          ! initialise identity matrix
          MatrixTemp = 0.0
          MatrixTemp(1,1) = 1.0
          MatrixTemp(2,2) = 1.0
          MatrixTemp(3,3) = 1.0
          
          ! overwrite (part of) identity matrix with CheckMatrix
          do i = 1, MatrixDimension
            do j = 1, MatrixDimension
              MatrixTemp(i,j) = CheckMatrix(i,j)
            end do
          end do
          
          a=0.0
          b=0.0
          c=0.0
          
          a = MatrixTemp(1,1) * (MatrixTemp(2,2) * MatrixTemp(3,3) - MatrixTemp(3,2) * MatrixTemp(2,3))
          b = MatrixTemp(1,2) * (MatrixTemp(2,1) * MatrixTemp(3,3) - MatrixTemp(3,1) * MatrixTemp(2,3))
          c = MatrixTemp(1,3) * (MatrixTemp(2,1) * MatrixTemp(3,2) - MatrixTemp(3,1) * MatrixTemp(2,2))

          CalculateDeterminant = a - b + c

        end function CalculateDeterminant
      
        
      logical function IS0ARR(A, N)
      !-----------------------------------------------
      !
      !  Function: check for zero array
      !
      !  A   I     R()    Real array
      !  N   I     I      Length of real array
          !
      !  Note:  This routine is based on the subroutine MatVec, from Bentley System Inc. 
      !         and is subject to the copyright above based on MIT licensing.   
      !-----------------------------------------------

      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: N
        real(REAL_TYPE), dimension(N), intent(in) :: A

        ! local variables
        integer(INTEGER_TYPE) :: I
        logical :: temp

        ! initialise temp
        temp = .true.
        
        ! check each entry to be zero
        do I = 1,N
          if (A(I) /= 0.0) then
            temp = .false.
            EXIT
          end if
        end do
        
        IS0ARR = temp
      end function IS0ARR
      
      
    end module ModMatrixMath

    
   
      Subroutine MatVec(xMat,IM,Vec,N,VecR) 
    !
    !  Note:  This routine is based on the subroutine MatVec, from Bentley System Inc. 
    !         and is subject to the copyright above based on MIT licensing.   
    ! ******************************************************************************** 
        
      use ModGlobalConstants
      implicit none
      integer(INTEGER_TYPE) :: N, IM, I, J ! N is the dimension
      real(REAL_TYPE) :: xMat(IM,*), Vec(*), VecR(*)
      real(REAL_TYPE) :: X

      Do I=1,N
        X=0.0d0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
    End

    ! ********************************************************************************
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
    !
    !  Note:  This routine is based on the subroutine MatVec, from Bentley System Inc. 
    !         and is subject to the copyright above based on MIT licensing. 
    !  
    ! ******************************************************************************** 
      use ModGlobalConstants
      implicit none 
      integer(INTEGER_TYPE) :: Id1, Id2, nR1, nC1, nC2, IdR ! nC1, nC2, nR1 are the dimensions
      real(REAL_TYPE) :: xMat1(Id1,*), xMat2(Id2,*), xMatR(IdR,*)
      integer(INTEGER_TYPE) :: I, J, K
      real(REAL_TYPE) :: X

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      End 




      subroutine Norm( NNodes, NDof, V1, NDofEx, VNorm )
      !**********************************************
      !
      !  Function: Sums norm of V1 for each node
      !
      !**********************************************
      use ModGlobalConstants
      implicit none
        integer(INTEGER_TYPE), intent(in) :: NNodes, NDof
        integer(INTEGER_TYPE), dimension(NNodes), intent(in) :: NDofEx
        real(REAL_TYPE), dimension(NNodes*NDof), intent(in) :: V1
        real(REAL_TYPE), intent(out) :: VNorm
        
        ! Local variables
        integer(INTEGER_TYPE) :: I, J, Index
        real(REAL_TYPE) :: VI, F, F2
      
        VNorm = 0.0
        do I = 1,NNodes
          VI = 0.0
          do J = 1,NDof
            Index = NDofEx(I) + J
            F = V1(Index)
            F2 = F * F
            VI = VI + F2
          end do
          VI = SQRT(VI)
          VNorm = VNorm + VI
        end do 

    end subroutine Norm
    
   