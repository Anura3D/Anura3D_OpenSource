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
    !   Copyright © 2017-2021 Bentley Systems, Incorporated. All rights reserved.
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
	  
	  ! Module GetPrinStress
      !**********************************************************************
      !
      !     $Revision: 8878 $
      !     $Date: 2020-09-17 03:59:14 -0400 (Thu, 17 Sep 2020) $
      !
      !**********************************************************************
    
      ! Note these modules are based on the third party modules presented in 
    
    
      subroutine PrnSig(IOpt, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions
      !            from cartesian stress vector 
      !
      !  Note:  This routine is based on the subroutine PrnSig, from Bentley System Inc. 
      !         and is subject to the copyright notice above based on MIT licensing. 
      !         
      !
      !  IOpt            I   I     flag to calculate principal direction (IOpt = 1)
      !  IntGlo          I   I     global ID of Gauss point or particle 
      !  S               I   R()   cartesian stress
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------
      use ModGlobalConstants
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IOpt
        real(REAL_TYPE), intent(in), dimension(NTENSOR) :: S ! size of stress tensor (6 for 3D, 4 for 2D)
        real(REAL_TYPE), intent(out), dimension(NPRINCIPAL) :: xN1, xN2, xN3
        real(REAL_TYPE), intent(out) :: S1, S2, S3, P, Q

        if (IOpt == 1) then
          call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! Calculate principal direction
        else
          call Eig_3a(0,S,S1,S2,S3,P,Q) ! Do not calculate principal direction
        end if
     
      end subroutine PrnSig
      
      
      subroutine Eig_3(iOpt, St, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions 
      !            from cartesian stress vector
      !
      !  Note:  This routine is based on the subroutine Eig_3, from Bentley System Inc. 
      !         and is subject to the copyright notice above based on MIT licensing. 
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      ! 
      !  IOpt            I   I     flag for output writing (IOpt = 1) 
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  xN1, xN2, xN3   O   R()   principal direction
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------
      use ModString
      use ModFeedback
      use ModGlobalConstants

      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IOpt
        real(REAL_TYPE), intent(in), dimension(NTENSOR) :: St
        real(REAL_TYPE), intent(out), dimension(NPRINCIPAL) :: xN1, xN2, xN3 ! is (3) for 2D and 3D
        real(REAL_TYPE), intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        real(REAL_TYPE), dimension(NPRINCIPAL, NPRINCIPAL) :: A, V ! is (3,3) for 2D and 3D
        real(REAL_TYPE) :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
        integer(INTEGER_TYPE) :: i, k, it, itmax, ip, iq, iS1, iS2, iS3

        ! Put cartesian stress vector into matrix A
        select case(NTENSOR)
          case(4)
              A = 0.0
              A(1,1) = St(1) ! xx
              A(1,2) = St(4) ! xy = yx
              A(2,1) = St(4) ! xy = yx
              A(2,2) = St(2) ! yy
              A(3,3) = St(3) ! zz
          
              ! Set V to unity matrix
              V = 0.0
              V(1,1) = 1
              V(2,2) = 1
              V(3,3) = 1
              
          case(6)
            A(1,1) = St(1) ! xx
            A(1,2) = St(4) ! xy = yx
            A(1,3) = St(6) ! zx = xz

            A(2,1) = St(4) ! xy = yx
            A(2,2) = St(2) ! yy
            A(2,3) = St(5) ! zy = yz

            A(3,1) = St(6) ! zx = xz
            A(3,2) = St(5) ! zy = yz
            A(3,3) = St(3) ! zz
          
           ! Set V to unity matrix
            V(1,1) = 1
            V(2,1) = 0
            V(3,1) = 0

            V(1,2) = 0
            V(2,2) = 1
            V(3,2) = 0

            V(1,3) = 0
            V(2,3) = 0
            V(3,3) = 1
        end select 
        
        ! get maximum value of cartesian stress vector
        abs_max_s = 0.0
        do i = 1, NTENSOR
          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-16 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it < itmax) .and.  &
                   (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
        
          it = it + 1
          do k = 1,3
                
            if (k == 1) then
              ip = 1
              iq = 2
            else if (k ==2) then
              ip = 2
              iq = 3
            else
              ip = 1
              iq = 3
            end if
            
            if (abs(A(ip,iq)) > 1d-50) then
                
              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
              if (tau >= 0.0) then
                sign_tau = 1.0
              else
                sign_tau = -1.0
              end if
              
              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
              c = 1.0 / sqrt(1.0 + t**2)
              s = t * c
              
              temp1 = c * A(1, ip) - s * A(1, iq)
              temp2 = c * A(2, ip) - s * A(2, iq)
              temp3 = c * A(3, ip) - s * A(3, iq)
              A(1, iq) = s * A(1, ip) + c * A(1, iq)
              A(2, iq) = s * A(2, ip) + c * A(2, iq)
              A(3, iq) = s * A(3, ip) + c * A(3, iq)
              A(1, ip) = temp1
              A(2, ip) = temp2
              A(3, ip) = temp3

              temp1 = c * V(1, ip) - s * V(1, iq)
              temp2 = c * V(2, ip) - s * V(2, iq)
              temp3 = c * V(3, ip) - s * V(3, iq)
              V(1, iq) = s * V(1, ip) + c * V(1, iq)
              V(2, iq) = s * V(2, ip) + c * V(2, iq)
              V(3, iq) = s * V(3, ip) + c * V(3, iq)
              V(1, ip) = temp1
              V(2, ip) = temp2
              V(3, ip) = temp3

              temp1 = c * A(ip, 1) - s * A(iq, 1)
              temp2 = c * A(ip, 2) - s * A(iq, 2)
              temp3 = c * A(ip, 3) - s * A(iq, 3)
              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
              A(ip, 1) = temp1
              A(ip, 2) = temp2
              A(ip, 3) = temp3
            end if ! A(ip,iq)<>0
            
          end do ! k
          
          ! optional output writing
          if (iOpt == 1) then
            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
          end if
          
        end do ! while
     
        ! get principal stresses from diagonal of A
        S1 = A(1, 1)
        S2 = A(2, 2)
        S3 = A(3, 3)
      
        ! derived invariants
        P = (S1 + S2 + S3) / 3.
        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )

        ! Sort eigenvalues S1 <= S2 <= S3
        iS1 = 1
        iS2 = 2
        iS3 = 3
        
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
          it  = iS2
          iS2 = iS1
          iS1 = it
        end if
        
        if (S2 > S3) then
          t   = S3
          S3  = S2
          S2  = t
          it  = iS3
          iS3 = iS2
          iS2 = it
        end if
        
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
          it  = iS2
          iS2 = iS1
          iS1 = it
        end if
        
        ! get corresponding principal directions from V
        do i = 1, NPRINCIPAL ! is 3 for 2D and 3D
          xN1(i) = V(i, is1)
          xN2(i) = V(i, is2)
          xN3(i) = V(i, is3)
        end do
      
      end subroutine Eig_3

      
      subroutine Eig_3a(iOpt, St, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses from cartesian stress vector
      !
      !  Note:  This routine is based on the subroutine Eig_3a, from Bentley System Inc. 
      !         and is subject to the copyright notice above based on MIT licensing. 
      !
      !  NB: Wim Bomhof 15/11/'01, adapted to principal stress calculation
      ! 
      !  IOpt            I   I     flag for output writing (IOpt = 1) 
      !  St              I   R()   cartesian stress (XX, YY, ZZ, XY, YZ, ZX)
      !  S1, S2, S3      O   R     principal stress
      !  P               O   R     isotropic stress (positive for tension)
      !  Q               O   R     deviatoric stress
      ! 
      !-------------------------------------------------------------------
      use ModString
      use ModFeedback
      use ModGlobalConstants
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: IOpt
        real(REAL_TYPE), intent(in), dimension(NTENSOR) :: St ! is (6) for 3D and (4) for 2D
        real(REAL_TYPE), intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        integer(INTEGER_TYPE) :: IDim
        real(REAL_TYPE), dimension(:,:), allocatable :: A ! (3,3) for 2D and 3D
        real(REAL_TYPE) :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
        integer(INTEGER_TYPE) :: i, k, it, itmax, ip, iq
        
        IDim = NPRINCIPAL
        allocate(A(IDim,IDim)) 
        
        select case (NTENSOR)       
        case(4) ! for 2D
          A = 0.0
          A(1,1) = St(1) ! xx
          A(1,2) = St(4) ! xy = yx
          A(2,1) = St(4) ! xy = yx
          A(2,2) = St(2) ! yy
          A(3,3) = St(3) ! zz
        case(6) ! for 3D
          A(1,1) = St(1) ! xx
          A(1,2) = St(4) ! xy = yx
          A(1,3) = St(6) ! zx = xz

          A(2,1) = St(4) ! xy = yx
          A(2,2) = St(2) ! yy
          A(2,3) = St(5) ! zy = yz

          A(3,1) = St(6) ! zx = xz
          A(3,2) = St(5) ! zy = yz
          A(3,3) = St(3) ! zz
        end select

        ! get maximum value of cartesian stress vector
        abs_max_s = 0.0
        do i = 1,NTENSOR
          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-20 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it < itmax) .and.  &
                   (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
        
          it = it + 1
          do k = 1,3
            if (k == 1) then
              ip = 1
              iq = 2
            else if (k ==2) then
              ip = 2
              iq = 3
            else
              ip = 1
              iq = 3
            end if
            
            if (abs(A(ip,iq)) > 1d-50) then
                
              tau = ( A(iq,iq) - A(ip,ip) ) / ( 2.0 * A(ip,iq) )
              if (tau >= 0.0) then
                sign_tau = 1.0
              else
                sign_tau = -1.0
              end if
              
              t = sign_tau / ( abs(tau) + sqrt(1.0 + tau**2) )
              c = 1.0 / sqrt(1.0 + t**2)
              s = t * c

              temp1 = c * A(1, ip) - s * A(1, iq)
              temp2 = c * A(2, ip) - s * A(2, iq)
              temp3 = c * A(3, ip) - s * A(3, iq)
              A(1, iq) = s * A(1, ip) + c * A(1, iq)
              A(2, iq) = s * A(2, ip) + c * A(2, iq)
              A(3, iq) = s * A(3, ip) + c * A(3, iq)
              A(1, ip) = temp1
              A(2, ip) = temp2
              A(3, ip) = temp3

              temp1 = c * A(ip, 1) - s * A(iq, 1)
              temp2 = c * A(ip, 2) - s * A(iq, 2)
              temp3 = c * A(ip, 3) - s * A(iq, 3)
              A(iq, 1) = s * A(ip, 1) + c * A(iq, 1)
              A(iq, 2) = s * A(ip, 2) + c * A(iq, 2)
              A(iq, 3) = s * A(ip, 3) + c * A(iq, 3)
              A(ip, 1) = temp1
              A(ip, 2) = temp2
              A(ip, 3) = temp3
  
            end if ! A(ip,iq)<>0
            
          end do ! k
          
          ! optional output writing
          if (iOpt == 1) then
            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
          end if
          
        end do ! while
     
        ! get principal stresses from diagonal of A
        S1 = A(1, 1)
        S2 = A(2, 2)
        S3 = A(3, 3)
      
        ! derived invariants
        P = (S1 + S2 + S3) / 3.
        Q = sqrt( ( (S1 - S2)**2 + (S2 - S3)**2 + (S3 - S1)**2 ) / 2. )

        ! Sort eigenvalues S1 <= S2 <= S3
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
        end if
        
        if (S2 > S3) then
          t   = S3
          S3  = S2
          S2  = t
        end if
        
        if (S1 > S2) then
          t   = S2
          S2  = S1
          S1  = t
        end if
        
    end subroutine Eig_3a
      
    
    
      subroutine CalculatePrincipalStresses(IntGlo,Sig,SigPrin)
      !**********************************************************************
      !
      ! Implemented in the frame of the MPM project.
      !
      !  Note: This subroutine calculates the principal stresses
      !
      !**********************************************************************

      use ModGlobalConstants
      implicit none
      
      !Local variables
      real(REAL_TYPE), dimension(3) :: xN1,xN2,xN3
      real(REAL_TYPE) :: Sig1,Sig2,Sig3,p,q
      !In Variables
      integer(INTEGER_TYPE), intent(in) :: IntGlo ! Global ID of Gauss point or particle
      real(REAL_TYPE), intent(in), dimension(6) :: Sig
      !Out Variables
      real(REAL_TYPE), intent(out), dimension(6) :: SigPrin

      call PrnSig(1,Sig,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

      If (Sig1 >= Sig2.and.Sig2 >= Sig3) then
        SigPrin(1) = Sig1
        SigPrin(2) = Sig2
        SigPrin(3) = Sig3
      else if (Sig1 >= Sig3.and.Sig3 >= Sig2) then
        SigPrin(1) = Sig1
        SigPrin(2) = Sig3
        SigPrin(3) = Sig2
      else if (Sig3 >= Sig1.and.Sig1 >= Sig2) then
        SigPrin(1) = Sig3
        SigPrin(2) = Sig1
        SigPrin(3) = Sig2
      else if (Sig3 >= Sig2.and.Sig2 >= Sig1) then
        SigPrin(1) = Sig3
        SigPrin(2) = Sig2
        SigPrin(3) = Sig1
      else if (Sig2 >= Sig1.and.Sig1 >= Sig3) then
        SigPrin(1) = Sig2
        SigPrin(2) = Sig1
        SigPrin(3) = Sig3
      else if (Sig2 >= Sig3.and.Sig3 >= Sig1) then
        SigPrin(1) = Sig2
        SigPrin(2) = Sig3
        SigPrin(3) = Sig1
      end if

        SigPrin(4) = 0.0d0
        SigPrin(5) = 0.0d0
        SigPrin(6) = 0.0d0

      end subroutine CalculatePrincipalStresses