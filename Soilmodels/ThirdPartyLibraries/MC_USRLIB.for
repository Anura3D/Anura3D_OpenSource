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
    !	Copyright (C) 2020  Members of the Anura3D MPM Research Community
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
      
      Subroutine MZEROR(R,K)
C
C***********************************************************************
C
C     Function: To make a real array R with dimension K to zero
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J) = 0.0D0
      End Do

      Return
      End


      Subroutine MZEROI(I,K)
C
C***********************************************************************
C
C     Function: To make an integre array I with Dimension K to zero
C
C***********************************************************************
C
      Dimension I(*)

      Do J=1,K
        I(J)=0
      End Do

      Return
      End

      Subroutine SETRVAL(R,K,V)
C
C***********************************************************************
C
C     Function: To fill a real array R with Dimension K with value V
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J)=V
      End Do

      Return
      End

      Subroutine SETIVAL(I,K,IV)
C
C***********************************************************************
C
C     Function: To fill an integer array I with Dimension K with value IV
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension I(*)

      Do J=1,K
        I(J)=IV
      End Do

      Return
      End

      Subroutine COPYIVEC(I1,I2,K)
C
C***********************************************************************
C
C     Function: To copy an integer array I1 with Dimension K to I2
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension I1(*),I2(*)

      Do  J=1,K
        I2(J)=I1(J)
      End Do

      Return
      End

      Subroutine COPYRVEC(R1,R2,K)
C
C***********************************************************************
C
C     Function: To copy a Double array R1 with Dimension K to R2
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R1(*),R2(*)

      Do J=1,K
        R2(J)=R1(J)
      End Do

      Return
      End


      Logical Function IS0ARR(A,N)
C
C***********************************************************************
C    Function :  To check whether a real array contains only zero values.
C                When an array contains only zero's is might not need to be
C                written to the XXX file.
C                exit Function when first non-zero value occured or when
C                all elements are checked and are zero.
C
C    Input:  A : array to be checked
C            N : number of elements in array that should be checked
C
C    Output : .TRUE.  when all elements are 0
C             .FALSE. when at least one element is not zero
C
C    Called by :  Subroutine TOBXX
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*)
      Is0Arr=.False.
      Do I=1,N
        If ( A(I) .Ne. 0 ) Return
      End Do
      Is0Arr=.True.
      Return
      End

      Logical Function IS0IARR(IARR,N)
C
C***********************************************************************
C    Function :  To check whether a integer array contains only zero values.
C                Similar to IS0ARR
C
C    Input:  IARR : array to be checked
C            N    : number of elements in array that should be checked
C
C    Output : .TRUE.  when all elements are 0
C             .FALSE. when at least one element is not zero
C
C    Called by :  Subroutine TOBXX
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension IARR(*)

      Is0IArr=.False.
      Do I=1,N
        If ( IARR(I) .Ne. 0 ) Return
      End Do
      Is0IArr=.True.
      Return
      End
C***********************************************************************
      Subroutine MulVec(V,F,K)
C***********************************************************************
C
C     Function: To multiply a real vector V with dimension K by F
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*)

      Do J=1,K
        V(J)=F*V(J)
      End Do

      Return
      End     ! Subroutine Mulvec
C***********************************************************************
      Subroutine MatVec(xMat,IM,Vec,N,VecR)
C***********************************************************************
C
C     Calculate VecR = xMat*Vec
C
C I   xMat  : (Square) Matrix (IM,*)
C I   Vec   : Vector
C I   N     : Number of rows/colums
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine MatVec

C***********************************************************************
      Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
C***********************************************************************
C
C     Calculate VecR() = R1*Vec1()+R2*Vec2()
C
C I   Vec1,
C I   Vec2  : Vectors
C I   R1,R2 : Multipliers
C I   N     : Number of rows
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec
C
C***********************************************************************
      Double Precision Function DInProd(A,B,N)
C***********************************************************************
C
C     Returns the Inproduct of two vectors
C
C I   A,B  : Two vectors
C I   N    : Used length of vectors
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*),B(*)
C***********************************************************************

      X = 0
      Do I=1,N
        X = X + A(I)*B(I)
      End Do
      DInProd = X
      Return
      End     ! Function DInProd
C
C***********************************************************************
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2
C
C I   xMat1 : Matrix (Id1,*)
C I   xMat2 : Matrix (Id2,*)
C I   nR1   : Number of rows in resulting matrix    (= No rows in xMat1)
C I   nC2   : Number of columns in resulting matrix (= No cols in xMat2)
C I   nC1   : Number of columns in matrix xMat1
C             = Number  rows    in matrix xMat2
C O   xMatR : Resulting matrix (IdR,*)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(Id1,*),xMat2(Id2,*),xMatR(IdR,*)
C**********************************************************************

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMat

C***********************************************************************
      Subroutine MatMatSq(n, xMat1, xMat2, xMatR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2 for square matrices, size n
C
C I   n     : Dimension of matrices
C I   xMat1 : Matrix (n,*)
C I   xMat2 : Matrix (n,*)
C O   xMatR : Resulting matrix (n,*)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(n,*),xMat2(n,*),xMatR(n,*)
C**********************************************************************

      Do I=1,n
        Do J=1,n
          X=0
          Do K=1,n
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMatSq

C***********************************************************************
      Subroutine WriVal ( io, C , V )
C***********************************************************************
C
C Write (Double) value to file unit io (when io>0)
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,V
    1 Format( A,3x, 1x,1p,e12.5)
      Return
      End
C***********************************************************************
      Subroutine WriIVl ( io, C , I )
C***********************************************************************
C
C Write (integer) value to file unit io (when io>0)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,I
    1 Format( A,3x, 1x,I6)
      Return
      End
C***********************************************************************
      Subroutine WriIVc ( io, C , iV , n )
C***********************************************************************
C
C Write (integer) vector to file unit io (when io>0)
C
C***********************************************************************
      Character C*(*)
      Dimension iV(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) (iv(i),i=1,n)
    1 Format( ( 2(3x,5i4) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec ( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      If (Len_Trim(C).Le.6) Then
        Write(io,2) C,( V(i),i=1,n)
      Else
        Write(io,*) C
        Write(io,1) ( V(i),i=1,n)
      End If
    1 Format( ( 2(1x, 3(1x,1p,e10.3) ) ) )
    2 Format( A, ( T7, 2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec5( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 5 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) ( V(i),i=1,n)
    1 Format( 5(1x,1p,e12.5) )
      Return
      End
C***********************************************************************
      Subroutine WriMat ( io, C , V , nd, nr, nc )
C***********************************************************************
C
C Write (Double) matrix to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(nd,*)

      If (io.Le.0) Return

      Write(io,*) C
      Do j=1,nr
        Write(io,1) j,( V(j,i),i=1,nc)
      End Do
    1 Format(i4, (  T7,2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************

      Subroutine MatInvPiv(Aorig,B,N)
      Implicit Double Precision (A-H,O-Z)
      Dimension Aorig(n,*), B(n,*),A(:,:)
      Allocatable :: A
      Allocate ( A(n,n) ) ! No error checking !!
      Call CopyRVec(AOrig, A, n*n )
      Call MZeroR(B,n*n)
      Do i=1,n
        B(i,i) = 1d0
      End Do
      Do I=1,n
        T=A(I,I)
        iPiv=i
        Do j=i+1,n
          If ( Abs(A(j,i)) .Gt. Abs(A(iPiv,i))  ) iPiv=j
        End Do
        If (iPiv.Ne.i) Then
          Do j=1,n
            x         = A( i  ,j)
            A( i  ,j) = A(iPiv,j)
            A(iPiv,j) = x
            x         = B( i  ,j)
            B( i  ,j) = B(iPiv,j)
            B(iPiv,j) = x
          End Do
          T=A(I,I)
        End If
        Do J=1,n
          A(I,J)=A(I,J)/T
          B(I,J)=B(I,J)/T
        End Do
        Do K=1,n
          If (K.Ne.I) Then
            T=A(K,I)
            Do J=1,n
              A(K,J)=A(K,J)-T*A(I,J)
              B(K,J)=B(K,J)-T*B(I,J)
            End Do
          End If
        End Do
      End Do
      DeAllocate ( A  )
      Return
      End ! MatinvPiv

C***********************************************************************  
            subroutine PrnSig(IOpt, ntens, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions
      !            from cartesian stress vector
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
      
      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt, ntens
        double precision, intent(in), dimension(ntens) :: S ! size of stress tensor (6 for 3D, 4 for 2D)
        double precision, intent(out), dimension(3) :: xN1, xN2, xN3
        double precision, intent(out) :: S1, S2, S3, P, Q

        if (IOpt == 1) then
          call Eig_3(0,ntens,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! Calculate principal direction
        else
          call Eig_3a(0,ntens,S,S1,S2,S3,P,Q) ! Do not calculate principal direction
        end if
     
      end subroutine PrnSig
C***********************************************************************
      
      subroutine Eig_3(iOpt, ntens, St, xN1, xN2, xN3, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses and directions 
      !            from cartesian stress vector
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

      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt, ntens
        double precision, intent(in), dimension(ntens) :: St
        double precision, intent(out), dimension(3) :: xN1, xN2, xN3 ! is (3) for 2D and 3D
        double precision, intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        double precision, dimension(3, 3) :: A, V ! is (3,3) for 2D and 3D
        double precision :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
        integer :: i, k, it, itmax, ip, iq, iS1, iS2, iS3

        ! Put cartesian stress vector into matrix A
        select case(ntens)
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
        do i = 1, ntens
          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-16 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it < itmax) .and.  
     &            (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
        
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
!            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
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
        do i = 1, 3 ! is 3 for 2D and 3D
          xN1(i) = V(i, is1)
          xN2(i) = V(i, is2)
          xN3(i) = V(i, is3)
        end do
      
      end subroutine Eig_3

      
      subroutine Eig_3a(iOpt, ntens, St, S1, S2, S3, P, Q)
      !-------------------------------------------------------------------
      !
      !  Function: calculate principal stresses from cartesian stress vector
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

      
      implicit none
      
        ! arguments
        integer, intent(in) :: IOpt, ntens
        double precision, intent(in), dimension(ntens) :: St ! is (6) for 3D and (4) for 2D
        double precision, intent(out) :: S1, S2, S3, P, Q
        
        ! local variables
        integer :: IDim
        double precision, dimension(:,:), allocatable :: A ! (3,3) for 2D and 3D
        double precision :: abs_max_s, tol, tau, sign_tau, t, c, s, temp1, temp2, temp3
        integer :: i, k, it, itmax, ip, iq
        
        IDim = 3
        allocate(A(IDim,IDim)) 
        
        select case (ntens)       
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
        do i = 1,ntens
          if (abs(St(i)) > abs_max_s) abs_max_s = abs(St(i))
        end do
      
        ! set tolerance
        tol = 1d-20 * abs_max_s
        
        ! get principal stresses and directions iteratively
        it = 0
        itmax = 50
        do while ( (it < itmax) .and.  
     &              (abs(A(1,2)) + abs(A(2,3)) + abs(A(1,3)) > tol) )
        
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
!            call GiveMessage('error: ' // trim(String(abs(A(1, 2)) + abs(A(2, 3)) + abs(A(1, 3)))))
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
      
C
C***********************************************************************
      Logical Function LEqual(A,B,Eps)
C***********************************************************************
C
C     Returns .TRUE.  when two real values are (almost) equal,
C             .FALSE. otherwise
C
C I   A,B  : Two real values to be compared
C I   Eps  : Toleration (Magnitude ~= 1E-5)
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
C***********************************************************************
      LEqual =.True.
      If (A .Eq. B) Return
      If (DAbs(A-B) .LT. 0.5D0*Eps*( DAbs(A) + DAbs(B) + Eps ) )Return
      LEqual =.False.
      Return
      End     ! function LEqual
C
C***********************************************************************
      Subroutine CrossProd(xN1,xN2,xN3)
C***********************************************************************
C
C     Returns cross product of xN1 and xN2
C
C I   xN1,xN2 : Two basic vectors
C O   xN3     : Resulting vector
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*)
C***********************************************************************

      xN3(1) = xN1(2)*xN2(3) - xN1(3)*xN2(2)
      xN3(2) = xN1(3)*xN2(1) - xN1(1)*xN2(3)
      xN3(3) = xN1(1)*xN2(2) - xN1(2)*xN2(1)

      Return
      End     ! Subroutine CrossProd
C
C***********************************************************************
      Double Precision Function ArcSin(X,ie)
C***********************************************************************
C
C     Returns the Arc Sine of X
C
C I   X : Input value
C
C     Note : In stead of using default routine DASIN we use this one
C            because �X� can be slightly beyond 1 and this will give
C            a RTE using DASIN(X)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
C***********************************************************************
      Ie=0
      S = (1-X*X)
!      If (S .Lt. -1E-10) Ie=1
!      If (S .Lt. -1E-10) Write(*,1) X,S
!      If (S .Lt. -1E-10) Write(2,1) X,S
    1 Format(' ArcSin(',1x,1p,e13.5e3,') , S =',1x,1p,e13.5e3)
      If (S.LT.0) S = 0
      S = DSQRT(S)
      ArcSin = DATan2(X,S)
      Return
      End     ! function ArcSin
C
C***********************************************************************
      subroutine CarSig(S1, S2, S3, xN1, xN2, xN3, ntens, Stress)
      !----------------------------------------------------------
      !
      !  Function:  Returns the Cartesian stresses using 
      !             principal stress and principal direction
      !
      !  S1, S2, S3      I   R     principal stress
      !  xN1, xN2, xN3   I   R()   principal direction
      !  Stress          O   R()   cartesian stress
      !
      !----------------------------------------------------------
      
      implicit none
      
        ! arguments
        double precision, intent(in) :: S1, S2, S3
        double precision, intent(in), dimension(3) :: xN1, xN2, xN3
        integer, intent(in):: ntens
        double precision, intent(out), dimension(ntens) :: Stress
        
        ! local variables
        integer :: I
        integer :: IDim
        double precision, dimension(:, :), allocatable :: SM, T, TT, STT
        
        IDim = 3
        allocate(SM(IDim, IDim), T(IDim, IDim), TT(IDim, IDim), STT(IDim, IDim))

        ! put principal directions in matrix T
        do I = 1,3
            T(I,1) = xN1(I)
            T(I,2) = xN2(I)
            T(I,3) = xN3(I)
            TT(1,I) = T(I,1)
            TT(2,I) = T(I,2)
            TT(3,I) = T(I,3)            
        end do
      


        ! put principal stresses in matrix SM
        SM = 0.0
        SM(1,1) = S1
        SM(2,2) = S2
        SM(3,3) = S3

        ! calculate matrix multiplication T.S.TT
        call MatMat(SM, IDim, TT,  IDim, IDim, IDim, IDim , STT, IDim)
        call MatMat(T,  IDim, STT, IDim, IDim, IDim, IDim , SM,  IDim)
       
        ! Extract cartesian stress vector from stress matrix
        do I = 1, IDim
          Stress(I) = SM(I, I)
        end do        
        
        ! Extract cartesian stress vector from stress matrix
        Stress(4) = SM(2, 1)  
        if (ntens == 6) then
          Stress(5) = SM(3, 2)
          Stress(6) = SM(3, 1)
        end if
        
      end subroutine CarSig
      subroutine setveclen(xn,n,xl)
C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN(*)
      x=0
      do i=1,n
        x=x+xn(i)**2
      end do
      if (x.Ne.0) Then
        f=xl/sqrt(x)
        do i=1,3
          xn(i)=xn(i)*f
        end do
      end if
      return
      end ! setveclen

C**********************************************************************
C End Of file
C**********************************************************************
