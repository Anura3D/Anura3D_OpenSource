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
	   
	   
	  module ModNURBS
      !**********************************************************************
      !
      !    Function:  This module provides specific routines for evaluating shape functions with NURBS
      !
      !     $Revision: 1 $
      !     $Date: 2022-25-05 07:58:40 -0400 (Wed, 25 May 2022) $
      !
      !**********************************************************************

      use ModGeometryMath
      use ModString
      use ModReadCalculationData
      use ModGlobalConstants

      implicit none

    contains ! routines of this module
    
    !
    !
    !
    !
    !subroutine InitialiseShapeFunctionsNURBS()
    !!**********************************************************************
    !!
    !!    Function: This function is the kernel for evaluating the NURBS shapefunctions. 
    !!              The user needs to define the order of the NURBS shape functions (NURBSOrder),
    !!               which tells you how many times you need to loop. 
    !!
    !!
    !!    @note : 2D element... maybe?
    !!
    !! I  :
    !! I  :
    !! O  :
    !! O  :
    !! O  :
    !!
    !!**********************************************************************
    !
    !! N_i,p(xi) = (xi - xi_i)   N_i,p-1(xi)  +  (xi_i+p+1 - xi)     N_i+1,p-1(xi)
    !!             _____________                ___________________
    !!             (xi_i+p - xi)                (xi_i+p+1 - xi_i+1)
    !!
    !!             left fraction                   right fraction
    !
    !implicit none 
    !
    !! Input 
    !integer(INTEGER_TYPE), intent(in) :: pp ! polynomial order(s)
    !integer(INTEGER_TYPE), intent(in) :: nn ! number of basis functions used to construct the B-spline 
    !!:: quadrature point location 
    !!:: element number 
    !!:: control net BB --> nodal coordinates NodeCoord
    !!:: knot vector 
    !!:: INC --> connectivity array 
    !!:: IEN --> connectivity array 
    !!:: n_en --> number of local shape functions
    !
    !
    !! Local     
    !integer(INTEGER_TYPE) :: ii ! knot index that we need to loop on
    !
    !
    !
    !! Output 
    !real(REAL_TYPE), intent(out), dimension(:) :: RR ! vector of local shape functions values 
    !real(REAL_TYPE), intent(out), dimension(:) :: dRR ! vector of derivatives of local shape functions 
    !real(REAL_TYPE), intent(out), dimension(:) :: JJ ! Jacobian determinant 
    !
    !! connectivity array links every local shape function number to a global shape function number 
    !! IEN arrays need to be built directly from knot vectors and polynomial orders 
    !! 
    !! - IGA Initialization 
    !! -- Define computational parameters here 
    !RR = 0
    !
    !NDIM = 2 !2D implementation for quad elements 
    !
    !
    !! Read NURBS data 
    !call ReadNURBSData(uKnot, vKnot, ControlPoints, uKnot_size, vKnot_size, uPP_Order, vPP_Order)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! - Build connectivities and allocate global arrays 
    !                    
    !nn = uKnot_size - uPP_Order - 1 !number of univariate basis functions for xi direction 
    !mm = vKnot_size - vPP_Order - 1 !number of univariate basis functions for eta direction 
    !
    !!ww = wKnot_size - wPP_Order - 1 ! need to include this for 3D implementation 
    !
    !call Build_INC_IEN_Array(uPP_Order, vPP_Order, nn, mm, & ! input
    !                         INN, IEN, nel, nnp, nen & !output 
    !    )
    !                    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !! gauss point details input here 
    !
    !xi_tilde_size = 1
    !eta_tilde_size = 1
    !
    !call NURBSGaussPoints(NoGP, xi_tilde, eta_tilde, xi_tilde_size, eta_tilde_size) ! note that the gauss point information here 
    !
    !
    !! - Evaluate basis functions and their derivatives 
    !do ii = 1,nel !loop over elements 
    !    do aa = 1,xi_tilde_size !xi gauss point
    !        do bb = 1,eta_tilde_size !eta gauss point
    !            
    !            
    !            ![RR, dRR_dx]
    !            
    !            call ShapeFunctionInitialization(&
    !                & !!!!!! INPUTS: 
    !            xi_tilde(aa), eta_tilde(bb), & 
    !            & ! quadrature points... if you were to plot basis functions then you have to use up all your points as quadrature points
    !            ii, & !element number
    !            uPP_Order, vPP_Order, & !polynomial order
    !            ControlPoints, & !B
    !            uKnot, vKnot, & !Knot vector 
    !            INN, & !%connectivity array (global to local)
    !            IEN, & !connectivity array (local to global)
    !            nen, & !number of local shape functions 
    !            NDIM, & !number of dimensions 
    !            nn, &
    !            mm, &
    !            & !!!!!! OUTPUTS:
    !            RR, & ! vector of local shape functions (2D tensor product)
    !            dR_dx, & !vector of shape function derivatives (2D tensor product) 
    !            JJ, xx_yy_coordinates & !jacobian for the element 
    !            )
    !            
    !
    !            
    !            
    !    
    !    
    !    
    !    
    !    
    !    
    !    
    !
    !        end
    !    end
    !end
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !!! 
    !
    !
    !
    !
    !!
    !
    !!do ! loop accross every element 
    !
    !
    !
    !!! Shape function evaluation 
    !!ii = nn+pp+1 ! this is the size of the knot vector that we need to loop on 
    !!
    !!
    !!do ii 
    !!
    !!left = xi - knot_xi 
    !!right = knot_xi - xi 
    !!
    !!! Shape function derivative evaluation 
    !!dR_dx 
    !
    !end subroutine InitialiseShapeFunctionsNURBS
    !
    !
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !subroutine NURBSGaussPoints(xi_tilde, eta_tilde, xi_tilde_size, eta_tilde_size)
    !implicit none 
    !! initialise variables 
    !integer(INTEGER_TYPE) :: xi_tilde_size, eta_tilde_size
    !integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: xi_tilde, eta_tilde, !zeta_tilde 
    !
    !! depends on how many gauss point you want -> hard coded one gauss point here 
    !xi_tilde_size = 1 !note that these should be GOM file inputs 
    !eta_tilde_size = 1 !note that these should be GOM file inputs 
    !
    !! allocate variables to the number of gauss points selected 
    !allocate(xi_tilde(xi_tilde_size), stat=IError)
    !allocate(eta_tilde(eta_tilde_size), stat=IError)
    !
    !!assume one gauss point for now point for now 
    !xi_tilde = 0.0 ![-1, 0, 1]; !gauss point location in xi direction (GOM file inputs hardcoded herein)
    !eta_tilde = 0.0 ![-1, 0, 1]; !gauss point location in eta direction (GOM file inputs hardcoded herein) 
    !
    !! here we are selecting the middle of the element 
    !
    !end subroutine NURBSGaussPoints
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    !!!!!!!!!!!!!!!! Subroutine revision 1 
    !!!!!!!!!!!!!!!!!!!! START: Reading NURBS data here 
    !subroutine ReadNURBSData(uKnot, vKnot, ControlPoints)
    !! initialize variables 
    !! Input: 
    !! Output: 
    !
    !integer(INTEGER_TYPE) :: uPP_Order, vPP_Order
    !integer(INTEGER_TYPE) :: uKnot_size, vKnot_size
    !integer(INTEGER_TYPE) :: stat, IError 
    !integer(INTEGER_TYPE) :: NumberOfControlPoints
    !real(REAL_TYPE), dimension(:), allocatable :: uKnot 
    !real(REAL_TYPE), dimension(:), allocatable :: vKnot
    !
    !real(REAL_TYPE), dimension(:,:), allocatable :: ControlPoints
    !uPP_Order = 2 ! order in the xi direction
    !vPP_Order = 2 ! order in the eta direction 
    !
    !! Knot size 
    !uKnot_size = 7
    !vKnot_size = 6 
    !
    !allocate(uKnot(1, uKnot_size), stat=IError)
    !allocate(vKnot(1, vKnot_size), stat=IError)
    !
    !! -- Define one dimensional knot vector here 
    !! -- Number of basis functions = # of items in knot vectors - order - 1 
    !uKnot = [0, 0, 0, 0.5, 1, 1, 1] ! nn = 7 - 2 - 1 = 4
    !vKnot = [0, 0, 0,   1, 1, 1] ! mm = 6 - 2 - 1 = 3
    !
    !! --Control Points are assumed to be the nodal coordinates 
    !! ControlPoints = NodalCoordinates
    !! Example from page 34 from Cottrell et al. (2009) 
    !
    !NumberOfControlPoints = 12
    !
    !allocate(ControlPoints(NumberOfControlPoints,1), stat=IError)
    !
    !ControlPoints = [0	0 ! i =1 ; j=1
    !                -1	0 ! i =1 ; j=2
    !                -2	0 ! i =1 ; j=3
    !                 0	1 ! i =2 ; j=1
    !                -1	2 ! i =2 ; j=2
    !                -2	2 ! i =2 ; j=3
    !                 1	1.5 ! i =3 ; j=1
    !                 1	4 ! i =3 ; j=2
    !                 1	5 ! i =3 ; j=3
    !                 3	1.5 ! i =4 ; j=1
    !                 3	4 ! i =4 ; j=2
    !                 3  5] ! i =4 ; j=3 %note that number of control points equal number of basis functions 
    ! 
    !                    
    !     
    !! no of columns in CDB matrix should equal to no. of itembers in knot
    !! vector minus one. As we go higher order, we leave out columns in the
    !! beginning of the matrix equal to the order of the basis function. For
    !! example, if order = 1, then we leave out the first column. Similarly, if
    !! order = 2, then we leave out the first two columns. 
    !
    !
    !end subroutine ReadNURBSData 
    !!!!!!!!!!!!!!!!!!!! END: Reading NURBS data here
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    !subroutine ShapeFunctionInitialization(&
    !                & !!!!!! INPUTS: 
    !            xi_tilde(aa), eta_tilde(bb), & 
    !            & ! quadrature points... if you were to plot basis functions then you have to use up all your points as quadrature points
    !            ii, & !element number
    !            uPP_Order, vPP_Order, & !polynomial order
    !            ControlPoints, & !B
    !            uKnot, vKnot, & !Knot vector 
    !            INN, & !%connectivity array (global to local)
    !            IEN, & !connectivity array (local to global)
    !            nen, & !number of local shape functions 
    !            NDIM, & !number of dimensions 
    !            nn, &
    !            mm, &
    !            & !!!!!! OUTPUTS:
    !            RR, & ! vector of local shape functions (2D tensor product)
    !            dR_dx, & !vector of shape function derivatives (2D tensor product) 
    !            JJ, xx_yy_coordinates & !jacobian for the element 
    !            )
    !
    !implicit none 
    !
    !! initialization 
    !integer(INTEGER_TYPE) :: ii, jj, aa, bb, cc, JJ !kk=0;
    !integer(INTEGER_TYPE) :: ni, nj !nk 
    !integer(INTEGER_TYPE), intent(in) :: nen, nel, nnp, NDIM,
    !integer(INTEGER_TYPE), intent(in) :: uPP_Order, vPP_Order
    !integer(INTEGER_TYPE), intent(in), dimension(nen,nel) :: IEN 
    !integer(INTEGER_TYPE), intent(in), dimension(nnp, NDIM) :: INN
    !
    !
    !allocate(RR, nen, 1) !basis functions 
    !allocate(dR_dx, nen, 1) !basis function derivatives 
    !
    !JJ = 0 !jacobian determinant 
    !
    !! local variable initialization 
    !ni = 0 !NURBS coordinates 
    !nj = 0 !NURBS coordinates 
    !!nk = 0 -> note 2D implementation 
    !
    !! array of univariate Bspline basis function
    !allocate(NN, pp+1)
    !allocate(MM, qq+1) 
    !!allocate(LL, rr+1) 
    !
    !
    !!Univariate Bspline function derivative w.r.t. appropriate parametric coordinates
    !allocate(dN_dxi, pp+1)
    !allocate(dN_deta, qq+1) 
    !!allocate(dN_dzeta, rr+1) 
    !
    !! Trivariate NURBS function derivatives w.r.t. parametric coordinates 
    !allocate(dR_dxi, nen, NDIM)
    !! Derivative of physical coordinates w.r.t. parametric coordinates 
    !allocate(dx_dxi, NDIM, NDIM)
    !! Inverse of dx_dxi
    !allocate(dxi_dx, NDIM, NDIM)
    !! Derivative of parametric coordinates w.r.t. parent element coordinates
    !allocate(dxi_dtildexi, NDIM, NDIM)
    !! Jacobian matrix 
    !allocate(J_mat, NDIM, NDIM)
    !
    !! Loop counters  
    !ii=0
    !jj=0
    !aa=0
    !bb=0 
    !cc=0 !kk=0; 
    !
    !! Local basis function counter 
    !loc_num = 0 
    !
    !! Dummy sums for calculating rational derivatives 
    !sum_xi = 0 
    !sum_eta = 0 
    !sum_tot = 0 
    !!sum_zeta = 0 
    !
    !! find the local gauss point in NURBS coordinate (xi, eta) 
    !
    !call InputLocalGaussPointToOutputNURBSGaussPoint(...) ! add in the argument here
    !
    !
    !! zeta = ( (KV_Zeta(nk+1) - KV_Zeta(nk) ) * eta_zeta ...
    !!         + (KV_Zeta(nk+1) + KV_Zeta(nk)) ) * 0.5;   
    !
    !
    !! - Calculate univariate Bspline functions using (2.1) and (2.2) and their
    !!   derivatives using (2.12)
    !! xi direction 
    !call Bspline_basis_and_deriv(ni, pp, KV_Xi, xi, & ! inputs 
    !                            NN, dN_dxi & ! outputs 
    !                                        )
    !! eta direction    
    !call Bspline_basis_and_deriv(nj, qq, KV_Eta, eta, & ! inputs
    !                        MM, dM_deta & ! outputs 
    !                                        )
    !
    !! zeta direction 
    !!call Bspline_basis_and_deriv(nk, rr, KV_Zeta, & ! inputs 
    !!                        LL, dL_dzeta & ! outputs
    !!                                        ) 
    !
    !
    !! note that there is a subroutine called Plot_NURBS_Surface_(ControlPoints, NN, MM, & !inputs
    !!                        xx_yy_coordinates & !output 
    !!)
    !
    !!Build numerator and denominator 
    !do jj = 0,qq 
    !    do ii = 0,pp
    !        loc_num = loc_num + 1 !local basis function number 
    !        RR(loc_num) = NN(pp+1-ii) * MM(qq+1-jj) !* LL(rr+1-kk)
    !        !ControlPoints(ni-ii, nj-jj); % Function numerator %, nk-kk)
    !        ! Note that we need to multiply by the weight in NURBS, but in
    !        ! Bsplines that is not necessary 
    !        
    !                    
    !        sum_tot = sum_tot + RR(loc_num) ! Function denominator 
    !        
    !        dR_dxi(loc_num,1) = dN_dxi(pp+1-ii) * MM(qq+1-jj) !* LL(rr+1-kk);!...
    !        !ControlPoints(ni-ii, nj-jj); % Derivative numerator %, nk-kk)
    !    
    !        sum_xi = sum_xi + dR_dxi(loc_num, 1) ! Derivative denom 
    !    
    !        dR_dxi(loc_num,2) = NN(pp+1-ii) * dM_deta(qq+1-jj) !* LL(rr+1-kk); 
    !   
    !        ! * ControlPoints(ni-ii, nj-jj); ! Derivative numerator 
    !    
    !        sum_eta = sum_eta + dR_dxi(loc_num,2) ! Derivative denominator 
    !     
    !
    !        !dR_dxi(loc_num, 3) = NN_dx(pp+1) * MM_dx(qq+1-jj); %dL_dzeta(rr+1-kk)
    !     
    !        ! * ControlPoints(ni-ii, nj-jj, nk-kk, 4) % Derivative numerator 
    !     
    !     
    !        !sum_zeta = sum_zeta + dR_dxi(loc_num,1); % Derivative denominator 
    !end
    !end 
    !
    !
    !
    !! - Divide by denominator to complete definitions of functions and derivatives w.r.t. parametric coordinates 
    !!do loc_num = 1,nen
    !!    RR(loc_num) = RR(loc_num)/sum_tot
    !!    
    !!    dR_dxi(loc_num,1) = (dR_dxi(loc_num,1) * sum_tot &
    !!        - RR(loc_num) * sum_xi)/ sum_tot^2
    !!    dR_dxi(loc_num,2) = (dR_dxi(loc_num,2) * sum_tot &
    !!        - RR(loc_num) * sum_eta)/ sum_tot^2
    !!    
    !!end do 
    !
    !
    !        
    !
    !
    !            
    !            
    !            end subroutine ShapeFunctionInitialization 
    !            
    !            
    !            
    !            
    !            
    !    subroutine Bspline_basis_and_deriv(ni, PP_Order, uKnot, xi, size_uKnot, & !input 
    !                                       NN, dN_dxi & !output 
    !                                       )
    !    ! Cox De Boor 1D equation
    !    
    !    implicit none
    !    
    !    real(REAL_TYPE), allocatable, dimension(:,:,:) :: BasisFunctionCDB
    !    real(REAL_TYPE), allocatable, dimension(:,:,:) :: BasisFunctionCDB_derivative
    !    
    !    !%This function evaluates basis functions and derivative at a gauss point xi
    !
    !    !%resolution
    !
    !    !%other inputs? : ni (NURBS coordinate),  
    !    !% output 
    !    !% 1- vector of p+1 function values corresponding to the p+1 functions that
    !    !%    are non-zero  on [xi, xi_i+1] -> [N_1, N_2, N_3, N_4,...]
    !    !% 2- 
    !    !% 3- 
    !    !% 4- 
    !
    !    !%pertinent variables calculated. 
    !    !% Unique_uKnot = unique(uKnot); % - Find unique elements in the knot vector 
    !    !% Max_uKnot = max(Unique_uKnot); 
    !    !% Min_uKnot = min(Unique_uKnot); 
    !            
    !    
    !    NoOfBasisFunctions_uKnot = size_uKnot - 1 !number of columns should be equal to this 
    !            
    !    nGP = 1 ! this should be implemented better here 
    !    
    !    !initialize 
    !    allocate(BasisFunctionCDB(nGP, NoOfBasisFunctions_uKnot, PP_Order+1))
    !    allocate(BasisFunctionCDB_derivative(nGP, NoOfBasisFunctions_uKnot, PP_Order+1))
    !    
    !    !Zero order basis functions 
    !    do ii = 1,NoOfBasisFunctions_uKnot !loop accross basis function
    !        uknot_left = uKnot(ii)
    !        uknot_right = uKnot(ii+1)
    !        do jj = 1,nGP !loop accross domain data points 
    !            if (uKnot_left<=xi(jj)) && (xi(jj)<uKnot_right) then 
    !                BasisFunctionCDB(jj,ii,1) = 1               
    !                !note that the derivative of zero order is zero 
    !            end if 
    !        end do 
    !    end do 
    !         
    !                          
    !    if (PP_Order > 0) then 
    !        do kk = 1, PP_Order then 
    !            do ii = 1,(NoOfBasisFunctions_uKnot-kk) then !loop accross number of basis functions (bandwidth is equal to pp+1)
    !               uknot_ii = uKnot(ii) ! -> I am giving ii here as an input ni 
    !               uknot_iiPlusPP = uKnot(ii+kk) ! -> k here is related to the order 
    !               uknot_iiPlusPPPlus1 = uKnot(ii+kk+1)
    !               uknot_iiPlus1 = uKnot(ii+1)
    !            
    !               do jj = 1,nGP !loop accross the domain data points 
    !               
    !                   LeftBasis = BasisFunctionCDB(jj,ii+kk-1,kk)
    !                   RightBasis = BasisFunctionCDB(jj,ii+kk,kk)
    !                
    !                   DenominatorLeft = (uknot_iiPlusPP - uknot_ii)
    !                   DenominatorRight = (uknot_iiPlusPPPlus1 - uknot_iiPlus1)
    !                
    !                   BasisFunctionCDBLeft = ( LeftBasis * (xi(jj)-uknot_ii)/DenominatorLeft)
    !                   BasisFunctionCDBLeft_Derivative = ( LeftBasis * kk/DenominatorLeft)
    !                
    !                   BasisFunctionCDBRight = ( RightBasis * (uknot_iiPlusPPPlus1-xi(jj))/DenominatorRight)
    !                   BasisFunctionCDBRight_Derivative = - ( RightBasis * kk/DenominatorRight)
    !            
    !               end do
    !               
    !               if (BasisFunctionCDBLeft == NaN) then 
    !                   BasisFunctionCDBLeft = 0 
    !                   BasisFunctionCDBLeft_Derivative = 0
    !               end if 
    !               
    !               if (BasisFunctionCDBRight == NaN) then 
    !                   BasisFunctionCDBRight = 0 
    !                   BasisFunctionCDBRight_Derivative = 0
    !               end if
    !               
    !               BasisFunctionCDB(jj,ii+kk,kk+1) =  BasisFunctionCDBLeft + BasisFunctionCDBRight
    !               BasisFunctionCDB_Derivative(jj,ii+kk,kk+1) =  BasisFunctionCDBLeft_Derivative + BasisFunctionCDBRight_Derivative
    !               
    !            
    !            end do 
    !        end do
    !        
    !    end if
    !    
    !    NN = BasisFunctionCDB(:, &
    !                        PP_Order+1:size(BasisFunctionCDB),&
    !                        PP_Order+1)
    !
    !    dN_dxi = BasisFunctionCDB_Derivative(:, &
    !                        PP_Order+1:size(BasisFunctionCDB), &
    !                        PP_Order+1)
    !    
    !    
    !    
    !    
    !    
    !    end subroutine 
    !            
    !            
    ! subroutine InputLocalGaussPointToOutputNURBSGaussPoint(ee, KV_Xi_size, KV_Eta_size, KV_Xi, KV_Eta, xi_tilde, eta_tilde, & ! input 
    !                                                        xi, eta & ! output 
    !                                                         )
    ! 
    ! implicit none 
    ! 
    ! integer(INTEGER_TYPE), intent(in) :: ee 
    ! integer(INTEGER_TYPE), intent(in) :: KV_Xi_size, KV_Eta_size
    ! 
    ! real(REAL_TYPE), intent(in), dimension(KV_Xi_size) :: KV_Xi 
    ! real(REAL_TYPE), intent(in), dimension(KV_Eta_size) :: KV_Eta
    ! 
    ! real(REAL_TYPE), intent(in) :: xi_tilde, eta_tilde
    ! 
    ! real(REAL_TYPE), intent(out) :: xi, eta      
    ! 
    ! integer(INTEGER_TYPE) :: ni, nj, 
    !
    ! 
    ! 
    !! - NURBS coordinates; convection consistent with Algorithm 7 in Cottrell et al. (2009)
    !ni = INN(IEN(1,ee),1)
    !nj = INN(IEN(2,ee),2)
    !!nk = INN(IEN(3,ee),3) ! implementation only in 2D 
    !
    !! - Calculate parametric coordinates from parent element coordinates 
    !!   Knot vectors KV_Xi, KV_Eta, and KV_Zeta, and 
    !!   parent element coordinates xi_tilde, eta_tilde, zeta_tilde are given as 
    !!   inputs 
    !xi = ( (KV_Xi(ni+1) - KV_Xi(ni) ) * xi_tilde &
    !    + (KV_Xi(ni+1) + KV_Xi(ni)) ) * 0.5
    !
    !eta = ( (KV_Eta(nj+1) - KV_Eta(nj) ) * eta_tilde &
    !    + (KV_Eta(nj+1) + KV_Eta(nj)) ) * 0.5   
    !
    !
    !
    ! 
    ! end subroutine InputLocalGaussPointToOutputNURBSGaussPoint
    !
    
    !!!!!!!!!!!!!!!! Subroutine revision 1 !!!!!!!!!!!!!!!!!!!!!!!
    subroutine Build_INC_IEN_Array(pp, qq, nn, mm, & ! input
                             INN, IEN, nel, nnp, nen & !output 
        )
    implicit none 
    ! Description here about inputs/outputs/what does this subroutine does 
    !% pp = 2;
    !% qq = 2;
    !% nn = 4; 
    !% mm = 3; 
    !%note that this is implemented in 2D here 

    !%inputs 
    !% 1- polynomial orders (p,q,r)
    !% 2- number of univariate basis functions (n, m, l)

    !%outputs 
    !% 1- total number of elements, nel
    !% 2- total number of global basis functions, nnp 
    !% 3- number of local basis functions, nen
    !% 4- INC: consumes a global basis function number and a parametric
    !%         direction number and returns the corresponding NURBS coordinate 
    !%         
    !% 5- IEN: 
    
    !Initialise variables 
    integer(INTEGER_TYPE) :: NDIM 
    integer(INTEGER_TYPE) :: ee, AA, BB, CC, ii, jj, iloc, jloc, stat, IError
    
    !input
    integer(INTEGER_TYPE), intent(in) :: pp, qq, nn, mm 
    
    !output
    integer(INTEGER_TYPE), intent(out) :: nnp, nen, nel
    integer(INTEGER_TYPE), intent(out), allocatable, array(:,:) :: IEN !connectivity array 
    integer(INTEGER_TYPE), intent(out), allocatable, array(:,:) :: INN !NURBS coordinate array (also called INC)
        
    NDIM = 2 !2D implementation  ! hardcoded 2 dimensional
    
    ! - global variable definitions and initializations: 
    nel = (nn-pp) * (mm-qq) !number of elements -> note 2D implementation = 2 elements in the example 
    ! nel = (4-2)*(3-2) = 2 elements
    !     element 1   element 2
    !     __________ __________
    !    |          |          |
    !    |          |          |
    !    |          |          |
    !    |          |          |
    !    |__________|__________|
    nnp = nn*mm !number of global basis functions (global here refers to its global domain within the 'super' element)
    ! nnp = 4*3 = 12 ... This is also equal to the number of control points  
    nen = (pp+1) * (qq+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
    ! nen = (2+1)*(2+1) = 9 local basis functions 
    
    allocate(INN(nnp, NDIM), stat=IError) ! INN has the size of number of control points(or global basis functions x NDIM )
    allocate(IEN(nen, nel), stat=IError)  ! IEN has the size of number of local basis functions x NDIM 
    
    INN = 0 !NURBS coordinate array (also called INC)
    IEN = 0 !connectivity array
    
    !local variable initialization 
    ee = 0 
    AA = 0
    BB = 0 
    CC = 0
    ii = 0
    jj = 0 
    ! kk = 0
    iloc = 0 
    jloc = 0
    ! kloc = 0
    
    do jj = 1,mm ! loop over the eta univariate basis function
        do ii = 1,nn ! loop over the xi univariate basis function
            
            AA=AA+1 !increment global function number (AA should have a max of mm*nn = 12 = number of global basis = number of control points)
            
            !assign NURBS coordinate 
            INN(AA, 1) = ii
            INN(AA, 2) = jj
            
            if ( (ii>=pp+1) .and. (jj>=qq+1) ) then 
                ee=ee+1 !increment element number 
                
                do jloc = 0,qq
                    do iloc = 0,pp
                        BB = AA - jloc*nn - iloc !global function number 
                        CC = (jloc*(pp+1)) + iloc + 1
                        IEN(CC,ee) = BB
                    end do 
                end do 
            end if 
        end do 
    end do 
    
        
    end subroutine Build_INC_IEN_Array
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !
    !subroutine 
    !
    !
    !end subroutine 
    
    
    
    
    
    
    !subroutine LinearShapeFunctionNURBS
    
    end module ModNURBS