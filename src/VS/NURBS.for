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
      !use ModMPMData
      !use ModMeshInfo

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
        subroutine Bspline_basis_and_deriv(NXiKnotOrder, NXiKnotEntries, nGP, Xi_ParametricDomain, XiKnotEntries, & !input 
                                           NN_IncludesZeroValues, dN_dxi_IncludesZeroValues & !output 
                                           )
        ! Cox De Boor 1D equation
        ! [NN, dN_dxi] = Bspline_basis_and_deriv(ni, NXiKnotOrder, XiKnotEntries, Xi_ParametricDomain) ! -> NN should be nn+pp in size
        
        !(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
        !                            NN_IncludesZeroValues, dN_dxi_IncludesZeroValues)

        implicit none
        

        
        ! local 
        integer(INTEGER_TYPE) :: ii, jj, kk
        real(REAL_TYPE) :: uknot_left, uknot_right
        real(REAL_TYPE) :: uknot_ii, uknot_iiPlusPP, uknot_iiPlusPPPlus1, uknot_iiPlus1
        real(REAL_TYPE) :: LeftBasis, RightBasis, DenominatorLeft, DenominatorRight 
        real(REAL_TYPE) :: BasisFunctionCDBLeft, BasisFunctionCDBLeft_Derivative 
        real(REAL_TYPE) :: BasisFunctionCDBRight, BasisFunctionCDBRight_Derivative
        integer(INTEGER_TYPE) :: stat,IError
        integer(INTEGER_TYPE) :: NoOfBasisFunctions_uKnot

        
        
        ! input 
        !integer(INTEGER_TYPE), intent(in) :: ni_NURBS
        !integer(INTEGER_TYPE), intent(in) :: nn_NURBS_NumberOfUnivariateXiKnots
        integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
        integer(INTEGER_TYPE), intent(in) :: nGP
        real(REAL_TYPE), intent(in), dimension(nGP) :: Xi_ParametricDomain 
        real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
        
        ! output 
        real(REAL_TYPE), allocatable, dimension(:, :, :), intent(inout) :: NN_IncludesZeroValues
        real(REAL_TYPE), allocatable, dimension(:, :, :), intent(inout) :: dN_dxi_IncludesZeroValues
        
        
        !real(REAL_TYPE) :: uknot_ii
        !real(REAL_TYPE) :: uknot_iiPlusPP
        !real(REAL_TYPE) :: uknot_iiPlusPPPlus1
        !real(REAL_TYPE) :: uknot_iiPlus1
        !real(REAL_TYPE) :: LeftBasis
        !real(REAL_TYPE) :: RightBasis
        !real(REAL_TYPE) :: DenominatorLeft
        !real(REAL_TYPE) :: DenominatorRight
        !real(REAL_TYPE) :: BasisFunctionCDBLeft
        !real(REAL_TYPE) :: BasisFunctionCDBLeft_Derivative
        !real(REAL_TYPE) :: BasisFunctionCDBRight
        !real(REAL_TYPE) :: BasisFunctionCDBRight_Derivative
        !
        !integer(INTEGER_TYPE) :: kk
        !integer(INTEGER_TYPE) :: ii
        !integer(INTEGER_TYPE) :: jj
        
        
        !real(REAL_TYPE), allocatable, dimension(:,:,:) :: BasisFunctionCDB
        !real(REAL_TYPE), allocatable, dimension(:,:,:) :: BasisFunctionCDB_derivative
        
        !%This function evaluates basis functions and derivative at a gauss point xi
    
        !%resolution
    
        !%other inputs? : ni (NURBS coordinate),  
        !% output 
        !% 1- vector of p+1 function values corresponding to the p+1 functions that
        !%    are non-zero  on [xi, xi_i+1] -> [N_1, N_2, N_3, N_4,...]
        !% 2- 
        !% 3- 
        !% 4- 
    
        !%pertinent variables calculated. 
        !% Unique_uKnot = unique(uKnot); % - Find unique elements in the knot vector 
        !% Max_uKnot = max(Unique_uKnot); 
        !% Min_uKnot = min(Unique_uKnot); 
                
        
                  
        NoOfBasisFunctions_uKnot = NXiKnotEntries-1 ! calculate number of knot spans !number of columns should be equal to this 
                
        !nGP = 1 ! this should be implemented better here 
        
        ! allocate memory to basis function and basis function derivative matrices 
        !allocate(BasisFunctionCDB(nGP, NoOfBasisFunctions_uKnot, PP_Order+1))
        !allocate(BasisFunctionCDB_derivative(nGP, NoOfBasisFunctions_uKnot, PP_Order+1))
        allocate(NN_IncludesZeroValues(nGP, NoOfBasisFunctions_uKnot, NXiKnotOrder+1), stat=IError)
        allocate(dN_dxi_IncludesZeroValues(nGP, NoOfBasisFunctions_uKnot, NXiKnotOrder+1), stat=IError)
        
        NN_IncludesZeroValues = 0
        dN_dxi_IncludesZeroValues = 0
        
        !Zero order basis functions 
        do ii = 1,NoOfBasisFunctions_uKnot !loop accross basis function
            uknot_left = XiKnotEntries(ii)
            uknot_right = XiKnotEntries(ii+1)
            do jj = 1,nGP !loop accross domain data points 
                if ( (uKnot_left<=Xi_ParametricDomain(jj)) .and. (Xi_ParametricDomain(jj)<uKnot_right) ) then 
                    NN_IncludesZeroValues(jj,ii,1) = 1  ! zero order values are in slot 1     
                    !note that the derivative of zero order is zero 
                end if 
            end do 
        end do 
             
        
        
                              
        if (NXiKnotOrder > 0) then
            do kk = 1, NXiKnotOrder 
                do ii = 1,(NoOfBasisFunctions_uKnot-kk)  !loop accross number of basis functions (bandwidth is equal to pp+1)
                   uknot_ii = XiKnotEntries(ii) ! -> I am giving ii here as an input ni 
                   uknot_iiPlusPP = XiKnotEntries(ii+kk) ! -> k here is related to the order 
                   uknot_iiPlusPPPlus1 = XiKnotEntries(ii+kk+1)
                   uknot_iiPlus1 = XiKnotEntries(ii+1)
                
                   do jj = 1,nGP !loop accross the domain data points 
                   
                       LeftBasis = NN_IncludesZeroValues(jj,ii+kk-1,kk)
                       RightBasis = NN_IncludesZeroValues(jj,ii+kk,kk)
                    
                       DenominatorLeft = (uknot_iiPlusPP - uknot_ii)
                       DenominatorRight = (uknot_iiPlusPPPlus1 - uknot_iiPlus1)
                    
                       
                       if (DenominatorLeft > 1e-15) then 
                           
                           BasisFunctionCDBLeft = ( LeftBasis * (Xi_ParametricDomain(jj)-uknot_ii)/DenominatorLeft)
                           BasisFunctionCDBLeft_Derivative = ( LeftBasis * kk/DenominatorLeft)
                       else 
                           BasisFunctionCDBLeft = 0.0
                           BasisFunctionCDBLeft_Derivative = 0.0
                           
                       end if 
                       
                       if (DenominatorRight > 1e-15) then                       
                          BasisFunctionCDBRight = ( RightBasis * (uknot_iiPlusPPPlus1-Xi_ParametricDomain(jj))/DenominatorRight)
                          BasisFunctionCDBRight_Derivative = - ( RightBasis * kk/DenominatorRight)
                       
                       else
                        BasisFunctionCDBRight = 0.0
                        BasisFunctionCDBRight_Derivative = 0.0
                       end if
                       
                       
                   
                       NN_IncludesZeroValues(jj,ii+kk,kk+1) =  BasisFunctionCDBLeft + BasisFunctionCDBRight
                       dN_dxi_IncludesZeroValues(jj,ii+kk,kk+1) =  BasisFunctionCDBLeft_Derivative + BasisFunctionCDBRight_Derivative
                   
                    
                
                   end do
                   
                   !if (isnan(BasisFunctionCDBLeft)) then 
                   !    BasisFunctionCDBLeft = 0 
                   !    BasisFunctionCDBLeft_Derivative = 0
                   !end if 
                   !
                   !if (isnan(BasisFunctionCDBRight)) then 
                   !    BasisFunctionCDBRight = 0 
                   !    BasisFunctionCDBRight_Derivative = 0
                   !end if
                   
                   
                
                end do 
            end do
            
        end if
        
        !NN = BasisFunctionCDB(:, &
        !                    PP_Order+1:size(BasisFunctionCDB),&
        !                    PP_Order+1)
        !
        !dN_dxi = BasisFunctionCDB_Derivative(:, &
        !                    PP_Order+1:size(BasisFunctionCDB), &
        !                    PP_Order+1)
        
        
        
        
        
                                           end subroutine Bspline_basis_and_deriv
                                           
                                           
                                           
                                           
                                           
    !--------------------------------------------------------------------------------------------------------------
    ! SINGLE PARTICLE 
    
    subroutine Bspline_basis_and_deriv_SINGLEPARTICLE(NXiKnotOrder, NXiKnotEntries, nGP, Xi_ParametricDomain, XiKnotEntries, & !input 
                                           NN_IncludesZeroValues, dN_dxi_IncludesZeroValues & !output 
                                           )
        ! Cox De Boor 1D equation
        ! [NN, dN_dxi] = Bspline_basis_and_deriv(ni, NXiKnotOrder, XiKnotEntries, Xi_ParametricDomain) ! -> NN should be nn+pp in size
        
        !(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
        !                            NN_IncludesZeroValues, dN_dxi_IncludesZeroValues)

        implicit none
        

        
        ! local 
        integer(INTEGER_TYPE) :: ii, jj, kk
        real(REAL_TYPE) :: uknot_left, uknot_right
        real(REAL_TYPE) :: uknot_ii, uknot_iiPlusPP, uknot_iiPlusPPPlus1, uknot_iiPlus1
        real(REAL_TYPE) :: LeftBasis, RightBasis, DenominatorLeft, DenominatorRight 
        real(REAL_TYPE) :: BasisFunctionCDBLeft, BasisFunctionCDBLeft_Derivative 
        real(REAL_TYPE) :: BasisFunctionCDBRight, BasisFunctionCDBRight_Derivative
        integer(INTEGER_TYPE) :: stat,IError
        integer(INTEGER_TYPE) :: NoOfBasisFunctions_uKnot

        
        
        ! input 
        !integer(INTEGER_TYPE), intent(in) :: ni_NURBS
        !integer(INTEGER_TYPE), intent(in) :: nn_NURBS_NumberOfUnivariateXiKnots
        integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
        integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
        integer(INTEGER_TYPE), intent(in) :: nGP
        real(REAL_TYPE), intent(in), dimension(1) :: Xi_ParametricDomain 
        real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
        
        ! output 
        real(REAL_TYPE), allocatable, dimension(:, :, :), intent(inout) :: NN_IncludesZeroValues
        real(REAL_TYPE), allocatable, dimension(:, :, :), intent(inout) :: dN_dxi_IncludesZeroValues
        
        
        !real(REAL_TYPE) :: uknot_ii
        !real(REAL_TYPE) :: uknot_iiPlusPP
        !real(REAL_TYPE) :: uknot_iiPlusPPPlus1
        !real(REAL_TYPE) :: uknot_iiPlus1
        !real(REAL_TYPE) :: LeftBasis
        !real(REAL_TYPE) :: RightBasis
        !real(REAL_TYPE) :: DenominatorLeft
        !real(REAL_TYPE) :: DenominatorRight
        !real(REAL_TYPE) :: BasisFunctionCDBLeft
        !real(REAL_TYPE) :: BasisFunctionCDBLeft_Derivative
        !real(REAL_TYPE) :: BasisFunctionCDBRight
        !real(REAL_TYPE) :: BasisFunctionCDBRight_Derivative
        !
        !integer(INTEGER_TYPE) :: kk
        !integer(INTEGER_TYPE) :: ii
        !integer(INTEGER_TYPE) :: jj
        
        
        !real(REAL_TYPE), allocatable, dimension(:,:,:) :: BasisFunctionCDB
        !real(REAL_TYPE), allocatable, dimension(:,:,:) :: BasisFunctionCDB_derivative
        
        !%This function evaluates basis functions and derivative at a gauss point xi
    
        !%resolution
    
        !%other inputs? : ni (NURBS coordinate),  
        !% output 
        !% 1- vector of p+1 function values corresponding to the p+1 functions that
        !%    are non-zero  on [xi, xi_i+1] -> [N_1, N_2, N_3, N_4,...]
        !% 2- 
        !% 3- 
        !% 4- 
    
        !%pertinent variables calculated. 
        !% Unique_uKnot = unique(uKnot); % - Find unique elements in the knot vector 
        !% Max_uKnot = max(Unique_uKnot); 
        !% Min_uKnot = min(Unique_uKnot); 
                
        
                  
        NoOfBasisFunctions_uKnot = NXiKnotEntries-1 ! calculate number of knot spans !number of columns should be equal to this 
                
        !nGP = 1 ! this should be implemented better here 
        
        ! allocate memory to basis function and basis function derivative matrices 
        !allocate(BasisFunctionCDB(nGP, NoOfBasisFunctions_uKnot, PP_Order+1))
        !allocate(BasisFunctionCDB_derivative(nGP, NoOfBasisFunctions_uKnot, PP_Order+1))
        allocate(NN_IncludesZeroValues(nGP, NoOfBasisFunctions_uKnot, NXiKnotOrder+1), stat=IError)
        allocate(dN_dxi_IncludesZeroValues(nGP, NoOfBasisFunctions_uKnot, NXiKnotOrder+1), stat=IError)
        
        NN_IncludesZeroValues = 0
        dN_dxi_IncludesZeroValues = 0
        
        !Zero order basis functions 
        do ii = 1,NoOfBasisFunctions_uKnot !loop accross basis function
            uknot_left = XiKnotEntries(ii)
            uknot_right = XiKnotEntries(ii+1)
            do jj = 1,nGP !loop accross domain data points 
                if ( (uKnot_left<=Xi_ParametricDomain(jj)) .and. (Xi_ParametricDomain(jj)<uKnot_right) ) then 
                    NN_IncludesZeroValues(jj,ii,1) = 1  ! zero order values are in slot 1     
                    !note that the derivative of zero order is zero 
                end if 
            end do 
        end do 
             
        
        
                              
        if (NXiKnotOrder > 0) then
            do kk = 1, NXiKnotOrder 
                do ii = 1,(NoOfBasisFunctions_uKnot-kk)  !loop accross number of basis functions (bandwidth is equal to pp+1)
                   uknot_ii = XiKnotEntries(ii) ! -> I am giving ii here as an input ni 
                   uknot_iiPlusPP = XiKnotEntries(ii+kk) ! -> k here is related to the order 
                   uknot_iiPlusPPPlus1 = XiKnotEntries(ii+kk+1)
                   uknot_iiPlus1 = XiKnotEntries(ii+1)
                
                   do jj = 1,nGP !loop accross the domain data points 
                   
                       LeftBasis = NN_IncludesZeroValues(jj,ii+kk-1,kk)
                       RightBasis = NN_IncludesZeroValues(jj,ii+kk,kk)
                    
                       DenominatorLeft = (uknot_iiPlusPP - uknot_ii)
                       DenominatorRight = (uknot_iiPlusPPPlus1 - uknot_iiPlus1)
                    
                       
                       if (DenominatorLeft > 1e-15) then 
                           
                           BasisFunctionCDBLeft = ( LeftBasis * (Xi_ParametricDomain(jj)-uknot_ii)/DenominatorLeft)
                           BasisFunctionCDBLeft_Derivative = ( LeftBasis * kk/DenominatorLeft)
                       else 
                           BasisFunctionCDBLeft = 0.0
                           BasisFunctionCDBLeft_Derivative = 0.0
                           
                       end if 
                       
                       if (DenominatorRight > 1e-15) then                       
                          BasisFunctionCDBRight = ( RightBasis * (uknot_iiPlusPPPlus1-Xi_ParametricDomain(jj))/DenominatorRight)
                          BasisFunctionCDBRight_Derivative = - ( RightBasis * kk/DenominatorRight)
                       
                       else
                        BasisFunctionCDBRight = 0.0
                        BasisFunctionCDBRight_Derivative = 0.0
                       end if
                       
                       
                   
                       NN_IncludesZeroValues(jj,ii+kk,kk+1) =  BasisFunctionCDBLeft + BasisFunctionCDBRight
                       dN_dxi_IncludesZeroValues(jj,ii+kk,kk+1) =  BasisFunctionCDBLeft_Derivative + BasisFunctionCDBRight_Derivative
                   
                    
                
                   end do
                   
                   !if (isnan(BasisFunctionCDBLeft)) then 
                   !    BasisFunctionCDBLeft = 0 
                   !    BasisFunctionCDBLeft_Derivative = 0
                   !end if 
                   !
                   !if (isnan(BasisFunctionCDBRight)) then 
                   !    BasisFunctionCDBRight = 0 
                   !    BasisFunctionCDBRight_Derivative = 0
                   !end if
                   
                   
                
                end do 
            end do
            
        end if
        
        !NN = BasisFunctionCDB(:, &
        !                    PP_Order+1:size(BasisFunctionCDB),&
        !                    PP_Order+1)
        !
        !dN_dxi = BasisFunctionCDB_Derivative(:, &
        !                    PP_Order+1:size(BasisFunctionCDB), &
        !                    PP_Order+1)
        
        
        
        
        
        end subroutine Bspline_basis_and_deriv_SINGLEPARTICLE
    
    
                                           
    ! SINGLE PARTICLE                                 
    !---------------------------------------------------------------------------------------------------------------
                                           
    
                                           
                                           
                                           
                                           
                                           
                                           
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
    subroutine Build_INC_IEN_Array()
    
    !pp, qq, nn, mm, & ! input
    !                         INN, IEN, nel, nnp, nen & !output 
    !    )
    !pp -> NXiKnotOrder 
    !qq -> NEtaKnotOrder 
    !nn -> NumberOfUnivariateXiKnots
    !mm -> NumberOfUnivariateEtaKnots 
    
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
    !integer(INTEGER_TYPE) :: NDIM 
    integer(INTEGER_TYPE) :: ee, AA, BB, CC, ii, jj, iloc, jloc, stat, IError
    
    !input
    !integer(INTEGER_TYPE), intent(in) :: pp, qq, nn, mm 
    
    !output
    !integer(INTEGER_TYPE), intent(out) :: nnp, nen, nel
    !integer(INTEGER_TYPE), intent(out), allocatable, dimension(:,:) :: IEN !connectivity array 
    !integer(INTEGER_TYPE), intent(out), allocatable, dimension(:,:) :: INN !NURBS coordinate array (also called INC)
        
    !NDIM = 2 !2D implementation  ! hardcoded 2 dimensional
    
    ! - global variable definitions and initializations: 
    nel_NURBS = (nn_NURBS_NumberOfUnivariateXiKnots-NXiKnotOrder) * (mm_NURBS_NumberOfUnivariateEtaKnots-NEtaKnotOrder) !number of elements -> note 2D implementation = 2 elements in the example 
    
    
    ! overwrite the number of elements by nel_NURBS --> this should be calculated by the code and not given as an input
    !nel_NURBS
    
    ! nel = (4-2)*(3-2) = 2 elements
    !     element 1   element 2
    !     __________ __________
    !    |          |          |
    !    |          |          |
    !    |          |          |
    !    |          |          |
    !    |__________|__________|
    nnp_NURBS = nn_NURBS_NumberOfUnivariateXiKnots*mm_NURBS_NumberOfUnivariateEtaKnots !number of global basis functions (global here refers to its global domain within the 'super' element)
    ! nnp = 4*3 = 12 ... This is also equal to the number of control points  
    nen_NURBS = (NXiKnotOrder+1) * (NEtaKnotOrder+1) !number of local basis functions (local here refers to a knot span i.e. accross one single element)
    ! nen = (2+1)*(2+1) = 9 local basis functions 
    
    allocate(INN(nnp_NURBS, NVECTOR), stat=IError) ! INN has the size of number of control points(or global basis functions x NDIM )
    allocate(IEN(nen_NURBS, nel_NURBS), stat=IError)  ! IEN has the size of number of local basis functions x NDIM 
    
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
    
    do jj = 1,mm_NURBS_NumberOfUnivariateEtaKnots ! loop over the eta univariate basis function
        do ii = 1,nn_NURBS_NumberOfUnivariateXiKnots ! loop over the xi univariate basis function
            
            AA=AA+1 !increment global function number (AA should have a max of mm*nn = 12 = number of global basis = number of control points)
            
            !assign NURBS coordinate 
            INN(AA, 1) = ii
            INN(AA, 2) = jj
            
            if ( (ii>=NXiKnotOrder+1) .and. (jj>=NEtaKnotOrder+1) ) then 
                ee=ee+1 !increment element number 
                
                do jloc = 0,NEtaKnotOrder
                    do iloc = 0,NXiKnotOrder
                        BB = AA - jloc*nn_NURBS_NumberOfUnivariateXiKnots - iloc !global function number 
                        CC = (jloc*(NXiKnotOrder+1)) + iloc + 1
                        IEN(CC,ee) = BB
                    end do 
                end do 
            end if 
        end do 
    end do 
    
    
    
    !call BuildKnotBezierMesh()
    
    !ElementConnectivities = IEN
    
        
    end subroutine Build_INC_IEN_Array
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !
    !subroutine 
    !
    !
    !end subroutine 
    
    
    
    
    
    
    !subroutine LinearShapeFunctionNURBS
    
    
    
    
    
    
    
    !subroutine BuildKnotBezierMesh()
    ! loops over the knots and build the Bezier mesh through unique size of the knot vectors 
    
    !XiKnotEntries
    
    !EtaKnotEntries
    
    
    
    
    
    
    
    !end subroutine 
        
    
    
    
    
    
        subroutine FindUniqueEnteriesInVector(OriginalVector, VectorSize, UniqueVector, counter)
        ! this rountine finds the unique items in a vector 
        ! this is used in the context of NURBS
        
        implicit none
        
        ! local 
        integer(INTEGER_TYPE) :: stat, IError ! used for error control
        integer(INTEGER_TYPE) :: ii
        
        ! inputs
        integer(INTEGER_TYPE), intent(in) :: VectorSize
        real(REAL_TYPE), dimension(:), intent(in) :: OriginalVector ! input originial vector 
        
        
        ! outputs 
        real(REAL_TYPE), allocatable, dimension(:), intent(out) :: UniqueVector ! outout unique vector
        integer(INTEGER_TYPE), intent(out) :: counter ! this tells you how big is your matrix 
        !integer(INTEGER_TYPE),  :: UniqueSize
        
        
        ii = 1
        
        ! assign first entry of the vector to the unique vector
        ! in this way the first entry of the vector is secured. 
        !UniqueVector(ii) = OriginalVector(ii)
        
        ! establish counter 
        counter = 1 
        
        
        do ii = 1, VectorSize-1 
            
           ! check if the second entry is similar to the first one. 
            if ( OriginalVector(ii) == OriginalVector(ii+1) ) then ! if the second entry is the same as the previous entry 
                
                !... do nothing 
                
            else 
                !... if it is different, then write into out unique vector and increase counter 
                !UniqueVector(ii) = OriginalVector(ii)
                
                ! increase counter 
                counter = counter + 1
                
            end if 
            
            
        end do 
        
            
            
        
        ! allocate unique vector with the number of unique elements
        
        allocate(UniqueVector(counter), stat = IError)
            
        
        UniqueVector = 0.0
            
        
        ! it looks like we need to allocate the matrix first before we do anything else when it is allocatable
        ii = 1
        
        
            
        ! assign first entry of the vector to the unique vector
        ! in this way the first entry of the vector is secured. 
        UniqueVector(1) = OriginalVector(1)
        
  
            
        ! establish counter 
        counter = 2
        
        
        
            
        do ii = 1, VectorSize-1 
            ! check if the second entry is similar to the first one. 
            if ( OriginalVector(ii) == OriginalVector(ii+1) ) then ! if the second entry is the same as the previous entry 
                
                
                    !... do nothing 
            
            else 
                
                    !... if it is different, then write into out unique vector and increase counter 
                
                    UniqueVector(counter) = OriginalVector(ii+1)
                
                
                    ! increase counter 
                counter = counter + 1
            end if 
            
            
        
            
        end do 
        
        counter = counter - 1
        
        
        end subroutine 
    
    
    
        
        
        subroutine InitialiseShapeFunctionsLINE2_NURBS (HS, dHS, Wt, &               !classic inout parameters
                                                        XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder) !NURBS related inputs
        ! This is 1D implementation        
        !(ElementNumber, NumberOfUnivariateXiKnots, NXiKnotOrder, NXiKnotEntries, XiKnotEntries, ee, nen, nel, nnp, NDIM, IEN, INN)
        !**********************************************************************
        !      !!!!------      Adapting this for NURBS      ------!!!!
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
        !               (-1)            (1)
        !**********************************************************************
        implicit none
      
          real(REAL_TYPE), dimension(:,:), intent(inout) :: HS !these should not be allocatables at this point
          real(REAL_TYPE), dimension(:,:,:), intent(inout) :: dHS !these should not be allocatables at this point 
          real(REAL_TYPE), dimension(:), intent(inout) :: Wt !these should not be allocatables at this point 
          
          integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
          integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
          real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries !dimension(NXiKnotEntries)
          
          real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Xi_ParametricDomain !dimension(NXiGaussPoints)
          
          ! Basis functions allocatables 
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: NN_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dN_dxi_IncludesZeroValues
          
                    
          integer(INTEGER_TYPE) :: IError, stat     
          integer(INTEGER_TYPE) :: ii, jj, kk

          
          !real(REAL_TYPE), dimension(NXiKnotOrder+1) :: NN_WithoutZeroValues
                    
              ! - evaluate each basis function value at the gauss point 
              call Bspline_basis_and_deriv(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
                                    NN_IncludesZeroValues, dN_dxi_IncludesZeroValues) !output 
           
              !allocate(HS(NXiKnotOrder+1), stat=IError)
              !allocate(dHS(NXiKnotOrder+1), stat=IError)
              !
              
              jj = NXiKnotOrder+1
              kk = (2*NXiKnotOrder)+1
              
              do ii = jj, kk
                  ! note that the first entry that says '1' should be looping accross gaussian points
                 HS(1,ii-1) = NN_IncludesZeroValues(NXiGaussPoints, ii, NXiKnotOrder+1)
                 dHS(1, ii-1, 1) = dN_dxi_IncludesZeroValues(NXiGaussPoints, ii, NXiKnotOrder+1)
              end do 
             
              Wt = 2
          
        
        end subroutine InitialiseShapeFunctionsLINE2_NURBS
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        ! local variables 
          !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Xi
          !integer(INTEGER_TYPE) :: nGP_Xi
          !integer(INTEGER_TYPE) :: ee_NURBS
          !integer(INTEGER_TYPE) :: ni_NURBS
          !integer(INTEGER_TYPE) :: nj_NURBS
          !integer(INTEGER_TYPE) :: nk_NURBS
          !          
          !integer(INTEGER_TYPE) :: ii, jj, kk
          !
          !real(REAL_TYPE), allocatable, dimension(:) :: xi_tilde
          !real(REAL_TYPE), allocatable, dimension(:) :: eta_tilde
          !real(REAL_TYPE), allocatable, dimension(:) :: zeta_tilde
          
             
          
          !
          !real(REAL_TYPE), allocatable, dimension(:) :: Eta_ParametricDomain
          !real(REAL_TYPE), allocatable, dimension(:) :: Zeta_ParametricDomain
          !
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: NN_IncludesZeroValues
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dN_dxi_IncludesZeroValues
          !
          !
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: MM_IncludesZeroValues
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dM_dxi_IncludesZeroValues
          !
          !
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: LL_IncludesZeroValues
          !real(REAL_TYPE), allocatable, dimension(:,:,:) :: dL_dxi_IncludesZeroValues     
        
        
        
        
        
        
        
        
        
        
        !
          !! initialize basis functions 
          !Number_of_Knot_Spans_Xi = nn_NURBS_NumberOfUnivariateXiKnots-1 ! calculate number of knot spans
          !Number_of_Knot_Spans_Eta = mm_NURBS_NumberOfUnivariateEtaKnots-1 ! calculate number of knot spans
          !!Number_of_Knot_Spans_Zeta = ll_NURBS_NumberOfUnivariateZetaKnots-1 ! calculate number of knot spans
          !
          !nGP_Xi = 1
          !nGP_Eta = 1
          !!nGP_Zeta = 1
          !
          !do ee_NURBS = 1, nel_NURBS !loop over elements 
          !    !need to loop over elements (assume 2D)
          !    ni_NURBS = INN(IEN(1,ee_NURBS), 1) !get ni which tells you the left knot of your knot span in the xi direction  
          !    nj_NURBS = INN(IEN(1,ee_NURBS), 2) !get nj which tells you the left knot of your knot span in the eta direction 
          !    !nk_NURBS = INN(IEN(1,ee_NURBS), 3) !get nk which tells you the left knot of your knot span in the eta direction
          !
          !
          !    allocate(xi_tilde(nGP_Xi), stat=IError) ! allocating one gauss point in each element (this should be between -1 and 1)
          !    allocate(Xi_ParametricDomain(nGP_Xi), stat=IError) ! allocating one gauss point in each element (this should be in between the knot span of element) 
          !    
          !    
          !    allocate(eta_tilde(nGP_Eta), stat=IError) ! allocating one gauss point in each element (this should be between -1 and 1)
          !    allocate(Eta_ParametricDomain(nGP_Eta), stat=IError) ! allocating one gauss point in each element (this should be in between the knot span of element) 
          !    
          !    !allocate(zeta_tilde(nGP_Zeta), stat=IError) ! allocating one gauss point in each element (this should be between -1 and 1)
          !    !allocate(Zeta_ParametricDomain(nGP_Zeta), stat=IError) ! allocating one gauss point in each element (this should be in between the knot span of element) 
          !               
          !    
          !    ! specify local gauss point in parent domain 
          !    xi_tilde = 0.0 !xi in local domain
          !    eta_tilde = 0.0 !eta in local domain 
          !    !zeta_tilde = 0.0 !zeta in local domain 
          !
          !    ! find this local gauss point in the parametric domain 
              !Xi_ParametricDomain = ( (XiKnotEntries(ni_NURBS+1) - XiKnotEntries(ni_NURBS) ) * xi_tilde &
              !                      + (XiKnotEntries(ni_NURBS+1) + XiKnotEntries(ni_NURBS)) ) * 0.5 ! this will give you your gaussian location in the parametric domain
              !
              !Eta_ParametricDomain = ( (EtaKnotEntries(nj_NURBS+1) - EtaKnotEntries(nj_NURBS) ) * eta_tilde &
              !                      + (EtaKnotEntries(nj_NURBS+1) + EtaKnotEntries(nj_NURBS)) ) * 0.5 ! this will give you your gaussian location in the parametric domain
              
              !Zeta_ParametricDomain = ( (ZetaKnotEntries(nk_NURBS+1) - ZetaKnotEntries(nk_NURBS) ) * zeta_tilde &
              !                      + (ZetaKnotEntries(nk_NURBS+1) + ZetaKnotEntries(nk_NURBS)) ) * 0.5 ! this will give you your gaussian location in the parametric domain
              !
              
        
        
        subroutine InitialiseShapeFunctionsQUAD4_NURBS(HS, dHS, Wt, & !classic inout parameters
                                                    HS_Xi, dHS_Xi, Wt_Xi, &
                                                    HS_Eta, dHS_Eta, Wt_Eta, &
                                                    XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
                                                    EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
                                                    ni, nj) !NURBS related inputs in the eta direction 
        !**********************************************************************
        !
        !    SUBROUTINE: InitialiseShapeFunctionsQUAD4_NURBS
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at  one Gaussian integration point for a 4-noded 2D quadrilateral element using NURBS.
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
        
        implicit none
        
          !!real(REAL_TYPE), dimension(:), intent(inout) :: LocPos
          !real(REAL_TYPE), dimension(:,:), intent(inout) :: HS
          !real(REAL_TYPE), dimension(:,:,:), intent(inout) :: dHS
          !real(REAL_TYPE), dimension(:), intent(inout) :: Wt
          !
          !! local variables
          !real(REAL_TYPE) :: Xi, Eta
          !integer(INTEGER_TYPE) :: int, I1, Nint1
          
          ! Note this is two dimensional 
        
          !real(REAL_TYPE), dimension(:, :), intent(inout) :: HS
          !real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS
          !real(REAL_TYPE), dimension(:), intent(inout) :: Wt
        
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt !these should not be allocatables at this point  
           
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS_Xi !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS_Xi !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt_Xi !these should not be allocatables at this point  
           
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS_Eta !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS_Eta !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt_Eta !these should not be allocatables at this point  
           
           integer(INTEGER_TYPE), intent(inout) :: ni, nj
          
          !NURBS related inputs in the xi direction 
          integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
          integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
          real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
          
          real(REAL_TYPE), intent(in), dimension(NXiGaussPoints)  :: Xi_ParametricDomain !, dimension(NXiGaussPoints) 
          
          !NURBS related inputs in the eta direction 
          integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
          integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
          real(REAL_TYPE), dimension(NEtaKnotEntries), intent(in) :: EtaKnotEntries
          
          real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Eta_ParametricDomain
          
          
          ! Basis functions allocatables 
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: NN_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dN_dxi_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:) :: RR
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dR_dxi
          real(REAL_TYPE), dimension(NXiKnotOrder+1) :: NN_WithoutZeroValues
          real(REAL_TYPE), dimension(NXiKnotOrder+1) :: dN_dxi_WithoutZeroValues
          
          
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: MM_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dM_deta_IncludesZeroValues         
          real(REAL_TYPE), dimension(NEtaKnotOrder+1) :: MM_WithoutZeroValues
          real(REAL_TYPE), dimension(NEtaKnotOrder+1) :: dM_deta_WithoutZeroValues

          
          ! local variables 
          integer(INTEGER_TYPE) :: counter, ww, kk
          integer(INTEGER_TYPE) :: ii, jj, loc_num 
          real(REAL_TYPE) :: sum_tot
          real(REAL_TYPE) :: sum_xi
          real(REAL_TYPE) :: sum_eta
          
          !integer(INTEGER_TYPE), dimension(ELEMENTNODES,NDIM) :: Indices_NURBS
          
          
          !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Xi
          !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Eta
          !integer(INTEGER_TYPE) :: nGP_xi, nGP_eta, nGP_zeta 
          !integer(INTEGER_TYPE) :: ee_NURBS
          !integer(INTEGER_TYPE) :: ni_NURBS
          !integer(INTEGER_TYPE) :: nj_NURBS
          !integer(INTEGER_TYPE) :: nk_NURBS
          !
          !
          !real(REAL_TYPE), allocatable, dimension(:) :: xi_tilde
          !real(REAL_TYPE), allocatable, dimension(:) :: eta_tilde
          !real(REAL_TYPE), allocatable, dimension(:) :: zeta_tilde
          
          integer(INTEGER_TYPE) :: IError, stat    
          
        
          
          
          
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: LL_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dL_dxi_IncludesZeroValues          
          
              ! - evaluate each basis function value at the gauss point 
              call Bspline_basis_and_deriv(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
                                    NN_IncludesZeroValues_Print, dN_dxi_IncludesZeroValues_Print) !output 
                                    
              call Bspline_basis_and_deriv(NEtaKnotOrder, NEtaKnotEntries, NEtaGaussPoints, Eta_ParametricDomain, EtaKnotEntries, & !input 
                                    MM_IncludesZeroValues_Print, dM_deta_IncludesZeroValues_Print) !output 
              
              !call Bspline_basis_and_deriv(nk_NURBS, ll_NURBS_NumberOfUnivariateEtaKnots, NZetaKnotOrder, NZetaKnotEntries, nGP_Zeta, Zeta_ParametricDomain, ZetaKnotEntries, & !input 
              !                      LL_IncludesZeroValues, dL_dxi_IncludesZeroValues) !output 
              
              
              ! do we need to update ni and nj here so that we can use a large courant number????????????
              !-loop over knot spans and find ni and nj 
              
              !! Xi
              !ii = 1
              !do 
              !    
              !    if    ( (XiKnotEntries(ii)<Xi_ParametricDomain(1)) .and. (Xi_ParametricDomain(1)<XiKnotEntries(ii+1)) )     then 
              !        exit 
              !    end if 
              !    ii = ii + 1
              !    
              !end do
              !ni = ii
              !
              !! Eta
              !ii = 1
              !do 
              !    
              !    if    ( (EtaKnotEntries(ii)<Eta_ParametricDomain(1)) .and. (Eta_ParametricDomain(1)<EtaKnotEntries(ii+1)) )     then 
              !        exit 
              !    end if 
              !    ii = ii + 1
              !    
              !end do
              !nj = ii 
              
              
              counter = 0
              ! Xi is analogous to the x-coordinate in the parametric domain 
              do jj = 1, NXiGaussPoints
                  do ii = ni, ni+NXiKnotOrder
                 counter = counter + 1
                 HS_Xi(jj,counter) = NN_IncludesZeroValues_Print(jj,ii,NXiKnotOrder+1)
                 dHS_Xi(jj,counter,1) = dN_dxi_IncludesZeroValues_Print(jj,ii,NXiKnotOrder+1) ! note that this is the derivative in the parameter space... might need to normalize this somehow and add that term to the jacobian 
                 Wt_Xi(jj) = 2.0/NXiGaussPoints ! this weight is wrong 
                  end do 
                  counter = 0
              end do 
              
              
              counter = 0
              ! Eta is analogous to the y-coordinate in the parametric domain 
              do jj = 1, NEtaGaussPoints
                  do ii = nj, nj+NEtaKnotOrder
                 counter = counter + 1
                 HS_Eta(jj,counter) = MM_IncludesZeroValues_Print(jj,ii,NEtaKnotOrder+1) 
                 dHS_Eta(jj,counter,1) = dM_deta_IncludesZeroValues_Print(jj,ii,NEtaKnotOrder+1) ! note that this is the derivative in the parameter space... might need to normalize this somehow and add that term to the jacobian  
                 Wt_Eta(jj) = 2.0/NEtaGaussPoints ! this weight is wrong 
                  end do
                  counter = 0
              end do 
              
                  
                  
              !    NXiKnotOrder+1, (2*NXiKnotOrder)+1
              !    counter = counter + 1
              !    ! picking out the non-zero terms for shape functions 
              !    HS_Xi(NXiGaussPoints,counter) = NN_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
              !    HS_Eta(NEtaGaussPoints,counter) = MM_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
              !    ! picking out the non-zero terms for shape function derivatives 
              !    dHS_Xi(NXiGaussPoints,counter,1) = dN_dxi_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
              !    dHS_Eta(NEtaGaussPoints,counter,1) = dM_deta_IncludesZeroValues(NEtaGaussPoints,ii,NEtaKnotOrder+1) 
              !    Wt_Xi(NXiGaussPoints) = 2.0
              !    Wt_Eta(NEtaGaussPoints) = 2.0
              !end do 
              
              allocate(RR    (NXiGaussPoints*NEtaGaussPoints, (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
              allocate(dR_dxi(NXiGaussPoints*NEtaGaussPoints, (NXiKnotOrder+1) * (NEtaKnotOrder+1), NDIM ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
          
              RR = 0.0
              dR_dxi = 0.0
              sum_tot = 0.0
              sum_xi = 0.0
              sum_eta = 0.0
            ! - need to include tensor product multiplication here between NN and MM 
            ! build numerator and denominators 
            
              loc_num = 0
              
              
              ! Indices to take into account when arranging the RR and dR_dxi matrices 
              
              !Indices_NURBS = reshape( (/  2, 2, &
              !                             1, 2, &
              !                             1, 1,  &
              !                             2, 1/), &
              !                          (/ 4, 2 /) )
              
              !Indices_NURBS = [2, 2,
              !                 1, 2,
              !                 1, 1,
              !                 2, 1]
              counter = 0
              do ww = 1, NEtaGaussPoints
                  do kk = 1, NXiGaussPoints
                        
                      loc_num=0
                      counter = counter + 1
                      
                      do jj = 0, NEtaKnotOrder!1, NEtaKnotOrder+1 !0, NEtaKnotOrder
                          do ii = 0, NXiKnotOrder !0, NXiKnotOrder !1, NXiKnotOrder+1
                      
                      
                      ! shape functions
                      !RR(NXiGaussPoints, loc_num) = HS_Xi(NXiGaussPoints,Indices_NURBS(ii,jj)) * HS_Eta(NXiGaussPoints,Indices_NURBS(ii,jj))
                      
                              loc_num = loc_num + 1
                              RR(counter, loc_num) = HS_Xi(kk,NXiKnotOrder+1-ii) * HS_Eta(ww,NEtaKnotOrder+1-jj)
                      
                      
                      ! shape function derivatives 
                      !dR_dxi(NXiGaussPoints,loc_num,1) = dHS_Xi(NXiGaussPoints,Indices_NURBS(ii,jj),1) * HS_Eta(NEtaGaussPoints,Indices_NURBS(ii,jj))
                      !dR_dxi(NXiGaussPoints,loc_num,2) = HS_Xi(NXiGaussPoints,Indices_NURBS(ii,jj)) * dHS_Eta(NEtaGaussPoints,Indices_NURBS(ii,jj),1)
                      
                      dR_dxi(counter,loc_num,1) = dHS_Xi(kk,NXiKnotOrder+1-ii,1) * HS_Eta(ww,NEtaKnotOrder+1-jj)
                      dR_dxi(counter,loc_num,2) = HS_Xi(kk,NXiKnotOrder+1-ii) * dHS_Eta(ww,NEtaKnotOrder+1-jj,1)
                      
                      ! these are required when we are using weights 
                      sum_tot = (sum_tot + RR(counter, loc_num) )/(NXiGaussPoints*NEtaGaussPoints)
                      sum_xi = sum_xi + dR_dxi(counter,loc_num,1)
                      sum_eta = sum_eta + dR_dxi(counter,loc_num,2)
                end do     
                      end do
                      
                      !counter = 0
                  end do
              end do
              
            
              !allocate(HS_Xi(NXiKnotOrder+1), stat=IError)
              !allocate(dHS_Xi(NXiKnotOrder+1), stat=IError)
              !
              !allocate(HS_Eta(NEtaKnotOrder+1), stat=IError)
              !allocate(dHS_Eta(NEtaKnotOrder+1), stat=IError)
           
              !allocate(HS_Zeta(NZetaKnotOrder+1), stat=IError)
              !allocate(dHS_Zeta(NZetaKnotOrder+1), stat=IError)
              
          
          !Nint1=1 !number of gauss points
          !Int = 0 !counter
          !
          !do I1 = 1, Nint1
          !
          !    Xi = 0.0 !local position in Xi (local) direction 
          !    Eta = 0.0 !local position in Eta (local) direction
          !    
          !    Int = Int+1
          !    
          !    Wt(Int) = 2.0 !1d0 / Nint1 * 0.5 !This should be =2... double check!!!! 
          ! 
          !    ! HS(i)
          !    HS(Int, 1) = (1.0 - Xi) * (1.0 - Eta) / 4.0 ! a=1
          !    HS(Int, 2) = (1.0 + Xi) * (1.0 - Eta) / 4.0 ! a=2
          !    HS(Int, 3) = (1.0 + Xi) * (1.0 + Eta) / 4.0 ! a=3
          !    HS(Int, 4) = (1.0 - Xi) * (1.0 + Eta) / 4.0 ! a=4
          !
          !    ! dHS(i,1) = dHS / dXi
          !    dHS(Int,1,1) =  - (1.0 - Eta) / 4.0 ! a=1
          !    dHS(Int,2,1) =    (1.0 - Eta) / 4.0 ! a=2
          !    dHS(Int,3,1) =    (1.0 + Eta) / 4.0 ! a=3
          !    dHS(Int,4,1) =  - (1.0 + Eta) / 4.0 ! a=4
          !
          !    ! dHS(i,2) = dHS / dEta
          !    dHS(Int,1,2) =  - (1.0 - Xi) / 4.0 ! a=1
          !    dHS(Int,2,2) =  - (1.0 + Xi) / 4.0 ! a=2
          !    dHS(Int,3,2) =    (1.0 + Xi) / 4.0 ! a=3
          !    dHS(Int,4,2) =    (1.0 - Xi) / 4.0 ! a=4
          !
          !end do
              
              
              !allocate(HS    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
              !allocate(dHS( (NXiKnotOrder+1) * (NEtaKnotOrder+1), 2 ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
              !allocate(Wt    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
              
              
              !allocate(NN_IncludesZeroValues(NXiKnotEntries_uKnot-1), stat = IError)
              !allocate(dN_dxi_IncludesZeroValues(NXiKnotEntries_uKnot-1), stat = IError)
              !
              !
              !allocate(MM_IncludesZeroValues(NEtaKnotEntries_uKnot-1), stat = IError)
              !allocate(dM_deta_IncludesZeroValues(NEtaKnotEntries_uKnot-1), stat = IError)
              

              ! -----------------------------------------------------------------------------------

              !! initialize 
              !NN_IncludesZeroValues_Print = 0.0
              !dN_dxi_IncludesZeroValues_Print = 0.0
              !
              !MM_IncludesZeroValues_Print = 0.0
              !dM_deta_IncludesZeroValues_Print = 0.0
              !
              !! write debug parameters 
              !NN_IncludesZeroValues_Print = NN_IncludesZeroValues(1,:, NXiKnotOrder+1)
              !dN_dxi_IncludesZeroValues_Print = dN_dxi_IncludesZeroValues(1,:, NXiKnotOrder+1)
              !
              !MM_IncludesZeroValues_Print = MM_IncludesZeroValues(1,:, NEtaKnotOrder+1)
              !dM_deta_IncludesZeroValues_Print = dM_deta_IncludesZeroValues(1,:, NEtaKnotOrder+1)
              
              ! -----------------------------------------------------------------------------------
              
              HS = RR
              dHS = dR_dxi
              
              counter = 0
              do ww = 1, NEtaGaussPoints
                  do kk = 1, NXiGaussPoints
                      
                      counter = counter + 1
                      Wt(counter) = Wt_Xi(kk) * Wt_Eta(ww)
                      
                  end do 
              end do 
              
          
        
                                                    end subroutine InitialiseShapeFunctionsQUAD4_NURBS
                                                    
                                                    
                                                    
        !-----------------------------------------------------------------------------------------------
        ! SINGLE PARTICLE 
        subroutine InitialiseShapeFunctionsQUAD4_NURBS_SINGLEPARTICLE(HS, dHS, Wt, & !classic inout parameters
                                                    HS_Xi, dHS_Xi, Wt_Xi, &
                                                    HS_Eta, dHS_Eta, Wt_Eta, &
                                                    XiKnotEntries, NXiKnotEntries, Xi_ParametricDomain, NXiKnotOrder, & !NURBS related inputs in the xi direction 
                                                    EtaKnotEntries, NEtaKnotEntries, Eta_ParametricDomain, NEtaKnotOrder, &
                                                    ni, nj) !NURBS related inputs in the eta direction 
        !**********************************************************************
        !
        !    SUBROUTINE: InitialiseShapeFunctionsQUAD4_NURBS
        !
        !    DESCRIPTION:
        !>   To calculate the values of shape functions and their
        !>   derivatives at  one Gaussian integration point for a 4-noded 2D quadrilateral element using NURBS.
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
        
        implicit none
        
          !!real(REAL_TYPE), dimension(:), intent(inout) :: LocPos
          !real(REAL_TYPE), dimension(:,:), intent(inout) :: HS
          !real(REAL_TYPE), dimension(:,:,:), intent(inout) :: dHS
          !real(REAL_TYPE), dimension(:), intent(inout) :: Wt
          !
          !! local variables
          !real(REAL_TYPE) :: Xi, Eta
          !integer(INTEGER_TYPE) :: int, I1, Nint1
          
          ! Note this is two dimensional 
        
          !real(REAL_TYPE), dimension(:, :), intent(inout) :: HS
          !real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS
          !real(REAL_TYPE), dimension(:), intent(inout) :: Wt
        
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt !these should not be allocatables at this point  
           
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS_Xi !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS_Xi !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt_Xi !these should not be allocatables at this point  
           
           real(REAL_TYPE), dimension(:, :), intent(inout) :: HS_Eta !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:, :, :), intent(inout) :: dHS_Eta !these should not be allocatables at this point 
           real(REAL_TYPE), dimension(:), intent(inout) :: Wt_Eta !these should not be allocatables at this point  
           
           integer(INTEGER_TYPE), intent(inout) :: ni, nj
          
          !NURBS related inputs in the xi direction 
          integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
          integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
          real(REAL_TYPE), dimension(1), intent(in) :: XiKnotEntries
          
          real(REAL_TYPE), intent(in), dimension(1)  :: Xi_ParametricDomain !, dimension(NXiGaussPoints) 
          
          !NURBS related inputs in the eta direction 
          integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
          integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
          real(REAL_TYPE), dimension(NEtaKnotEntries), intent(in) :: EtaKnotEntries
          
          real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Eta_ParametricDomain
          
          
          ! Basis functions allocatables 
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: NN_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dN_dxi_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:) :: RR
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dR_dxi
          real(REAL_TYPE), dimension(NXiKnotOrder+1) :: NN_WithoutZeroValues
          real(REAL_TYPE), dimension(NXiKnotOrder+1) :: dN_dxi_WithoutZeroValues
          
          
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: MM_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dM_deta_IncludesZeroValues         
          real(REAL_TYPE), dimension(NEtaKnotOrder+1) :: MM_WithoutZeroValues
          real(REAL_TYPE), dimension(NEtaKnotOrder+1) :: dM_deta_WithoutZeroValues

          
          ! local variables 
          integer(INTEGER_TYPE) :: counter, ww, kk
          integer(INTEGER_TYPE) :: ii, jj, loc_num 
          real(REAL_TYPE) :: sum_tot
          real(REAL_TYPE) :: sum_xi
          real(REAL_TYPE) :: sum_eta
          
          !integer(INTEGER_TYPE), dimension(ELEMENTNODES,NDIM) :: Indices_NURBS
          
          
          !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Xi
          !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Eta
          !integer(INTEGER_TYPE) :: nGP_xi, nGP_eta, nGP_zeta 
          !integer(INTEGER_TYPE) :: ee_NURBS
          !integer(INTEGER_TYPE) :: ni_NURBS
          !integer(INTEGER_TYPE) :: nj_NURBS
          !integer(INTEGER_TYPE) :: nk_NURBS
          !
          !
          !real(REAL_TYPE), allocatable, dimension(:) :: xi_tilde
          !real(REAL_TYPE), allocatable, dimension(:) :: eta_tilde
          !real(REAL_TYPE), allocatable, dimension(:) :: zeta_tilde
          
          integer(INTEGER_TYPE) :: IError, stat    
          
        
          
          
          
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: LL_IncludesZeroValues
          real(REAL_TYPE), allocatable, dimension(:,:,:) :: dL_dxi_IncludesZeroValues          
          
              ! - evaluate each basis function value at the gauss point 
              call Bspline_basis_and_deriv_SINGLEPARTICLE(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
                                    NN_IncludesZeroValues_Print, dN_dxi_IncludesZeroValues_Print) !output 
                                    
              call Bspline_basis_and_deriv_SINGLEPARTICLE(NEtaKnotOrder, NEtaKnotEntries, NEtaGaussPoints, Eta_ParametricDomain, EtaKnotEntries, & !input 
                                    MM_IncludesZeroValues_Print, dM_deta_IncludesZeroValues_Print) !output 
              
              !call Bspline_basis_and_deriv(nk_NURBS, ll_NURBS_NumberOfUnivariateEtaKnots, NZetaKnotOrder, NZetaKnotEntries, nGP_Zeta, Zeta_ParametricDomain, ZetaKnotEntries, & !input 
              !                      LL_IncludesZeroValues, dL_dxi_IncludesZeroValues) !output 
              
              
              ! do we need to update ni and nj here so that we can use a large courant number????????????
              !-loop over knot spans and find ni and nj 
              
              !! Xi
              !ii = 1
              !do 
              !    
              !    if    ( (XiKnotEntries(ii)<Xi_ParametricDomain(1)) .and. (Xi_ParametricDomain(1)<XiKnotEntries(ii+1)) )     then 
              !        exit 
              !    end if 
              !    ii = ii + 1
              !    
              !end do
              !ni = ii
              !
              !! Eta
              !ii = 1
              !do 
              !    
              !    if    ( (EtaKnotEntries(ii)<Eta_ParametricDomain(1)) .and. (Eta_ParametricDomain(1)<EtaKnotEntries(ii+1)) )     then 
              !        exit 
              !    end if 
              !    ii = ii + 1
              !    
              !end do
              !nj = ii 
              
              
              counter = 0
              ! Xi is analogous to the x-coordinate in the parametric domain 
              do jj = 1, NXiGaussPoints
                  do ii = ni, ni+NXiKnotOrder
                 counter = counter + 1
                 HS_Xi(jj,counter) = NN_IncludesZeroValues_Print(jj,ii,NXiKnotOrder+1)
                 dHS_Xi(jj,counter,1) = dN_dxi_IncludesZeroValues_Print(jj,ii,NXiKnotOrder+1) ! note that this is the derivative in the parameter space... might need to normalize this somehow and add that term to the jacobian 
                 Wt_Xi(jj) = 2.0/NXiGaussPoints ! this weight is wrong 
                  end do 
                  counter = 0
              end do 
              
              
              counter = 0
              ! Eta is analogous to the y-coordinate in the parametric domain 
              do jj = 1, NEtaGaussPoints
                  do ii = nj, nj+NEtaKnotOrder
                 counter = counter + 1
                 HS_Eta(jj,counter) = MM_IncludesZeroValues_Print(jj,ii,NEtaKnotOrder+1) 
                 dHS_Eta(jj,counter,1) = dM_deta_IncludesZeroValues_Print(jj,ii,NEtaKnotOrder+1) ! note that this is the derivative in the parameter space... might need to normalize this somehow and add that term to the jacobian  
                 Wt_Eta(jj) = 2.0/NEtaGaussPoints ! this weight is wrong 
                  end do
                  counter = 0
              end do 
              
                  
                  
              !    NXiKnotOrder+1, (2*NXiKnotOrder)+1
              !    counter = counter + 1
              !    ! picking out the non-zero terms for shape functions 
              !    HS_Xi(NXiGaussPoints,counter) = NN_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
              !    HS_Eta(NEtaGaussPoints,counter) = MM_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
              !    ! picking out the non-zero terms for shape function derivatives 
              !    dHS_Xi(NXiGaussPoints,counter,1) = dN_dxi_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
              !    dHS_Eta(NEtaGaussPoints,counter,1) = dM_deta_IncludesZeroValues(NEtaGaussPoints,ii,NEtaKnotOrder+1) 
              !    Wt_Xi(NXiGaussPoints) = 2.0
              !    Wt_Eta(NEtaGaussPoints) = 2.0
              !end do 
              
              allocate(RR    (NXiGaussPoints*NEtaGaussPoints, (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
              allocate(dR_dxi(NXiGaussPoints*NEtaGaussPoints, (NXiKnotOrder+1) * (NEtaKnotOrder+1), NDIM ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
          
              RR = 0.0
              dR_dxi = 0.0
              sum_tot = 0.0
              sum_xi = 0.0
              sum_eta = 0.0
            ! - need to include tensor product multiplication here between NN and MM 
            ! build numerator and denominators 
            
              loc_num = 0
              
              
              ! Indices to take into account when arranging the RR and dR_dxi matrices 
              
              !Indices_NURBS = reshape( (/  2, 2, &
              !                             1, 2, &
              !                             1, 1,  &
              !                             2, 1/), &
              !                          (/ 4, 2 /) )
              
              !Indices_NURBS = [2, 2,
              !                 1, 2,
              !                 1, 1,
              !                 2, 1]
              counter = 0
              do ww = 1, NEtaGaussPoints
                  do kk = 1, NXiGaussPoints
                        
                      loc_num=0
                      counter = counter + 1
                      
                      do jj = 0, NEtaKnotOrder!1, NEtaKnotOrder+1 !0, NEtaKnotOrder
                          do ii = 0, NXiKnotOrder !0, NXiKnotOrder !1, NXiKnotOrder+1
                      
                      
                      ! shape functions
                      !RR(NXiGaussPoints, loc_num) = HS_Xi(NXiGaussPoints,Indices_NURBS(ii,jj)) * HS_Eta(NXiGaussPoints,Indices_NURBS(ii,jj))
                      
                              loc_num = loc_num + 1
                              RR(counter, loc_num) = HS_Xi(kk,NXiKnotOrder+1-ii) * HS_Eta(ww,NEtaKnotOrder+1-jj)
                      
                      
                      ! shape function derivatives 
                      !dR_dxi(NXiGaussPoints,loc_num,1) = dHS_Xi(NXiGaussPoints,Indices_NURBS(ii,jj),1) * HS_Eta(NEtaGaussPoints,Indices_NURBS(ii,jj))
                      !dR_dxi(NXiGaussPoints,loc_num,2) = HS_Xi(NXiGaussPoints,Indices_NURBS(ii,jj)) * dHS_Eta(NEtaGaussPoints,Indices_NURBS(ii,jj),1)
                      
                      dR_dxi(counter,loc_num,1) = dHS_Xi(kk,NXiKnotOrder+1-ii,1) * HS_Eta(ww,NEtaKnotOrder+1-jj)
                      dR_dxi(counter,loc_num,2) = HS_Xi(kk,NXiKnotOrder+1-ii) * dHS_Eta(ww,NEtaKnotOrder+1-jj,1)
                      
                      ! these are required when we are using weights 
                      sum_tot = (sum_tot + RR(counter, loc_num) )/(NXiGaussPoints*NEtaGaussPoints)
                      sum_xi = sum_xi + dR_dxi(counter,loc_num,1)
                      sum_eta = sum_eta + dR_dxi(counter,loc_num,2)
                end do     
                      end do
                      
                      !counter = 0
                  end do
              end do
              
            
              !allocate(HS_Xi(NXiKnotOrder+1), stat=IError)
              !allocate(dHS_Xi(NXiKnotOrder+1), stat=IError)
              !
              !allocate(HS_Eta(NEtaKnotOrder+1), stat=IError)
              !allocate(dHS_Eta(NEtaKnotOrder+1), stat=IError)
           
              !allocate(HS_Zeta(NZetaKnotOrder+1), stat=IError)
              !allocate(dHS_Zeta(NZetaKnotOrder+1), stat=IError)
              
          
          !Nint1=1 !number of gauss points
          !Int = 0 !counter
          !
          !do I1 = 1, Nint1
          !
          !    Xi = 0.0 !local position in Xi (local) direction 
          !    Eta = 0.0 !local position in Eta (local) direction
          !    
          !    Int = Int+1
          !    
          !    Wt(Int) = 2.0 !1d0 / Nint1 * 0.5 !This should be =2... double check!!!! 
          ! 
          !    ! HS(i)
          !    HS(Int, 1) = (1.0 - Xi) * (1.0 - Eta) / 4.0 ! a=1
          !    HS(Int, 2) = (1.0 + Xi) * (1.0 - Eta) / 4.0 ! a=2
          !    HS(Int, 3) = (1.0 + Xi) * (1.0 + Eta) / 4.0 ! a=3
          !    HS(Int, 4) = (1.0 - Xi) * (1.0 + Eta) / 4.0 ! a=4
          !
          !    ! dHS(i,1) = dHS / dXi
          !    dHS(Int,1,1) =  - (1.0 - Eta) / 4.0 ! a=1
          !    dHS(Int,2,1) =    (1.0 - Eta) / 4.0 ! a=2
          !    dHS(Int,3,1) =    (1.0 + Eta) / 4.0 ! a=3
          !    dHS(Int,4,1) =  - (1.0 + Eta) / 4.0 ! a=4
          !
          !    ! dHS(i,2) = dHS / dEta
          !    dHS(Int,1,2) =  - (1.0 - Xi) / 4.0 ! a=1
          !    dHS(Int,2,2) =  - (1.0 + Xi) / 4.0 ! a=2
          !    dHS(Int,3,2) =    (1.0 + Xi) / 4.0 ! a=3
          !    dHS(Int,4,2) =    (1.0 - Xi) / 4.0 ! a=4
          !
          !end do
              
              
              !allocate(HS    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
              !allocate(dHS( (NXiKnotOrder+1) * (NEtaKnotOrder+1), 2 ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
              !allocate(Wt    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
              
              
              !allocate(NN_IncludesZeroValues(NXiKnotEntries_uKnot-1), stat = IError)
              !allocate(dN_dxi_IncludesZeroValues(NXiKnotEntries_uKnot-1), stat = IError)
              !
              !
              !allocate(MM_IncludesZeroValues(NEtaKnotEntries_uKnot-1), stat = IError)
              !allocate(dM_deta_IncludesZeroValues(NEtaKnotEntries_uKnot-1), stat = IError)
              

              ! -----------------------------------------------------------------------------------

              !! initialize 
              !NN_IncludesZeroValues_Print = 0.0
              !dN_dxi_IncludesZeroValues_Print = 0.0
              !
              !MM_IncludesZeroValues_Print = 0.0
              !dM_deta_IncludesZeroValues_Print = 0.0
              !
              !! write debug parameters 
              !NN_IncludesZeroValues_Print = NN_IncludesZeroValues(1,:, NXiKnotOrder+1)
              !dN_dxi_IncludesZeroValues_Print = dN_dxi_IncludesZeroValues(1,:, NXiKnotOrder+1)
              !
              !MM_IncludesZeroValues_Print = MM_IncludesZeroValues(1,:, NEtaKnotOrder+1)
              !dM_deta_IncludesZeroValues_Print = dM_deta_IncludesZeroValues(1,:, NEtaKnotOrder+1)
              
              ! -----------------------------------------------------------------------------------
              
              HS = RR
              dHS = dR_dxi
              
              counter = 0
              do ww = 1, NEtaGaussPoints
                  do kk = 1, NXiGaussPoints
                      
                      counter = counter + 1
                      Wt(counter) = Wt_Xi(ww) * Wt_Eta(kk)
                      
                  end do 
              end do 
              
          
        
            end subroutine InitialiseShapeFunctionsQUAD4_NURBS_SINGLEPARTICLE
                                                    
                                                    
        ! SINGLE PARTICLE 
        !-----------------------------------------------------------------------------------------------
                                                    
                                                    
                                                    
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP1(NElementParticles, ParticleStatus)
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
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

          ParticleStatus(1) = .true.

        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP1
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP4(NElementParticles, ParticleStatus)
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
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NElementParticles
          logical, dimension(NElementParticles), intent(inout) :: ParticleStatus

          ParticleStatus(1) = .true.
          ParticleStatus(2) = .true.
          ParticleStatus(3) = .true.
          ParticleStatus(4) = .true.

        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP4
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP9(NElementParticles, ParticleStatus)
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

        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP9
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP16(NElementParticles, ParticleStatus)
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


        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP16
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP25(NElementParticles, ParticleStatus)
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
          

        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP25
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP36(NElementParticles, ParticleStatus)
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
          

        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP36
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP49(NElementParticles, ParticleStatus)
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


        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP49
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP64(NElementParticles, ParticleStatus)
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



        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP64
        
        
        
        
        subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP81(NElementParticles, ParticleStatus)
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
          



        end subroutine DetermineAdjacentParticlesQUAD4_NURBS_MP81
        
        
        !subroutine GradientOfMappingFromParameterSpaceToPhysicalSpace & !NURBS
        !( dR_dxi, NodalCorrdinates, ni, nj, dR_dxi, NEtaKnotOrder, NXiKnotOrder !inputs 
        !dx_dxi ) ! output 
        !
        !! This code is part of the NURBS implementation
        !! 
        !
        !implicit none 
        !
        !!input 
        !integer(INTEGER_TYPE), intent(in) :: ni, nj 
        !integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
        !integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
        !real(REAL_TYPE), intent(in) :: dR_dxi
        !real(REAL_TYPE), intent(in), dimension(:,:) :: NodalCoordinates ! these are control points
        !
        !!output 
        !real(REAL_TYPE), intent(out), dimension(:,:) :: dx_dxi
        !
        !!local 
        !integer(INTEGER_TYPE) :: ii, jj
        !integer(INTEGER_TYPE) :: loc_num, NumberOfDimensions
        !integer(INTEGER_TYPE) :: aa, bb
        !
        !
        !
        !NumberOfDimensions = 2 
        !
        !loc_num = 0 
        !
        !do jj = 0, NEtaKnotOrder
        !    do ii = 0, NXiKnotOrder
        !        loc_num = loc_num + 1 
        !        
        !        do aa = 1, NumberOfDimensions
        !            do bb = 1, NumberOfDimensions
        !               dx_dxi(aa,bb) = dx_dxi(aa,bb) + &
        !                            NodalCoordinates(ni - ii, nj - jj) * dR_dxi(loc_num, bb) ! control points * dR_dxi
        !                
        !            end do 
        !        end do     
        !            
        !        
        !        
        !        
        !        
        !        
        !    end do 
        !end do 
        !
        !
        !
        !end subroutine GradientOfMappingFromParameterSpaceToPhysicalSpace
                                                    
                                                    
                                                    
                                                    
                                                    
                                                    
                                                    
                                                    
        !                                            
        !                                            
        !                                            
        !                                            
        !
        !!**********************************************************************
        !!
        !!    SUBROUTINE: ShapeLocPosQUAD4
        !!
        !!    DESCRIPTION:
        !!>   To calculate the values of shape functions and their
        !!>   derivatives at LocPos for a 4-noded 2D quadrilateral element.
        !!
        !!>   @note : 2D element
        !!>   @note : https://ses.library.usyd.edu.au/bitstream/2123/709/8/adt-NU20060210.15574814appendixD.pdf
        !!>   @note : R. K. Livesley, Finite Elements: An Introduction for Engineers, CUP Archive 1983
        !!
        !!>   @param[in] LocPos : Local coordinates of a point inside an element
        !!
        !!>   @param[out] HS(i) : Value of shape function i at LocPos
        !!>   @param[out] dHS(i,j) : Value of derivative of shape function i at LocPos with respect to direction j
        !!
        !!             4) (-1,1)   ^ Eta    3) (1,1)
        !!                 4       |
        !!                +---------------+ 3
        !!                |        |      |
        !!                |        |      |
        !!                |        |      |
        !!                |        -------|---> Xi
        !!                |               |
        !!                |               |
        !!                |1              | 2
        !!                +---------------+-
        !!             1) (-1,-1)           2) (-1,1)
        !!**********************************************************************
        !subroutine ShapeLocPosQUAD4_NURBS(LocPos, HS, dHS)
        !
        !implicit none
        !
        !
        !
        !
        !
        !
        !   real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: HS 
        !   real(REAL_TYPE), allocatable, dimension(:, :), intent(inout) :: dHS 
        !   real(REAL_TYPE), allocatable, dimension(:), intent(inout) :: Wt ! 
        !  
        !  
        !  !NURBS related inputs in the xi direction 
        !  integer(INTEGER_TYPE), intent(in) :: NXiKnotOrder
        !  integer(INTEGER_TYPE), intent(in) :: NXiKnotEntries
        !  real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: XiKnotEntries
        !  
        !  real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Xi_ParametricDomain
        !  
        !  !NURBS related inputs in the eta direction 
        !  integer(INTEGER_TYPE), intent(in) :: NEtaKnotOrder
        !  integer(INTEGER_TYPE), intent(in) :: NEtaKnotEntries
        !  real(REAL_TYPE), dimension(NXiKnotEntries), intent(in) :: EtaKnotEntries
        !  
        !  real(REAL_TYPE), intent(in), dimension(NXiGaussPoints) :: Eta_ParametricDomain
        !  
        !  
        !  ! Basis functions allocatables 
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: NN_IncludesZeroValues
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dN_dxi_IncludesZeroValues
        !  real(REAL_TYPE), allocatable, dimension(:) :: RR
        !  real(REAL_TYPE), allocatable, dimension(:,:) :: dR_dxi
        !  real(REAL_TYPE), dimension(NXiKnotOrder+1) :: NN_WithoutZeroValues
        !  real(REAL_TYPE), dimension(NXiKnotOrder+1) :: dN_dxi_WithoutZeroValues
        !  
        !  
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: MM_IncludesZeroValues
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dM_deta_IncludesZeroValues         
        !  real(REAL_TYPE), dimension(NEtaKnotOrder+1) :: MM_WithoutZeroValues
        !  real(REAL_TYPE), dimension(NEtaKnotOrder+1) :: dM_deta_WithoutZeroValues
        !
        !  
        !  ! local variables 
        !  integer(INTEGER_TYPE) :: counter
        !  integer(INTEGER_TYPE) :: ii, jj, loc_num 
        !  real(REAL_TYPE) :: sum_tot
        !  real(REAL_TYPE) :: sum_xi
        !  real(REAL_TYPE) :: sum_eta
        !  
        !  
        !  
        !  !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Xi
        !  !integer(INTEGER_TYPE) :: Number_of_Knot_Spans_Eta
        !  !integer(INTEGER_TYPE) :: nGP_xi, nGP_eta, nGP_zeta 
        !  !integer(INTEGER_TYPE) :: ee_NURBS
        !  !integer(INTEGER_TYPE) :: ni_NURBS
        !  !integer(INTEGER_TYPE) :: nj_NURBS
        !  !integer(INTEGER_TYPE) :: nk_NURBS
        !  !
        !  !
        !  !real(REAL_TYPE), allocatable, dimension(:) :: xi_tilde
        !  !real(REAL_TYPE), allocatable, dimension(:) :: eta_tilde
        !  !real(REAL_TYPE), allocatable, dimension(:) :: zeta_tilde
        !  
        !  integer(INTEGER_TYPE) :: IError, stat    
        !  
        !
        !  
        !  
        !  
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: LL_IncludesZeroValues
        !  real(REAL_TYPE), allocatable, dimension(:,:,:) :: dL_dxi_IncludesZeroValues          
        !  
        !      ! - evaluate each basis function value at the gauss point 
        !      call Bspline_basis_and_deriv(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
        !                            NN_IncludesZeroValues, dN_dxi_IncludesZeroValues) !output 
        !                            
        !      call Bspline_basis_and_deriv(NEtaKnotOrder, NEtaKnotEntries, NEtaGaussPoints, Eta_ParametricDomain, EtaKnotEntries, & !input 
        !                            MM_IncludesZeroValues, dM_deta_IncludesZeroValues) !output 
        !      
        !      !call Bspline_basis_and_deriv(nk_NURBS, ll_NURBS_NumberOfUnivariateEtaKnots, NZetaKnotOrder, NZetaKnotEntries, nGP_Zeta, Zeta_ParametricDomain, ZetaKnotEntries, & !input 
        !      !                      LL_IncludesZeroValues, dL_dxi_IncludesZeroValues) !output 
        !      counter = 0
        !      
        !      do ii = NXiKnotOrder+1, (2*NXiKnotOrder)+1
        !          counter = counter + 1
        !          ! picking out the non-zero terms for shape functions 
        !          NN_WithoutZeroValues(counter) = NN_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
        !          MM_WithoutZeroValues(counter) = MM_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
        !          ! picking out the non-zero terms for shape function derivatives 
        !          dN_dxi_WithoutZeroValues(counter) = dN_dxi_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
        !          dM_deta_WithoutZeroValues(counter) = dM_deta_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1) 
        !      end do 
        !      
        !      allocate(RR    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
        !      allocate(dR_dxi( (NXiKnotOrder+1) * (NEtaKnotOrder+1), 2 ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
        !  
        !    ! - need to include tensor product multiplication here between NN and MM 
        !    ! build numerator and denominators 
        !    
        !      do jj = 0, NEtaKnotOrder
        !          do ii = 0, NXiKnotOrder
        !              loc_num = loc_num + 1
        !              ! shape functions
        !              RR(loc_num) = NN_WithoutZeroValues(NXiKnotOrder+1-ii) * MM_WithoutZeroValues(NEtaKnotOrder+1-jj)
        !              ! shape function derivatives 
        !              dR_dxi(loc_num,1) = dN_dxi_WithoutZeroValues(NXiKnotOrder+1-ii) * MM_WithoutZeroValues(NEtaKnotOrder+1-jj)
        !              dR_dxi(loc_num,2) = NN_WithoutZeroValues(NXiKnotOrder+1-ii) * dM_deta_WithoutZeroValues(NEtaKnotOrder+1-jj)
        !              ! these are required when we are using weights 
        !              sum_tot = sum_tot + RR(loc_num) 
        !              sum_xi = sum_xi + dR_dxi(loc_num,1)
        !              sum_eta = sum_eta + dR_dxi(loc_num,2)
        !        end do     
        !      end do
        !    
        !      !allocate(HS_Xi(NXiKnotOrder+1), stat=IError)
        !      !allocate(dHS_Xi(NXiKnotOrder+1), stat=IError)
        !      !
        !      !allocate(HS_Eta(NEtaKnotOrder+1), stat=IError)
        !      !allocate(dHS_Eta(NEtaKnotOrder+1), stat=IError)
        !   
        !      !allocate(HS_Zeta(NZetaKnotOrder+1), stat=IError)
        !      !allocate(dHS_Zeta(NZetaKnotOrder+1), stat=IError)
        !      
        !  
        !  !Nint1=1 !number of gauss points
        !  !Int = 0 !counter
        !  !
        !  !do I1 = 1, Nint1
        !  !
        !  !    Xi = 0.0 !local position in Xi (local) direction 
        !  !    Eta = 0.0 !local position in Eta (local) direction
        !  !    
        !  !    Int = Int+1
        !  !    
        !  !    Wt(Int) = 2.0 !1d0 / Nint1 * 0.5 !This should be =2... double check!!!! 
        !  ! 
        !  !    ! HS(i)
        !  !    HS(Int, 1) = (1.0 - Xi) * (1.0 - Eta) / 4.0 ! a=1
        !  !    HS(Int, 2) = (1.0 + Xi) * (1.0 - Eta) / 4.0 ! a=2
        !  !    HS(Int, 3) = (1.0 + Xi) * (1.0 + Eta) / 4.0 ! a=3
        !  !    HS(Int, 4) = (1.0 - Xi) * (1.0 + Eta) / 4.0 ! a=4
        !  !
        !  !    ! dHS(i,1) = dHS / dXi
        !  !    dHS(Int,1,1) =  - (1.0 - Eta) / 4.0 ! a=1
        !  !    dHS(Int,2,1) =    (1.0 - Eta) / 4.0 ! a=2
        !  !    dHS(Int,3,1) =    (1.0 + Eta) / 4.0 ! a=3
        !  !    dHS(Int,4,1) =  - (1.0 + Eta) / 4.0 ! a=4
        !  !
        !  !    ! dHS(i,2) = dHS / dEta
        !  !    dHS(Int,1,2) =  - (1.0 - Xi) / 4.0 ! a=1
        !  !    dHS(Int,2,2) =  - (1.0 + Xi) / 4.0 ! a=2
        !  !    dHS(Int,3,2) =    (1.0 + Xi) / 4.0 ! a=3
        !  !    dHS(Int,4,2) =    (1.0 - Xi) / 4.0 ! a=4
        !  !
        !  !end do
        !      
        !      
        !      allocate(HS    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
        !      allocate(dHS( (NXiKnotOrder+1) * (NEtaKnotOrder+1), 2 ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
        !      allocate(Wt    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
        !      
        !      HS = RR
        !      dHS = dR_dxi
        !      Wt = 1.0
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
        !
        !! ------------- ORIGNAL IMPLEMENTAION
        !
        !  real(REAL_TYPE), dimension(:), intent(in) :: LocPos
        !  real(REAL_TYPE), dimension(:), intent(out) :: HS
        !  real(REAL_TYPE), dimension(:, :), intent(out) :: dHS
        !
        !  ! local variables
        !  real(REAL_TYPE) :: Xi, Eta
        !
        !  ! LocPos is an input here which allows the input of the local position of the material point 
        !  Xi_ParametricDomain = LocPos(1) !local position in Xi (local) direction 
        !  Eta_ParametricDomain = LocPos(2) !local position in Eta (local) direction
        !  
        !  
        !  ! HS(i) and dHS(i) 
        !  call Bspline_basis_and_deriv(NXiKnotOrder, NXiKnotEntries, NXiGaussPoints, Xi_ParametricDomain, XiKnotEntries, & !input 
        !                            NN_IncludesZeroValues, dN_dxi_IncludesZeroValues) !output 
        !                            
        !      
        !  call Bspline_basis_and_deriv(NEtaKnotOrder, NEtaKnotEntries, NEtaGaussPoints, Eta_ParametricDomain, EtaKnotEntries, & !input 
        !                            MM_IncludesZeroValues, dM_deta_IncludesZeroValues) !output 
        !      
        !  ! some rearrangement to pick out the basis functions that we need 
        !  do ii = NXiKnotOrder+1, (2*NXiKnotOrder)+1
        !      counter = counter + 1
        !      ! picking out the non-zero terms for shape functions 
        !      NN_WithoutZeroValues(counter) = NN_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
        !      MM_WithoutZeroValues(counter) = MM_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
        !      ! picking out the non-zero terms for shape function derivatives 
        !      dN_dxi_WithoutZeroValues(counter) = dN_dxi_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)
        !      dM_deta_WithoutZeroValues(counter) = dM_deta_IncludesZeroValues(NXiGaussPoints,ii,NXiKnotOrder+1)     
        !  end do 
        !      
        !  ! allocate variables 
        !  allocate(RR    ( (NXiKnotOrder+1) * (NEtaKnotOrder+1)), stat=IError) ! no of rows = 4, no of columns = 1 for linear element 
        !  allocate(dR_dxi( (NXiKnotOrder+1) * (NEtaKnotOrder+1), 2 ), stat=IError) ! no of rows = 4, no of columns = 2 for linear element 
        !  
        !  ! - need to include tensor product multiplication here between NN and MM 
        !  ! build numerator and denominators 
        !  do jj = 0, NEtaKnotOrder
        !      do ii = 0, NXiKnotOrder
        !          loc_num = loc_num + 1
        !          
        !          ! shape functions
        !          RR(loc_num) = NN_WithoutZeroValues(NXiKnotOrder+1-ii) * MM_WithoutZeroValues(NEtaKnotOrder+1-jj)
        !          
        !          ! shape function derivatives 
        !          dR_dxi(loc_num,1) = dN_dxi_WithoutZeroValues(NXiKnotOrder+1-ii) * MM_WithoutZeroValues(NEtaKnotOrder+1-jj)
        !          dR_dxi(loc_num,2) = NN_WithoutZeroValues(NXiKnotOrder+1-ii) * dM_deta_WithoutZeroValues(NEtaKnotOrder+1-jj)
        !          
        !          ! these are required when we are using weights 
        !          sum_tot = sum_tot + RR(loc_num) 
        !          sum_xi = sum_xi + dR_dxi(loc_num,1)
        !          sum_eta = sum_eta + dR_dxi(loc_num,2)
        !      end do         
        !  end do
        !  
        !  
        !  ! HS(i)
        !  HS(1) = RR(1)          
        !  HS(1) = RR(1)
        !  HS(1) = RR(1)          
        !  HS(1) = RR(1)          
        !  
        !  
        !  ! dHS(i,1) = dHS / dXi
        !  ! CHECK THESE INDEX NUMBERINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !  dHS(1,1) =  dR_dxi(1,1) !- (1.0 - Eta) / 4.0 ! a=1
        !  dHS(2,1) =  dR_dxi(2,1) !(1.0 - Eta) / 4.0 ! a=2
        !  dHS(3,1) =  dR_dxi(1,2) !(1.0 + Eta) / 4.0 ! a=3
        !  dHS(4,1) =  dR_dxi(1,2) !- (1.0 + Eta) / 4.0 ! a=4
        !
        !  !! HS(i)
        !  !HS(1) = (1.0 - Xi) * (1.0 - Eta) / 4.0 ! a=1
        !  !HS(2) = (1.0 + Xi) * (1.0 - Eta) / 4.0 ! a=2
        !  !HS(3) = (1.0 + Xi) * (1.0 + Eta) / 4.0 ! a=3
        !  !HS(4) = (1.0 - Xi) * (1.0 + Eta) / 4.0 ! a=4
        !  !
        !  !! dHS(i,1) = dHS / dXi
        !  !dHS(1,1) =  - (1.0 - Eta) / 4.0 ! a=1
        !  !dHS(2,1) =    (1.0 - Eta) / 4.0 ! a=2
        !  !dHS(3,1) =    (1.0 + Eta) / 4.0 ! a=3
        !  !dHS(4,1) =  - (1.0 + Eta) / 4.0 ! a=4
        !  !
        !  !! dHS(i,2) = dHS / dEta
        !  !dHS(1,2) =  - (1.0 - Xi) / 4.0 ! a=1
        !  !dHS(2,2) =  - (1.0 + Xi) / 4.0 ! a=2
        !  !dHS(3,2) =    (1.0 + Xi) / 4.0 ! a=3
        !  !dHS(4,2) =    (1.0 - Xi) / 4.0 ! a=4
        !
        !end subroutine ShapeLocPosQUAD4_NURBS
        !
        !                                            
        !                                            
        !                                            
                                                    
        
        
    subroutine FindControlPointNumbersForTractionApplication( NumberOfControlPointsForNURBSTraction, NURBSTractionNodes, ILoadCon, NumberOfTractionElements )
    
    !use ModCounters
    !  use ModElementEvaluation
    !  use ModGlobalConstants
    !  use ModMeshAdjacencies
    !    use ModMeshInfo
    
    
    implicit none
    
    ! initialize variables 
    ! input 
    integer(INTEGER_TYPE), intent(inout) :: NumberOfControlPointsForNURBSTraction
    integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: NURBSTractionNodes
    integer(INTEGER_TYPE), dimension(:), intent(in) :: ILoadCon
    ! output 
    integer(INTEGER_TYPE), intent(inout) :: NumberOfTractionElements
    ! local 
    integer(INTEGER_TYPE) :: MinimumControlPointNumber, MaximumControlPointNumber
    integer(INTEGER_TYPE) :: ii
    
    
    MinimumControlPointNumber = ILoadCon(1)
    MaximumControlPointNumber = ILoadCon(2)
    
    NumberOfControlPointsForNURBSTraction = MaximumControlPointNumber - MinimumControlPointNumber + 1
    NumberOfTractionElements = 1 !NumberOfControlPointsForNURBSTraction - (NXiKnotOrder+1) + 1 ! hardcoded
    allocate(NURBSTractionNodes ( NumberOfControlPointsForNURBSTraction ) )

    !NURBSTractionNodes(1) = MinimumControlPointNumber
    
    do ii = 1, NumberOfControlPointsForNURBSTraction ! start from 1 
        
        NURBSTractionNodes(ii) = MinimumControlPointNumber
        
        MinimumControlPointNumber = MinimumControlPointNumber + 1
 
    end do
    
        
        
    !    NumberOfControlPointsForNURBSTraction =  NumberOfControlPointsForNURBSTraction + 1
    !   
    !   if (NumberOfControlPointsForNURBSTraction>=MaximumControlPointNumber) then 
    !       
    !       exit 
    !       
    !       end if
    !
    !!end do 
    
    
    
    
    
    end subroutine FindControlPointNumbersForTractionApplication  
    
                                                    
                                                    
        
        
    
    end module ModNURBS