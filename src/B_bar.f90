module ModB_Bar
    
    use ModElementEvaluation, only : BMatrix
    use ModGlobalConstants  , only : INTEGER_TYPE, REAL_TYPE 
    
    implicit none

contains
    
    !*************************************************************************************
    !    SUBROUTINE: eval_B_bar_i
    ! 
    !    DESCRIPTION:
    !>   Evaluates Equation (4.5.17) of Hughes (1987) The Finite Element Method
    !    Eqn: $\bar{B}_{i}(\xi) = \sum_{\tilde{a} = 1}^{\tilde{n}_{int}} \tilde{N}_{\tilde{a}} (\xi) B_{i \tilde{a}}$
    !
    !>   @note : Notes
    ! 
    !    @param[in] reduced_gauss_order : Reduced Gauss integration order for the evaluation of /B_{i} - pg 234 Hughes (1987)
    !    @param[in] IElTyp              : Number of nodes per element
    !    @param[in] NEl                 : Total number of elements
    !    @param[in] NodTot              : Total number of nodes
    !    @param[in] IDim                : Number of dimensions
    !    @param[in] IElement            : ID of the element
    !    @param[in] ICon                : Element connectivities ICon(I, J): global node number of local node I in element J
    !    @param[in] NodeCoord           : Global nodal coordinates Co(I, J): J-coordinate of node I
    !>   @param[out] B_Bar_i : First three
    !
    !*************************************************************************************
    pure subroutine eval_B_bar_i(reduced_gauss_order, IElTyp, NEL, NodTot, IDim, IElement, ICon, NodeCoord, B_bar_i, HS)
        integer(INTEGER_TYPE), intent(in) :: reduced_gauss_order, IElTyp, NEL, NodTot, IDim, IElement
                                             
        real(REAL_TYPE), dimension(IElTyp, NEl) , intent(in)  :: ICon
        real(REAL_TYPE), dimension(NodTot, IDim), intent(in)  :: NodeCoord
        real(REAL_TYPE), dimension(IDim, IElTyp), intent(out) :: B_bar_i
        real(REAL_TYPE), dimension(IElTyp)      , intent(out) :: HS
        

        ! Local variables
        real(REAL_TYPE), dimension(:, :), allocatable :: gauss_point_locations 
        real(REAL_TYPE)                               :: DetJac
        integer(INTEGER_TYPE)                         :: one_point_gauss
        real(REAL_TYPE), dimension(IElTyp, IDim)      :: dHS ! Derivatives of the shape function
        ! Set local variables

        one_point_gauss   = 1 !: 2D, 3D: 1 gauss point at the element center
        
        ! Not implemented - Run the test and implement if necessary
        ! two_point_gauss   = 2 !: 2D: 2x2 points, 3D: 2x2 points
        ! three_point_gauss = 3 !: 2D: 3x3 points, 3D: 3x3 points

        select case(reduced_gauss_order)
            
        case(one_point_gauss)
            ! Only one gauss point needs to evaluated
            allocate(gauss_point_locations(model_dim, one_point_gauss))
            
            ! Set the location to zero 
            gauss_point_locations = 0.0
            
            ! Determine $\bar{B}_{i} (\xi) = B_{i}(0) i = 1, n_{sd}
            ! /B_{1}, /B_{2}, /B_{3} is in each column for each node in the element
            call BMatrix(gauss_point_locations(:, 1), NEl, NodTot, IDim, IElement, ICon, Co, RJac, InvRJac, DetJac)   
            
            call ShapeFunctionData(gauss_point_locations(:, 1), IElTyp, HS, dHS)
        case default
            ! TODO: Add error message here
            print *, "Add error message here about higher order gauss not being implemented"

    end subroutine eval_B_bar_i


    !*************************************************************************************
    !    subroutine:     Get_Strain
    !
    !    DESCRIPTION:
    !>   Calculates the strains at an integration point
    !
    !>   @param[in] IEl        : Element number
    !    @param[in] IMP        : Material point number
    !    @param[in] Icon       : Element conectivities
    !    @param[in] B          : B matrix
    !    @param[in] B_bar      : B bar matrix from Hughes (1987)
    !    @param[in] B_bar_flag : Logical that controls whether B_bar should be usedd
    !    @param[in] disp       : Displacements of global dof
    !    @param[in] NDofEx     : Reduced dof
    ! 
    !>   @param[out] Eps       : Strain
    !
    !*************************************************************************************
    subroutine Get_Strain(IEl, IPoint, ICon, B, B_bar_flag, B_bar, Disp, NDofEx, Eps )

        use ModCounters
        use ModMeshInfo
        use ModMPMData
     
        implicit none
     
        integer(INTEGER_TYPE)                                        ,intent(in)  :: IEl, IPoint
        integer(INTEGER_TYPE) ,dimension(Counters%NodTot+1)          ,intent(in)  :: NDofEx
        integer(INTEGER_TYPE) ,dimension(ELEMENTNODES, Counters%Nel) ,intent(in)  :: ICon
        real(REAL_TYPE)       ,dimension(NVECTOR, ELEMENTNODES)      ,intent(in)  :: B
        real(REAL_TYPE)       ,dimension(Counters%N)                 ,intent(in)  :: Disp
        real(REAL_TYPE)       ,dimension(NTENSOR)                    ,intent(out) :: Eps
     
        ! local variables
        integer(INTEGER_TYPE) :: J, NN, ParticleIndex
        real(REAL_TYPE)       :: Position
        real(REAL_TYPE) ,dimension(NVECTOR)      :: U
        real(REAL_TYPE) ,dimension(ELEMENTNODES) :: ShapeValues
        real(REAL_TYPE), dimension(ELEMENTNODES) :: B_0
     
        Eps = 0.0
     
        select case(NDIM)
     
         case(2)
            ! B_0 = N/r Eqn 4.5.25 Hughes (1987)- radial strain only required for 2D axi-symmetric
            B_0 = 0.0

           if ( ISAXISYMMETRIC ) then ! 2D-Axisymmetrix flag
              if ( IsParticleIntegration(IEl) ) then ! MP-integration
                 ParticleIndex = GetParticleIndex(IPoint, IEl) ! MP global ID
                 Position = GlobPosArray(ParticleIndex, 1) ! index 1 is radial direction
                 ShapeValues(:) = ShapeValuesArray(ParticleIndex, :)
              else ! GP-integration
                 Position = GPGlobalPositionElement(1, IPoint, IEl) ! index 1 is radial direction
                 ShapeValues(:) = GPShapeFunction(IPoint, :)
              end if
              
              ! Set values - helps modify B matrix 
              B_0 = ShapeValues/ Position
                 
              ! TODO: Need to evaluate /B_0} -  evaluate N_{a} using eqn. 4.5.17 Hughes (1987)
           end if
            
           if (B_bar_flag) then
            ! modify the B-matrix 
           end if
           
           !IMPORTANT - /B_{a} has dimensions 4x3 instead of the standard B_{a} having dimensions 3x3
           do J = 1, ELEMENTNODES
              NN = Icon(J, Iel)
              U(1) = Disp(NDofEx(NN)+1)
              U(2) = Disp(NDofEx(NN)+2)
              Eps(1) = Eps(1) + B(1,J) * U(1)                  ! Eps_XX = dUx/dX
              Eps(2) = Eps(2) + B(2,J) * U(2)                  ! Eps_YY = dUy/dY
     
              if ( ISAXISYMMETRIC ) then
                 Eps(3) = Eps(3) + U(1) * B_0(J)            ! Eps_tt = u_r*N/r
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
           !TODO: Add condtion that if B-Bar then update the B matrix
        end select
     
     end subroutine Get_Strain 
end module ModB_Bar