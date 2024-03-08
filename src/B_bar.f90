module ModB_Bar
    
    use ModElementEvaluation, only : BMatrix
    use ModGlobalConstants  , only : INTEGER_TYPE, REAL_TYPE 
    use ModFeedback         , only : WriteInLogFile
    use ModString           , only : String
    use ModElementEvaluation, only : ShapeFunctionData
    implicit none

contains
    !**********************************************************************
   !
   !  SUBROUTINE : DetJacob_shape
   !             Determination of the Jacobian matrix and the determinant of
   !             the Jacobian matrix given the values of the shape functions
   !             evaluated at a point
   !
   !  I  LocPos : Local coordinates of the considered point inside an element
   !  I  NEl : Total number of elements
   !  I  NodTot : Total number of nodes
   !  I  IDim : Number of dimensions
   !  I  IElement : ID of the element
   !  I  ICon : Element connectivities ICon(I, J): global node number of local node I in element J
   !  I  Co : Global nodal coordinates Co(I, J): j-coordinate of node I
   !
   !  O  RJac : Jacobian matrix
   !  O  InvRJac : Inverse of the Jacobian matrix
   !  O  DetJac : Determinate of the Jacobian matrix
   !
   !**********************************************************************
   subroutine DetJacob_given_shape(dHS, NEl, NodTot, IDim, IElement, ICon, Co, RJac, InvRJac, DetJac)

      implicit none

      integer(INTEGER_TYPE)                  ,intent(in)    :: NEl, NodTot, IDim, IElement
      integer(INTEGER_TYPE), dimension(:, :) ,intent(in)    :: ICon
      real(REAL_TYPE)      , dimension(:, :) ,intent(in)    :: dHS 
      real(REAL_TYPE)      , dimension(:, :) ,intent(in)    :: Co
      real(REAL_TYPE)      , dimension(:, :) ,intent(inout) :: RJac
      real(REAL_TYPE)      , dimension(:, :) ,intent(inout) :: InvRJac
      real(REAL_TYPE)      , intent(inout)                  :: DetJac

      ! local variables
      integer(INTEGER_TYPE)                      :: I, J, INode, NodeID, NNodes
      real(REAL_TYPE) :: Det1

      NNodes = size(ICon,1)

      ! Determine the Jacobian matrix RJac
      RJac = 0.0
      do INode = 1, NNodes ! loop nodes of each element
         NodeID = ICon(INode, IElement)
         do I = 1, IDim
            do J = 1, IDim
               RJac(I, J) = RJac(I, J) + dHS(INode, I) * Co(NodeID, J)
            end do
         end do
      end do

      ! Calculate inverse and determinant of Jacobian matrix
      call RJacInv(IDim, RJac, InvRJac, DetJac, Det1)

      if ( (DetJac < 0.0) .or. (Det1 < 0.0) ) then
         call WriteInLogFile('Negative determinant '// trim(String(DetJac)) // trim(String(Det1)) // ' element ' // trim(String(IElement)))
      end if

   end subroutine DetJacob_given_shape
   

   !--------------------------------------------------------------------------------------------------
   !--------------------------------------------------------------------------------------------------
   !--------------------------------------------------------------------------------------------------

   subroutine B_bar_matrix(dHS, IElTyp, NEl, NodTot, IDim, IElement, ICon, NodeCoord, B_bar)
  
      integer(INTEGER_TYPE)                          ,intent(in)  :: IElTyp, NEl, NodTot, IDim, IElement
      integer(INTEGER_TYPE) ,dimension(IElTyp, NEl)  ,intent(in)  :: ICon
      real(REAL_TYPE)       ,dimension(:, :)         ,intent(in)  :: NodeCoord
      real(REAL_TYPE)       ,dimension(:, :)         ,intent(in)  :: dHS ! Derivatives of shape functions
      real(REAL_TYPE)       ,dimension(IDim, IElTyp) ,intent(out) :: B_bar

      ! Local variables
      integer(INTEGER_TYPE)                  :: I, J, K
      real(REAL_TYPE)                        :: DetJac
      real(REAL_TYPE), dimension(IDim, IDim) :: RJac, RJacInv ! Jacobian matrix, inverse of Jacobian matrix

      ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix RJacInv
      call DetJacob_given_shape(dHS, NEl, NodTot, IDim, IElement, ICon, NodeCoord, RJac, RJacInv, DetJac)

      ! Assemble B matrix (cartesian derivatives)
      B_bar = 0.0
      do J = 1, IDim
         do I = 1, IElTyp
            do K = 1, IDim
               B_bar(K, I) = B_bar(K,I) + RJacInv(K, J) * dHS(I, J)
            end do
         end do
      end do

   end subroutine B_bar_matrix

    
   ! Evaluate the shape and deriv shape functions at the local element center
   subroutine eval_local_elem_center(model_dimension, NNodes, HS_center, dHS_center)
   
      integer(INTEGER_TYPE)                         :: model_dimension, NNodes
      real(REAL_TYPE), dimension(:)   , intent(out) :: HS_center
      real(REAL_TYPE), dimension(:, :), intent(out) :: dHS_center

      ! Local variables
      real(REAL_TYPE), dimension(model_dimension) :: zero_vector

      zero_vector(:) = 0.0
      call ShapeFunctionData(zero_vector, NNodes, HS_center, dHS_center)

   end subroutine

end module ModB_Bar