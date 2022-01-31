      module ModSolver
      !**********************************************************************
      !
      !    Function:  This module contains all routines related to solvers.
      !
      !               At present the Intel MKL Pardiso solver is considered:
      !               - direct sparse solver
      !               - shared-memory multiprocessing parallelised
      !               - row-major upper triangular storage format
      !
      !               The DSS interface to the Pardiso solver is used.
      !               The routines XXX are called outside this module for
      !               initialising, for updating the equation system, for solving
      !               it, for destroying the data structure.
      !               The code can be easily extended to incorporate other solvers.
      !               No code has as yet been added for this purpose.
      !
      ! Implemented in the frame of the Deltares project
      ! 'Extension of the Deltares MPM for quasi-static 2-phase analyses'.
      ! Based on source code written by Lars Beuth.
      !
      !**********************************************************************

      use mkl_dss
      use ModCounters
      use ModMeshInfo
      use ModMeshAdjacencies
      use ModGlobalConstants
      
      implicit none

      type TypeEquations
        
        type(MKL_DSS_HANDLE) :: HandleDSS
        integer(INTEGER_TYPE) :: DSSStatus = MKL_DSS_SUCCESS
        integer(INTEGER_TYPE) :: NNonZeros = 0
        integer(INTEGER_TYPE) :: NRHS = 1
        logical :: IsSymmetric = .true.
        integer(INTEGER_TYPE) :: NReducedDoF = 0

        real(REAL_TYPE), dimension(:), allocatable :: RHS
        real(REAL_TYPE), dimension(:), allocatable :: Matrix
        real(REAL_TYPE), dimension(:), allocatable :: Solution

        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDoFMapping
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDoFMappingInverse
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDoFMappingInverseNodes
        integer(INTEGER_TYPE), dimension(:), allocatable :: ReducedDoFMappingMode 

        integer(INTEGER_TYPE), dimension(:), allocatable :: ColumnsDSS
        integer(INTEGER_TYPE), dimension(:), allocatable :: RowIndexDSS

      end type TypeEquations

      type(TypeEquations), private, save :: Equations

      contains ! Routines of this module

      subroutine ResetGlobalStiffnessMatrix()

      implicit none

      Equations%Matrix = 0.0

      end subroutine ResetGlobalStiffnessMatrix
      
      subroutine ResetGlobalStiffnessRHS()

      implicit none

      Equations%RHS = 0.0

      end subroutine ResetGlobalStiffnessRHS
      
      logical function IsSymmetricStiffnessMatrix()

      implicit none

      IsSymmetricStiffnessMatrix = Equations%IsSymmetric

      end function IsSymmetricStiffnessMatrix

      subroutine InitialiseReducedSolution()
      !**********************************************************************
      !
      !  Function:  Initialises data related to the system of equations to be solved.
      !
      !**********************************************************************

      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IError

      Equations%IsSymmetric = (.not.CalParams%ApplyObjectiveStress)

      call CreateProfileDSS()

      allocate(Equations%RHS(Equations%NReducedDoF), stat = IError)
      Equations%RHS = 0.0

      allocate(Equations%Solution(Equations%NReducedDoF), stat = IError)
      Equations%Solution = 0.0

      end subroutine InitialiseReducedSolution

      subroutine CreateProfileDSS()
      !**********************************************************************
      !
      !  Function:  Creates the profile of non-zeros for the stiffness matrix
      !             as required by the Intel MKL DSS:
      !             row-major, in case of symmetric matrix for the upper
      !             triangular half of the matrix.
      !             Computes NNonZeros, initialises ColumnsDSS and RowIndexDSS,
      !             allocates Matrix, filling it with 0.0.
      !
      !**********************************************************************

      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IError, I, J, NodeID, NNodes, IDof, JDof
      integer(INTEGER_TYPE) :: NFreeDoFI, NFreeDoFJ, INodeID, INode, IReducedDoF, JReducedDoF
      integer(INTEGER_TYPE) :: Opt, Counter
      integer(INTEGER_TYPE), dimension(:), allocatable :: PermutationVector
      logical :: FirstNonZeroDoF

      allocate(Equations%ReducedDoFMapping(Counters%NodTot * NDOFLQS), stat = IError)  
      Equations%ReducedDoFMapping = 0
      allocate(Equations%ReducedDoFMappingMode(Counters%NodTot * NDOFLQS), stat = IError)
      Equations%ReducedDoFMappingMode = 0

      Equations%NReducedDoF = 0
      do NodeID = 1, Counters%NodTot
        if (IsElementCornerNode(NodeID)) then
          do I = 1, NVECTOR
            IDoF = ReducedDoFQuasiStatic(NodeID) + I
            if (PBoundaryQuasiStatic(IDoF)/=0) then
              Equations%NReducedDoF = Equations%NReducedDoF + 1
              Equations%ReducedDoFMapping(IDoF) = Equations%NReducedDoF
              Equations%ReducedDoFMappingMode(IDoF) = I
            end if
          end do
          IDoF = ReducedDoFQuasiStatic(NodeID) + I
          if (PBoundaryQuasiStatic(IDoF)/=0) then
            Equations%NReducedDoF = Equations%NReducedDoF + 1
            Equations%ReducedDoFMapping(IDoF) = Equations%NReducedDoF
            Equations%ReducedDoFMappingMode(IDoF) = I
          end if
        end if
      end do

      allocate(Equations%ReducedDoFMappingInverse(Equations%NReducedDoF), stat = IError)
      Equations%ReducedDoFMappingInverse = 0
      allocate(Equations%ReducedDoFMappingInverseNodes(Equations%NReducedDoF), stat = IError)
      Equations%ReducedDoFMappingInverseNodes = 0

      Counter = 0
      do NodeID = 1, Counters%NodTot
        if (IsElementCornerNode(NodeID)) then
          do I = 1, NDOFLQS 
            IDoF = ReducedDoFQuasiStatic(NodeID) + I
            if (PBoundaryQuasiStatic(IDoF)/=0) then
              Counter = Counter + 1
              Equations%ReducedDoFMappingInverse(Counter) = IDoF 
              Equations%ReducedDoFMappingInverseNodes(Counter) = NodeID
            end if
          end do
        end if
      end do

      Equations%NNonZeros = 0
      do NodeID = 1, Counters%NodTot
        if (IsElementCornerNode(NodeID)) then
          NFreeDoFI = GetNDoFOfNode(NodeID)
          if (NFreeDoFI/=0) then
            NNodes = GetNNodesOfNode(NodeID)
            NFreeDoFJ = 0
            do I = 1, NNodes
              NFreeDoFJ = NFreeDoFJ + GetNDoFOfNode(GetNodeIOfNode(NodeID, I))
            end do
            Equations%NNonZeros = Equations%NNonZeros + NFreeDoFJ * NFreeDoFI
          end if
        end if
      end do

      if (Equations%IsSymmetric) then
        Equations%NNonZeros =  (Equations%NNonZeros - Equations%NReducedDoF) / 2 + Equations%NReducedDoF
      end if

      allocate(Equations%ColumnsDSS(Equations%NNonZeros), stat = IError)
      Equations%ColumnsDSS = 0
      allocate(Equations%RowIndexDSS(Equations%NReducedDoF + 1), stat = IError)
      Equations%RowIndexDSS = 0

      Counter = 0
      do NodeID = 1, Counters%NodTot
        if (IsElementCornerNode(NodeID)) then
          NNodes = GetNNodesOfNode(NodeID)
          do I = 1, NDOFLQS
            IDoF = ReducedDoFQuasiStatic(NodeID) + I
            if (PBoundaryQuasiStatic(IDoF)/=0) then
              IReducedDoF = Equations%ReducedDoFMapping(IDoF)
              FirstNonZeroDoF = .false.
              do INode = 1, NNodes
                INodeID = GetNodeIOfNode(NodeID, INode)
                do J = 1, NDOFLQS
                  JDoF = ReducedDoFQuasiStatic(INodeID) + J
                  if (PBoundaryQuasiStatic(JDoF)/=0) then
                    JReducedDoF = Equations%ReducedDoFMapping(JDoF)
                    if ((.not.Equations%IsSymmetric).or.(JReducedDoF>=IReducedDoF)) then

                      Counter = Counter + 1
                      Equations%ColumnsDSS(Counter) = JReducedDoF

                      if (.not.FirstNonZeroDoF) then
                        FirstNonZeroDoF = .true.
                        Equations%RowIndexDSS(IReducedDoF) = Counter
                      end if
                    end if
                  end if
                end do
              end do
            end if
          end do
        end if
      end do

      Equations%RowIndexDSS(Equations%NReducedDoF + 1) = Equations%NNonZeros + 1

      allocate(Equations%Matrix(Equations%NNonZeros), stat = IError)
      Equations%Matrix = 0.0

      Opt = MKL_DSS_DEFAULTS
      Equations%DSSStatus = DSS_CREATE(Equations%HandleDSS, Opt)

      if (Equations%DSSStatus==MKL_DSS_SUCCESS) then
        if (Equations%IsSymmetric) then
          Opt = MKL_DSS_SYMMETRIC
        else
        !Opt = MKL_DSS_SYMMETRIC_STRUCTURE ! In case of large deformations?
          Opt = MKL_DSS_NON_SYMMETRIC
        end if
        Equations%DSSStatus = DSS_DEFINE_STRUCTURE(Equations%HandleDSS, Opt, Equations%RowIndexDSS, Equations%NReducedDoF, Equations%NReducedDoF, Equations%ColumnsDSS, Equations%NNonZeros)

        if (Equations%DSSStatus==MKL_DSS_SUCCESS) then
          Opt = MKL_DSS_AUTO_ORDER
          allocate(PermutationVector(Equations%NReducedDoF), stat = IError)
          Equations%DSSStatus = DSS_REORDER(Equations%HandleDSS, Opt, PermutationVector)
          if (Equations%DSSStatus/=MKL_DSS_SUCCESS) then
            write(OUTunit, *) "DSS_REORDER returned error code  ", Equations%DSSStatus, "."
          end if
          deallocate(PermutationVector, stat = IError)
        else
          write(OUTunit, *) "DSS_DEFINE_STRUCTURE returned error code ", Equations%DSSStatus, "."
        end if
      else
        write(OUTunit, *) "DSS_CREATE returned error code ", Equations%DSSStatus, "."
      end if

      end subroutine CreateProfileDSS

      integer(INTEGER_TYPE) function GetNDoFOfNode(NodeID)

      implicit none

      integer(INTEGER_TYPE), intent(in) :: NodeID

      GetNDoFOfNode = PBoundaryQuasiStatic(ReducedDoFQuasiStatic(NodeID) + 1) + PBoundaryQuasiStatic(ReducedDoFQuasiStatic(NodeID) + 2) + PBoundaryQuasiStatic(ReducedDoFQuasiStatic(NodeID) + 3) + PBoundaryQuasiStatic(ReducedDoFQuasiStatic(NodeID) + 4)

      end function GetNDoFOfNode

      integer(INTEGER_TYPE) function GetMatrixIndexOfIDoF(IDoFRow, IDoFColumn)
      !**********************************************************************
      !
      !  Function:  Returns the index of (IDoFRow, IDoFColumn) in Equations%Matrix.
      !             Returns -1 if the matrix is symmetric and IDoFRow is smaller
      !             than IDoFColumn or (IDoFRow, IDoFColumn) marks a zero entry.
      !
      !**********************************************************************

      implicit none

      integer(INTEGER_TYPE), intent(in) :: IDoFRow, IDoFColumn
      ! Local variables
      integer(INTEGER_TYPE) :: I, N, FirstNonZeroInRow, Column

      GetMatrixIndexOfIDoF = -1

      if ((IDoFColumn>=IDoFRow).or.(.not.Equations%IsSymmetric)) then
        FirstNonZeroInRow = Equations%RowIndexDSS(IDoFRow)
        N = Equations%RowIndexDSS(IDoFRow + 1) - FirstNonZeroInRow

        do I = 0, N - 1 ! Loop over non-zeros of row
          Column = Equations%ColumnsDSS(FirstNonZeroInRow + I)
          if (Column==IDoFColumn) then
            GetMatrixIndexOfIDoF = FirstNonZeroInRow + I
            EXIT
          end if
        end do
      end if

      end function GetMatrixIndexOfIDoF

      subroutine AddElementStiffnessMatrix(ElementID, ElementStiffnessMatrix) 
      !**********************************************************************
      !
      !  Function:  Assembles the element stiffness matrix for ElementID and
      !             adds it to Equations%Matrix.
      !
      !**********************************************************************

      implicit none

      integer(INTEGER_TYPE), intent(in) :: ElementID
      real(REAL_TYPE),dimension(ELEMENTNODES * NDOFLQS, ELEMENTNODES * NDOFLQS) ::ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: NodeLocalRow, NodeIDRow, IDimRow, IDoFRow, MatrixIndex
      integer(INTEGER_TYPE) :: NodeLocalColumn, NodeIDColumn, IDimColumn, IDoFColumn
      integer(INTEGER_TYPE) :: IReducedDoFRow, IReducedDoFColumn

      integer(INTEGER_TYPE) :: i,j,row,dof,node,column,mode
      real(REAL_TYPE) :: value
      integer(INTEGER_TYPE), parameter :: SOLVERFILE = 11
      
      ! Add to matrix
      do NodeLocalRow = 1, ELEMENTNODES
        NodeIDRow = ElementConnectivities(NodeLocalRow, ElementID)
        do IDimRow = 1, NDOFLQS
          IDoFRow = ReducedDoFQuasiStatic(NodeIDRow) + IDimRow
          if (PBoundaryQuasiStatic(IDoFRow)/=0) then
            IReducedDoFRow = Equations%ReducedDoFMapping(IDoFRow)
            do NodeLocalColumn = 1, ELEMENTNODES
              NodeIDColumn = ElementConnectivities(NodeLocalColumn, ElementID)
              do IDimColumn = 1, NDOFLQS
                IDoFColumn = ReducedDoFQuasiStatic(NodeIDColumn) + IDimColumn
                if (PBoundaryQuasiStatic(IDoFColumn)/=0) then
                  IReducedDoFColumn = Equations%ReducedDoFMapping(IDoFColumn)
                  MatrixIndex = GetMatrixIndexOfIDoF(IReducedDoFRow, IReducedDoFColumn)
                  if (MatrixIndex/=-1) then
                    Equations%Matrix(MatrixIndex) = Equations%Matrix(MatrixIndex) + ElementStiffnessMatrix((NodeLocalRow - 1) * NDOFLQS + IDimRow, (NodeLocalColumn - 1) * NDOFLQS + IDimColumn)
                  end if
                else if (IDimColumn<=NDIM) then
                  Equations%RHS(IReducedDoFRow) = Equations%RHS(IReducedDoFRow) 
                else
                  Equations%RHS(IReducedDoFRow) = Equations%RHS(IReducedDoFRow) - ElementStiffnessMatrix((NodeLocalRow - 1) * NDOFLQS + IDimRow, (NodeLocalColumn - 1) * NDOFLQS + IDimColumn) * SubIncrementalPressure(NodeIDColumn)
                end if
              end do
            end do
          end if
        end do
      end do

      if (CalParams%OutputDebugData) then
        open(SOLVERFILE,file='d:\tmp\solvermatrix.txt')
        if (ElementID==1) then
          write(SOLVERFILE,*) 'ElementID',ElementID
          write(SOLVERFILE,*)
          write(SOLVERFILE,*) 'ReducedDof'
          do NodeLocalRow = 1, ELEMENTNODES
            NodeIDRow = ElementConnectivities(NodeLocalRow, ElementID)
            do IDimRow = 1, NDOFLQS
              IDoFRow = ReducedDoFQuasiStatic(NodeIDRow) + IDimRow
              write(SOLVERFILE,*) NodeLocalRow, NodeIDRow, IDimRow, IDoFRow
            end do
          end do
        end if
        write(SOLVERFILE,*) 'Equations%ColumnsDSS', Equations%NNonZeros
        do i = 1, Equations%NNonZeros
          write(SOLVERFILE,*) i, Equations%ColumnsDSS(i)
        end do
        write(SOLVERFILE,*) 'Equations%RowIndexDSS', Equations%NReducedDoF+1
        do i = 1, Equations%NReducedDoF+1
          write(SOLVERFILE,*) i, Equations%RowIndexDSS(i)
        end do
        do i = 1, Equations%NReducedDoF
          row = i;
          dof = Equations%ReducedDoFMappingInverse(row)
          mode = Equations%ReducedDoFMappingMode(dof)
          node = Equations%ReducedDoFMappingInverseNodes(row)
          write(SOLVERFILE,*) 'row', row, 'dof', dof, 'mode', mode, 'node', node
          do j = Equations%RowIndexDSS(i), Equations%RowIndexDSS(i+1) - 1
            column = Equations%ColumnsDSS(j)
            dof = Equations%ReducedDoFMappingInverse(column)
            mode = Equations%ReducedDoFMappingMode(dof)
            node = Equations%ReducedDoFMappingInverseNodes(column)
            value = Equations%Matrix(j)
            write(SOLVERFILE,102) ' column ', column, ' dof ', dof, ' mode ', mode, ' node ', node, ' value ', value
102     format(A, I, A, I, A, I, A, I, A, E)
          end do
          write(SOLVERFILE,*)
        end do
        close(SOLVERFILE)
      end if

      end subroutine AddElementStiffnessMatrix

      subroutine AdjustGlobalStiffnessMatrixForInactiveElements()

      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: INode, IDim, IDoF, IReducedDoF, MatrixIndex

      if (.not.IsMPMComputation()) RETURN

      do INode = 1, Counters%NodTot
        if (IsElementCornerNode(INode)) then
          if (.not.ActiveNode(INode)) then
            do IDim = 1, NDOFLQS
              IDoF = ReducedDoFQuasiStatic(INode) + IDim
              if (PBoundaryQuasiStatic(IDoF)/=0) then
                IReducedDoF = Equations%ReducedDoFMapping(IDoF)
                MatrixIndex = GetMatrixIndexOfIDoF(IReducedDoF, IReducedDoF)
                Equations%Matrix(MatrixIndex) = 1.0
              end if
            end do
          end if
        end if
      end do

      end subroutine AdjustGlobalStiffnessMatrixForInactiveElements

      subroutine SolveEquations(RateOfMomentum, RateOffFlow, SubIncrementalDisplacement, SubIncrementalPressure)
      !**********************************************************************
      !
      !  Function:  Solves the system of equations provided by Equations.
      !
      !**********************************************************************

      implicit none

      real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: RateOfMomentum
      real(REAL_TYPE), dimension(Counters%N), intent(inout) :: SubIncrementalDisplacement
      real(REAL_TYPE), dimension(Counters%NodTot), intent(in) :: RateOffFlow
      real(REAL_TYPE), dimension(Counters%NodTot), intent(inout) :: SubIncrementalPressure
      
      ! Local variables
      integer(INTEGER_TYPE) :: Opt, IReducedDoF, IDoF, INode, I, IMode, FirstNonZeroInRow, N, Row, Column, J
      real(REAL_TYPE) :: value
      integer(INTEGER_TYPE), parameter :: SOLVERFILE2 = 12
      real(REAL_TYPE), dimension(:), allocatable :: ErrorVector
      integer(INTEGER_TYPE) :: IError
      allocate(ErrorVector(Equations%NReducedDoF), stat = IError)  
      ErrorVector = 0
      
      if (CalParams%OutputDebugData) then
        open(SOLVERFILE2,file='d:\tmp\struct.txt')
        write(SOLVERFILE2,*) 'PBoundary', Counters%N
        do I = 1, Counters%N
          write(SOLVERFILE2,*) I, PBoundary(I)
        end do
        write(SOLVERFILE2,*) 'ReducedDoF', Counters%NodTot+1
        do I = 1, Counters%NodTot+1
          write(SOLVERFILE2,*) I, ReducedDoF(I)
        end do
        write(SOLVERFILE2,*) 'ReducedDofQuasiStatic', Counters%NodTot+1
        do I = 1, Counters%NodTot+1
          write(SOLVERFILE2,*) I, ReducedDofQuasiStatic(I)
        end do
        write(SOLVERFILE2,*) 'PBoundaryQuasiStatic', Counters%NodTot*4
        do I = 1, Counters%NodTot*NDOFLQS
          write(SOLVERFILE2,*) I, PBoundaryQuasiStatic(I)
        end do
        write(SOLVERFILE2,*) 'PBoundaryWater', Counters%N
        do I = 1, Counters%N
          write(SOLVERFILE2,*) I, PBoundaryWater(I)
        end do
        write(SOLVERFILE2,*) 'ReducedDoFMapping ReducedDoFMappingMode', Counters%NodTot * NDOFLQS
        do I = 1, Counters%NodTot * NDOFLQS
          write(SOLVERFILE2,*) I, Equations%ReducedDoFMapping(I), Equations%ReducedDoFMappingMode(I)
        end do
        write(SOLVERFILE2,*) 'ReducedDoFInverse ReducedDoFMappingInverseNodes', Equations%NReducedDoF
        do I = 1, Equations%NReducedDoF
          write(SOLVERFILE2,*) I, Equations%ReducedDoFMappingInverse(I), Equations%ReducedDoFMappingInverseNodes(I)
        end do
        close(SOLVERFILE2)
      end if 
      
      do IReducedDoF = 1, Equations%NReducedDoF
        IDoF = Equations%ReducedDoFMappingInverse(IReducedDoF)
        IMode = Equations%ReducedDoFMappingMode(IDoF)
        INode = Equations%ReducedDoFMappingInverseNodes(IReducedDoF)
        IDof = ReducedDof(INode)
        if (IMode<=NDIM) then
           Equations%RHS(IReducedDoF) = Equations%RHS(IReducedDoF) + RateOfMomentum(IDof+IMode, 1)
        else
           Equations%RHS(IReducedDoF) = Equations%RHS(IReducedDoF) + RateOffFlow(INode)
        end if
      end do

      Opt = MKL_DSS_INDEFINITE !MKL_DSS_DEFAULTS !MKL_DSS_POSITIVE_DEFINITE
      Equations%DSSStatus = DSS_FACTOR_REAL(Equations%HandleDSS, Opt, Equations%Matrix)

      Opt = MKL_DSS_DEFAULTS
      Equations%DSSStatus = DSS_SOLVE_REAL_D(Equations%HandleDSS, Opt, Equations%RHS, 1, Equations%Solution)

      do IReducedDoF = 1, Equations%NReducedDoF
        IDoF = Equations%ReducedDoFMappingInverse(IReducedDoF)
        IMode = Equations%ReducedDoFMappingMode(IDoF)
        INode = Equations%ReducedDoFMappingInverseNodes(IReducedDoF)
        IDof = ReducedDof(INode)
        if (IMode<=NDIM) then
           SubIncrementalDisplacement(IDoF+IMode) = Equations%Solution(IReducedDoF)
        else
           SubIncrementalPressure(INode) = Equations%Solution(IReducedDoF)
        end if
      end do
      
      if (CalParams%OutputDebugData) then
        open(SOLVERFILE2,file='d:\tmp\solvervector.txt')
        write(SOLVERFILE2,*) 'rate of momentum'
        do I = 1, Counters%N
          write(SOLVERFILE2,*) I, RateOfMomentum(I, 1)
        end do
        write(SOLVERFILE2,*) 'steady flow'
        do I = 1, Counters%NodTot
          write(SOLVERFILE2,*) I, RateOffFlow(I)
        end do
        write(SOLVERFILE2,*) 'RHS'
        do IReducedDoF = 1, Equations%NReducedDoF
          IDoF = Equations%ReducedDoFMappingInverse(IReducedDoF)
          INode = Equations%ReducedDoFMappingInverseNodes(IReducedDoF)
          Value = Equations%RHS(IReducedDoF)
          write(SOLVERFILE2,102) ' IReducedDof ', IReducedDof, ' IDof ', IDof, ' INode ', INode, ' Value ', Value
        end do
102     format(A, I, A, I, A, I, A, E)
        write(SOLVERFILE2,*) 'matrix', Equations%NReducedDoF
        do Row = 1, Equations%NReducedDoF
          FirstNonZeroInRow = Equations%RowIndexDSS(Row)
          N = Equations%RowIndexDSS(Row + 1) - FirstNonZeroInRow
          write(SOLVERFILE2,*) Row
          do J = 0, N - 1 
            write(SOLVERFILE2,*) J, Equations%Matrix(FirstNonZeroInRow + J)
          end do
        end do
        do Row = 1, Equations%NReducedDoF
          FirstNonZeroInRow = Equations%RowIndexDSS(Row)
          N = Equations%RowIndexDSS(Row + 1) - FirstNonZeroInRow
          ErrorVector(Row) = Equations%RHS(Row)
          do J = 0, N - 1 
            Column = Equations%ColumnsDSS(FirstNonZeroInRow + J)
            ErrorVector(Row) = ErrorVector(Row) - Equations%Matrix(FirstNonZeroInRow + J) * Equations%Solution(Column)
          end do
        end do
        write(SOLVERFILE2,*) 'solution', Equations%NReducedDoF
        do I = 1, Equations%NReducedDoF
          write(SOLVERFILE2,*) I, Equations%Solution(I)
        end do
        write(SOLVERFILE2,*) 'RHSerror', Equations%NReducedDoF
        do I = 1, Equations%NReducedDoF
          write(SOLVERFILE2,*) I, ErrorVector(I)
        end do
        write(SOLVERFILE2,*) 'displacement', Counters%N
        do I = 1, Counters%N
          write(SOLVERFILE2,*) I, SubIncrementalDisplacement(I)
        end do
        write(SOLVERFILE2,*) 'pressure', Counters%NodTot
        do I = 1, Counters%NodTot
          write(SOLVERFILE2,*) I, SubIncrementalPressure(I)
        end do
        close(SOLVERFILE2)
      end if
      
      end subroutine SolveEquations

      subroutine DestroyEquations()
      !**********************************************************************
      !
      !  Function:  Deallocates data of Equations.
      !
      !**********************************************************************

      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IError, Opt

      Equations%NNonZeros = 0
      Equations%NRHS = 0
      Equations%IsSymmetric = .true.
      Equations%NReducedDoF = 0

      if (allocated(Equations%Matrix)) then
        Opt = MKL_DSS_DEFAULTS
        Equations%DSSStatus = DSS_DELETE(Equations%HandleDSS, Opt)
      end if

      if (allocated(Equations%RHS)) then
        deallocate(Equations%RHS, stat = IError)
      end if

      if (allocated(Equations%Matrix)) then
        deallocate(Equations%Matrix, stat = IError)
      end if

      if (allocated(Equations%ColumnsDSS)) then
        deallocate(Equations%ColumnsDSS, stat = IError)
      end if

      if (allocated(Equations%RowIndexDSS)) then
        deallocate(Equations%RowIndexDSS, stat = IError)
      end if

      if (allocated(Equations%ReducedDoFMapping)) then
        deallocate(Equations%ReducedDoFMapping, stat = IError)
      end if

      if (allocated(Equations%ReducedDoFMappingMode)) then
        deallocate(Equations%ReducedDoFMappingMode, stat = IError)
      end if

      if (allocated(Equations%ReducedDoFMappingInverse)) then
        deallocate(Equations%ReducedDoFMappingInverse, stat = IError)
      end if

      if (allocated(Equations%ReducedDoFMappingInverseNodes)) then
        deallocate(Equations%ReducedDoFMappingInverseNodes, stat = IError)
      end if

      end subroutine DestroyEquations

      end module ModSolver