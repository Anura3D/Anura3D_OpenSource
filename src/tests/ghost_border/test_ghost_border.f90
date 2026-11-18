! Unit harness for ModGhostBorder (src/GhostBorderDetection.FOR): tests T1-T7 of
! notes/FN-T-002.md (implementation note 6, predictions P1-P7) on a structured
! 6 x 6 TRI3 mesh, plus the Q-010 fill check and small edge cases.
! Build and run: tests/ghost_border/run_tests.bat (ifx, same real kind as the kernel).
! Exit code 0 if every check passes, 1 otherwise.

      module ModTestMesh
      use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
      implicit none

        integer(INTEGER_TYPE), parameter :: NC = 6, NR = 6            ! squares per row and column
        integer(INTEGER_TYPE), parameter :: NEL = 2 * NC * NR, NNOD = (NC + 1) * (NR + 1)
        real(REAL_TYPE), parameter :: L = 0.1d0, RHO = 2.0d0, AREA = 0.5d0 * L * L
        integer(INTEGER_TYPE) :: Conn(3, NEL), Adj(NEL, 3)
        real(REAL_TYPE) :: XY(2, NNOD)
        integer(INTEGER_TYPE) :: NFail = 0, NPass = 0

      contains

        integer(INTEGER_TYPE) function NodeId(i, j)
        ! node at (i L, j L), i = 0..NC, j = 0..NR
          integer(INTEGER_TYPE), intent(in) :: i, j
          NodeId = i * (NR + 1) + j + 1
        end function NodeId

        integer(INTEGER_TYPE) function ElLower(c, r)
        ! lower triangle (BL, BR, TR) of square (c, r) = [c-1, c] x [r-1, r] (FN-T-002 note 6)
          integer(INTEGER_TYPE), intent(in) :: c, r
          ElLower = 2 * ((c - 1) * NR + (r - 1)) + 1
        end function ElLower

        integer(INTEGER_TYPE) function ElUpper(c, r)
        ! upper triangle (BL, TR, TL) of square (c, r)
          integer(INTEGER_TYPE), intent(in) :: c, r
          ElUpper = 2 * ((c - 1) * NR + (r - 1)) + 2
        end function ElUpper

        subroutine BuildMesh()
        ! connectivity and side adjacency with the side convention of BuildGhostFaces (1-2, 3-1, 2-3)
          integer(INTEGER_TYPE) :: c, r, i, j, IEl, JEl, si, sj, a1, a2, b1, b2
          integer(INTEGER_TYPE), dimension(2, 3) :: SV
          SV = reshape((/ 1, 2, 3, 1, 2, 3 /), (/ 2, 3 /))
          do i = 0, NC
            do j = 0, NR
              XY(1, NodeId(i, j)) = i * L
              XY(2, NodeId(i, j)) = j * L
            end do
          end do
          do c = 1, NC
            do r = 1, NR
              Conn(:, ElLower(c, r)) = (/ NodeId(c-1, r-1), NodeId(c, r-1), NodeId(c, r) /)
              Conn(:, ElUpper(c, r)) = (/ NodeId(c-1, r-1), NodeId(c, r), NodeId(c-1, r) /)
            end do
          end do
          Adj = 0
          do IEl = 1, NEL
            do si = 1, 3
              a1 = Conn(SV(1, si), IEl); a2 = Conn(SV(2, si), IEl)
              do JEl = 1, NEL
                if (JEl == IEl) CYCLE
                do sj = 1, 3
                  b1 = Conn(SV(1, sj), JEl); b2 = Conn(SV(2, sj), JEl)
                  if ((a1 == b1 .and. a2 == b2) .or. (a1 == b2 .and. a2 == b1)) Adj(IEl, si) = JEl
                end do
              end do
            end do
          end do
        end subroutine BuildMesh

        subroutine CheckInt(Name, Got, Expected)
          character(len=*), intent(in) :: Name
          integer(INTEGER_TYPE), intent(in) :: Got, Expected
          if (Got == Expected) then
            NPass = NPass + 1
            write(*, '(A, A, A, I6, A, I6)') 'PASS  ', Name, ': ', Got, ' expected ', Expected
          else
            NFail = NFail + 1
            write(*, '(A, A, A, I6, A, I6)') 'FAIL  ', Name, ': ', Got, ' expected ', Expected
          end if
        end subroutine CheckInt

        subroutine CheckReal(Name, Got, Expected, Tol)
          character(len=*), intent(in) :: Name
          real(REAL_TYPE), intent(in) :: Got, Expected, Tol
          if (abs(Got - Expected) <= Tol) then
            NPass = NPass + 1
            write(*, '(A, A, A, F12.6, A, F12.6, A, ES9.2)') 'PASS  ', Name, ': ', Got, ' expected ', Expected, ' tol ', Tol
          else
            NFail = NFail + 1
            write(*, '(A, A, A, F12.6, A, F12.6, A, ES9.2)') 'FAIL  ', Name, ': ', Got, ' expected ', Expected, ' tol ', Tol
          end if
        end subroutine CheckReal

        subroutine CheckTrue(Name, Cond)
          character(len=*), intent(in) :: Name
          logical, intent(in) :: Cond
          if (Cond) then
            NPass = NPass + 1
            write(*, '(A, A)') 'PASS  ', Name
          else
            NFail = NFail + 1
            write(*, '(A, A)') 'FAIL  ', Name
          end if
        end subroutine CheckTrue

      end module ModTestMesh


      program test_ghost_border
      use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
      use ModGhostBorder
      use ModTestMesh
      implicit none

        logical :: Candidate(NEL), Cut(NEL), Touch(NNOD), Expected(NEL), InPatch(NNOD)
        real(REAL_TYPE) :: ElDensity(NEL), ElSpace(NEL), SumW(NEL), NodalMass(NNOD), MassRef(NNOD)
        real(REAL_TYPE) :: MinMu, Mu
        integer(INTEGER_TYPE) :: c, r, i, NBorder, NAdded, NFaces, Sparse(6), Node32
        logical :: Ok

        call BuildMesh()
        ElSpace = AREA
        ElDensity = RHO

        ! ---------------- T1: flat surface on y = 3L (P1) ----------------
        call Fill(3, NC)
        call Border(NBorder, NFaces)
        call CheckInt('T1 node-based border elements', NBorder, 12)
        call CheckInt('T1 node-based ghost faces', NFaces, 17)
        Expected = .false.
        do c = 1, NC
          Expected(ElLower(c, 3)) = .true.
          Expected(ElUpper(c, 3)) = .true.
        end do
        call CheckTrue('T1 border set = both triangles of row 3', all(Cut .eqv. Expected))
        call FillCriterion(NBorder, NFaces)
        call CheckInt('T1 fill border elements', NBorder, 0)
        call CheckInt('T1 fill ghost faces', NFaces, 0)

        ! ---------------- T2: corner, c <= 3, r <= 3 (P2) ----------------
        Candidate = .false.
        SumW = 0.0
        do c = 1, 3
          do r = 1, 3
            Candidate(ElLower(c, r)) = .true.; Candidate(ElUpper(c, r)) = .true.
            SumW(ElLower(c, r)) = AREA; SumW(ElUpper(c, r)) = AREA
          end do
        end do
        call Border(NBorder, NFaces)
        call CheckInt('T2 node-based border elements', NBorder, 10)
        call CheckInt('T2 node-based ghost faces', NFaces, 13)
        Expected = .false.
        do r = 1, 3
          Expected(ElLower(3, r)) = .true.; Expected(ElUpper(3, r)) = .true.
        end do
        do c = 1, 2
          Expected(ElLower(c, 3)) = .true.; Expected(ElUpper(c, 3)) = .true.
        end do
        call CheckTrue('T2 border set = both triangles of (3,1) (3,2) (3,3) (1,3) (2,3)', all(Cut .eqv. Expected))
        call FillCriterion(NBorder, NFaces)
        call CheckInt('T2 fill border elements', NBorder, 0)

        ! ---------------- T3a (topological part): loaded nodes on y = 3L of T1 are ghost-patch nodes (P4, S2a) ----
        call Fill(3, NC)
        call Border(NBorder, NFaces, InPatch)
        Ok = .true.
        do i = 0, NC
          Ok = Ok .and. InPatch(NodeId(i, 3))
        end do
        call CheckTrue('T3a all 7 surface nodes y = 3L lie in ghost patches (100 %)', Ok)

        ! ---------------- T4: sparse isolated element, upper (3,2) (P6) ----------------
        call Fill(3, NC)
        SumW(ElUpper(3, 2)) = AREA / 3.0
        call UniformMass()
        ! one MP of mass rho |K| / 3 at the centroid instead of the full element mass rho |K|
        do i = 1, 3
          NodalMass(Conn(i, ElUpper(3, 2))) = NodalMass(Conn(i, ElUpper(3, 2))) - RHO * AREA / 3.0 + RHO * AREA / 9.0
        end do
        call Border(NBorder, NFaces)
        call CheckInt('T4 node-based border elements (unchanged)', NBorder, 12)
        call CheckTrue('T4 sparse element not flagged by P0', .not.Cut(ElUpper(3, 2)))
        call FillCriterion(NBorder, NFaces)
        call CheckInt('T4 fill border elements', NBorder, 1)
        call CheckInt('T4 fill ghost faces', NFaces, 3)
        call GhostNodalMassFlag(NEL, NNOD, Conn, Candidate, ElDensity, ElSpace, NodalMass, 0.25d0, MassRef, Cut, NAdded, MinMu)
        Mu = huge(Mu)
        do i = 1, 3
          Mu = min(Mu, NodalMass(Conn(i, ElUpper(3, 2))) / MassRef(Conn(i, ElUpper(3, 2))))
        end do
        call CheckReal('T4 min mu over the sparse element nodes', Mu, (5.0d0/3.0d0 + 1.0d0/9.0d0) / 2.0d0, 1.0d-12)

        ! ---------------- T5: sparse region at node (3L, 2L), delta = 3 (P7) ----------------
        call Fill(5, NC)
        call UniformMass()
        Node32 = NodeId(3, 2)
        Sparse = (/ ElLower(3, 2), ElUpper(3, 2), ElUpper(4, 2), ElLower(3, 3), ElLower(4, 3), ElUpper(4, 3) /)
        do i = 1, 6
          call SparseElement(Sparse(i), Node32)
        end do
        call Border(NBorder, NFaces, InPatch)
        call CheckInt('T5 node-based border elements', NBorder, 12)
        call CheckInt('T5 node-based ghost faces', NFaces, 17)
        call CheckTrue('T5 node (3L, 2L) in no ghost patch under P0', .not.InPatch(Node32))
        ! FN-T-002 self-check 13: the horizontal faces y = 4L add the vertex (c-1, 3) of upper (c, 4), c = 1..NC
        Ok = .not.InPatch(NodeId(NC, 3))
        do i = 0, NC
          if (i < NC) Ok = Ok .and. InPatch(NodeId(i, 3))
          Ok = Ok .and. .not.InPatch(NodeId(i, 2)) .and. .not.InPatch(NodeId(i, 1)) .and. .not.InPatch(NodeId(i, 0))
        end do
        call CheckTrue('T5 lowest ghost-patch nodes on y = 3L: (c-1, 3), c = 1..6; none below', Ok)
        call FillCriterion(NBorder, NFaces)
        call CheckInt('T5 fill border elements', NBorder, 0)
        call CheckInt('T5 fill ghost faces', NFaces, 0)
        call Border(NBorder, NFaces)
        call GhostNodalMassFlag(NEL, NNOD, Conn, Candidate, ElDensity, ElSpace, NodalMass, 0.25d0, MassRef, Cut, NAdded, MinMu)
        call CheckReal('T5 mu at (3L, 2L)', NodalMass(Node32) / MassRef(Node32), 0.15d0, 1.0d-12)
        call CheckReal('T5 min mu (S4a diagnostic)', MinMu, 0.15d0, 1.0d-12)
        Ok = .true.
        do i = 1, 6 ! the 6 neighbours of the sparse node
          Ok = Ok .and. abs(NodalMass(Nbr(i)) / MassRef(Nbr(i)) - (4.0d0/3.0d0 + 2.0d0 * 0.475d0) / 2.0d0) <= 1.0d-12
        end do
        call CheckTrue('T5 mu = 1.1417 at the 6 neighbours of (3L, 2L)', Ok)
        call CheckInt('T5 S4a added elements (alpha = 0.25)', NAdded, 6)
        Ok = .true.
        do i = 1, 6
          Ok = Ok .and. Cut(Sparse(i))
        end do
        call CheckTrue('T5 S4a adds exactly the 6 elements around (3L, 2L)', Ok .and. count(Cut) == 18)
        call CheckInt('T5 node-based + S4a ghost faces', CountFaces(), 29)

        ! ---------------- T6: gas row 3 (S1, P5) ----------------
        call Fill(3, NC)
        do c = 1, NC
          Candidate(ElLower(c, 3)) = .false.; Candidate(ElUpper(c, 3)) = .false. ! gas-only: active, not candidate
        end do
        call Border(NBorder, NFaces)
        call CheckInt('T6 S1 border elements (row 2)', NBorder, 12)
        call CheckInt('T6 S1 ghost faces', NFaces, 17)

        ! ---------------- T7: full mesh (P3) ----------------
        call Fill(NR, NC)
        call Border(NBorder, NFaces)
        call CheckInt('T7 border elements', NBorder, 0)
        call CheckInt('T7 ghost faces', NFaces, 0)

        ! ---------------- E3: isolated element and a vertex-only pair ----------------
        Candidate = .false.
        Candidate(ElUpper(2, 2)) = .true.
        call Border(NBorder, NFaces)
        call CheckInt('E3 isolated element flagged', NBorder, 1)
        call CheckInt('E3 isolated element ghost faces', NFaces, 0)
        Candidate(ElLower(3, 3)) = .true. ! shares only the vertex (2L, 2L) with upper (2,2)
        call Border(NBorder, NFaces)
        call CheckInt('E3 vertex-only pair flagged', NBorder, 2)
        call CheckInt('E3 vertex-only pair ghost faces', NFaces, 0)

        ! ---------------- Q-010: axisymmetric fill is radius independent ----------------
        ! element of area A at mean radius rbar: ElementSpace = A rbar, sum w_p = phi A rbar
        Ok = .true.
        do i = 1, 3
          Mu = 10.0d0**(i - 2) ! rbar = 0.1, 1, 10 m
          Ok = Ok .and. GhostFillIsCut(0.5d0 * AREA * Mu, AREA * Mu, 0.6d0)
          Ok = Ok .and. .not.GhostFillIsCut(0.7d0 * AREA * Mu, AREA * Mu, 0.6d0)
        end do
        call CheckTrue('Q-010 fill criterion phi = sum w / ElementSpace independent of radius', Ok)
        call CheckTrue('Q-010 old axisymmetric fill sum(w/r)/ElementSpace depends on radius (reference)', &
                       (0.7d0 * AREA / 0.1d0) / AREA >= 0.6d0 .and. (0.7d0 * AREA / 10.0d0) / AREA < 0.6d0)

        write(*, '(A, I4, A, I4)') 'SUMMARY passed ', NPass, ' failed ', NFail
        if (NFail > 0) call exit(1)

      contains

        subroutine Fill(Rows, Cols)
        ! squares with c <= Cols, r <= Rows are candidates with full fill
          integer(INTEGER_TYPE), intent(in) :: Rows, Cols
          integer(INTEGER_TYPE) :: cc, rr
          Candidate = .false.
          SumW = 0.0
          do cc = 1, Cols
            do rr = 1, Rows
              Candidate(ElLower(cc, rr)) = .true.; Candidate(ElUpper(cc, rr)) = .true.
              SumW(ElLower(cc, rr)) = AREA; SumW(ElUpper(cc, rr)) = AREA
            end do
          end do
        end subroutine Fill

        subroutine UniformMass()
        ! lumped mass of uniformly filled candidates: rho |K| / 3 per vertex (exact for linear N)
          integer(INTEGER_TYPE) :: e, k
          NodalMass = 0.0
          do e = 1, NEL
            if (.not.Candidate(e)) CYCLE
            do k = 1, 3
              NodalMass(Conn(k, e)) = NodalMass(Conn(k, e)) + RHO * AREA / 3.0
            end do
          end do
        end subroutine UniformMass

        subroutine SparseElement(e, Node)
        ! replace the uniform mass of e by one MP of mass rho |K| at barycentric 0.05 (Node), 0.475 (others)
          integer(INTEGER_TYPE), intent(in) :: e, Node
          integer(INTEGER_TYPE) :: k
          do k = 1, 3
            NodalMass(Conn(k, e)) = NodalMass(Conn(k, e)) - RHO * AREA / 3.0
            if (Conn(k, e) == Node) then
              NodalMass(Conn(k, e)) = NodalMass(Conn(k, e)) + 0.05d0 * RHO * AREA
            else
              NodalMass(Conn(k, e)) = NodalMass(Conn(k, e)) + 0.475d0 * RHO * AREA
            end if
          end do
          SumW(e) = AREA
        end subroutine SparseElement

        integer(INTEGER_TYPE) function Nbr(k)
        ! the 6 mesh neighbours of node (3L, 2L) in the diagonal mesh
          integer(INTEGER_TYPE), intent(in) :: k
          integer(INTEGER_TYPE), dimension(2, 6) :: D
          D = reshape((/ 1, 0, -1, 0, 0, 1, 0, -1, 1, 1, -1, -1 /), (/ 2, 6 /))
          Nbr = NodeId(3 + D(1, k), 2 + D(2, k))
        end function Nbr

        integer(INTEGER_TYPE) function CountFaces(Patch)
        ! faces passing IsGhostFacePair for the current Candidate and Cut; optional patch-node marks
          logical, dimension(:), intent(out), optional :: Patch
          integer(INTEGER_TYPE) :: e, s, f, k
          CountFaces = 0
          if (present(Patch)) Patch = .false.
          do e = 1, NEL
            if (.not.Candidate(e)) CYCLE
            do s = 1, 3
              f = Adj(e, s)
              if (.not.IsGhostFacePair(e, f, Candidate, Cut)) CYCLE
              CountFaces = CountFaces + 1
              if (present(Patch)) then
                do k = 1, 3
                  Patch(Conn(k, e)) = .true.
                  Patch(Conn(k, f)) = .true.
                end do
              end if
            end do
          end do
        end function CountFaces

        subroutine Border(NB, NF, Patch)
        ! P0 node-based with S1 on the current Candidate set
          integer(INTEGER_TYPE), intent(out) :: NB, NF
          logical, dimension(:), intent(out), optional :: Patch
          call GhostNodeBasedBorder(NEL, NNOD, Conn, Candidate, Touch, Cut, NB)
          if (present(Patch)) then
            NF = CountFaces(Patch)
          else
            NF = CountFaces()
          end if
        end subroutine Border

        subroutine FillCriterion(NB, NF)
        ! coded fill criterion with phi_G = 0.6 on the current SumW
          integer(INTEGER_TYPE), intent(out) :: NB, NF
          integer(INTEGER_TYPE) :: e
          do e = 1, NEL
            Cut(e) = Candidate(e) .and. GhostFillIsCut(SumW(e), ElSpace(e), 0.6d0)
          end do
          NB = count(Cut)
          NF = CountFaces()
        end subroutine FillCriterion

      end program test_ghost_border


      subroutine GiveError(Message)
      ! stub of the kernel error handler, referenced by ModGlobalConstants (SetElementType)
        character(len=*), intent(in) :: Message
        write(*, '(A, A)') 'GiveError: ', trim(Message)
        call exit(3)
      end subroutine GiveError
