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
	!   and soil-water-structure interaction using the material point method (MPM)
	!
	!	Copyright (C) 2025  Members of the Anura3D MPM Research Community
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


      module ModGhostBorder
      !**********************************************************************
      !
      !    Function:  Border ("cut") element detection of the ghost
      !               stabilisation, on plain arrays (no mesh or MP modules),
      !               so that the same code runs in the kernel and in the
      !               unit test harness (tests/ghost_border).
      !
      !               Design reference: notes/FN-T-002.md (P0 node-based,
      !               form (ii) of D2, S1 candidate activity, S4a nodal-mass
      !               safety net), decision D-002.
      !
      !**********************************************************************

      use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE

      implicit none

        ! $$GHOST_STABILISATION field 7: border detection
        integer(INTEGER_TYPE), parameter :: GHOST_CUT_FILL = 0        ! fill fraction below threshold (coded criterion)
        integer(INTEGER_TYPE), parameter :: GHOST_CUT_NODE_BASED = 1  ! P0 node-based, activity = candidate flag (S1)
        ! $$GHOST_STABILISATION field 10: ghost-mass scheme
        integer(INTEGER_TYPE), parameter :: GHOST_MASS_CODED = 0      ! (M_L + M^g) a = R on ghost DOFs (coded)
        integer(INTEGER_TYPE), parameter :: GHOST_MASS_F1 = 1         ! symmetric operator G_s, FN-T-014 D7 (F1)
        ! $$GHOST_STABILISATION field 11: diagnostics level
        integer(INTEGER_TYPE), parameter :: GHOST_DIAG_BASIC = 0      ! one log line per load step
        integer(INTEGER_TYPE), parameter :: GHOST_DIAG_EXTENDED = 1   ! + border, S4a, A6, contact (S2a) counts
        integer(INTEGER_TYPE), parameter :: GHOST_DIAG_FILE = 2       ! + border elements and ghost faces to file
        integer(INTEGER_TYPE), parameter :: GHOST_DIAG_OPERATOR = 3   ! + F1 operator dump (unit checks FN-T-014 D7.7)

      contains


        logical function GhostFillIsCut(FillVolume, ElSpace, Threshold)
        !**********************************************************************
        !    GhostFillIsCut(FillVolume, ElSpace, Threshold)
        !
        !    Return the coded fill criterion phi_K < phi_G, phi_K = FillVolume / ElSpace.
        !
        !    FillVolume : sum of the MP integration weights of the element, in the
        !                 measure of ElSpace (m^2 plane strain; r-weighted m^3 per
        !                 radian axisymmetric, Q-010)
        !    ElSpace    : element measure ElementSpace
        !    Threshold  : fill threshold phi_G (-)
        !
        !    Example: GhostFillIsCut(0.1d0, 0.5d0, 0.6d0) is .true. (phi_K = 0.2)
        !**********************************************************************
        implicit none

          real(REAL_TYPE), intent(in) :: FillVolume, ElSpace, Threshold

          GhostFillIsCut = .false.
          if (ElSpace > 0.0) GhostFillIsCut = (FillVolume / ElSpace) < Threshold

        end function GhostFillIsCut


        subroutine GhostNodeBasedBorder(NEl, NNod, Conn, Candidate, TouchesNonCandidate, Cut, NBorder)
        !**********************************************************************
        !    GhostNodeBasedBorder(NEl, NNod, Conn, Candidate, TouchesNonCandidate, Cut, NBorder)
        !
        !    Flag the node-based border elements of FN-T-002 D1.4 with the
        !    candidate set as activity (S1): a candidate element is a border
        !    element if one of its vertices is also a vertex of a non-candidate
        !    element (inactive, gas-only, hole-filled, rigid-body or zero-weight).
        !    Form (ii) of D2: two passes over the elements, no adjacency look-up.
        !    Mesh-boundary sides do not flag an element (A4).
        !
        !    Conn                : TRI3 vertices, Conn(1:3, IEl) (ElementConnectivities)
        !    Candidate           : candidate flag of FN-T-001 D3 per element
        !    TouchesNonCandidate : work array (NNod), mutated
        !    Cut                 : border flag per element, overwritten
        !    NBorder             : number of border elements
        !
        !    Example: two triangles sharing a side, only one a candidate: that one
        !    is flagged, NBorder = 1 (tests/ghost_border, case E6).
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NEl, NNod
          integer(INTEGER_TYPE), dimension(:,:), intent(in) :: Conn
          logical, dimension(:), intent(in) :: Candidate
          logical, dimension(:), intent(inout) :: TouchesNonCandidate
          logical, dimension(:), intent(inout) :: Cut
          integer(INTEGER_TYPE), intent(out) :: NBorder

          integer(INTEGER_TYPE) :: IEl, I

          TouchesNonCandidate(1:NNod) = .false.
          do IEl = 1, NEl
            if (Candidate(IEl)) CYCLE
            do I = 1, 3
              TouchesNonCandidate(Conn(I, IEl)) = .true.
            end do
          end do

          NBorder = 0
          do IEl = 1, NEl
            Cut(IEl) = .false.
            if (.not.Candidate(IEl)) CYCLE
            do I = 1, 3
              if (TouchesNonCandidate(Conn(I, IEl))) then
                Cut(IEl) = .true.
                EXIT
              end if
            end do
            if (Cut(IEl)) NBorder = NBorder + 1
          end do

        end subroutine GhostNodeBasedBorder


        subroutine GhostNodalMassFlag(NEl, NNod, Conn, Candidate, ElDensity, ElSpace, NodalMass, Alpha, &
                                      MassRef, Cut, NAdded, MinMu)
        !**********************************************************************
        !    GhostNodalMassFlag(NEl, NNod, Conn, Candidate, ElDensity, ElSpace, NodalMass, Alpha,
        !                       MassRef, Cut, NAdded, MinMu)
        !
        !    Add to Cut every candidate element with a vertex a of mass ratio
        !    mu_a = m_a / mhat_a < Alpha (S4a of FN-T-002), with the reference
        !    mass over the candidate elements of the node only,
        !    mhat_a = sum_{K in E(a), K candidate} rho_K |K| / 3.
        !
        !    ElDensity : rho_K = sum m_p / sum w_p over the non-gas MPs (t/m^3)
        !    ElSpace   : |K|, ElementSpace (same measure as the MP weights)
        !    NodalMass : soil lumped mass m_a per node (LumpedMassDry, t per m)
        !    Alpha     : threshold on mu_a (-)
        !    MassRef   : work array (NNod), mutated; mhat_a on return
        !    Cut       : mutated, flags are only added
        !    NAdded    : number of elements flagged here and not before
        !    MinMu     : min mu_a over nodes with mhat_a > 0 (huge if none)
        !
        !    Example: test T5 of FN-T-002 (mu = 0.15 at the sparse node, Alpha =
        !    0.25) adds the 6 elements around that node (tests/ghost_border).
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NEl, NNod
          integer(INTEGER_TYPE), dimension(:,:), intent(in) :: Conn
          logical, dimension(:), intent(in) :: Candidate
          real(REAL_TYPE), dimension(:), intent(in) :: ElDensity, ElSpace, NodalMass
          real(REAL_TYPE), intent(in) :: Alpha
          real(REAL_TYPE), dimension(:), intent(inout) :: MassRef
          logical, dimension(:), intent(inout) :: Cut
          integer(INTEGER_TYPE), intent(out) :: NAdded
          real(REAL_TYPE), intent(out) :: MinMu

          integer(INTEGER_TYPE) :: IEl, I, NodeID
          real(REAL_TYPE) :: Mu

          MassRef(1:NNod) = 0.0
          do IEl = 1, NEl
            if (.not.Candidate(IEl)) CYCLE
            do I = 1, 3
              NodeID = Conn(I, IEl)
              MassRef(NodeID) = MassRef(NodeID) + ElDensity(IEl) * ElSpace(IEl) / 3.0 ! int_K N_a dA = |K|/3 (TRI3)
            end do
          end do

          MinMu = huge(MinMu)
          do NodeID = 1, NNod
            if (MassRef(NodeID) > 0.0) MinMu = min(MinMu, NodalMass(NodeID) / MassRef(NodeID))
          end do

          NAdded = 0
          do IEl = 1, NEl
            if (.not.Candidate(IEl)) CYCLE
            if (Cut(IEl)) CYCLE
            do I = 1, 3
              NodeID = Conn(I, IEl)
              if (MassRef(NodeID) <= 0.0) CYCLE
              Mu = NodalMass(NodeID) / MassRef(NodeID)
              if (Mu < Alpha) then
                Cut(IEl) = .true.
                NAdded = NAdded + 1
                EXIT
              end if
            end do
          end do

        end subroutine GhostNodalMassFlag


        logical function IsGhostFacePair(IEl, JEl, Candidate, Cut)
        !**********************************************************************
        !    IsGhostFacePair(IEl, JEl, Candidate, Cut)
        !
        !    Return the topological face rules D5.3.1-4 of FN-T-001 for the side
        !    of IEl shared with JEl = ElementAdjacencies(IEl, ISide): interior
        !    side (JEl > 0), counted once (JEl > IEl), both elements candidates,
        !    at least one cut. The geometric and contact rules stay in
        !    BuildGhostFaces.
        !
        !    Example: IsGhostFacePair(1, 0, C, X) is .false. (mesh boundary).
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl, JEl
          logical, dimension(:), intent(in) :: Candidate, Cut

          IsGhostFacePair = .false.
          if (JEl <= 0) RETURN
          if (JEl <= IEl) RETURN
          if (.not.Candidate(IEl)) RETURN
          if (.not.Candidate(JEl)) RETURN
          IsGhostFacePair = Cut(IEl) .or. Cut(JEl)

        end function IsGhostFacePair

      end module ModGhostBorder
