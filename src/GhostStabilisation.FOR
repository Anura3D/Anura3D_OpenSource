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


      module ModGhostStabilisation
      !**********************************************************************
      !
      !    Function:  Ghost stabilisation (Burman ghost penalty, Coombs 2023)
      !               adapted to the explicit dynamic MPM of Anura3D.
      !
      !               Design reference: docs/Ghost_Stabilisation_Anura3D.md
      !
      !               Scope of this implementation:
      !               - 2D plane strain and 2D axisymmetric, TRI3 elements
      !               - NFORMULATION = 1, single phase (dry)
      !               - acts on the soil entity (SOFT_ENTITY) only
      !
      !               For linear elements only the jump of the first normal
      !               derivative is penalised. Per ghost face F:
      !                 g_a  = (grad N_a|K+ - grad N_a|K-) . n_F
      !                 K^g  = gammaK * M* * h_F * w_F * |F| * g g^T   (per direction)
      !                 M^g  = gammaM * rho * h_F^3 * w_F * |F| * g g^T
      !               with w_F = 1 (plane strain) or r_F (axisymmetric, per radian,
      !               consistent with the r-scaled MP weights).
      !
      !               Explicit adaptation:
      !               - ghost force  f^g = dt * K^g * v   (ghost viscosity),
      !                 subtracted from the soil rate of momentum
      !               - ghost mass   (M_L + M^g) a = R, solved by Jacobi sweeps
      !                 (or diagonal only), on ghost-patch DOFs
      !               - critical time step: dt <= sqrt(2 m~_a / sum_b |K^g_ab|)
      !
      !               Sum_a g_a = 0, so f^g conserves linear momentum exactly and
      !               the rows of M^g sum to zero.
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModCounters
      use ModReadCalculationData
      use ModReadMaterialData
      use ModMeshInfo
      use ModMeshAdjacencies
      use ModMPMData
      use ModParticle
      use ModMPMDynContact
      use ModConvectivePhase
      use ModString
      use ModFeedback
      use ModGasification

      implicit none

        integer(INTEGER_TYPE), parameter, private :: GHOST_PATCH_NODES = 4 ! TRI3: 3 + 3 nodes, 2 shared

        integer(INTEGER_TYPE), private :: NGhostFaces = 0
        integer(INTEGER_TYPE), dimension(:,:), allocatable, private :: GhostFaceNodes   ! (GHOST_PATCH_NODES, MaxFaces) global node IDs
        real(REAL_TYPE), dimension(:,:), allocatable, private :: GhostFaceJump          ! (GHOST_PATCH_NODES, MaxFaces) normal-gradient jumps g_a
        real(REAL_TYPE), dimension(:), allocatable, private :: GhostFaceCoefK            ! gammaK * M* * h * w * |F|
        real(REAL_TYPE), dimension(:), allocatable, private :: GhostFaceCoefM            ! gammaM * rho * h^3 * w * |F|
        real(REAL_TYPE), dimension(:), allocatable, private :: GhostMassDiagonal         ! (Counters%N) diagonal of M^g
        real(REAL_TYPE), dimension(:), allocatable, private :: GhostStiffRowSum          ! (Counters%N) sum_b |K^g_ab| (without dt)
        real(REAL_TYPE), dimension(:), allocatable, private :: GhostMassOffDiagAbs       ! (Counters%N) sum_{b/=a} |M^g_ab| (Gershgorin, Jacobi damping)
        logical, dimension(:), allocatable, private :: IsGhostDof                        ! (Counters%N) DOF belongs to a ghost patch
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: GhostNodeParent     ! (Counters%NodTot) union-find parent, 0 = not a ghost node
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: GhostNodeCluster    ! (Counters%NodTot) cluster index of ghost nodes, 0 otherwise
        integer(INTEGER_TYPE), private :: NGhostClusters = 0
        logical, private :: GhostIsAvailable = .false.  ! element type / formulation supported
        logical, private :: GhostChecked = .false.      ! support check done
        real(REAL_TYPE), private :: GhostLastMomentumError = 0.0

      contains


        logical function IsGhostActive()
        !**********************************************************************
        !    Function:  True if ghost stabilisation is requested and supported
        !**********************************************************************
        implicit none

          IsGhostActive = .false.
          if (.not.CalParams%ApplyGhostStabilisation) RETURN

          if (.not.GhostChecked) then
            GhostChecked = .true.
            GhostIsAvailable = (NDIM == 2) .and. (ELEMENTNODES == 3) .and. (NFORMULATION == 1) &
                               .and. (CalParams%NumberOfPhases == 1)
            if (.not.GhostIsAvailable) then
              call WriteInLogFile('Ghost stabilisation requested but only implemented for 2D TRI3, ' // &
                                  'single-phase, 1-layer formulation. It is switched off.')
            else
              call WriteInLogFile('Ghost stabilisation active: gammaK = ' // trim(String(CalParams%GhostGammaK)) // &
                                  ', gammaM = ' // trim(String(CalParams%GhostGammaM)) // &
                                  ', fill threshold = ' // trim(String(CalParams%GhostFillThreshold)) // &
                                  ', mass solver = ' // trim(String(CalParams%GhostMassSolver)) // &
                                  ', Jacobi sweeps = ' // trim(String(CalParams%GhostNJacobi)))
            end if
          end if

          IsGhostActive = GhostIsAvailable

        end function IsGhostActive


        subroutine AllocateGhostArrays()
        !**********************************************************************
        !    Function:  (Re)allocate the ghost arrays for the current mesh size
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: MaxFaces, IError

          MaxFaces = max(1, (Counters%NEl * ELEMENTSIDES) / 2 + 1)

          if (allocated(GhostFaceNodes)) then
            if (size(GhostFaceNodes, 2) >= MaxFaces .and. size(IsGhostDof) == Counters%N &
                .and. size(GhostNodeParent) == Counters%NodTot) RETURN
            call DestroyGhostArrays()
          end if

          allocate(GhostFaceNodes(GHOST_PATCH_NODES, MaxFaces), stat = IError)
          allocate(GhostFaceJump(GHOST_PATCH_NODES, MaxFaces), stat = IError)
          allocate(GhostFaceCoefK(MaxFaces), stat = IError)
          allocate(GhostFaceCoefM(MaxFaces), stat = IError)
          allocate(GhostMassDiagonal(Counters%N), stat = IError)
          allocate(GhostStiffRowSum(Counters%N), stat = IError)
          allocate(GhostMassOffDiagAbs(Counters%N), stat = IError)
          allocate(IsGhostDof(Counters%N), stat = IError)
          allocate(GhostNodeParent(Counters%NodTot), stat = IError)
          allocate(GhostNodeCluster(Counters%NodTot), stat = IError)

          GhostNodeParent = 0
          GhostNodeCluster = 0
          NGhostClusters = 0
          GhostFaceNodes = 0
          GhostFaceJump = 0.0
          GhostFaceCoefK = 0.0
          GhostFaceCoefM = 0.0
          GhostMassDiagonal = 0.0
          GhostStiffRowSum = 0.0
          GhostMassOffDiagAbs = 0.0
          IsGhostDof = .false.
          NGhostFaces = 0

        end subroutine AllocateGhostArrays


        subroutine DestroyGhostArrays()
        !**********************************************************************
        !    Function:  Deallocate the ghost arrays
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IError

          if (allocated(GhostFaceNodes)) deallocate(GhostFaceNodes, stat = IError)
          if (allocated(GhostFaceJump)) deallocate(GhostFaceJump, stat = IError)
          if (allocated(GhostFaceCoefK)) deallocate(GhostFaceCoefK, stat = IError)
          if (allocated(GhostFaceCoefM)) deallocate(GhostFaceCoefM, stat = IError)
          if (allocated(GhostMassDiagonal)) deallocate(GhostMassDiagonal, stat = IError)
          if (allocated(GhostStiffRowSum)) deallocate(GhostStiffRowSum, stat = IError)
          if (allocated(GhostMassOffDiagAbs)) deallocate(GhostMassOffDiagAbs, stat = IError)
          if (allocated(IsGhostDof)) deallocate(IsGhostDof, stat = IError)
          if (allocated(GhostNodeParent)) deallocate(GhostNodeParent, stat = IError)
          if (allocated(GhostNodeCluster)) deallocate(GhostNodeCluster, stat = IError)
          NGhostFaces = 0
          NGhostClusters = 0

        end subroutine DestroyGhostArrays


        subroutine GetElementGradientsTRI3(IEl, GradN, Area2)
        !**********************************************************************
        !    Function:  Constant shape function gradients of a TRI3 element
        !               from the current nodal coordinates
        !
        !    GradN(I, J) : dN_J / dx_I
        !    Area2       : twice the signed element area
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl
          real(REAL_TYPE), dimension(2, 3), intent(out) :: GradN
          real(REAL_TYPE), intent(out) :: Area2

          real(REAL_TYPE) :: x1, y1, x2, y2, x3, y3

          x1 = NodalCoordinatesUpd(ElementConnectivities(1, IEl), 1)
          y1 = NodalCoordinatesUpd(ElementConnectivities(1, IEl), 2)
          x2 = NodalCoordinatesUpd(ElementConnectivities(2, IEl), 1)
          y2 = NodalCoordinatesUpd(ElementConnectivities(2, IEl), 2)
          x3 = NodalCoordinatesUpd(ElementConnectivities(3, IEl), 1)
          y3 = NodalCoordinatesUpd(ElementConnectivities(3, IEl), 2)

          Area2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
          GradN = 0.0
          if (abs(Area2) <= 1.0d-300) RETURN

          GradN(1, 1) = (y2 - y3) / Area2
          GradN(2, 1) = (x3 - x2) / Area2
          GradN(1, 2) = (y3 - y1) / Area2
          GradN(2, 2) = (x1 - x3) / Area2
          GradN(1, 3) = (y1 - y2) / Area2
          GradN(2, 3) = (x2 - x1) / Area2

        end subroutine GetElementGradientsTRI3


        real(REAL_TYPE) function GetElementPWaveModulus(IEl)
        !**********************************************************************
        !    Function:  max rho*c^2 over the MPs of element IEl, using the same
        !               wave speed as the critical time step
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl
          integer(INTEGER_TYPE) :: IPart, ParticleIndex
          real(REAL_TYPE) :: WaveSpeed, Weight

          GetElementPWaveModulus = 0.0
          do IPart = 1, NPartEle(IEl)
            ParticleIndex = GetParticleIndex(IPart, IEl)
            Weight = Particles(ParticleIndex)%IntegrationWeight
            if (Weight <= 0.0) CYCLE
            WaveSpeed = 0.0
            call GetWaveSpeedForTimeStepSize(ParticleIndex, IEl, WaveSpeed)
            if (WaveSpeed > 0.0) then
              GetElementPWaveModulus = max(GetElementPWaveModulus, MassArray(ParticleIndex) / Weight * WaveSpeed**2)
            end if
          end do

        end function GetElementPWaveModulus


        subroutine GetElementGhostData(IEl, IsCandidate, IsCut, Density)
        !**********************************************************************
        !    Function:  Classify an element for ghost stabilisation
        !
        !    IsCandidate  : active, holds MPs, and all MPs are soil entity
        !    IsCut        : partially filled (particle integration, or fill below threshold)
        !    Density      : mean MP density in the element (mass / weight)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl
          logical, intent(out) :: IsCandidate, IsCut
          real(REAL_TYPE), intent(out) :: Density

          integer(INTEGER_TYPE) :: IPart, ParticleIndex, MaterialIndex
          real(REAL_TYPE) :: SumMass, SumWeight, SumVolume, Weight, Radius

          IsCandidate = .false.
          IsCut = .false.
          Density = 0.0

          if (.not.IsActiveElement(IEl)) RETURN
          if (NPartEle(IEl) <= 0) RETURN

          SumMass = 0.0
          SumWeight = 0.0
          SumVolume = 0.0

          do IPart = 1, NPartEle(IEl)
            ParticleIndex = GetParticleIndex(IPart, IEl)

            if (CalParams%ApplyContactAlgorithm) then
              if (EntityIDArray(ParticleIndex) /= SOFT_ENTITY) RETURN ! element touches another entity (e.g. rigid body)
            end if
            MaterialIndex = MaterialIDArray(ParticleIndex)
            if (MatParams(MaterialIndex)%MaterialModel == ESM_RIGID_BODY) RETURN

            if (IsGasParticle(ParticleIndex)) CYCLE ! gasification: gas MPs carry no stiffness and do not fill

            Weight = Particles(ParticleIndex)%IntegrationWeight
            if (Weight <= 0.0) CYCLE

            SumMass = SumMass + MassArray(ParticleIndex)
            SumWeight = SumWeight + Weight

            ! geometric volume of the MP (the axisymmetric weight carries the radius)
            if (ISAXISYMMETRIC) then
              Radius = GlobPosArray(ParticleIndex, 1)
              if (Radius > 0.0) then
                SumVolume = SumVolume + Weight / Radius
              end if
            else
              SumVolume = SumVolume + Weight
            end if
          end do

          if (SumWeight <= 0.0) RETURN

          IsCandidate = .true.
          Density = SumMass / SumWeight

          ! Cut = partially filled. Only the MP fill is used: IsParticleIntegration follows the
          ! degree of filling (0.9 by default), and the fill of interior elements fluctuates
          ! around that value as MPs migrate, which would put ghost faces everywhere.
          if (ElementSpace(IEl) > 0.0) then
            IsCut = (SumVolume / ElementSpace(IEl)) < CalParams%GhostFillThreshold
          end if

        end subroutine GetElementGhostData


        integer(INTEGER_TYPE) function GhostFindRoot(NodeID)
        !**********************************************************************
        !    Function:  Union-find root of a ghost node (path halving)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodeID
          integer(INTEGER_TYPE) :: N

          N = NodeID
          do while (GhostNodeParent(N) /= N)
            GhostNodeParent(N) = GhostNodeParent(GhostNodeParent(N))
            N = GhostNodeParent(N)
          end do
          GhostFindRoot = N

        end function GhostFindRoot


        subroutine GhostUnion(NodeA, NodeB)
        !**********************************************************************
        !    Function:  Merge the ghost clusters of two nodes
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodeA, NodeB
          integer(INTEGER_TYPE) :: RootA, RootB

          if (GhostNodeParent(NodeA) == 0) GhostNodeParent(NodeA) = NodeA
          if (GhostNodeParent(NodeB) == 0) GhostNodeParent(NodeB) = NodeB
          RootA = GhostFindRoot(NodeA)
          RootB = GhostFindRoot(NodeB)
          if (RootA /= RootB) GhostNodeParent(max(RootA, RootB)) = min(RootA, RootB)

        end subroutine GhostUnion


        logical function IsContactNode(NodeID)
        !**********************************************************************
        !    Function:  True if the node takes part in the Bardenhagen contact
        !               or carries mass of an entity other than the soil
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodeID
          integer(INTEGER_TYPE) :: IEntity, IDof

          IsContactNode = .false.
          if (.not.CalParams%ApplyContactAlgorithm) RETURN

          if (allocated(InterfaceNodes)) then
            if (InterfaceNodes(NodeID)) then
              IsContactNode = .true.
              RETURN
            end if
          end if

          IDof = ReducedDof(NodeID) + 1
          do IEntity = 1, Counters%nEntity
            if (IEntity == SOFT_ENTITY) CYCLE
            if (LumpedMassDry(IDof, IEntity) /= 0.0) then
              IsContactNode = .true.
              RETURN
            end if
          end do

        end function IsContactNode


        subroutine BuildGhostFaces()
        !**********************************************************************
        !    Function:  Build the list of ghost faces for the current MP
        !               configuration (called once per time step, after the
        !               nodal masses have been mapped)
        !
        !    A face is a ghost face if:
        !     - it is shared by two candidate elements (active, soil only)
        !     - at least one of the two elements is cut (partially filled)
        !     - none of its patch nodes is a contact / other-entity node
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IEl, JEl, ISide, I, J, K, NodeID, IDof
          integer(INTEGER_TYPE) :: SideNode1, SideNode2
          integer(INTEGER_TYPE), dimension(2, 3) :: SideVertices
          integer(INTEGER_TYPE), dimension(GHOST_PATCH_NODES) :: PatchNodes
          integer(INTEGER_TYPE) :: NPatch
          logical, dimension(:), allocatable :: Candidate, Cut
          real(REAL_TYPE), dimension(:), allocatable :: ElModulus, ElDensity
          real(REAL_TYPE), dimension(2, 3) :: GradI, GradJ
          real(REAL_TYPE), dimension(2) :: Normal, GradJumpI, GradJumpJ
          real(REAL_TYPE), dimension(GHOST_PATCH_NODES) :: Jump
          real(REAL_TYPE) :: AreaI, AreaJ, FaceLength, H, Weight, Modulus, Rho, dx, dy
          logical :: IsValid
          integer(INTEGER_TYPE) :: IError

          NGhostFaces = 0
          if (.not.IsGhostActive()) RETURN
          if (.not.allocated(ElementAdjacencies)) RETURN

          call AllocateGhostArrays()
          GhostMassDiagonal = 0.0
          GhostStiffRowSum = 0.0
          GhostMassOffDiagAbs = 0.0
          IsGhostDof = .false.
          GhostNodeParent = 0
          GhostNodeCluster = 0
          NGhostClusters = 0

          ! local vertices of the TRI3 sides, same convention as CheckEdgeNodeAdjacencies
          SideVertices = reshape((/ 1, 2, 3, 1, 2, 3 /), (/ 2, 3 /))

          allocate(Candidate(Counters%NEl), Cut(Counters%NEl), ElModulus(Counters%NEl), ElDensity(Counters%NEl), stat = IError)
          Candidate = .false.
          Cut = .false.
          ElModulus = -1.0 ! not yet evaluated (evaluated only for elements of ghost faces)
          ElDensity = 0.0

          do IEl = 1, Counters%NEl
            if (.not.IsActiveElement(IEl)) CYCLE
            call GetElementGhostData(IEl, Candidate(IEl), Cut(IEl), ElDensity(IEl))
          end do

          do IEl = 1, Counters%NEl
            if (.not.Candidate(IEl)) CYCLE

            do ISide = 1, ELEMENTSIDES
              JEl = ElementAdjacencies(IEl, ISide)
              if (JEl <= 0) CYCLE              ! mesh boundary
              if (JEl <= IEl) CYCLE            ! each face once
              if (.not.Candidate(JEl)) CYCLE
              if (.not.(Cut(IEl).or.Cut(JEl))) CYCLE

              ! patch nodes: nodes of IEl followed by the vertex of JEl not shared with IEl
              PatchNodes = 0
              NPatch = 0
              do I = 1, 3
                NPatch = NPatch + 1
                PatchNodes(NPatch) = ElementConnectivities(I, IEl)
              end do
              do I = 1, 3
                NodeID = ElementConnectivities(I, JEl)
                if (all(PatchNodes(1:3) /= NodeID)) then
                  NPatch = NPatch + 1
                  if (NPatch > GHOST_PATCH_NODES) EXIT
                  PatchNodes(NPatch) = NodeID
                end if
              end do
              if (NPatch /= GHOST_PATCH_NODES) CYCLE ! not a proper side neighbour

              IsValid = .true.
              do I = 1, GHOST_PATCH_NODES
                if (IsContactNode(PatchNodes(I))) then
                  IsValid = .false.
                  EXIT
                end if
              end do
              if (.not.IsValid) CYCLE

              ! face geometry
              SideNode1 = ElementConnectivities(SideVertices(1, ISide), IEl)
              SideNode2 = ElementConnectivities(SideVertices(2, ISide), IEl)
              dx = NodalCoordinatesUpd(SideNode2, 1) - NodalCoordinatesUpd(SideNode1, 1)
              dy = NodalCoordinatesUpd(SideNode2, 2) - NodalCoordinatesUpd(SideNode1, 2)
              FaceLength = sqrt(dx*dx + dy*dy)
              if (FaceLength <= 1.0d-300) CYCLE
              Normal(1) =  dy / FaceLength
              Normal(2) = -dx / FaceLength

              if (ISAXISYMMETRIC) then
                Weight = 0.5 * (NodalCoordinatesUpd(SideNode1, 1) + NodalCoordinatesUpd(SideNode2, 1))
                if (Weight <= 0.0) CYCLE
              else
                Weight = 1.0
              end if

              ! normal-gradient jumps of the patch shape functions
              call GetElementGradientsTRI3(IEl, GradI, AreaI)
              call GetElementGradientsTRI3(JEl, GradJ, AreaJ)
              if ((abs(AreaI) <= 1.0d-300) .or. (abs(AreaJ) <= 1.0d-300)) CYCLE

              Jump = 0.0
              do I = 1, GHOST_PATCH_NODES
                GradJumpI = 0.0
                GradJumpJ = 0.0
                do J = 1, 3
                  if (ElementConnectivities(J, IEl) == PatchNodes(I)) GradJumpI = GradI(:, J)
                  if (ElementConnectivities(J, JEl) == PatchNodes(I)) GradJumpJ = GradJ(:, J)
                end do
                Jump(I) = (GradJumpI(1) - GradJumpJ(1)) * Normal(1) + (GradJumpI(2) - GradJumpJ(2)) * Normal(2)
              end do

              if (ElModulus(IEl) < 0.0) ElModulus(IEl) = GetElementPWaveModulus(IEl)
              if (ElModulus(JEl) < 0.0) ElModulus(JEl) = GetElementPWaveModulus(JEl)
              Modulus = max(ElModulus(IEl), ElModulus(JEl))
              Rho = 0.5 * (ElDensity(IEl) + ElDensity(JEl))
              H = min(ElementLMin(IEl), ElementLMin(JEl))
              if ((Modulus <= 0.0) .or. (Rho <= 0.0) .or. (H <= 0.0)) CYCLE

              NGhostFaces = NGhostFaces + 1
              GhostFaceNodes(:, NGhostFaces) = PatchNodes
              GhostFaceJump(:, NGhostFaces) = Jump
              GhostFaceCoefK(NGhostFaces) = CalParams%GhostGammaK * Modulus * H * Weight * FaceLength
              GhostFaceCoefM(NGhostFaces) = CalParams%GhostGammaM * Rho * H**3 * Weight * FaceLength

              ! connected ghost clusters (used for the momentum correction of the ghost mass)
              do I = 2, GHOST_PATCH_NODES
                call GhostUnion(PatchNodes(1), PatchNodes(I))
              end do

              ! diagonal ghost mass, Gershgorin row sum of K^g, ghost DOF flags
              do I = 1, GHOST_PATCH_NODES
                do K = 1, NVECTOR
                  IDof = ReducedDof(PatchNodes(I)) + K
                  IsGhostDof(IDof) = .true.
                  GhostMassDiagonal(IDof) = GhostMassDiagonal(IDof) + GhostFaceCoefM(NGhostFaces) * Jump(I)**2
                  GhostStiffRowSum(IDof) = GhostStiffRowSum(IDof) + &
                                           GhostFaceCoefK(NGhostFaces) * abs(Jump(I)) * sum(abs(Jump))
                  GhostMassOffDiagAbs(IDof) = GhostMassOffDiagAbs(IDof) + &
                                           GhostFaceCoefM(NGhostFaces) * abs(Jump(I)) * (sum(abs(Jump)) - abs(Jump(I)))
                end do
              end do

            end do ! sides
          end do ! elements

          deallocate(Candidate, Cut, ElModulus, ElDensity, stat = IError)

          ! number the clusters 1..NGhostClusters (roots first, then all members)
          do NodeID = 1, Counters%NodTot
            if (GhostNodeParent(NodeID) == 0) CYCLE
            if (GhostFindRoot(NodeID) == NodeID) then
              NGhostClusters = NGhostClusters + 1
              GhostNodeCluster(NodeID) = NGhostClusters
            end if
          end do
          do NodeID = 1, Counters%NodTot
            if (GhostNodeParent(NodeID) == 0) CYCLE
            GhostNodeCluster(NodeID) = GhostNodeCluster(GhostFindRoot(NodeID))
          end do

        end subroutine BuildGhostFaces


        subroutine ApplyGhostForce(RateOfMomentumSoil)
        !**********************************************************************
        !    Function:  Subtract the ghost (viscous) force f^g = dt K^g v from
        !               the soil rate of momentum. The force is divided by the
        !               mass scaling factor, as all forces in the rate of momentum.
        !
        !    RateOfMomentumSoil : rate of momentum of the soil entity (Counters%N)
        !**********************************************************************
        implicit none

          real(REAL_TYPE), dimension(:), intent(inout) :: RateOfMomentumSoil

          integer(INTEGER_TYPE) :: IFace, I, K, IDof
          real(REAL_TYPE) :: JumpVel, Coef

          if (.not.IsGhostActive()) RETURN

          call BuildGhostFaces()
          if (NGhostFaces == 0) RETURN

          do IFace = 1, NGhostFaces
            Coef = GhostFaceCoefK(IFace) * CalParams%TimeIncrement / CalParams%ScalingMassFactor
            do K = 1, NVECTOR
              ! jump of the normal velocity gradient across the face, direction K
              JumpVel = 0.0
              do I = 1, GHOST_PATCH_NODES
                IDof = ReducedDof(GhostFaceNodes(I, IFace)) + K
                JumpVel = JumpVel + GhostFaceJump(I, IFace) * TotalVelocitySoil(IDof, SOFT_ENTITY)
              end do
              do I = 1, GHOST_PATCH_NODES
                IDof = ReducedDof(GhostFaceNodes(I, IFace)) + K
                RateOfMomentumSoil(IDof) = RateOfMomentumSoil(IDof) - Coef * GhostFaceJump(I, IFace) * JumpVel
              end do
            end do
          end do

        end subroutine ApplyGhostForce


        subroutine GhostMassTimesVector(Vec, Res)
        !**********************************************************************
        !    Function:  Res = (M^g - diag(M^g)) * Vec, assembled face by face
        !**********************************************************************
        implicit none

          real(REAL_TYPE), dimension(:), intent(in) :: Vec
          real(REAL_TYPE), dimension(:), intent(out) :: Res

          integer(INTEGER_TYPE) :: IFace, I, K, IDof
          real(REAL_TYPE) :: JumpVec

          Res = 0.0
          do IFace = 1, NGhostFaces
            do K = 1, NVECTOR
              JumpVec = 0.0
              do I = 1, GHOST_PATCH_NODES
                IDof = ReducedDof(GhostFaceNodes(I, IFace)) + K
                JumpVec = JumpVec + GhostFaceJump(I, IFace) * Vec(IDof)
              end do
              do I = 1, GHOST_PATCH_NODES
                IDof = ReducedDof(GhostFaceNodes(I, IFace)) + K
                Res(IDof) = Res(IDof) + GhostFaceCoefM(IFace) * GhostFaceJump(I, IFace) * JumpVec
              end do
            end do
          end do

          ! remove the diagonal part
          do IDof = 1, Counters%N
            if (IsGhostDof(IDof)) Res(IDof) = Res(IDof) - GhostMassDiagonal(IDof) * Vec(IDof)
          end do

        end subroutine GhostMassTimesVector


        subroutine ApplyGhostMass()
        !**********************************************************************
        !    Function:  Recompute the soil acceleration on ghost-patch DOFs with
        !               the ghost mass: (M_L + M^g) a = R.
        !               MassSolver = 0 : diagonal only, a = R / (m + M^g_aa)
        !               MassSolver = 1 : Jacobi sweeps on the consistent system
        !               followed by a per-cluster uniform momentum correction, so
        !               that linear momentum is conserved exactly.
        !               The system acceleration of these DOFs is set equal to the
        !               soil acceleration (ghost patches never touch contact nodes,
        !               so the soil is the only entity there).
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IDof, ISweep, IError, NodeID, IC, K
          real(REAL_TYPE), dimension(:), allocatable :: Acc, OffDiag
          real(REAL_TYPE), dimension(:,:), allocatable :: ClusterDP, ClusterMass
          real(REAL_TYPE) :: Diag, SumF, SumMA, MaxF

          if (.not.IsGhostActive()) RETURN
          if (NGhostFaces == 0) RETURN
          if (CalParams%GhostGammaM <= 0.0) RETURN

          allocate(Acc(Counters%N), OffDiag(Counters%N), stat = IError)
          Acc = AccelerationSoil(:, SOFT_ENTITY)

          ! initial guess: diagonal solve
          do IDof = 1, Counters%N
            if (.not.IsGhostDof(IDof)) CYCLE
            Diag = LumpedMassDry(IDof, SOFT_ENTITY) + GhostMassDiagonal(IDof)
            if (Diag > 0.0) then
              Acc(IDof) = RateofMomentum(IDof, SOFT_ENTITY) / Diag * PBoundary(IDof)
            else
              Acc(IDof) = 0.0
            end if
          end do

          if (CalParams%GhostMassSolver == 1) then
            ! Preconditioned Richardson (damped Jacobi) iteration on A a = R, A = M_L + M^g:
            !   a <- a + (R - A a) / P,   P_i = D_i + sum_{j/=i} |M^g_ij|
            ! 2P - A is diagonally dominant with positive diagonal, hence SPD, so the
            ! iteration converges for any mass distribution; at convergence
            ! sum_a m_a a_a = sum_a R_a (rows of M^g sum to zero).
            do ISweep = 1, CalParams%GhostNJacobi
              call GhostMassTimesVector(Acc, OffDiag)
              do IDof = 1, Counters%N
                if (.not.IsGhostDof(IDof)) CYCLE
                Diag = LumpedMassDry(IDof, SOFT_ENTITY) + GhostMassDiagonal(IDof)
                if (Diag > 0.0) then
                  Acc(IDof) = Acc(IDof) + (RateofMomentum(IDof, SOFT_ENTITY) - Diag * Acc(IDof) - OffDiag(IDof)) &
                                        / (Diag + GhostMassOffDiagAbs(IDof))
                  Acc(IDof) = Acc(IDof) * PBoundary(IDof)
                end if
              end do
            end do
          end if

          ! Momentum correction per connected ghost cluster and direction:
          !   a_i <- a_i + (sum_c R - sum_c m a) / sum_c m   on the free DOFs of cluster c
          ! A uniform acceleration is affine, so it produces no ghost jump, and it
          ! restores sum_c m a = sum_c R exactly (also for the diagonal solver and for
          ! modes the few Jacobi sweeps have not converged).
          allocate(ClusterDP(NGhostClusters, NVECTOR), ClusterMass(NGhostClusters, NVECTOR), stat = IError)
          ClusterDP = 0.0
          ClusterMass = 0.0
          do NodeID = 1, Counters%NodTot
            IC = GhostNodeCluster(NodeID)
            if (IC <= 0) CYCLE
            do K = 1, NVECTOR
              IDof = ReducedDof(NodeID) + K
              if (PBoundary(IDof) == 0.0) CYCLE
              ClusterDP(IC, K) = ClusterDP(IC, K) + RateofMomentum(IDof, SOFT_ENTITY) - LumpedMassDry(IDof, SOFT_ENTITY) * Acc(IDof)
              ClusterMass(IC, K) = ClusterMass(IC, K) + LumpedMassDry(IDof, SOFT_ENTITY)
            end do
          end do
          do NodeID = 1, Counters%NodTot
            IC = GhostNodeCluster(NodeID)
            if (IC <= 0) CYCLE
            do K = 1, NVECTOR
              IDof = ReducedDof(NodeID) + K
              if (PBoundary(IDof) == 0.0) CYCLE
              if (ClusterMass(IC, K) > 0.0) Acc(IDof) = Acc(IDof) + ClusterDP(IC, K) / ClusterMass(IC, K)
            end do
          end do
          deallocate(ClusterDP, ClusterMass, stat = IError)

          ! momentum error on the ghost DOFs after correction: |sum(m a) - sum(R)| / sum(|R|)
          SumF = 0.0
          SumMA = 0.0
          MaxF = 0.0
          do IDof = 1, Counters%N
            if (.not.IsGhostDof(IDof)) CYCLE
            if (PBoundary(IDof) == 0.0) CYCLE
            SumF = SumF + RateofMomentum(IDof, SOFT_ENTITY)
            SumMA = SumMA + LumpedMassDry(IDof, SOFT_ENTITY) * Acc(IDof)
            MaxF = MaxF + abs(RateofMomentum(IDof, SOFT_ENTITY))
          end do
          if (MaxF > 0.0) then
            GhostLastMomentumError = abs(SumMA - SumF) / MaxF
          else
            GhostLastMomentumError = 0.0
          end if

          do IDof = 1, Counters%N
            if (.not.IsGhostDof(IDof)) CYCLE
            AccelerationSoil(IDof, SOFT_ENTITY) = Acc(IDof)
            if (CalParams%ApplyContactAlgorithm) then
              if (allocated(AccelerationSys)) AccelerationSys(IDof) = Acc(IDof)
            end if
          end do

          deallocate(Acc, OffDiag, stat = IError)

        end subroutine ApplyGhostMass


        real(REAL_TYPE) function GhostCriticalTimeStep()
        !**********************************************************************
        !    Function:  Time step bound of the explicit ghost viscosity
        !               (Gershgorin): dt^2 * sum_b |K^g_ab| / m~_a <= 2
        !               Returns huge() if ghost stabilisation is not active.
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IDof
          real(REAL_TYPE) :: MassTilde

          GhostCriticalTimeStep = huge(GhostCriticalTimeStep)
          if (.not.IsGhostActive()) RETURN
          if (NGhostFaces == 0) RETURN

          do IDof = 1, Counters%N
            if (.not.IsGhostDof(IDof)) CYCLE
            if (GhostStiffRowSum(IDof) <= 0.0) CYCLE
            ! (m + M^g) a = R - (dt K^g / s) v  =>  dt^2 sum|K^g| / (s (m + M^g)) <= 2
            MassTilde = (LumpedMassDry(IDof, SOFT_ENTITY) + GhostMassDiagonal(IDof)) * CalParams%ScalingMassFactor
            if (MassTilde <= 0.0) CYCLE
            GhostCriticalTimeStep = min(GhostCriticalTimeStep, sqrt(2.0 * MassTilde / GhostStiffRowSum(IDof)))
          end do

        end function GhostCriticalTimeStep


        subroutine WriteGhostDiagnostics()
        !**********************************************************************
        !    Function:  Write ghost statistics to the log file
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IDof, NDofs
          real(REAL_TYPE) :: MaxRatio, MinMassRatio, MeanMass, SumMass

          if (.not.IsGhostActive()) RETURN

          NDofs = 0
          MaxRatio = 0.0
          SumMass = 0.0
          do IDof = 1, Counters%N
            if (LumpedMassDry(IDof, SOFT_ENTITY) > 0.0) then
              NDofs = NDofs + 1
              SumMass = SumMass + LumpedMassDry(IDof, SOFT_ENTITY)
            end if
          end do
          MeanMass = 0.0
          if (NDofs > 0) MeanMass = SumMass / NDofs

          MinMassRatio = huge(MinMassRatio)
          do IDof = 1, Counters%N
            if (LumpedMassDry(IDof, SOFT_ENTITY) <= 0.0) CYCLE
            if (MeanMass > 0.0) MinMassRatio = min(MinMassRatio, LumpedMassDry(IDof, SOFT_ENTITY) / MeanMass)
            if (IsGhostDof(IDof)) then
              MaxRatio = max(MaxRatio, GhostMassDiagonal(IDof) / LumpedMassDry(IDof, SOFT_ENTITY))
            end if
          end do
          if (MinMassRatio > 1.0e30) MinMassRatio = 0.0

          call WriteInLogFile('Ghost: faces = ' // trim(String(NGhostFaces)) // &
                              ', max Mg/m = ' // trim(String(MaxRatio)) // &
                              ', min m/mean m = ' // trim(String(MinMassRatio)) // &
                              ', ghost dt = ' // trim(String(GhostCriticalTimeStep())) // &
                              ', momentum error = ' // trim(String(GhostLastMomentumError)))

        end subroutine WriteGhostDiagnostics


        logical function IsGhostElementForOutput(IEl)
        !**********************************************************************
        !    Function:  True if element IEl is part of a ghost face (output only)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl
          integer(INTEGER_TYPE) :: I

          IsGhostElementForOutput = .false.
          if (.not.CalParams%ApplyGhostStabilisation) RETURN
          if (.not.allocated(IsGhostDof)) RETURN
          if (NGhostFaces == 0) RETURN

          do I = 1, ELEMENTNODES
            if (IsGhostDof(ReducedDof(ElementConnectivities(I, IEl)) + 1)) then
              IsGhostElementForOutput = .true.
              RETURN
            end if
          end do

        end function IsGhostElementForOutput

      end module ModGhostStabilisation
