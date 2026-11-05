    !*****************************************************************************
    !                                       ____  _____
    !           /\                         |___ \|  __ \
    !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
    !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
    !        / ____ \| | | | |_| | | | (_| |___) | |__| |
    !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/
    !
    !
    !   Anura3D - Numerical modelling and simulation of large deformations
    !   and soil-water-structure interaction using the material point method (MPM)
    !
    !   Copyright (C) 2025  Members of the Anura3D MPM Research Community
    !   (See Contributors file "Contributors.txt")
    !
    !   This program is free software: you can redistribute it and/or modify
    !   it under the terms of the GNU Lesser General Public License as published by
    !   the Free Software Foundation, either version 3 of the License, or
    !   (at your option) any later version.
    !
    !   This program is distributed in the hope that it will be useful,
    !   but WITHOUT ANY WARRANTY; without even the implied warranty of
    !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !   GNU Lesser General Public License for more details.
    !
    !   You should have received a copy of the GNU Lesser General Public License
    !   along with this program.  If not, see <https://www.gnu.org/licenses/>.
    !
    !*****************************************************************************


      module ModGasification
      !**********************************************************************
      !
      !    Function:  Gasification of stress-free, loosened material points
      !               (see docs/Gasification_Anura3D.md).
      !
      !               Each soil material point (MP) is either CONTINUUM or GAS.
      !               A gas MP keeps its mass, momentum and volume and still maps
      !               them to the grid (so total momentum is conserved and the
      !               Bardenhagen contact still sees it), but it carries no stress,
      !               the constitutive model is not called, it adds no internal
      !               force and no wave-speed time-step limit. A kinematic limit
      !               dt <= h / |v| replaces the wave-speed limit.
      !
      !               Continuum -> gas, when for NHold consecutive steps:
      !                 element not touching the rigid body or the contact
      !                 interface, and either
      !                 (a) stress-free (|p|, q <= tol*pRef, or p > 0 tension),
      !                     loosened (rho/rho0 <= chiOff) or sparse (fill <= fillOff),
      !                     and not compressing (dEps_v >= 0)
      !                 (b) isolated: the element and all its edge neighbours are
      !                     sparse (fill <= fillOff). A fragment smaller than an
      !                     element is traction free on all sides, so its mean
      !                     stress is physically zero whatever the MP stress says.
      !                 (c) loose: loosened (rho/rho0 <= chiOff) and in a sparse
      !                     element (fill <= fillOff), whatever the MP stress says.
      !               The element fill is measured by material content,
      !               fill = sum V0_p / |K| (V0 = reference volume), because the
      !               current volume of loose flying MPs inflates many times.
      !
      !               Gas -> continuum, when for NHold consecutive steps:
      !                 rho/rho0 >= chiOn, or
      !                 (fill >= fillOn, element holds continuum MPs, dEps_v < 0)
      !               The MP restarts stress-free with the current volume as its
      !               new reference, and its strain increment passed to the
      !               constitutive model is ramped over NRamp steps.
      !
      !               Scope: explicit dynamic, NFORMULATION = 1, single phase.
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
      use ModString
      use ModFeedback
      use ModFileIO

      implicit none

        integer(INTEGER_TYPE), parameter :: GAS_STATE_CONTINUUM = 0
        integer(INTEGER_TYPE), parameter :: GAS_STATE_GAS = 1

        integer(INTEGER_TYPE), dimension(:), allocatable, private :: GasState   ! (NParticles) GAS_STATE_*
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: GasHold    ! (NParticles) consecutive steps the switch condition held
        integer(INTEGER_TYPE), dimension(:), allocatable, private :: GasRamp    ! (NParticles) steps since reactivation (0 = no ramp)
        real(REAL_TYPE), dimension(:), allocatable, private :: GasV0            ! (NParticles) reference volume for rho/rho0
        real(REAL_TYPE), dimension(:,:), allocatable, private :: GasStatev0     ! (NParticles, NSTATEVAR) initial ESM state variables

        real(REAL_TYPE), dimension(:), allocatable, private :: GasElemFill      ! (NEl) soil fill of the element (gas + continuum), sum V0 / |K|
        logical, dimension(:), allocatable, private :: GasElemHasContinuum      ! (NEl) element holds continuum soil MPs
        logical, dimension(:), allocatable, private :: GasElemBlocked           ! (NEl) element touches rigid body / contact

        logical, private :: GasChecked = .false.
        logical, private :: GasIsAvailable = .false.
        logical, private :: GasStateLoaded = .false.  ! state read from the BRF file of the previous stage
        integer(INTEGER_TYPE), private :: NGasifiedSinceReport = 0
        integer(INTEGER_TYPE), private :: NReactivatedSinceReport = 0
        ! candidates of the last UpdateGasState call (tuning diagnostics)
        integer(INTEGER_TYPE), private :: NCandIsolated = 0   ! continuum MPs in isolated elements
        integer(INTEGER_TYPE), private :: NCandStressFree = 0 ! continuum MPs stress-free and loosened or sparse
        integer(INTEGER_TYPE), private :: NCandExpanding = 0  ! continuum MPs meeting a gasify condition this step
        integer(INTEGER_TYPE), private :: NCandLoose = 0      ! continuum MPs loosened and in a sparse element

      contains


        logical function IsGasificationActive()
        !**********************************************************************
        !    Function:  True if gasification is requested and supported
        !**********************************************************************
        implicit none

          IsGasificationActive = .false.
          if (.not.CalParams%ApplyGasification) RETURN

          if (.not.GasChecked) then
            GasChecked = .true.
            GasIsAvailable = (NFORMULATION == 1) .and. (CalParams%NumberOfPhases == 1)
            if (.not.GasIsAvailable) then
              call WriteInLogFile('Gasification requested but only implemented for the single-phase, ' // &
                                  '1-layer formulation. It is switched off.')
            else
              call WriteInLogFile('Gasification active: stress tol = ' // trim(String(CalParams%GasStressTol)) // &
                                  ', pRef = ' // trim(String(CalParams%GasPRef)) // &
                                  ', chi off/on = ' // trim(String(CalParams%GasChiOff)) // ' / ' // trim(String(CalParams%GasChiOn)) // &
                                  ', fill off/on = ' // trim(String(CalParams%GasFillOff)) // ' / ' // trim(String(CalParams%GasFillOn)) // &
                                  ', hold steps = ' // trim(String(CalParams%GasNHold)) // &
                                  ', ramp steps = ' // trim(String(CalParams%GasNRamp)) // &
                                  ', gas FLIP factor = ' // trim(String(CalParams%GasFLIPFactor)))
            end if
          end if

          IsGasificationActive = GasIsAvailable

        end function IsGasificationActive


        subroutine AllocateGasArrays()
        !**********************************************************************
        !    Function:  Allocate the per-MP state, keep it if the sizes match
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IError

          if (allocated(GasState)) then
            if (size(GasState) == Counters%NParticles) RETURN
            call WriteInLogFile('Gasification: number of material points changed, gas state is reset.')
            deallocate(GasState, GasHold, GasRamp, GasV0, stat = IError)
            if (allocated(GasStatev0)) deallocate(GasStatev0, stat = IError)
            GasStateLoaded = .false.
          end if

          allocate(GasState(Counters%NParticles), stat = IError)
          allocate(GasHold(Counters%NParticles), stat = IError)
          allocate(GasRamp(Counters%NParticles), stat = IError)
          allocate(GasV0(Counters%NParticles), stat = IError)
          GasState = GAS_STATE_CONTINUUM
          GasHold = 0
          GasRamp = 0
          GasV0 = -1.0 ! not set

        end subroutine AllocateGasArrays


        subroutine InitialiseGas()
        !**********************************************************************
        !    Function:  Make sure the gas state exists. MPs without a reference
        !               volume take their current volume as reference.
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: I, IError

          call AllocateGasArrays()

          do I = 1, Counters%NParticles
            if (GasV0(I) <= 0.0) GasV0(I) = Particles(I)%IntegrationWeight
          end do

          if (.not.allocated(GasStatev0)) then
            if (allocated(ESMstatevArray)) then
              if (size(ESMstatevArray, 1) == Counters%NParticles .and. size(ESMstatevArray, 2) > 0) then
                allocate(GasStatev0(Counters%NParticles, size(ESMstatevArray, 2)), stat = IError)
                GasStatev0 = ESMstatevArray
              end if
            end if
          end if

        end subroutine InitialiseGas


        logical function IsGasParticle(ParticleIndex)
        !**********************************************************************
        !    Function:  True if the MP is in the gas state
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex

          IsGasParticle = .false.
          if (.not.CalParams%ApplyGasification) RETURN
          if (.not.allocated(GasState)) RETURN
          if (ParticleIndex < 1 .or. ParticleIndex > size(GasState)) RETURN
          IsGasParticle = (GasState(ParticleIndex) == GAS_STATE_GAS)

        end function IsGasParticle


        real(REAL_TYPE) function GasRampFactor(ParticleIndex)
        !**********************************************************************
        !    Function:  Scaling of the strain increment passed to the
        !               constitutive model after reactivation (1 = no ramp)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex

          GasRampFactor = 1.0
          if (.not.CalParams%ApplyGasification) RETURN
          if (.not.allocated(GasRamp)) RETURN
          if (ParticleIndex < 1 .or. ParticleIndex > size(GasRamp)) RETURN
          if (GasRamp(ParticleIndex) > 0) then
            GasRampFactor = min(1.0d0, real(GasRamp(ParticleIndex), REAL_TYPE) / real(max(1, CalParams%GasNRamp), REAL_TYPE))
          end if

        end function GasRampFactor


        integer(INTEGER_TYPE) function GasStateForOutput(ParticleIndex)
        !**********************************************************************
        !    Function:  0 = continuum, 1 = gas, 2 = continuum, stiffness ramp
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex

          GasStateForOutput = 0
          if (.not.allocated(GasState)) RETURN
          if (ParticleIndex < 1 .or. ParticleIndex > size(GasState)) RETURN
          if (GasState(ParticleIndex) == GAS_STATE_GAS) then
            GasStateForOutput = 1
          else if (GasRamp(ParticleIndex) > 0) then
            GasStateForOutput = 2
          end if

        end function GasStateForOutput


        logical function IsGasEligibleParticle(ParticleIndex)
        !**********************************************************************
        !    Function:  Soil MP whose material can gasify (not a rigid body,
        !               not another contact entity, solid or mixture MP)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex
          integer(INTEGER_TYPE) :: MaterialIndex

          IsGasEligibleParticle = .false.
          if (CalParams%ApplyContactAlgorithm) then
            if (EntityIDArray(ParticleIndex) /= SOFT_ENTITY) RETURN
          end if
          MaterialIndex = MaterialIDArray(ParticleIndex)
          if (MatParams(MaterialIndex)%MaterialModel == ESM_RIGID_BODY) RETURN
          if (.not.((MaterialPointTypeArray(ParticleIndex) == MaterialPointTypeMixture) .or. &
                    (MaterialPointTypeArray(ParticleIndex) == MaterialPointTypeSolid))) RETURN
          IsGasEligibleParticle = .true.

        end function IsGasEligibleParticle


        logical function IsGasContactNode(NodeID)
        !**********************************************************************
        !    Function:  True if the node takes part in the Bardenhagen contact
        !               or carries mass of an entity other than the soil
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: NodeID
          integer(INTEGER_TYPE) :: IEntity, IDof

          IsGasContactNode = .false.
          if (.not.CalParams%ApplyContactAlgorithm) RETURN

          if (allocated(InterfaceNodes)) then
            if (InterfaceNodes(NodeID)) then
              IsGasContactNode = .true.
              RETURN
            end if
          end if

          if (.not.allocated(LumpedMassDry)) RETURN
          IDof = ReducedDof(NodeID) + 1
          do IEntity = 1, Counters%nEntity
            if (IEntity == SOFT_ENTITY) CYCLE
            if (LumpedMassDry(IDof, IEntity) /= 0.0) then
              IsGasContactNode = .true.
              RETURN
            end if
          end do

        end function IsGasContactNode


        subroutine ClassifyGasElements()
        !**********************************************************************
        !    Function:  Per active element: soil fill (gas + continuum), whether
        !               it holds continuum soil MPs, whether it is blocked
        !               (rigid body / other entity MPs or contact nodes)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IAEl, IEl, IPart, ParticleIndex, INode, IError
          real(REAL_TYPE) :: SumVolume, Weight, Radius

          if (allocated(GasElemFill)) then
            if (size(GasElemFill) /= Counters%NEl) then
              deallocate(GasElemFill, GasElemHasContinuum, GasElemBlocked, stat = IError)
            end if
          end if
          if (.not.allocated(GasElemFill)) then
            allocate(GasElemFill(Counters%NEl), stat = IError)
            allocate(GasElemHasContinuum(Counters%NEl), stat = IError)
            allocate(GasElemBlocked(Counters%NEl), stat = IError)
          end if

          GasElemFill = 0.0
          GasElemHasContinuum = .false.
          GasElemBlocked = .false.

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            SumVolume = 0.0

            do IPart = 1, NPartEle(IEl)
              ParticleIndex = GetParticleIndex(IPart, IEl)

              if (.not.IsGasEligibleParticle(ParticleIndex)) then
                GasElemBlocked(IEl) = .true.
                CYCLE
              end if

              if (GasState(ParticleIndex) /= GAS_STATE_GAS) GasElemHasContinuum(IEl) = .true.

              ! fill by material content: reference volume V0 (= m / rho0), not the current
              ! volume, which inflates for loose, flying MPs and would hide sparse elements
              Weight = GasV0(ParticleIndex)
              if (Weight <= 0.0) Weight = Particles(ParticleIndex)%IntegrationWeight
              if (Weight <= 0.0) CYCLE
              if (ISAXISYMMETRIC) then ! the axisymmetric weight carries the radius
                Radius = GlobPosArray(ParticleIndex, 1)
                if (Radius > 0.0) SumVolume = SumVolume + Weight / Radius
              else
                SumVolume = SumVolume + Weight
              end if
            end do

            if (ElementSpace(IEl) > 0.0) GasElemFill(IEl) = SumVolume / ElementSpace(IEl)

            if (.not.GasElemBlocked(IEl)) then
              do INode = 1, ELEMENTNODES
                if (IsGasContactNode(ElementConnectivities(INode, IEl))) then
                  GasElemBlocked(IEl) = .true.
                  EXIT
                end if
              end do
            end if
          end do

        end subroutine ClassifyGasElements


        logical function IsIsolatedGasElement(IEl)
        !**********************************************************************
        !    Function:  The element and all its edge neighbours are sparse
        !               (inactive or fill <= fillOff)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl
          integer(INTEGER_TYPE) :: ISide, JEl

          IsIsolatedGasElement = .false.
          if (GasElemFill(IEl) > CalParams%GasFillOff) RETURN
          if (.not.allocated(ElementAdjacencies)) RETURN

          do ISide = 1, ELEMENTSIDES
            JEl = ElementAdjacencies(IEl, ISide)
            if (JEl <= 0) CYCLE
            if (.not.IsActiveElement(JEl)) CYCLE
            if (GasElemBlocked(JEl)) RETURN ! never next to the rigid body
            if (GasElemFill(JEl) > CalParams%GasFillOff) RETURN
          end do

          IsIsolatedGasElement = .true.

        end function IsIsolatedGasElement


        subroutine GetMeanAndDeviatoricStress(ParticleIndex, P, Q)
        !**********************************************************************
        !    Function:  p = tr(sigma')/3 (tension positive), q = sqrt(3/2 s:s)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex
          real(REAL_TYPE), intent(out) :: P, Q
          real(REAL_TYPE) :: SS
          integer(INTEGER_TYPE) :: I

          P = (SigmaEffArray(ParticleIndex, 1) + SigmaEffArray(ParticleIndex, 2) + SigmaEffArray(ParticleIndex, 3)) / 3.0
          SS = 0.0
          do I = 1, 3
            SS = SS + (SigmaEffArray(ParticleIndex, I) - P)**2
          end do
          do I = 4, NTENSOR
            SS = SS + 2.0 * SigmaEffArray(ParticleIndex, I)**2
          end do
          Q = sqrt(1.5 * SS)

        end subroutine GetMeanAndDeviatoricStress


        subroutine GasifyParticle(ParticleIndex)
        !**********************************************************************
        !    Function:  Continuum -> gas. Mass, velocity, position and volume
        !               are kept (momentum unchanged), the stress is removed.
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex

          GasState(ParticleIndex) = GAS_STATE_GAS
          GasHold(ParticleIndex) = 0
          GasRamp(ParticleIndex) = 0
          SigmaEffArray(ParticleIndex, :) = 0.0
          Particles(ParticleIndex)%DBulkViscousPressure = 0.0
          Particles(ParticleIndex)%IsBoundaryParticle = .true.
          NGasifiedSinceReport = NGasifiedSinceReport + 1

        end subroutine GasifyParticle


        subroutine ReactivateParticle(ParticleIndex)
        !**********************************************************************
        !    Function:  Gas -> continuum: stress-free fresh deposit at the
        !               current density, momentum unchanged, stiffness ramp
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ParticleIndex

          GasState(ParticleIndex) = GAS_STATE_CONTINUUM
          GasHold(ParticleIndex) = 0
          GasRamp(ParticleIndex) = 1
          SigmaEffArray(ParticleIndex, :) = 0.0
          Particles(ParticleIndex)%DBulkViscousPressure = 0.0
          Particles(ParticleIndex)%IsBoundaryParticle = .true.
          GasV0(ParticleIndex) = Particles(ParticleIndex)%IntegrationWeight ! new reference density

          ! external models: restart from the initial state variables
          if (allocated(GasStatev0) .and. allocated(ESMstatevArray)) then
            if (size(GasStatev0, 2) == size(ESMstatevArray, 2)) then
              ESMstatevArray(ParticleIndex, :) = GasStatev0(ParticleIndex, :)
            end if
          end if

          NReactivatedSinceReport = NReactivatedSinceReport + 1

        end subroutine ReactivateParticle


        subroutine UpdateGasState()
        !**********************************************************************
        !    Function:  Evaluate the switching conditions once per time step,
        !               after the stress and volume update
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IAEl, IEl, IPart, ParticleIndex, I
          real(REAL_TYPE) :: EpsVol, Chi, P, Q, Fill, StressLimit
          logical :: StressFree, Condition, Isolated, Loose

          if (.not.IsGasificationActive()) RETURN

          call InitialiseGas()
          call ClassifyGasElements()

          StressLimit = CalParams%GasStressTol * CalParams%GasPRef
          NCandIsolated = 0
          NCandStressFree = 0
          NCandExpanding = 0
          NCandLoose = 0

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            Fill = GasElemFill(IEl)

            do IPart = 1, NPartEle(IEl)
              ParticleIndex = GetParticleIndex(IPart, IEl)
              if (.not.IsGasEligibleParticle(ParticleIndex)) CYCLE

              EpsVol = 0.0
              do I = 1, NVECTOR ! hoop strain included in axisymmetry
                EpsVol = EpsVol + GetEpsStepI(Particles(ParticleIndex), I)
              end do
              Chi = 1.0
              if (Particles(ParticleIndex)%IntegrationWeight > 0.0) then
                Chi = GasV0(ParticleIndex) / Particles(ParticleIndex)%IntegrationWeight
              end if

              if (GasState(ParticleIndex) == GAS_STATE_GAS) then
                ! ---- gas -> continuum
                Condition = (Chi >= CalParams%GasChiOn) .or. &
                            ((Fill >= CalParams%GasFillOn) .and. GasElemHasContinuum(IEl) .and. (EpsVol < 0.0))
                if (Condition) then
                  GasHold(ParticleIndex) = GasHold(ParticleIndex) + 1
                else
                  GasHold(ParticleIndex) = 0
                end if
                if (GasHold(ParticleIndex) >= CalParams%GasNHold) call ReactivateParticle(ParticleIndex)

              else
                ! ---- continuum: advance the stiffness ramp
                if (GasRamp(ParticleIndex) > 0) then
                  GasRamp(ParticleIndex) = GasRamp(ParticleIndex) + 1
                  if (GasRamp(ParticleIndex) > CalParams%GasNRamp) GasRamp(ParticleIndex) = 0
                end if

                ! ---- continuum -> gas
                Condition = .false.
                if (.not.GasElemBlocked(IEl)) then
                  call GetMeanAndDeviatoricStress(ParticleIndex, P, Q)
                  StressFree = ((abs(P) <= StressLimit) .or. (P > 0.0)) .and. (Q <= StressLimit)
                  StressFree = StressFree .and. ((Chi <= CalParams%GasChiOff) .or. (Fill <= CalParams%GasFillOff))
                  Isolated = IsIsolatedGasElement(IEl)
                  if (Isolated) NCandIsolated = NCandIsolated + 1
                  if (StressFree) NCandStressFree = NCandStressFree + 1
                  ! an isolated fragment is traction free: its MP stress and strain rate are
                  ! spurious, so only the stress-free branch requires that it is not compressed
                  ! loosened and sparse at the same time: cannot carry confinement either
                  Loose = (Chi <= CalParams%GasChiOff) .and. (Fill <= CalParams%GasFillOff)
                  if (Loose) NCandLoose = NCandLoose + 1
                  Condition = (StressFree .and. (EpsVol >= 0.0)) .or. Isolated .or. Loose
                end if
                if (Condition) NCandExpanding = NCandExpanding + 1
                if (Condition) then
                  GasHold(ParticleIndex) = GasHold(ParticleIndex) + 1
                else
                  GasHold(ParticleIndex) = 0
                end if
                if (GasHold(ParticleIndex) >= CalParams%GasNHold) call GasifyParticle(ParticleIndex)
              end if

            end do
          end do

        end subroutine UpdateGasState


        subroutine ForceParticleIntegrationGasElements()
        !**********************************************************************
        !    Function:  Elements holding gas MPs are integrated with the MPs,
        !               so a Gauss point never spreads stress onto gas mass.
        !               Call after CheckFillingOfElements.
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IAEl, IEl, IPart

          if (.not.IsGasificationActive()) RETURN
          if (.not.allocated(GasState)) RETURN
          if (.not.allocated(IsParticleIntegration)) RETURN

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            if (IsParticleIntegration(IEl)) CYCLE
            do IPart = 1, NPartEle(IEl)
              if (IsGasParticle(GetParticleIndex(IPart, IEl))) then
                IsParticleIntegration(IEl) = .true.
                EXIT
              end if
            end do
          end do

        end subroutine ForceParticleIntegrationGasElements


        logical function IsGasOnlyElement(IEl)
        !**********************************************************************
        !    Function:  Active element whose MPs are all gas
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IEl
          integer(INTEGER_TYPE) :: IPart

          IsGasOnlyElement = .false.
          if (.not.CalParams%ApplyGasification) RETURN
          if (.not.allocated(GasState)) RETURN
          if (NPartEle(IEl) <= 0) RETURN
          do IPart = 1, NPartEle(IEl)
            if (.not.IsGasParticle(GetParticleIndex(IPart, IEl))) RETURN
          end do
          IsGasOnlyElement = .true.

        end function IsGasOnlyElement


        real(REAL_TYPE) function GasCriticalTimeStep()
        !**********************************************************************
        !    Function:  Kinematic limit for gas MPs, dt <= h / |v|, so a gas
        !               MP does not cross more than one element per step
        !               (the Courant number is applied by the caller)
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: IAEl, IEl, IPart, ParticleIndex
          real(REAL_TYPE) :: Speed

          GasCriticalTimeStep = huge(GasCriticalTimeStep)
          if (.not.IsGasificationActive()) RETURN
          if (.not.allocated(GasState)) RETURN

          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            do IPart = 1, NPartEle(IEl)
              ParticleIndex = GetParticleIndex(IPart, IEl)
              if (.not.IsGasParticle(ParticleIndex)) CYCLE
              Speed = sqrt(sum(VelocityArray(ParticleIndex, 1:NVECTOR)**2))
              if (Speed > 1.0d-12) then
                GasCriticalTimeStep = min(GasCriticalTimeStep, ElementLMin(IEl) / Speed)
              end if
            end do
          end do

        end function GasCriticalTimeStep


        subroutine WriteGasDiagnostics()
        !**********************************************************************
        !    Function:  Log line: gas MPs, switches since the last report,
        !               gas mass fraction, gas-only elements, kinematic dt
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE) :: I, IAEl, NGas, NRamp, NGasOnly
          real(REAL_TYPE) :: GasMass, TotalMass

          if (.not.IsGasificationActive()) RETURN
          if (.not.allocated(GasState)) RETURN

          NGas = 0
          NRamp = 0
          GasMass = 0.0
          TotalMass = 0.0
          do I = 1, Counters%NParticles
            if (.not.IsGasEligibleParticle(I)) CYCLE
            TotalMass = TotalMass + MassArray(I)
            if (GasState(I) == GAS_STATE_GAS) then
              NGas = NGas + 1
              GasMass = GasMass + MassArray(I)
            else if (GasRamp(I) > 0) then
              NRamp = NRamp + 1
            end if
          end do

          NGasOnly = 0
          do IAEl = 1, Counters%NAEl
            if (IsGasOnlyElement(ActiveElement(IAEl))) NGasOnly = NGasOnly + 1
          end do

          if (TotalMass > 0.0) GasMass = GasMass / TotalMass

          call WriteInLogFile('Gas: points = ' // trim(String(NGas)) // &
                              ', gasified = ' // trim(String(NGasifiedSinceReport)) // &
                              ', reactivated = ' // trim(String(NReactivatedSinceReport)) // &
                              ', ramping = ' // trim(String(NRamp)) // &
                              ', candidates isolated/loose/stress-free/all = ' // trim(String(NCandIsolated)) // ' / ' // &
                              trim(String(NCandLoose)) // ' / ' // trim(String(NCandStressFree)) // ' / ' // &
                              trim(String(NCandExpanding)) // &
                              ', mass fraction = ' // trim(String(GasMass)) // &
                              ', gas-only elements = ' // trim(String(NGasOnly)) // &
                              ', kinematic dt = ' // trim(String(GasCriticalTimeStep())))

          NGasifiedSinceReport = 0
          NReactivatedSinceReport = 0

        end subroutine WriteGasDiagnostics


        subroutine WriteBlockGASSTATE(FileUnit)
        !**********************************************************************
        !    Function:  Write the gas state to the binary result file (BRF),
        !               4 reals per MP: state, hold counter, ramp counter, V0
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE) :: NByts, I
          real(REAL_TYPE), dimension(:), allocatable :: Buffer

          if (.not.CalParams%ApplyGasification) RETURN
          if (.not.allocated(GasState)) RETURN
          if (size(GasState) /= Counters%NParticles) RETURN

          allocate(Buffer(4 * Counters%NParticles))
          do I = 1, Counters%NParticles
            Buffer(4*I - 3) = real(GasState(I), REAL_TYPE)
            Buffer(4*I - 2) = real(GasHold(I), REAL_TYPE)
            Buffer(4*I - 1) = real(GasRamp(I), REAL_TYPE)
            Buffer(4*I)     = GasV0(I)
          end do

          NByts = 4 * Counters%NParticles * 8
          write(FileUnit) '$$GASSTATE$$', NByts
          write(FileUnit) (Buffer(I), I = 1, 4 * Counters%NParticles)

          deallocate(Buffer)

        end subroutine WriteBlockGASSTATE


        subroutine ReadBlockGASSTATE(FileUnit, NBytes)
        !**********************************************************************
        !    Function:  Read the gas state of the previous stage from the BRF
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: FileUnit
          integer(INTEGER_TYPE), intent(in) :: NBytes
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE), dimension(:), allocatable :: Buffer

          if (NBytes /= 4 * Counters%NParticles * 8) then
            call Skip(FileUnit, NBytes)
            RETURN
          end if

          allocate(Buffer(4 * Counters%NParticles))
          read(FileUnit) (Buffer(I), I = 1, 4 * Counters%NParticles)

          call AllocateGasArrays()
          do I = 1, Counters%NParticles
            GasState(I) = nint(Buffer(4*I - 3))
            GasHold(I)  = nint(Buffer(4*I - 2))
            GasRamp(I)  = nint(Buffer(4*I - 1))
            GasV0(I)    = Buffer(4*I)
          end do
          GasStateLoaded = .true.

          deallocate(Buffer)

        end subroutine ReadBlockGASSTATE


      end module ModGasification
