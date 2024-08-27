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
    !	Copyright (C) 2023  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModQuasiStaticImplicit
      !**********************************************************************
      !
      !    Function:  This module contains all routines specific for the quasi-static
      !               MPM using implicit time integration.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModReadCalculationData
      use ModReadMaterialData
      use ModWriteTestData
      use ModMPMData
      use ModMeshInfo
      use ModRotBoundCond
      use ModMPMMeshAdjustment
      use ModReadMPMData
      use ModElementEvaluation
      use ModDYNConvectivePhase
      use ModMPMDynContact
      use ModCounters
      use ModEmptyElements
      use ModParticle
      use ModWriteResultData
      use ModSolver
      use ModGlobalConstants
      
      implicit none

      contains ! Routines of this module

      subroutine InitialiseQuasiStaticImplicit()
      !**********************************************************************
      !
      !  Function:  initialise quasistatic implicit data
      !             
      !
      !**********************************************************************
      implicit none

      if (.not.CalParams%ApplyImplicitQuasiStatic) RETURN

      call InitialiseReducedSolution()

      call InitialiseMultipliers()
      call GetNodalExtForces() ! Rotated to local c.s

      end subroutine InitialiseQuasiStaticImplicit

      subroutine DestroyQuasiStaticImplicit()
      !**********************************************************************
      !
      !  Function:  destroy quasistatic implicit data
      !             
      !
      !**********************************************************************
      implicit none

      if (.not.CalParams%ApplyImplicitQuasiStatic) RETURN

      call DestroyEquations()

      end subroutine DestroyQuasiStaticImplicit

      subroutine RunImplicitQuasiStaticLoadStep()
      !**********************************************************************
      !
      !  Function:  Routine called from the main routine for performing a
      !             quasi-static load step using implicit time integration.
      !
      !**********************************************************************

      implicit none

      integer(INTEGER_TYPE) I, J
      real(REAL_TYPE) :: DT
      integer(INTEGER_TYPE), parameter :: TESTFILE3 = 31

      DT = CalParams%TotalTime / CalParams%MaxTimeSteps
         
      TotalDisplacementSoil = 0
      TotalPressure = 0 
      
      do I = 1, CalParams%MaxTimeSteps
        
         CalParams%TimeStep = I
         call GetTimeStep(I, CalParams%FileNames%TimeStepExt)
         call UpdateImplicitMultipliersForLoadStep()
         call QuasiStaticMPMLagrangianPhase(DT)
         
         call GiveMessage('time step ' // trim(String(I)))
         
         if (CalParams%OutputDebugData) then
          open(TESTFILE3,file='d:\tmp\pressure.txt')
          write(TESTFILE3,*) 'time', CalParams%OverallRealTime + I * DT
          do J = 1, Counters%NodTot
             write(TESTFILE3,*) J,TotalPressure(J)
          end do
          close(TESTFILE3);
         end if
         

         call QuasiStaticMPMConvectivePhase()
         call MaterialPointOutput()
         !call DefineVariableData()
         !call WriteFEMNodeData()
        
         if (CalParams%OutputNumberOfTimeSteps>0) then 
            do J = 1, CalParams%OutputNumberOfTimeSteps
              if (I==CalParams%OutputTimeStepID(J)) then
                call GiveMessage('VTK ' // trim(String(J)))
                call WriteVTKOutput()
              end if
            end do
          end if         
         
      end do
      CalParams%OverallRealTime = CalParams%OverallRealTime + CalParams%TotalTime
             
      end subroutine RunImplicitQuasiStaticLoadStep

      subroutine QuasiStaticMPMLagrangianPhase(DT)
      !**********************************************************************
      !
      !  Function:  Lagrangian phase of the implicit cycle. 
      !             This is a displacement-pressure (u-p) formulation. 
      !             
      !
      !**********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: DT
      
      ! Local variables
      integer(INTEGER_TYPE) :: I, J
      integer(INTEGER_TYPE), parameter :: TESTFILE3 = 31
      
      IncrementalDisplacementSoil = 0.0
      IncrementalPressure = 0.0
      SubIncrementalDisplacement = 0.0 
      SubIncrementalPressure = 0.0
      
      call GetNodalIntForces() 
      call GetNodalIntFlow(DT)
      if ((IsMPMComputation().or.IsULFEMComputation()).or.(CalParams%IStep==1)) then
        call UpdateGlobalStiffnessMatrix(DT) 
      end if

      if ( CalParams%IStep > 1) then
        call CheckForZLS() 
      end if
      call GetNodalExtForces() 
      
      RateOfMomentum = ExtLoad + GravityLoad - IntLoad
      RateOfMomentum(1:Counters%N, 1) = RateOfMomentum(1:Counters%N, 1) * PBoundary(1:Counters%N) 
      RateOfFlux = ExtFlow - IntFlow
      do J = 1, Counters%NodTot
         RateOfFlux(J) = RateOfFlux(J) * PBoundaryQuasiStatic(J*(NVECTOR+1))
      end do

      call SetZeroParticleIncStrain()

      CalParams%ConvergenceCheck%DoesConverge = .false.
      CalParams%ImplicitIntegration%Iteration = 0
      CalParams%ImplicitIntegration%MaximumJacobianRatio = 0.0
      do while(.not.CalParams%ConvergenceCheck%DoesConverge)

        CalParams%ImplicitIntegration%Iteration = CalParams%ImplicitIntegration%Iteration + 1

        call SolveEquations(RateOfMomentum, RateOfFlux, SubIncrementalDisplacement, SubIncrementalPressure) 

        if (IS3DCYLINDRIC) then
          call RotVec(IRotation, NRotNodes, IRotMat, ReducedDof, SubIncrementalDisplacement, SubIncrementalDisplacement)
        end if

        IncrementalDisplacementSoil(1:Counters%N, 1) = IncrementalDisplacementSoil(1:Counters%N, 1) + SubIncrementalDisplacement(1:Counters%N) 
        IncrementalPressure = IncrementalPressure + SubIncrementalPressure 
        
        call UpdateNodesImplicit()

        if (.not.CalParams%ApplyStrainSmoothing) then
          call UpdateParticleStrains()
        end if
        call UpdateParticlePressure()
        if (.not.IsMPMComputation()) then
           call UpdateParticleDisplacement()
        end if
        
        call MPMDYNGetSig()
        call ComputeParticleVelocityWater()
        call GetNodalIntForces()
        call GetNodalIntFlow(DT)
        RateofMomentum = ExtLoad + GravityLoad - IntLoad
        RateOfMomentum(1:Counters%N, 1) = RateOfMomentum(1:Counters%N, 1) * PBoundary(1:Counters%N)
        RateofFlux = ExtFlow - IntFlow
        do J = 1, Counters%NodTot
           RateOfFlux(J) = RateOfFlux(J) * PBoundaryQuasiStatic(J*(NVECTOR+1))
        end do
        
        CalParams%ConvergenceCheck%DoesConverge = QuasiStaticConvergenceCheck()

      end do

      call Scaling()

      TotalDisplacementSoil(1:Counters%N) = TotalDisplacementSoil(1:Counters%N) + IncrementalDisplacementSoil(1:Counters%N, 1) 
      TotalPressure = TotalPressure + IncrementalPressure 
      
      call IncreaseParticleTotalStrains()
      call IncreaseParticleTotalPressure() 

      if (IsULFEMComputation()) then
        do I = 1, Counters%NEl
          call CoordLocalToGlobal(I, NodalCoordinatesUpd)
        end do
      end if

      if (.not.IsMPMComputation()) then
        call SetInitialParticleStressForNextLoadStep()
      end if

      end subroutine QuasiStaticMPMLagrangianPhase

      subroutine QuasiStaticMPMConvectivePhase()
      !**********************************************************************
      !
      !  Function:  Mapping data from nodes to particles and from first particle
      !             in an element to all other particles in the elements           
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IError
      real(REAL_TYPE), dimension(:), allocatable :: NodalIncDisplacementMeshAdjust

      if (.not.IsMPMComputation()) RETURN

      if (CalParams%ImplicitIntegration%IsZeroLoadStep) then
        call SetInitialStressForNextLoadStep()
        RETURN
      end if

      if (IsMPMWithMixedIntegration()) then
        call AssignStressesToParticles()            
        call AssignStateParametersToParticles()
      end if

      call UpdateParticleWeights()

      if (CalParams%ApplyMeshSmoothing) then

        allocate(NodalIncDisplacementMeshAdjust(Counters%N), stat = IError)
        NodalIncDisplacementMeshAdjust = 0.0
        NodalIncDisplacementMeshAdjust(1:Counters%N) = IncrementalDisplacementSoil(1:Counters%N, 1)

        call UpdateNodesMeshAdjust(NodalIncDisplacementMeshAdjust)

        call MeshAdjustment()

        deallocate(NodalIncDisplacementMeshAdjust, stat = IError)

        call UpdateMeshAdjacencyInformation(NodalCoordinates)
        
      end if

      call UpdateParticlePos()
      call UpdateParticleHouseKeeping()

      NodalCoordinatesUpd = NodalCoordinates

      call SetActiveElement()
      call setParticleIndex()

      if (CalParams%ApplyEmptyElements.or.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix) then
        call CheckEmptyElements()
      end if
      if (CalParams%ApplyEmptyElements) then
        call AdjustParticleDiscretisation()
      end if

      call CheckFillingOfElements()

      if (IsMPMWithMixedIntegration()) then
        call StressAndPorePressureSmoothening()
        call StateParametersSmoothening()
      else
        call SetInitialStressForNextLoadStep()
      end if

      end subroutine QuasiStaticMPMConvectivePhase

      subroutine UpdateNodesImplicit()
      !**********************************************************************
      !
      !  Function:  This can be used with small strain FEM to update nodal coordinates after the lagrangian phase. 
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: I, J, IDoF

      if (IsSmallDeformationFEMComputation().or.CalParams%ImplicitIntegration%IsZeroLoadStep) RETURN

      do I = 1, Counters%NodTot
        IDoF = ReducedDoF(I)
        do J = 1, NVECTOR
          NodalCoordinatesUpd(I, J) = NodalCoordinates(I, J) + TotalDisplacementSoil(IDoF + J)
          if (IsULFEMComputation()) then
            NodalCoordinatesUpd(I, J) = NodalCoordinatesUpd(I, J) + IncrementalDisplacementSoil(IDoF + J, 1)
          end if
        end do
      end do

      end subroutine UpdateNodesImplicit

      logical function QuasiStaticConvergenceCheck()
      !**********************************************************************
      !
      !  Function:  Convergence check for the implicit computational cycle. 
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      real(REAL_TYPE) :: OOBNorm, LoadNorm, Error, LocalError
      integer(INTEGER_TYPE) :: NAllowedInaccuratePlasticPoints

      call Norm(Counters%NodTot, NDOFL, RateOfMomentum, ReducedDof, OOBNorm)
      call Norm(Counters%NodTot, NDOFL, ExtLoad + GravityLoad, ReducedDof, LoadNorm)

      Error = OOBNorm / LoadNorm

      if (CalParams%ConvergenceCheck%SumIntegrationPointWeights<1.0) then
        CalParams%ConvergenceCheck%SumIntegrationPointWeights = 1.0
      end if
      LocalError = CalParams%ConvergenceCheck%SumLocalError / CalParams%ConvergenceCheck%SumIntegrationPointWeights

      NAllowedInaccuratePlasticPoints = nint(N_PERCENT_TOLERATED_INACCURATE_PLASTIC_POINTS * CalParams%IntegrationPointData%NPlasticPoints + N_ADDITIONAL_TOLERATED_INACCURATE_PLASTIC_POINTS)

      QuasiStaticConvergenceCheck = .true.

      if (Error>CalParams%ToleratedErrorForce) QuasiStaticConvergenceCheck = .false.
      if (LocalError>CalParams%ToleratedErrorForce) QuasiStaticConvergenceCheck = .false.
      if (CalParams%ImplicitIntegration%Iteration<2) QuasiStaticConvergenceCheck = .false.
      if (CalParams%ConvergenceCheck%NInaccuratePlasticPoints>NAllowedInaccuratePlasticPoints) then
        QuasiStaticConvergenceCheck = .false.
        if ((Error<0.01 * CalParams%ToleratedErrorForce).and.(CalParams%ImplicitIntegration%Iteration>1)) then
          QuasiStaticConvergenceCheck = .true.
        end if
      end if

      if (CalParams%ImplicitIntegration%Iteration==CalParams%ImplicitIntegration%MaxIterations) then
        QuasiStaticConvergenceCheck = .true.
      end if

      end function QuasiStaticConvergenceCheck

      subroutine CheckForZLS()
      !**********************************************************************
      !
      !  Function:  Check for zero load step - If error is still higher than tolerance then load step repeated.  
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      real(REAL_TYPE) :: OOBNorm, LoadNorm, Error

      if ((.not.IsMPMComputation()).or.(.not.CalParams%ImplicitIntegration%DoUseZLS)) RETURN

      RateOfMomentum = ExtLoad + GravityLoad - IntLoad
      RateOfMomentum(1:Counters%N, 1) = RateOfMomentum(1:Counters%N, 1) * PBoundary(1:Counters%N)
      call Norm(Counters%NodTot, NDOFL, RateOfMomentum, ReducedDof, OOBNorm)
      call Norm(Counters%NodTot, NDOFL, ExtLoad + GravityLoad, ReducedDof, LoadNorm)
      Error = OOBNorm / LoadNorm

      if (CalParams%ImplicitIntegration%IsZeroLoadStep.or.(CalParams%IStep==1)) then
        CalParams%ImplicitIntegration%IsZeroLoadStep = .false.
      else
        CalParams%ImplicitIntegration%IsZeroLoadStep = (Error>CalParams%ToleratedErrorForce)
        if (CalParams%ImplicitIntegration%IsZeroLoadStep) then
          CalParams%IStep = CalParams%IStep - 1
          call GetStepExt(CalParams%IStep, CalParams%FileNames%LoadStepExt)
        end if
      end if

      end subroutine CheckForZLS

      subroutine InitialiseMultipliers()
      !**********************************************************************
      !
      !  Function:  Initialise multipliers used in the quasistatic implicit scheme.  
      !             
      !
      !**********************************************************************
      implicit none

      if ((CalParams%IStep==1).or.(.not.CalParams%ImplicitIntegration%DoUseAutomaticLoadStepping)) then
        CalParams%ImplicitIntegration%InitialMultiplier = dble(CalParams%IStep - CalParams%PreviouslyRealisedLoadStep) / dble(CalParams%NLoadSteps - CalParams%PreviouslyRealisedLoadStep)
        ! TODO: NLoadSteps only purpose is to determine the initial step size.
        !       Take initial load multiplier of 50 per cent and use load step resetting to determine
        !       reasonable initial multiplier instead.
        CalParams%Multipliers%SolidAIncrement = CalParams%ImplicitIntegration%InitialMultiplier * (CalParams%Multipliers%SolidAFinal - CalParams%Multipliers%SolidARealised)
        CalParams%Multipliers%WaterAIncrement = CalParams%ImplicitIntegration%InitialMultiplier * (CalParams%Multipliers%WaterAFinal - CalParams%Multipliers%WaterARealised)
        CalParams%Multipliers%GravityIncrement = CalParams%ImplicitIntegration%InitialMultiplier * (CalParams%Multipliers%GravityFinal - CalParams%Multipliers%GravityRealised)
      end if

      CalParams%Multipliers%SolidACurrent = CalParams%Multipliers%SolidARealised
      CalParams%Multipliers%WaterACurrent = CalParams%Multipliers%WaterARealised
      CalParams%Multipliers%GravityCurrent = CalParams%Multipliers%GravityRealised

      end subroutine InitialiseMultipliers

      subroutine UpdateImplicitMultipliersForLoadStep()
      !**********************************************************************
      !
      !  Function:  Update load and gravity multipliers used in the quasistatic implicit scheme.  
      !             
      !
      !**********************************************************************
      implicit none

      CalParams%Multipliers%SolidAIncrement(1) = CalParams%ImplicitIntegration%ScaleFactor * CalParams%Multipliers%SolidAIncrement(1)
      CalParams%Multipliers%WaterAIncrement(1) = CalParams%ImplicitIntegration%ScaleFactor * CalParams%Multipliers%WaterAIncrement(1)
      CalParams%Multipliers%GravityIncrement = CalParams%ImplicitIntegration%ScaleFactor * CalParams%Multipliers%GravityIncrement

      if (CalParams%Multipliers%SolidAIncrement(1)>(CalParams%Multipliers%SolidAFinal(1) - CalParams%Multipliers%SolidACurrent(1))) then
        CalParams%Multipliers%SolidAIncrement(1) = CalParams%Multipliers%SolidAFinal(1)- CalParams%Multipliers%SolidACurrent(1)
      end if
      if (CalParams%Multipliers%WaterAIncrement(1)>(CalParams%Multipliers%WaterAFinal(1) - CalParams%Multipliers%WaterACurrent(1))) then
        CalParams%Multipliers%WaterAIncrement(1) = CalParams%Multipliers%WaterAFinal(1) - CalParams%Multipliers%WaterACurrent(1)
      end if
      if (CalParams%Multipliers%GravityIncrement>(CalParams%Multipliers%GravityFinal - CalParams%Multipliers%GravityCurrent)) then
        CalParams%Multipliers%GravityIncrement = CalParams%Multipliers%GravityFinal - CalParams%Multipliers%GravityCurrent
      end if

      if (.not.CalParams%ImplicitIntegration%IsZeroLoadStep) then
        CalParams%Multipliers%SolidACurrent(1) = CalParams%Multipliers%SolidACurrent(1)+ CalParams%Multipliers%SolidAIncrement(1)
        CalParams%Multipliers%WaterACurrent(1) = CalParams%Multipliers%WaterACurrent(1) + CalParams%Multipliers%WaterAIncrement(1)
        CalParams%Multipliers%GravityCurrent = CalParams%Multipliers%GravityCurrent + CalParams%Multipliers%GravityIncrement
      end if

      end subroutine UpdateImplicitMultipliersForLoadStep

      subroutine Scaling()
      !**********************************************************************
      !
      !  Function: Computes the scale factor to adjust the multipliers 
      !            according to the rate of convergence. Scale is under 1.0 if
      !            many iterations are needed for convergence or  
      !            greater than 1.0 if few iterations are needed for convergence.
      !**********************************************************************
      implicit none

      ! Local variables
      logical :: IsSevereMeshDistortion

      CalParams%ImplicitIntegration%ScaleFactor = 1.0

      if (.not.CalParams%ImplicitIntegration%DoUseAutomaticLoadStepping) RETURN
      if (CalParams%ImplicitIntegration%IsZeroLoadStep) RETURN

      if (CalParams%ImplicitIntegration%DoCheckChangeJacobian) then
        IsSevereMeshDistortion = CheckSevereMeshDistortion()
      end if

      if ((CalParams%ImplicitIntegration%Iteration>CalParams%ImplicitIntegration%MaxDesiredIterations).or.IsSevereMeshDistortion) then
        CalParams%ImplicitIntegration%ScaleFactor = CalParams%ImplicitIntegration%DownScaleFactor
      else if ((CalParams%ImplicitIntegration%Iteration<CalParams%ImplicitIntegration%MinDesiredIterations).and.(.not.IsSevereMeshDistortion)) then
        CalParams%ImplicitIntegration%ScaleFactor = CalParams%ImplicitIntegration%UpScaleFactor
      end if

      end subroutine Scaling

      subroutine IncreaseParticleTotalStrains()
      !**********************************************************************
      !
      !  Function:  update material point strain. 
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: I
      real(REAL_TYPE), dimension(NTENSOR) :: Eps

      do I = 1, Counters%NParticles

        Eps = Particles(I)%EpsStep

        call IncreaseEps(Particles(I), Eps)
      end do

      end subroutine IncreaseParticleTotalStrains

      subroutine CheckBoundaryLoadedParticles()
      !**********************************************************************
      !
      !  Function:  Check Boundary Loaded Particles
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IElement, IAdjacentElement, AdjacentElement, IParticle, ParticleIndex, IAEl
      logical :: IsBoundaryElement

      do IAEl = 1, Counters%NAEl
        IElement = ActiveElement(IAEl)

        if (.not.IsParticleIntegration(IElement)) then
          IsBoundaryElement = .false.
          do IAdjacentElement = 1, GetNElmOfElm(IElement)
            AdjacentElement = GetElmIOfElm(IElement, IAdjacentElement)
            IsBoundaryElement = IsBoundaryElement.or.(.not.IsActiveElement(IAdjacentElement))
          end do

          if (.not.IsBoundaryElement) then
            do IParticle = 1, NPartEle(IElement)
              ParticleIndex = GetParticleIndex(IParticle, IElement)
              Particles(ParticleIndex)%IsBoundaryParticle = .false.

            end do
          end if

          ! TODO : Check for interface element?
          ! Check whether element surrounded by activated elements?

        end if
      end do

      end subroutine CheckBoundaryLoadedParticles

      subroutine UpdateParticleWeights()
      !**********************************************************************
      !
      !  Function:  computes new particle weight given DetJac
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJac
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: InvRJac
      integer(INTEGER_TYPE) :: IElement, IAEl, NElemPart, IParticle, ParticleIndex
      real(REAL_TYPE) :: DetJacDeformed, DetJacUndeformed

      do IAEl = 1, Counters%NAEl ! Loop over all active elements for computation of stresses
        IElement = ActiveElement(IAEl)

        NElemPart = NPartEle(IElement)
        do IParticle = 1, NElemPart ! Loop over all particles of the element
          ParticleIndex = GetParticleIndex(IParticle, IElement)
          call DetJacob(Particles(ParticleIndex)%LocPos, Counters%NEl, Counters%NodTot, NVECTOR, IElement, ElementConnectivities, NodalCoordinates, RJac, InvRJac, DetJacUndeformed)
          call DetJacob(Particles(ParticleIndex)%LocPos, Counters%NEl, Counters%NodTot, NVECTOR, IElement, ElementConnectivities, NodalCoordinatesUpd, RJac, InvRJac, DetJacDeformed)
          if (dabs(DetJacUndeformed)>1d-10) then
            Particles(ParticleIndex)%IntegrationWeight = Particles(ParticleIndex)%IntegrationWeight * DetJacDeformed / DetJacUndeformed
          else
            write(5, *) 'UpdateParticleWeights encountered zero ', 'Jacobian for undeformed mesh: ', DetJacUndeformed
          end if

        end do
      end do

      end subroutine UpdateParticleWeights

      logical function CheckSevereMeshDistortion()
      !**********************************************************************
      !
      !  Function:  Checks for severe mesh distortion through comparison of 
      !             DeformedDetJac and UndeformedDetJac
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IElement, IAEl
      real(REAL_TYPE) :: UndeformedDetJac, DeformedDetJac, Ratio, WeiGP
      real(REAL_TYPE), dimension(NVECTOR) :: PosGP

      CheckSevereMeshDistortion = .false.

      do IAEl = 1, Counters%NAEl
        IElement = ActiveElement(IAEl)
        call GaussPointLocalCoordinates(1, WeiGP, PosGP)
        call ScaledDetJacob(PosGP, IElement, NodalCoordinatesUpd, DeformedDetJac)
        call ScaledDetJacob(PosGP, IElement, NodalCoordinates, UndeformedDetJac)

        if (DeformedDetJac>UndeformedDetJac) then ! Volume increase
          Ratio = 1.0 - dabs(UndeformedDetJac / DeformedDetJac)
        else ! No change or volume decrease
          Ratio = 1.0 - dabs(DeformedDetJac / UndeformedDetJac)
        end if

        if (Ratio>CalParams%ImplicitIntegration%MaximumJacobianRatio) then
          CalParams%ImplicitIntegration%MaximumJacobianRatio = Ratio
        end if

        if (Ratio>CalParams%ImplicitIntegration%LimitJacobianChange) then
          CheckSevereMeshDistortion = .true.
          EXIT
        end if

      end do

      end function CheckSevereMeshDistortion

      subroutine ScaledDetJacob(LocPos, IElement, NodalCoordinates, DetJac)
      !**********************************************************************
      !
      !  Function:  Scale the DetJac when severe mesh distortion happens. 
      !             Active only when $$APPLY_JACOBIAN_CHECK is true.
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: LocPos
      integer(INTEGER_TYPE), intent(in) :: IElement
      real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: NodalCoordinates
      real(REAL_TYPE), intent(out) :: DetJac
      ! Local variables
      integer(INTEGER_TYPE) :: I, J, K
      real(REAL_TYPE), dimension(ELEMENTNODES) :: HS
      real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: DHS
      real(REAL_TYPE) :: Det1
      real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: ScaledGlobCoord
      real(REAL_TYPE) :: ScaleFactor
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJac
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: InvRJac
      real(REAL_TYPE), dimension(NVECTOR) :: Minimum, Maximum, Delta

      Minimum = 1E200
      Maximum = -1E200

      ScaledGlobCoord = 0.0
      do I = 1, ELEMENTNODES
        do J = 1, NVECTOR
          ScaledGlobCoord(I, J) = NodalCoordinates(ElementConnectivities(I, IElement), J)
        end do
        call CheckMinMax(ScaledGlobCoord(I, :), Minimum, Maximum)
      end do

      do I = 1, NVECTOR
        Delta(I) = Maximum(I) - Minimum(I)
      end do  
      ScaleFactor = maxval(Delta)

      do I = 1, ELEMENTNODES
        do J = 1, NVECTOR  
          ScaledGlobCoord(I, J) = (ScaledGlobCoord(I, J) - Minimum(J)) / ScaleFactor
        end do  
      end do

      call ShapeFunctionData(LocPos, ELEMENTNODES, HS, DHS)

      RJac = 0.0
      do K = 1, ELEMENTNODES
        do I = 1, NVECTOR
          do J = 1, NVECTOR
            RJac(I, J) = RJac(I, J) + DHS(K, I) * ScaledGlobCoord(K, J)
          end do
        end do
      end do

      call RJacInv(NVECTOR, RJac, InvRJac, DetJac, Det1)

      end subroutine ScaledDetJacob

      subroutine UpdateGlobalStiffnessMatrix(DT)
      !**********************************************************************
      !
      !  Function:  Initialises or update the Global Stiffnes Matrix
      !             
      !
      !**********************************************************************
      implicit none
      real(REAL_TYPE), intent(in) :: DT
      
      ! Local variables
      integer(INTEGER_TYPE) :: IActiveElement, ElementID
      real(REAL_TYPE), dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)) ::ElementStiffnessMatrix 
      real(REAL_TYPE), dimension(ELEMENTNODES * (NVECTOR + 1)) ::ElementStiffnessVector
      logical :: IsEnhancedStiffnessMatrix

      call ResetGlobalStiffnessMatrix()
      call ResetGlobalStiffnessRHS()
      ExtFlow = 0
      
      do IActiveElement = 1, Counters%NAEl
        ElementID = ActiveElement(IActiveElement)

        IsEnhancedStiffnessMatrix = .false.
        call ComputeElementStiffnessMatrix(ElementID, ElementStiffnessMatrix, IsEnhancedStiffnessMatrix, ElementStiffnessVector, DT)
        call AddElementStiffnessMatrix(ElementID, ElementStiffnessMatrix)
        call AddExternalFlow(ElementID, ElementStiffnessVector)
        if (CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix.and.IsParticleIntegration(ElementID)) then
          IsEnhancedStiffnessMatrix = .true.
          call ComputeElementStiffnessMatrix(ElementID, ElementStiffnessMatrix, IsEnhancedStiffnessMatrix, ElementStiffnessVector, DT)
          call AddElementStiffnessMatrix(ElementID, ElementStiffnessMatrix)
          call AddExternalFlow(ElementID, ElementStiffnessVector)
        end if
      ! rotation
      end do

      if (IsMPMComputation()) then
        call AdjustGlobalStiffnessMatrixForInactiveElements()
      end if

      end subroutine UpdateGlobalStiffnessMatrix

      subroutine ComputeElementStiffnessMatrix(ElementID, ElementStiffnessMatrix, IsEnhancedStiffnessMatrix, ElementStiffnessVector, DT) 
      !**********************************************************************
      !
      !  Function:  Computes Element Stiffness Matrix
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: ElementID
      real(REAL_TYPE), intent(in) ::  DT
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)),intent(out) :: ElementStiffnessMatrix
      real(REAL_TYPE), dimension(ELEMENTNODES * (NVECTOR + 1)), intent(out) :: ElementStiffnessVector
      logical :: IsEnhancedStiffnessMatrix
      
      ! Local variables
      integer(INTEGER_TYPE) :: Int, NElemPart, ParticleIndex, ISet, NNodes, DoFI, DoFJ, NodeI, NodeJ, NodeID
      real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
      real(REAL_TYPE) :: Weight, Det, XNu, BFac, SPhi, SPsi, PC, KM, KMG
      real(REAL_TYPE) :: GG, Cohec, Tens

      real(REAL_TYPE) :: Fac, D1, D2, D3
      real(REAL_TYPE), dimension(NVECTOR) :: BI, BJ
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: NodalStiffnessTerms
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: NodalStiffnessTermsK
      real(REAL_TYPE), dimension(NVECTOR, 1) :: NodalStiffnessTermsQ
      real(REAL_TYPE), dimension(NVECTOR, 1) :: NodalStiffnessTermsR
      real(REAL_TYPE), dimension(1, 1) :: NodalStiffnessTermsS
      real(REAL_TYPE), dimension(1, 1) :: NodalStiffnessTermsH

      integer(INTEGER_TYPE) i,j,n
      real(REAL_TYPE) :: NI, NJ, MuL, KappaL, NL, GammaL, KL
      integer(INTEGER_TYPE), parameter :: TESTFILE = 10

      ElementStiffnessMatrix = 0.0
      ElementStiffnessVector = 0.0

      call FormB3(1, ElementID, ElementConnectivities, NodalCoordinates, B, Det, Weight)

      if (CalParams%OutputDebugData) then
        open(TESTFILE,file='d:\tmp\testfileA.txt')
        if (ElementID==1) then 
          write(TESTFILE,*) 'ElementID',ElementID
          write(TESTFILE,*)
          write(TESTFILE,*) 'NodalCoordinates'
          do i = 1,ELEMENTNODES
            n = ElementConnectivities(i,ElementID)
            write(TESTFILE,901) n, (NodalCoordinates(n,j),j=1,NVECTOR)
901         format(i5,3(1x,f10.5))
          end do
        end if
        if (ElementID==1) then
          write(TESTFILE,*) 'NVECTOR',NVECTOR
          write(TESTFILE,*) 'ELEMENTNODES',ELEMENTNODES
          write(TESTFILE,*)
          write(TESTFILE,*) 'B'
          do i = 1, NVECTOR
            write(TESTFILE,904) i,(B(i,j),j=1,ELEMENTNODES)
904         format(i5,10(1x,f10.5))
          end do
          write(TESTFILE,*)
        end if
        close(TESTFILE)
      end if
      
      NElemPart = NumberOfIntegrationPoints(ElementID)
      
      do Int = 1, 1 !NElemPart mixed mpm integration over gauss points
        ParticleIndex = GetParticleIndex(Int, ElementID)
        if (IsParticleIntegration(ElementID).and.(.not.IsEnhancedStiffnessMatrix)) then
          if ( ISAXISYMMETRIC ) then
            Weight = Weight * GPGlobalPositionElement(1, Int, ElementID)
          else
            Weight = Particles(ParticleIndex)%IntegrationWeight
          end if
        end if

        call GetMaterialData(ParticleIndex, ISet, &
                             XNu, BFac, SPhi,                &
                             SPsi, GG, Cohec, Tens)       
        
        if (CalParams%ApplyEffectiveStressAnalysis) then
          XNu = MatParams(MaterialIDArray(ParticleIndex))%UndrainedPoissonRatio
        end if
        
        MuL = MatParams(ISet)%ViscosityLiquid ! add to particle properties
        NL = MatParams(ISet)%InitialPorosity 
        KappaL = MatParams(ISet)%IntrinsicPermeabilityLiquid
        KL = MatParams(ISet)%BulkModulusLiquid
        GammaL = MatParams(ISet)%WeightLiquid
        Fac = 2.0 * GG / (1.0 - 2.0 * XNu)
        D1 = Fac * (1.0 - XNu)
        D2 = Fac * XNu
        D3 = D1 * (XNu/(1.0 - XNu))

        PC = NL / KL
        KM = KappaL / MuL * DT 
        KMG = KappaL * GammaL / MuL * DT
        
        do NodeI = 1, ELEMENTNODES
          DoFI = (NVECTOR + 1) * (NodeI - 1)
          BI = B(:, NodeI) * Weight
          NI = GPShapeFunction(Int, NodeI) * Weight

          if (IsSymmetricStiffnessMatrix()) then
            NNodes = NodeI
          else
            NNodes = ELEMENTNODES
          end if

          do NodeJ = 1, NNodes
            DoFJ = (NVECTOR + 1) * (NodeJ - 1)
            BJ = B(:, NodeJ)
            NJ = GPShapeFunction(Int, NodeJ)
            NodeID = ElementConnectivities(NodeJ, ElementID)
            
            call ComputeNodalStiffnessTermsK(BI, BJ, D1, D2, D3, GG, NodalStiffnessTermsK)
            call ComputeNodalStiffnessTermsQ(BI, NJ, NodalStiffnessTermsQ)
            call ComputeNodalStiffnessTermsR(NI, BJ, NodalStiffnessTermsR)
            call ComputeNodalStiffnessTermsS(NI, NJ, PC, NodalStiffnessTermsS)
            call ComputeNodalStiffnessTermsH(BI, BJ, KM, NodalStiffnessTermsH)

            ! call AddLargeDeformationStiffnessTerms(ParticleIndex, BI, BJ, NodalStiffnessTerms, IsEnhancedStiffnessMatrix)

            call AssembleElementStiffnessMatrixK(DoFI, DoFJ, NodalStiffnessTermsK, ElementStiffnessMatrix)
            call AssembleElementStiffnessMatrixQ(DoFI, DoFJ, NodalStiffnessTermsQ, ElementStiffnessMatrix)
            call AssembleElementStiffnessMatrixR(DoFI, DoFJ, NodalStiffnessTermsR, ElementStiffnessMatrix)
            call AssembleElementStiffnessMatrixS(DoFI, DoFJ, NodalStiffnessTermsS, ElementStiffnessMatrix)
            call AssembleElementStiffnessMatrixH(DoFI, DoFJ, NodalStiffnessTermsH, ElementStiffnessMatrix)
            ! call AssembleElementStiffnessVectorH(DoFI, NodeID, NodalStiffnessTermsH, ElementStiffnessVector) !replaced by intflow
            
          end do
          
          !call ComputeNodalStiffnessTermsV(BI, KMG, NodalStiffnessTermsV)
          !call AssembleElementStiffnessVectorV(DoFI, NodalStiffnessTermsV, ElementStiffnessVector)
          
        end do

      end do

      call CompleteSymmetricElementStiffnessMatrix(ElementStiffnessMatrix) !off

      call RotateElementStiffnessMatrix(ElementID, ElementStiffnessMatrix) !off

      call ComputeEnhancedStiffnessMatrix(ElementStiffnessMatrix, IsEnhancedStiffnessMatrix) !off

      call CheckStiffnessNegativeDiagonalTerms(ElementID, ElementStiffnessMatrix) !no flow

      call FixDanglingElement(ElementID, ElementStiffnessMatrix, IsEnhancedStiffnessMatrix) !off

      if (CalParams%OutputDebugData) then
        open(TESTFILE,file='d:\tmp\testfileB.txt')
        if (ElementID==2) then
          write(TESTFILE,*) 'Int', Int
          write(TESTFILE,*)
          write(TESTFILE,*) 'ElementStiffnessMatrix'
          do i = 1, (ELEMENTNODES * (NVECTOR + 1))
            write(TESTFILE,902) i,(ElementStiffnessMatrix(i,j),j=1,(ELEMENTNODES * (NVECTOR + 1)))
902         format(i5,16(1x,1p,e15.5e3)) 
          end do
          write(TESTFILE,*)
          write(TESTFILE,*) 'ElementStiffnessVector'
          write(TESTFILE,902) 1,(ElementStiffnessVector(j),j=1,(ELEMENTNODES * (NVECTOR + 1)))
          write(TESTFILE,*)
        end if    
        close(TESTFILE)
      end if

      end subroutine ComputeElementStiffnessMatrix 

      subroutine AddExternalFlow(ElementID, ElementStiffnessVector)
      !**********************************************************************
      !
      !  Function:  Adds the ExtFlow vector to the external water load imposed 
      !             through liquid pressure loading. 
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: ElementID
      real(REAL_TYPE), dimension(ELEMENTNODES * (NVECTOR + 1)), intent(in) :: ElementStiffnessVector
      
      !local vairiables
      integer(INTEGER_TYPE) :: INode, I, NodeID
      integer(INTEGER_TYPE), parameter :: SOLVERFILE = 11
      
      do INode = 1, ELEMENTNODES
        I = (INode-1) * (NVECTOR+1) + (NVECTOR+1)
        NodeID = ElementConnectivities(INode, ElementID)
        ExtFlow(NodeID) = ExtFlow(NodeID) + ElementStiffnessVector(I)
      end do
      
      if (CalParams%OutputDebugData) then
        open(SOLVERFILE,file='d:\tmp\steady flow.txt')
        write(SOLVERFILE,*) 'FlowLoad', Counters%NodTot
        do i = 1, Counters%NodTot
          write(SOLVERFILE,*) i, ExtFlow(i)
        end do
        close(SOLVERFILE)
      end if
      
      end subroutine AddExternalFlow
      
      subroutine RotateElementStiffnessMatrix(IElement, ElementStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Rotation of Element Stiffness Matrix
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IElement
      real(REAL_TYPE), dimension(ELEMENTNODES * NVECTOR, ELEMENTNODES * NVECTOR), intent(inout) :: ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: INode, I, J, K, NodeID, I0
      real(REAL_TYPE), dimension(NVECTOR) :: X

      if (.not.IS3DCYLINDRIC) RETURN

      do INode = 1, ELEMENTNODES
        NodeID = ElementConnectivities(INode, IElement)
        if (IRotation(NodeID)/=0) then
          K = IRotation(NodeID)
          I0 = (INode - 1) * NVECTOR

          do I = 1, ELEMENTNODES * NVECTOR
            do J = 1, NVECTOR
              X(J) = ElementStiffnessMatrix(I, I0 + 1) * IRotMat(1, J, K) + ElementStiffnessMatrix(I, I0 + 2) * IRotMat(2, J, K) + ElementStiffnessMatrix(I, I0 + 3) * IRotMat(3, J, K)
            end do
            do J = 1, NVECTOR
              ElementStiffnessMatrix(I, I0 + J) = X(J)
            end do  
          end do

          do J = 1, ELEMENTNODES * NVECTOR
            do I = 1, NVECTOR
              X(I) = IRotMat(1, I, K) * ElementStiffnessMatrix(I0 + 1, J) + IRotMat(2, I, K) * ElementStiffnessMatrix(I0 + 2, J) + IRotMat(3, I, K) * ElementStiffnessMatrix(I0 + 3, J)
            end do
            do I = 1, NVECTOR
              ElementStiffnessMatrix(I0 + I, J) = X(I)
            end do  
          end do

        end if
      end do

      end subroutine RotateElementStiffnessMatrix

      subroutine FixDanglingElement(ElementID, ElementStiffnessMatrix, IsEnhancedStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Determines the sum of the stiffness terms on the matrix diagonal
      !             Determines nodes that are not connected to two elements and can thereby move freely
      !             Increase stiffness of dof's belonging to the free hanging nodes of the dangling element
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: ElementID
      real(REAL_TYPE),dimension(ELEMENTNODES * NVECTOR, ELEMENTNODES * NVECTOR),intent(inout) :: ElementStiffnessMatrix
      logical, intent(in) :: IsEnhancedStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I, INode, IDim, IDoF
      logical, dimension(ELEMENTNODES) :: IsFreeNode
      real(REAL_TYPE) :: SumDiagonalStiffness

      if (.not.CalParams%ImplicitIntegration%DoUseEnhancedStiffnessMatrix.or.IsEnhancedStiffnessMatrix) RETURN

      if (.not.IsDanglingElement(ElementID)) RETURN

      ! Determine the sum of the stiffness terms on the matrix diagonal
      SumDiagonalStiffness = 0.0
      do I = 1, NDOFL * ELEMENTNODES
        SumDiagonalStiffness = SumDiagonalStiffness + ElementStiffnessMatrix(I, I)
      end do

      ! Determine nodes that are not connected to two elements and can thereby move freely
      call DetermineFreeMovingNodes(ElementID, IsFreeNode)

      ! Increase stiffness of dof's belonging to the free hanging nodes of the dangling element
      do INode = 1, ELEMENTNODES
        if (IsFreeNode(INode)) then
          do IDim = 1, NVECTOR
            IDoF = (INode - 1) * 3 + IDim
            ElementStiffnessMatrix(IDoF, IDoF) = ElementStiffnessMatrix(IDoF, IDoF) + CalParams%ImplicitIntegration%DanglingElementFactor * SumDiagonalStiffness
          end do
        end if
      end do

      end subroutine FixDanglingElement

      subroutine DetermineFreeMovingNodes(DanglingElement, IsFreeNode)
      !**********************************************************************
      !
      !  Function:  Determines Free Moving Nodes
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DanglingElement
      logical, dimension(ELEMENTNODES), intent(out) :: IsFreeNode
      ! Local variables
      integer(INTEGER_TYPE) :: INode, NodeID, NConnectedElements, IElement,ElementID

      IsFreeNode = .true. ! Mark as free hanging node

      do INode = 1, ELEMENTNODES ! Loop over all nodes of the dangling element
        NodeID = ElementConnectivities(INode, DanglingElement)
        NConnectedElements = GetNElmOfNode(NodeID)
        if (NConnectedElements>1) then ! There is more than one element connected to NodeID

          do IElement = 1, NConnectedElements
            ElementID = GetElmIOfNode(NodeID, IElement)
            if ( (ElementID/=DanglingElement).and.IsActiveElement(ElementID)) then ! More than one active element is connected to NodeID
              IsFreeNode(INode) = .false. ! Mark as hinge node
              EXIT
            end if
          end do

        end if
      end do

      end subroutine DetermineFreeMovingNodes

      subroutine CheckStiffnessNegativeDiagonalTerms(IElement, ElementStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Checks for Negative diagonal term in stiffness matrix of element
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IElement
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR+1), ELEMENTNODES * (NVECTOR+1)),intent(in) :: ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I

      do I = 1, (NVECTOR+1) * ELEMENTNODES
        if (ElementStiffnessMatrix(I, I)<0.0) then
          write(5, *) ' Negative diagonal term in ','stiffness matrix of element ',IElement, ' at dof ', I, ': ',ElementStiffnessMatrix(I, I)
        end if
      end do

      end subroutine CheckStiffnessNegativeDiagonalTerms

      subroutine ComputeEnhancedStiffnessMatrix(ElementStiffnessMatrix, IsEnhancedStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Computes Enhanced Stiffness Matrix
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE),dimension(ELEMENTNODES * NVECTOR, ELEMENTNODES * NVECTOR),intent(inout) :: ElementStiffnessMatrix
      logical, intent(in) :: IsEnhancedStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I, J

      if (.not.IsEnhancedStiffnessMatrix) RETURN

      do I = 1, ELEMENTNODES * NVECTOR
        do J = 1, ELEMENTNODES * NVECTOR
          ElementStiffnessMatrix(I, J) = ElementStiffnessMatrix(I, J) * CalParams%ImplicitIntegration%StiffnessIncreaseFactor
        end do
      end do

      end subroutine ComputeEnhancedStiffnessMatrix

      subroutine CompleteSymmetricElementStiffnessMatrix(ElementStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Complete Symmetric Element Stiffness Matrix
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE),dimension(ELEMENTNODES * NVECTOR, ELEMENTNODES * NVECTOR),intent(inout) :: ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I, J

      if (.not.IsSymmetricStiffnessMatrix()) RETURN

      do I = 1, ELEMENTNODES * NVECTOR - 1
        do J = I + 1, ELEMENTNODES * NVECTOR
          ElementStiffnessMatrix(I, J) = ElementStiffnessMatrix(J, I)
        end do
      end do

      end subroutine CompleteSymmetricElementStiffnessMatrix

      subroutine AssembleElementStiffnessMatrixK(DoFI, DoFJ, NodalStiffnessTerms, ElementStiffnessMatrix)  
      !**********************************************************************
      !
      !  Function:  Assemble the element stiffness matrix
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI, DoFJ
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I, J

      do I = 1, NVECTOR
        do J = 1, NVECTOR
          ElementStiffnessMatrix(DoFI + I, DoFJ + J) = ElementStiffnessMatrix(DoFI + I, DoFJ + J) + NodalStiffnessTerms(I, J)
        end do
      end do

      end subroutine AssembleElementStiffnessMatrixK

      subroutine AssembleElementStiffnessMatrixQ(DoFI, DoFJ, NodalStiffnessTerms, ElementStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Assemble the element coupling matrix Q
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI, DoFJ
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I

      do I = 1, NVECTOR
          ElementStiffnessMatrix(DoFI + I, DoFJ + NVECTOR + 1) = ElementStiffnessMatrix(DoFI + I, DoFJ + NVECTOR + 1) + NodalStiffnessTerms(I, 1)
      end do
      
      end subroutine AssembleElementStiffnessMatrixQ

      subroutine AssembleElementStiffnessMatrixR(DoFI, DoFJ, NodalStiffnessTerms, ElementStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Assemble the element transpose of the coupling matrix Q
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI, DoFJ
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessMatrix
      ! Local variables
      integer(INTEGER_TYPE) :: I

      do I = 1, NVECTOR
          ElementStiffnessMatrix(DoFI + NVECTOR + 1, DoFJ + I) = ElementStiffnessMatrix(DoFI + NVECTOR + 1, DoFJ + I) - NodalStiffnessTerms(I, 1)
      end do
      
      end subroutine AssembleElementStiffnessMatrixR
      
      subroutine AssembleElementStiffnessMatrixS(DoFI, DoFJ, NodalStiffnessTerms, ElementStiffnessMatrix) 
      !**********************************************************************
      !
      !  Function:  Assemble the element compressibility matrix S
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI, DoFJ
      real(REAL_TYPE), dimension(1, 1), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessMatrix
      ! Local variables

          ElementStiffnessMatrix(DoFI + NVECTOR + 1, DoFJ + NVECTOR + 1) = ElementStiffnessMatrix(DoFI + NVECTOR + 1, DoFJ + NVECTOR + 1) + NodalStiffnessTerms(1, 1)

      end subroutine AssembleElementStiffnessMatrixS
      
      subroutine AssembleElementStiffnessMatrixH(DoFI, DoFJ, NodalStiffnessTerms, ElementStiffnessMatrix)  
      !**********************************************************************
      !
      !  Function:  Assemble the element conductivity matrix H
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI, DoFJ
      real(REAL_TYPE), dimension(1, 1), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1), ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessMatrix

          ElementStiffnessMatrix(DoFI + NVECTOR + 1, DoFJ + NVECTOR + 1) = ElementStiffnessMatrix(DoFI + NVECTOR + 1, DoFJ + NVECTOR + 1) + NodalStiffnessTerms(1, 1)

      end subroutine AssembleElementStiffnessMatrixH
      
      subroutine AssembleElementStiffnessVectorV(DoFI, NodalStiffnessTerms, ElementStiffnessVector)   
      !**********************************************************************
      !
      !  Function:  Assemble the element hydraulic load V
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessVector
      ! Local variables

          ElementStiffnessVector(DoFI + NVECTOR + 1) = ElementStiffnessVector(DoFI + NVECTOR + 1) + NodalStiffnessTerms(1, 1)
      
      end subroutine AssembleElementStiffnessVectorV
      
      subroutine AssembleElementStiffnessVectorH(DoFI, NodeID, NodalStiffnessTerms, ElementStiffnessVector)   
      !**********************************************************************
      !
      !  Function:  Compute the element pressure load from the total pressure
      !             on the nodes of the elements (IntFlow).             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: DoFI
      integer(INTEGER_TYPE), intent(in) :: NodeID
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(in) :: NodalStiffnessTerms
      real(REAL_TYPE),dimension(ELEMENTNODES * (NVECTOR + 1)),intent(inout) :: ElementStiffnessVector
      ! Local variables

          ElementStiffnessVector(DoFI + NVECTOR + 1) = ElementStiffnessVector(DoFI + NVECTOR + 1) - NodalStiffnessTerms(1, 1) * TotalPressure(NodeID) 
      
      end subroutine AssembleElementStiffnessVectorH
      
      subroutine ComputeNodalStiffnessTermsK(BI, BJ, D1, D2, D3, GG, NodalStiffnessTerms)  
      !**********************************************************************
      !
      !  Function:  Compute the nodal stiffness matrix using the elastic parameters
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: BI, BJ
      real(REAL_TYPE) :: D1, D2, D3, GG
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR), intent(inout) :: NodalStiffnessTerms
      
      if (NVECTOR == 3) then ! 3D
      NodalStiffnessTerms(1,1) = BI(1) * D1 * BJ(1) + GG * (BI(2) * BJ(2) + BI(3) * BJ(3))
      NodalStiffnessTerms(1,2) = BI(1) * D2 * BJ(2) + BI(2) * GG * BJ(1)
      NodalStiffnessTerms(1,3) = BI(1) * D2 * BJ(3) + BI(3) * GG * BJ(1)
      NodalStiffnessTerms(2,1) = BI(2) * D2 * BJ(1) + BI(1) * GG * BJ(2)
      NodalStiffnessTerms(2,2) = BI(2) * D1 * BJ(2) + GG * (BI(1) * BJ(1) + BI(3) * BJ(3))   
      NodalStiffnessTerms(2,3) = BI(2) * D2 * BJ(3) + BI(3) * GG * BJ(2)
      NodalStiffnessTerms(3,1) = BI(3) * D2 * BJ(1) + BI(1) * GG * BJ(3)
      NodalStiffnessTerms(3,2) = BI(3) * D2 * BJ(2) + BI(2) * GG * BJ(3)
      NodalStiffnessTerms(3,3) = BI(3) * D1 * BJ(3) + GG * (BI(2) * BJ(2) + BI(1) * BJ(1))
      elseif (NVECTOR == 2) then !2D
        ! add in later sprint
      end if    
      
      end subroutine ComputeNodalStiffnessTermsK

      subroutine ComputeNodalStiffnessTermsQ(BI, NJ, NodalStiffnessTerms)  
      !**********************************************************************
      !
      !  Function:  Compute the nodal terms of the coupling matrix Q
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE) :: NJ
      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: BI
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(inout) :: NodalStiffnessTerms
      integer(INTEGER_TYPE) :: I

      do I = 1, NVECTOR
        NodalStiffnessTerms(I,1) = BI(I) * NJ 
      end do
            
      end subroutine ComputeNodalStiffnessTermsQ
      
      subroutine ComputeNodalStiffnessTermsR(NI, BJ, NodalStiffnessTerms)  
      !**********************************************************************
      !
      !  Function:  Compute the nodal terms of the transpose of the coupling matrix Q
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE) :: NI
      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: BJ
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(inout) :: NodalStiffnessTerms
      integer(INTEGER_TYPE) :: I

      do I = 1, NVECTOR
        NodalStiffnessTerms(I,1) = NI * BJ(I) 
      end do
            
      end subroutine ComputeNodalStiffnessTermsR
      
      subroutine ComputeNodalStiffnessTermsS(NI, NJ, PC, NodalStiffnessTerms) 
      !**********************************************************************
      !
      !  Function:  Compute the nodal terms of the compressibility matrix S
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE) :: NI, NJ, PC
      real(REAL_TYPE), dimension(1, 1), intent(inout) :: NodalStiffnessTerms
      ! Local variables

      NodalStiffnessTerms(1,1) = PC * NI * NJ 

      end subroutine ComputeNodalStiffnessTermsS
      
      subroutine ComputeNodalStiffnessTermsH(BI, BJ, KM, NodalStiffnessTerms)  
      !**********************************************************************
      !
      !  Function:  Compute the nodal terms of the conductivity matrix H
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE), intent(in) :: KM
      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: BI, BJ
      real(REAL_TYPE), dimension(1, 1), intent(inout) :: NodalStiffnessTerms
      ! Local variables

      if (NVECTOR == 3) then ! 3D
        NodalStiffnessTerms(1,1) = KM * BI(1) * BJ(1) + KM * BI(2) * BJ(2) + KM * BI(3) * BJ(3)
      elseif (NVECTOR == 2) then ! 2D
        ! add in later sprint
      end if

      end subroutine ComputeNodalStiffnessTermsH
      
      subroutine ComputeNodalStiffnessTermsV(BI, KMG, NodalStiffnessTerms) 
      !**********************************************************************
      !
      !  Function:  Compute the nodal terms of the hydraulic load V
      !             
      !
      !**********************************************************************
      implicit none

      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: BI
      real(REAL_TYPE) :: KMG
      real(REAL_TYPE), dimension(NVECTOR, 1), intent(inout) :: NodalStiffnessTerms
      ! Local variables

      NodalStiffnessTerms(1,1) = KMG * BI(2) 

      end subroutine ComputeNodalStiffnessTermsV
      
      subroutine AddLargeDeformationStiffnessTerms(IntGlo, BI, BJ, NodalStiffnessTerms, IsEnhancedStiffnessMatrix)
      !**********************************************************************
      !
      !  Function:  Add Large Deformation Stiffness Terms
      !             
      !
      !**********************************************************************
      implicit none

      integer(INTEGER_TYPE), intent(in) :: IntGlo
      real(REAL_TYPE), dimension(NVECTOR), intent(in) :: BI, BJ
      real(REAL_TYPE), dimension(NVECTOR, NVECTOR), intent(inout) :: NodalStiffnessTerms
      logical, intent(in) :: IsEnhancedStiffnessMatrix
      ! Local variables
      real(REAL_TYPE) :: SXX0, SYY0, SZZ0, SXY0, SYZ0, SZX0
      real(REAL_TYPE) :: S11, S12, S13, S21, S22, S23, S31, S32, S33

      if (.not.CalParams%ApplyObjectiveStress.or.IsEnhancedStiffnessMatrix) RETURN

      if (NVECTOR == 3) then ! 3D
      SXX0 = SigmaEff0Array(intGlo,1)
      SYY0 = SigmaEff0Array(intGlo,2)
      SZZ0 = SigmaEff0Array(intGlo,3)
      SXY0 = SigmaEff0Array(intGlo,4)
      SYZ0 = SigmaEff0Array(intGlo,5)
      SZX0 = SigmaEff0Array(intGlo,6)

      S11 = 0
      S11 = S11 - (2 * BI(1) * BJ(1) + BI(2) * BJ(2) + BI(3) * BJ(3)) * SXX0 / 2
      S11 = S11 + (BI(2) * BJ(2)) * SYY0 / 2
      S11 = S11 + (BI(3) * BJ(3)) * SZZ0 / 2
      S11 = S11 + (BI(2) * BJ(3) + BI(3) * BJ(2)) * SYZ0 / 2
      S12 = 0
      S12 = S12 - (BI(2) * BJ(1)) * SXX0 / 2
      S12 = S12 - (BI(2) * BJ(1)) * SYY0 / 2
      S12 = S12 - (2 * BI(1) * BJ(1) + 2 * BI(2) * BJ(2) + BI(3) * BJ(3)) * SXY0 / 2
      S12 = S12 - (BI(3) * BJ(1)) *SYZ0/2
      S12 = S12 - (BI(2) * BJ(3)) *SZX0/2
      S13 = 0
      S13 = S13 - (BI(3) * BJ(1)) * SXX0 / 2
      S13 = S13 - (BI(3) * BJ(1)) * SZZ0 / 2
      S13 = S13 - (BI(3) * BJ(2)) * SXY0 / 2
      S13 = S13 - (BI(2) * BJ(1)) * SYZ0 / 2
      S13 = S13 - (2 * BI(1) * BJ(1) + BI(2) * BJ(2) + 2 * BI(3) * BJ(3)) * SZX0 / 2
      S21 = 0
      S21 = S21 - (BI(1) * BJ(2)) * SXX0 / 2
      S21 = S21 - (BI(1) * BJ(2)) * SYY0 / 2
      S21 = S21 - (2 * BI(1) * BJ(1) + 2 * BI(2) * BJ(2) + BI(3) * BJ(3)) * SXY0 / 2
      S21 = S21 - (BI(1) * BJ(3)) * SYZ0 / 2
      S21 = S21 - (BI(3) * BJ(2)) * SZX0 / 2
      S22 = 0
      S22 = S22 + (BI(1) * BJ(1)) * SXX0 / 2
      S22 = S22 - (BI(1) * BJ(1) + 2 * BI(2) * BJ(2) + BI(3) * BJ(3)) * SYY0 / 2
      S22 = S22 + (BI(3) * BJ(3)) * SZZ0 / 2
      S22 = S22 + (BI(1) * BJ(3) + BI(3) * BJ(1)) * SZX0 / 2
      S23 = 0
      S23 = S23 - (BI(3) * BJ(2)) * SYY0 / 2
      S23 = S23 - (BI(3) * BJ(2)) * SZZ0 / 2
      S23 = S23 - (BI(3) * BJ(1)) * SXY0 / 2
      S23 = S23 - (BI(1) * BJ(1) + 2 * BI(2) * BJ(2) + 2 * BI(3) * BJ(3)) * SYZ0 / 2
      S23 = S23 - (BI(1) * BJ(2)) * SZX0 / 2
      S31 = 0
      S31 = S31 - (BI(1) * BJ(3)) * SXX0 / 2
      S31 = S31 - (BI(1) * BJ(3)) * SZZ0 / 2
      S31 = S31 - (BI(2) * BJ(3)) * SXY0 / 2
      S31 = S31 - (BI(1) * BJ(2)) * SYZ0 / 2
      S31 = S31 - (2 * BI(1) * BJ(1) + BI(2) * BJ(2) + 2 * BI(3) * BJ(3)) * SZX0 / 2
      S32 = 0
      S32 = S32 - (BI(2) * BJ(3)) * SYY0 / 2
      S32 = S32 - (BI(2) * BJ(3)) * SZZ0 / 2
      S32 = S32 - (BI(1) * BJ(3)) * SXY0 / 2
      S32 = S32 - (BI(1) * BJ(1) + 2 * BI(2) * BJ(2) + 2 * BI(3) * BJ(3)) * SYZ0 / 2
      S32 = S32 - (BI(2) * BJ(1)) * SZX0 / 2
      S33 = 0
      S33 = S33 + (BI(1) * BJ(1)) * SXX0 / 2
      S33 = S33 - (BI(1) * BJ(1) + BI(2) * BJ(2) + 2 * BI(3) * BJ(3)) * SZZ0 / 2
      S33 = S33 + (BI(2) * BJ(2)) * SYY0 / 2
      S33 = S33 + (BI(1) * BJ(2) + BI(2) * BJ(1)) * SXY0 / 2

      NodalStiffnessTerms(1, 1) = NodalStiffnessTerms(1, 1) + S11
      NodalStiffnessTerms(1, 2) = NodalStiffnessTerms(1, 2) + S12
      NodalStiffnessTerms(1, 3) = NodalStiffnessTerms(1, 3) + S13
      NodalStiffnessTerms(2, 1) = NodalStiffnessTerms(2, 1) + S21
      NodalStiffnessTerms(2, 2) = NodalStiffnessTerms(2, 2) + S22
      NodalStiffnessTerms(2, 3) = NodalStiffnessTerms(2, 3) + S23
      NodalStiffnessTerms(3, 1) = NodalStiffnessTerms(3, 1) + S31
      NodalStiffnessTerms(3, 2) = NodalStiffnessTerms(3, 2) + S32
      NodalStiffnessTerms(3, 3) = NodalStiffnessTerms(3, 3) + S33
      elseif (NVECTOR == 2) then ! 2D
          ! add in later sprint
      end if    

      end subroutine AddLargeDeformationStiffnessTerms

      subroutine SetZeroParticleIncStrain()
      !**********************************************************************
      !
      !  Function:  Reset particle strains to zero
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: IParticle

      do IParticle = 1, Counters%NParticles
        Particles(IParticle)%EpsStep = 0.0
      end do

      end subroutine SetZeroParticleIncStrain
      
      subroutine GetNodalIntFlow(DT)
      !**********************************************************************
      !
      !  Function: computes the internal flow vector 
      !             
      !
      !**********************************************************************
          implicit none
                 
          real(REAL_TYPE), intent(in) ::  DT
          
          ! Local variables
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
          integer(INTEGER_TYPE) :: IntGlo, IEl, Int, IAEl, NElemPart, NodeI, ISet, NodeID, I, L
          real(REAL_TYPE), dimension(NVECTOR) :: FAC
          real(REAL_TYPE) :: Det, Weight, XNu, BFac, SPhi, SPsi, GG, Cohec, Tens, &
                              PC, NL, KL, NI
                   
          IntGlo = 0
          IntFlow = 0.0
          
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
          
            if (IsParticleIntegration(IEl)) then
              NElemPart = NPartEle(IEl)
            else
              NElemPart = ELEMENTGAUSSPOINTS
            end if
             
            call FormB3(1, IEl, ElementConnectivities, NodalCoordinates, B, Det, Weight)

            do Int = 1, 1 !mixed
              
              IntGlo = GetParticleIndex(Int, IEl)
              
              call GetMaterialData(IntGlo, ISet, XNu, BFac, SPhi, SPsi, GG, Cohec, Tens)

              NL = MatParams(ISet)%InitialPorosity 
              KL = MatParams(ISet)%BulkModulusLiquid
              PC = NL / KL
        
              if (IsParticleIntegration(IEl)) then 
                Weight = Particles(IntGlo)%IntegrationWeight
              end if
              
              do I = 1, NVECTOR
                 FAC(I) = DT * Weight * VelocityWaterArray(IntGlo,I)
              end do
              
              DO NodeI=1,ELEMENTNODES
                NodeID = ElementConnectivities(NodeI, IEl)
                NI = GPShapeFunction(Int, NodeI) * Weight
                IntFlow(NodeID) = IntFlow(NodeID) - NI * (Particles(IntGlo)%EpsStep(1) + Particles(IntGlo)%EpsStep(2) + Particles(IntGlo)%EpsStep(3) - PC * (Particles(IntGlo)%WaterPressure - Particles(IntGlo)%WaterPressure0))
                do L = 1, NDIM
                 IntFlow(NodeID) = IntFlow(NodeID) + FAC(L) * B(L,NodeI)
                enddo
              end do
            end do
          
          end do

        end subroutine 

      subroutine UpdateParticlePressure()
      !**********************************************************************
      !
      !  Function: Update particle water pressure with increment 
      !             
      !
      !**********************************************************************     
          implicit none
        
          !Local variables
          real(REAL_TYPE) :: DWP
          integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, ParticleIndex, INode, NodeID
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape         
          
          do IElement = 1, Counters%NEl   
            if (IsActiveElement(IElement)) then 
              NElemPart = NPartEle(IElement)
              do IParticle = 1, NElemPart
                ParticleIndex = GetParticleIndex(IParticle, IElement) 
                ParticleShape = ShapeValuesArray(ParticleIndex,:)
                DWP = 0.0
                do INode = 1, ELEMENTNODES  
                     NodeID = ElementConnectivities(INode, IElement)
                     DWP = DWP + ParticleShape(INode) * IncrementalPressure(NodeID)
                end do
                Particles(ParticleIndex)%WaterPressure = Particles(ParticleIndex)%WaterPressure0 + DWP
              end do
            end if
          end do
          
      end subroutine UpdateParticlePressure

      subroutine UpdateParticleDisplacement()
      !**********************************************************************
      !
      !  Function:  Update Particle Displacement
      !             
      !
      !**********************************************************************     
          implicit none
        
          !Local variables
          real(REAL_TYPE) :: DU(NVECTOR)
          integer(INTEGER_TYPE) :: IElement, NElemPart, IParticle, ParticleIndex, INode, NodeID, IDof, IDim
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape         
          
          do IElement = 1, Counters%NEl   
            if (IsActiveElement(IElement)) then 
              NElemPart = NPartEle(IElement)
              do IParticle = 1, NElemPart
                ParticleIndex = GetParticleIndex(IParticle, IElement) 
                ParticleShape = ShapeValuesArray(ParticleIndex,:)
                DU = 0.0
                do INode = 1, ELEMENTNODES  
                   NodeID = ElementConnectivities(INode, IElement)
                   IDof = ReducedDof(NodeID)
                   do IDim = 1, NVECTOR
                     DU(IDim) = DU(IDim) + ParticleShape(INode) * IncrementalDisplacementSoil(IDof+IDim, 1)
                   end do
                end do
                UArray(ParticleIndex,:) = UArray(ParticleIndex,:) + DU
              end do
            end if
          end do
          
      end subroutine UpdateParticleDisplacement
      
      subroutine IncreaseParticleTotalPressure()
      !**********************************************************************
      !
      !  Function:  Increase Water Pressure 
      !             
      !
      !**********************************************************************
      implicit none

      ! Local variables
      integer(INTEGER_TYPE) :: I
      
      do I = 1, Counters%NParticles
        Particles(I)%WaterPressure0 = Particles(I)%WaterPressure
      end do
         
      end subroutine IncreaseParticleTotalPressure
      
      subroutine ComputeParticleVelocityWater()
      !**********************************************************************
      !
      !  Function:  ! computes pressure part of Darcy flux
      !             
      !
      !**********************************************************************
      
        implicit none
      
        !Local variables
        real(REAL_TYPE) :: DU(NVECTOR), MuL, KappaL, KM, Det, Weight
        integer(INTEGER_TYPE) :: IActiveElement, NElemPart, IParticle, ParticleIndex, INode, NodeID, IDim, ISet, ElementID
        real(REAL_TYPE), dimension(ELEMENTNODES) :: ParticleShape
        real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
        real(REAL_TYPE), dimension(Counters%NodTot) :: CurrentPressure
          
        ISet = 1 ! material data has to be added to material point properties
        MuL = MatParams(ISet)%ViscosityLiquid 
        KappaL = MatParams(ISet)%IntrinsicPermeabilityLiquid
        KM = KappaL / MuL
        
        CurrentPressure = TotalPressure + IncrementalPressure
        
        do IActiveElement = 1, Counters%NAEl
          ElementID = ActiveElement(IActiveElement)
          NElemPart = NPartEle(IActiveElement)
          call FormB3(1, ElementID, ElementConnectivities, NodalCoordinates, B, Det, Weight)
          do IParticle = 1, NElemPart
            ParticleIndex = GetParticleIndex(IParticle, IActiveElement) 
            ParticleShape = ShapeValuesArray(ParticleIndex,:)
            DU = 0.0
            do INode = 1, ELEMENTNODES  
              NodeID = ElementConnectivities(INode, IActiveElement)
              do IDim = 1, NVECTOR
                DU(IDim) = DU(IDim) + KM * B(IDim,INode) * CurrentPressure(NodeID)
              end do
            end do
            VelocityWaterArray(ParticleIndex,:) = DU
          end do
        end do
      
      end subroutine ComputeParticleVelocityWater
      
      end module ModQuasiStaticImplicit