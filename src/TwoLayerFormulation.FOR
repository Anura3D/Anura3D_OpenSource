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
	  
	  
	  module ModTwoLayerFormulation
      !**********************************************************************
      !
      !   Function : Contains routines for the two-layer formulation
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModCounters
      use ModReadCalculationData
      use ModMPMData
      use ModMeshAdjacencies
      use ModParticle
      
      implicit none
      
        type TwoLayerElementType(vsize)
        
          integer(INTEGER_TYPE), LEN :: vsize
          integer(INTEGER_TYPE) :: ContainedMaterialTypes = ContainedMaterialTypeUNDEFINED ! solid, liquid, or both
          logical :: IsBoundaryOfSolidElements, IsBoundaryOfLiquidElements, IsSolidMPwithPhaseLiquidInElement
          real(REAL_TYPE) :: ConcentrationRatioLiquidL = 0.0
          real(REAL_TYPE) :: ConcentrationRatioSolidL = 0.0
          real(REAL_TYPE) :: ConcentrationRatioLiquidS = 0.0
          real(REAL_TYPE) :: ConcentrationRatioSolidS = 0.0
          real(REAL_TYPE), dimension(vsize) :: GradientPorosityL 
        
        end type TwoLayerElementType
      
        type TwoLayerNodeType

          real(REAL_TYPE) :: DensitySolidL = 0.0 ! Solid density field in active elements which contain SOLID+LIQUID and LIQUID material points
          real(REAL_TYPE) :: DensityLiquidL = 0.0 ! Water density field in active elements which contain SOLID+LIQUID and LIQUID material points
          real(REAL_TYPE) :: DensitySolidS = 0.0 ! Solid density field in active elements which contain SOLID+LIQUID and SOLID material points
          real(REAL_TYPE) :: DensityLiquidS = 0.0 ! Water density field in active elements which contain SOLID+LIQUID and SOLID material points
          logical :: IsNodeOfBoundSolidDomain
          real(REAL_TYPE) :: FillingRatio = 0.0 ! Filling ratio field in all elements of liquid domain

        end type TwoLayerNodeType

        type TwoLayerType
        
            type(TwoLayerElementType(:)), dimension(:), allocatable :: Elements ! List of element objects
            type(TwoLayerNodeType), dimension(:), allocatable :: Nodes ! List of node objects
            real(REAL_TYPE), dimension(:,:), allocatable :: InteractionForceLiquid
            real(REAL_TYPE), dimension(:,:), allocatable :: InteractionForceSolid

          contains
          
            procedure :: Initialise => Initialise
            procedure :: Destroy => Destroy
            procedure :: DetermineContainedMaterialTypes => DetermineContainedMaterialTypes
            procedure :: DetermineConcentrationRatios => DetermineConcentrationRatios
            procedure :: DetermineConcentrationRatioElementFND => DetermineConcentrationRatioElementFND
            procedure :: DetermineInteractionForces => DetermineInteractionForces
            procedure :: DetermineTwoLayerStatus => DetermineTwoLayerStatus
            
        end type TwoLayerType
        
        type(TwoLayerType) :: TwoLayerData
        
      contains ! Routines of this module

      
        subroutine Initialise(This)
        !**********************************************************************
        !
        !    Function: initalise double point (2 layers of material points)   
        !
        !**********************************************************************
          implicit none
          
            class(TwoLayerType), intent(inout) :: This
            
            ! local variables
            integer(INTEGER_TYPE) :: IError
        
            if (.not.(NFORMULATION==2)) RETURN
            
            if (allocated(This%Elements)) then
              deallocate(This%Elements, stat = IError)
            end if
            allocate(TwoLayerElementType(NVECTOR)::This%Elements(Counters%NEl), stat = IError)
            
            if (allocated(This%Nodes)) then
              deallocate(This%Nodes, stat = IError)
            end if
            allocate(This%Nodes(Counters%NodTot), stat = IError)
            
            if (allocated(This%InteractionForceLiquid)) then
              deallocate(This%InteractionForceLiquid, stat = IError)
            end if
            allocate(This%InteractionForceLiquid(Counters%N, Counters%NEntity), stat = IError)
            This%InteractionForceLiquid = 0.0

            if (allocated(This%InteractionForceSolid)) then
              deallocate(This%InteractionForceSolid, stat = IError)
            end if
            allocate(This%InteractionForceSolid(Counters%N, Counters%NEntity), stat = IError)
            This%InteractionForceSolid = 0.0

        end subroutine Initialise
 
        
        subroutine Destroy(This)
        !**********************************************************************
        !
        !    Function: destroy double point (2 layers of material points)   
        !
        !**********************************************************************          
        
            class(TwoLayerType), intent(inout) :: This
            
            ! local variables
            integer(INTEGER_TYPE) :: IError
          
            if (.not.(NFORMULATION==2)) RETURN
            
            deallocate(This%Elements, stat = IError)
            deallocate(This%Nodes, stat = IError)
            deallocate(This%InteractionForceLiquid, stat = IError)
            deallocate(This%InteractionForceSolid, stat = IError)
            
        end subroutine Destroy

        
        subroutine DetermineConcentrationRatio(SolidNodalDensity, LiquidNodalDensity, MaterialPointType, OnlyElementsWLiquidMPs)
        !**********************************************************************
        !
        ! Function: Computes ConcentrationRatio and Porosity for each 
        !     material point within each element
        !
        !**********************************************************************
          implicit none
        
            real(REAL_TYPE), dimension(:), intent(in) :: SolidNodalDensity, LiquidNodalDensity
            integer(INTEGER_TYPE), intent(in) :: MaterialPointType
            logical, intent(in) :: OnlyElementsWLiquidMPs
            
            ! local variables
            integer(INTEGER_TYPE) :: IElement, IParticle, ParticleID, I, J
            integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
            real(REAL_TYPE), dimension(ELEMENTNODES) :: ElementNodalDensitySolid, ElementNodalDensityLiquid, ParticleShape
            real(REAL_TYPE) :: ConcentrationRatioSolid, ConcentrationRatioLiquid
            real(REAL_TYPE) :: ConstDensityLiquid, ConstDensitySolid
            
              !--- Calculate  density  ----
              ConstDensityLiquid = 0.0
              ConstDensitySolid = 0.0
                
              do J = 1, Counters%NLayers ! loop over all material points in element
                if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)    
                end if
              end do 
              ! ---------------------------------------------------
            
            do IElement = 1, Counters%NEl ! loop over all elements
                
              if(OnlyElementsWLiquidMPs) then
                if(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) CYCLE
              else
                if(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) CYCLE
              end if

              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              ElementNodalDensitySolid = SolidNodalDensity(NodeIDs)
              ElementNodalDensityLiquid = LiquidNodalDensity(NodeIDs)
              
              do IParticle = 1, NPartEle(IElement) ! loop over all material points in element
                ParticleID = GetParticleIndex(IParticle, IElement)
                
                if (MaterialPointTypeArray(ParticleID)/=MaterialPointType) CYCLE
                
                ParticleShape = ShapeValuesArray(ParticleID,:)
                ConcentrationRatioSolid = 0.0
                ConcentrationRatioLiquid = 0.0
                
                do I = 1, ELEMENTNODES
                  ConcentrationRatioSolid = ConcentrationRatioSolid + ParticleShape(I) * ElementNodalDensitySolid(I)
                  ConcentrationRatioLiquid = ConcentrationRatioLiquid + ParticleShape(I) * ElementNodalDensityLiquid(I)
                end do
                
                if(ConstDensityLiquid/=0.0) then
                    ConcentrationRatioLiquid = ConcentrationRatioLiquid / ConstDensityLiquid
                else
                    ConcentrationRatioLiquid = 0.0
                end if
                
                if(ConstDensitySolid/=0.0) then
                    ConcentrationRatioSolid = ConcentrationRatioSolid / ConstDensitySolid
                else
                    ConcentrationRatioSolid = 0.0
                end if
                
                if(OnlyElementsWLiquidMPs) then
                  Particles(ParticleID)%ConcentrationRatioLiquidL = ConcentrationRatioLiquid
                  Particles(ParticleID)%ConcentrationRatioSolidL = ConcentrationRatioSolid
                  Particles(ParticleID)%PorosityL = 1.0 - ConcentrationRatioSolid 
                else 
                  Particles(ParticleID)%ConcentrationRatioLiquidS = ConcentrationRatioLiquid
                end if
                 
              end do ! material point loop
              
            end do ! element loop
            
        end subroutine DetermineConcentrationRatio

        
        subroutine DetermineConcentrationRatioElementFND(This)
        !**********************************************************************
        !
        ! Function: Computes the smoothed ConcentrationRatio within each element
        !
        !**********************************************************************
        
          implicit none
 
            class(TwoLayerType), intent(inout) :: This
            
            ! local variables
            integer(INTEGER_TYPE) :: IElement, II, J
            integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
            real(REAL_TYPE), dimension(ELEMENTNODES) :: ElementNodalDensitySolidL, ElementNodalDensityLiquidL
            real(REAL_TYPE), dimension(ELEMENTNODES) :: ElementNodalDensitySolidS, ElementNodalDensityLiquidS
            real(REAL_TYPE) :: ConcentrationRatioSolidL, ConcentrationRatioLiquidL
            real(REAL_TYPE) :: ConcentrationRatioSolidS, ConcentrationRatioLiquidS
            real(REAL_TYPE) :: ConstDensityLiquid, ConstDensitySolid, LiquidViscosity, IntrinsicPermeability
            
              !--- Calculate  density  ----
              ConstDensityLiquid = 0.0
              ConstDensitySolid = 0.0
              LiquidViscosity = 0.0
              IntrinsicPermeability = 0.0
                
              do J = 1, Counters%NLayers ! loop over all material points in element
                if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
                    LiquidViscosity = MatParams(J)%ViscosityLiquid
                    IntrinsicPermeability = MatParams(J)%IntrinsicPermeabilityLiquid
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                    LiquidViscosity = MatParams(J)%ViscosityLiquid 
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000) 
                    IntrinsicPermeability = MatParams(J)%IntrinsicPermeabilityLiquid
                end if
              end do 
              ! ---------------------------------------------------
            
            do IElement = 1, Counters%NEl ! loop over all elements
                
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              ElementNodalDensitySolidL = This%Nodes(NodeIDs)%DensitySolidL
              ElementNodalDensityLiquidL = This%Nodes(NodeIDs)%DensityLiquidL
              ElementNodalDensitySolidS = This%Nodes(NodeIDs)%DensitySolidS
              ElementNodalDensityLiquidS = This%Nodes(NodeIDs)%DensityLiquidS
              
              ConcentrationRatioSolidL = 0.0
              ConcentrationRatioLiquidL= 0.0
              ConcentrationRatioSolidS = 0.0
              ConcentrationRatioLiquidS = 0.0
              
              do II = 1, ELEMENTNODES
                ConcentrationRatioSolidL = ConcentrationRatioSolidL + ElementNodalDensitySolidL(II)
                ConcentrationRatioLiquidL = ConcentrationRatioLiquidL + ElementNodalDensityLiquidL(II)
                ConcentrationRatioSolidS = ConcentrationRatioSolidS + ElementNodalDensitySolidS(II)
                ConcentrationRatioLiquidS = ConcentrationRatioLiquidS + ElementNodalDensityLiquidS(II)
              end do
              
              ConcentrationRatioSolidL = ConcentrationRatioSolidL / dble(ELEMENTNODES)
              ConcentrationRatioLiquidL = ConcentrationRatioLiquidL / dble(ELEMENTNODES)
              ConcentrationRatioSolidS = ConcentrationRatioSolidS / dble(ELEMENTNODES)
              ConcentrationRatioLiquidS = ConcentrationRatioLiquidS / dble(ELEMENTNODES)
              
              if(ConstDensityLiquid>0.0) then
                  ConcentrationRatioLiquidL = ConcentrationRatioLiquidL / ConstDensityLiquid
                  ConcentrationRatioLiquidS = ConcentrationRatioLiquidS / ConstDensityLiquid
              else
                  ConcentrationRatioLiquidL = 0.0
                  ConcentrationRatioLiquidS = 0.0
              end if
              
              if(ConstDensitySolid>0.0) then
                  ConcentrationRatioSolidL = ConcentrationRatioSolidL / ConstDensitySolid 
                  ConcentrationRatioSolidS = ConcentrationRatioSolidS / ConstDensitySolid
              else
                  ConcentrationRatioSolidL = 0.0
                  ConcentrationRatioSolidS = 0.0
              end if
              
              TwoLayerData%Elements(IElement)%ConcentrationRatioLiquidL = ConcentrationRatioLiquidL
              TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL = ConcentrationRatioSolidL
              
              TwoLayerData%Elements(IElement)%ConcentrationRatioLiquidS = ConcentrationRatioLiquidS
              TwoLayerData%Elements(IElement)%ConcentrationRatioSolidS = ConcentrationRatioSolidS
              
            end do ! element loop
            
        end subroutine DetermineConcentrationRatioElementFND
        

        subroutine DetermineFillingRatioLiquid(DoublePointData)
        !**********************************************************************
        !
        ! Function: Computes ConcentrationRatio and Porosity for each 
        !     material point within each element
        !
        !**********************************************************************
          implicit none
        
          type(TwoLayerType), intent(in):: DoublePointData

            ! local variables
            integer(INTEGER_TYPE) :: IElement, IParticle, ParticleID, I
            integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
            real(REAL_TYPE), dimension(ELEMENTNODES) :: ElementNodalFillingRatio, ParticleShape
            real(REAL_TYPE) :: FillingRatio
              
            do IElement = 1, Counters%NEl ! loop over all elements
                
              if(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) CYCLE

              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              ElementNodalFillingRatio = DoublePointData%Nodes(NodeIDs)%FillingRatio
              
              do IParticle = 1, NPartEle(IElement) ! loop over all material points in element
                ParticleID = GetParticleIndex(IParticle, IElement)
                
                if (MaterialPointTypeArray(ParticleID)/=MaterialPointTypeLiquid) CYCLE
                
                ParticleShape = ShapeValuesArray(ParticleID,:)
                FillingRatio = 0.0
                
                do I = 1, ELEMENTNODES
                  FillingRatio = FillingRatio + ParticleShape(I) * ElementNodalFillingRatio(I)
                end do
                
                Particles(ParticleID)%FillingRatio = FillingRatio

              end do ! material point loop
              
            end do ! element loop
            
      end subroutine DetermineFillingRatioLiquid 
        
      
        subroutine DetermineContainedMaterialTypes(This)
        !**********************************************************************
        !
        !    Function: Determines material point types in the element   
        !
        !**********************************************************************        
          implicit none
          
            class(TwoLayerType), intent(inout) :: This
            
            ! local variables
            integer(INTEGER_TYPE) :: I, J, ParticleIndex
            
            if (.not.(NFORMULATION==2)) RETURN

            do I = 1, Counters%NEl ! loop over all elements
              This%Elements(I)%ContainedMaterialTypes = ContainedMaterialTypeUNDEFINED
                              
              do J = 1, NPartEle(I) ! loop over all material points in element
                ParticleIndex = GetParticleIndex(J, I)
              
                if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeUNDEFINED) then
                  if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then
                    This%Elements(I)%ContainedMaterialTypes = ContainedMaterialTypeSOLID
                  else if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                    This%Elements(I)%ContainedMaterialTypes = ContainedMaterialTypeLIQUID
                  end if
                else if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) then
                  if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                    This%Elements(I)%ContainedMaterialTypes = ContainedMaterialTypeSOLIDLIQUID
                  end if
                else if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
                  if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then
                    This%Elements(I)%ContainedMaterialTypes = ContainedMaterialTypeSOLIDLIQUID
                  end if
                end if
              end do ! material point loop

              do J = 1, NPartEle(I) ! loop over all material points in element
                ParticleIndex = GetParticleIndex(J, I)
                if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) then
                  Particles(ParticleIndex)%ContainedMaterialTypes = ContainedMaterialTypeSOLID
                else if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
                  Particles(ParticleIndex)%ContainedMaterialTypes = ContainedMaterialTypeLIQUID
                else if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID) then
                  Particles(ParticleIndex)%ContainedMaterialTypes = ContainedMaterialTypeSOLIDLIQUID
                end if
              end do ! material point loop
              
            end do ! element loop
            
              
        end subroutine DetermineContainedMaterialTypes
 
        
        subroutine DetermineConcentrationRatios(This)
        !**********************************************************************
        !
        !    Function: determine density field, concentration ratios, boundary elements, and filling ratio   
        !
        !**********************************************************************        
          implicit none
        
            class(TwoLayerType), intent(inout) :: This
            
            if (.not.(NFORMULATION==2)) RETURN  
            if (.not.(CalParams%NumberOfPhases==2)) RETURN
            
            call This%DetermineContainedMaterialTypes()
            
            if (CalParams%IStep==1) then
              call DetermineDensityField(This%Nodes%DensityLiquidL, MaterialPointTypeLiquid, .false., .true.)         
              call DetermineDensityField(This%Nodes%DensitySolidL, MaterialPointTypeSolid, .false., .true.)
              call This%DetermineConcentrationRatioElementFND()
              call DetermineConcentrationRatio(This%Nodes%DensitySolidL, This%Nodes%DensityLiquidL,2,.true.)
              call AdjustMassLiquidMP()
            end if
            
            call DetermineDensityField(This%Nodes%DensityLiquidL, MaterialPointTypeLiquid, .false., .true.)         
            call DetermineDensityField(This%Nodes%DensitySolidL, MaterialPointTypeSolid, .false., .true.)
            call DetermineDensityField(This%Nodes%DensityLiquidS, MaterialPointTypeLiquid, .false., .false.)         
            call DetermineDensityField(This%Nodes%DensitySolidS, MaterialPointTypeSolid, .false., .false.)
            
            call DetermineConcentrationRatio(This%Nodes%DensitySolidL, This%Nodes%DensityLiquidL,MaterialPointTypeLiquid,.true.)
            call DetermineConcentrationRatio(This%Nodes%DensitySolidL, This%Nodes%DensityLiquidL,MaterialPointTypeSolid,.true.)
            call DetermineConcentrationRatio(This%Nodes%DensitySolidS, This%Nodes%DensityLiquidS, MaterialPointTypeSolid,.false.)
            
            call This%DetermineConcentrationRatioElementFND()
            
            call DetermineBoundaryElementSolidDomain2LayerForm()
            call DetermineBoundaryElementLiquidDomain2LayerForm()
            
            call DetermineFillingRatioField()
            call DetermineFillingRatioLiquid(This)
            
        end subroutine DetermineConcentrationRatios
  
        
        subroutine DetermineInteractionForces(This)
        !**********************************************************************
        !
        !    Function: determine interaction forces between solid and liquid  
        !
        !**********************************************************************
          implicit none
        
            class(TwoLayerType), intent(inout) :: This
            
            ! local variables
            integer(INTEGER_TYPE) :: I, J, K, ParticleIndex, IElement,II
            integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs, Nix, Niy, Niz
            !integer(INTEGER_TYPE) :: NumberPartLiquid, NumberPartSolid
            real(REAL_TYPE) :: WeiGP, Weighting, DetJac, VolumeMP
            real(REAL_TYPE), dimension(NVECTOR) :: PosGP, Term, GradientNf, GradientNs
            real(REAL_TYPE), dimension(NTENSOR) :: AverageSigma
            real(REAL_TYPE), dimension(ELEMENTNODES) :: HS, ElementNodalDensityLiquidL, ElementNodalDensitySolidL
            real(REAL_TYPE), dimension(ELEMENTNODES) :: ForcesX, ForcesY, ForcesZ, ParticleShape
            real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJac, InvRJac
            real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: DHS, CartD
            real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: CartDtemp
            real(REAL_TYPE) :: ConstDensitySolid, ConstDensityLiquid
            real(REAL_TYPE) :: TotIntWeightLiqMP
            
            if (.not.(NFORMULATION==2)) RETURN

            This%InteractionForceLiquid = 0.0
            This%InteractionForceSolid = 0.0
            
              !--- Calculate  density  ----
              ConstDensityLiquid = 0.0
              ConstDensitySolid = 0.0
                
              do J = 1, Counters%NLayers ! loop over all material points in element
                if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                    ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
                else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-solid'.or.MatParams(j)%MaterialPhases=='1-phase-solid') then
                    ConstDensitySolid = (MatParams(J)%DensitySolid/1000)    
                end if
              end do 
              ! ---------------------------------------------------
            

            do I = 1, Counters%NAEl ! loop over all active elements
            
              IElement = ActiveElement(I)
 
              if((This%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID).or. &
                   (This%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID)) then
                  

                NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
                
                  Nix = ReducedDof(NodeIDs) + 1
                  Niy = Nix + 1
                if (NVECTOR == 3) then ! 3D case  
                  Niz = Nix + 2
                end if   
                ElementNodalDensityLiquidL = This%Nodes(NodeIDs)%DensityLiquidL
                ElementNodalDensitySolidL = This%Nodes(NodeIDs)%DensitySolidL

                call GaussPointLocalCoordinates(1, WeiGP, PosGP)
                call ShapeFunctionData(PosGP, ELEMENTNODES, HS, DHS) !shape function derivatives dHS for LocPos
                call DetJacob(PosGP, Counters%NEl, Counters%NodTot, NVECTOR,  &
                              IElement, ElementConnectivities, NodalCoordinates, RJac, InvRJac, DetJac)  ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix InvRJac
                WeiGP = WeiGP * DetJac 
                
                ! Calculate the cartesian derivatives and transpose
                CartDtemp = 0.0
                do J = 1, NVECTOR
                  do II = 1, ELEMENTNODES
                    do K = 1, NVECTOR
                      CartDtemp(K, II) = CartDtemp(K,II) + InvRJac(K, J) * DHS(II, J)
                    end do
                  end do
                end do
                CartD = TRANSPOSE(CartDtemp)
                
                GradientNf = 0.0
                GradientNs = 0.0
                do J = 1, ELEMENTNODES
                  do K = 1, NVECTOR
                      GradientNf(K) = GradientNf(K) + CartD(J, K) * ElementNodalDensityLiquidL(J)
                      GradientNs(K) = GradientNs(K) + CartD(J, K) * ElementNodalDensitySolidL(J)
                  end do
                end do
                if(ConstDensityLiquid>0.0) then
                   GradientNf = GradientNf / ConstDensityLiquid
                else
                   GradientNf = 0.0
                end if
                
                if(ConstDensitySolid>0.0) then
                    GradientNs = GradientNs / ConstDensitySolid
                else 
                    GradientNs = 0.0
                end if
                
                This%Elements(IElement)%GradientPorosityL = -GradientNs
                  
                  Weighting = 0.0
                  AverageSigma = 0.0
                  VolumeMP = 0.0
                  TotIntWeightLiqMP = 0.0
                  do J = 1, NPartEle(IElement)
                      ParticleIndex = GetParticleIndex(J, IElement)
                      if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                          TotIntWeightLiqMP = TotIntWeightLiqMP + Particles(ParticleIndex)%IntegrationWeight
                      end if
                  end do
                  
                  do J = 1, NPartEle(IElement)
                    ParticleIndex = GetParticleIndex(J, IElement)
                    if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                      AverageSigma = SigmaEffArray(ParticleIndex,:)

                      VolumeMP = Particles(ParticleIndex)%IntegrationWeight
                      if(TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL>0.0) then                        
                         VolumeMP = ElementSpace(IElement) / TotIntWeightLiqMP * Particles(ParticleIndex)%IntegrationWeight
                         if(TwoLayerData%Elements(IElement)%IsBoundaryOfLiquidElements) then ! Boundary of Liquid domain
                          if((TotIntWeightLiqMP/ElementSpace(IElement))<CalParams%RequiredDegreeOfFilling) then
                              VolumeMP = Particles(ParticleIndex)%IntegrationWeight 
                          end if
                         end if
                      end if
                      ParticleShape = HS !use values of the ShapeFunction at GP

                    if (NVECTOR == 3) then ! 3D case
                      ! liquid interaction forces
                      Term(1) = GradientNf(1) * AverageSigma(1) + GradientNf(2) * AverageSigma(4) + GradientNf(3) * AverageSigma(6)
                      Term(2) = GradientNf(1) * AverageSigma(4) + GradientNf(2) * AverageSigma(2) + GradientNf(3) * AverageSigma(5)
                      Term(3) = GradientNf(1) * AverageSigma(6) + GradientNf(2) * AverageSigma(5) + GradientNf(3) * AverageSigma(3)
                  
                      ForcesX(:) = ParticleShape(:) * Term(1) * VolumeMP !WeiGP
                      ForcesY(:) = ParticleShape(:) * Term(2) * VolumeMP !WeiGP
                      ForcesZ(:) = ParticleShape(:) * Term(3) * VolumeMP !WeiGP
     
                      This%InteractionForceLiquid(Nix, 1) = This%InteractionForceLiquid(Nix, 1) - ForcesX
                      This%InteractionForceLiquid(Niy, 1) = This%InteractionForceLiquid(Niy, 1) - ForcesY
                      This%InteractionForceLiquid(Niz, 1) = This%InteractionForceLiquid(Niz, 1) - ForcesZ
                  
                      ! solid interaction forces
                      Term(1) = GradientNs(1) * AverageSigma(1) + GradientNs(2) * AverageSigma(4) + GradientNs(3) * AverageSigma(6)
                      Term(2) = GradientNs(1) * AverageSigma(4) + GradientNs(2) * AverageSigma(2) + GradientNs(3) * AverageSigma(5)
                      Term(3) = GradientNs(1) * AverageSigma(6) + GradientNs(2) * AverageSigma(5) + GradientNs(3) * AverageSigma(3)
                  
                      ForcesX(:) = ParticleShape(:) * Term(1) * VolumeMP !WeiGP
                      ForcesY(:) = ParticleShape(:) * Term(2) * VolumeMP !WeiGP
                      ForcesZ(:) = ParticleShape(:) * Term(3) * VolumeMP !WeiGP
                  
                      This%InteractionForceSolid(Nix, 1) = This%InteractionForceSolid(Nix, 1) - ForcesX
                      This%InteractionForceSolid(Niy, 1) = This%InteractionForceSolid(Niy, 1) - ForcesY
                      This%InteractionForceSolid(Niz, 1) = This%InteractionForceSolid(Niz, 1) - ForcesZ
                     
                    elseif (NVECTOR == 2) then ! 2D case plane strain?
                      ! liquid interaction forces
                      Term(1) = GradientNf(1) * AverageSigma(1) + GradientNf(2) * AverageSigma(4)
                      Term(2) = GradientNf(1) * AverageSigma(4) + GradientNf(2) * AverageSigma(2) 
                  
                      ForcesX(:) = ParticleShape(:) * Term(1) * VolumeMP !WeiGP
                      ForcesY(:) = ParticleShape(:) * Term(2) * VolumeMP !WeiGP
     
                      This%InteractionForceLiquid(Nix, 1) = This%InteractionForceLiquid(Nix, 1) - ForcesX
                      This%InteractionForceLiquid(Niy, 1) = This%InteractionForceLiquid(Niy, 1) - ForcesY
                  
                      ! solid interaction forces
                      Term(1) = GradientNs(1) * AverageSigma(1) + GradientNs(2) * AverageSigma(4)
                      Term(2) = GradientNs(1) * AverageSigma(4) + GradientNs(2) * AverageSigma(2)
                  
                      ForcesX(:) = ParticleShape(:) * Term(1) * VolumeMP !WeiGP
                      ForcesY(:) = ParticleShape(:) * Term(2) * VolumeMP !WeiGP

                      This%InteractionForceSolid(Nix, 1) = This%InteractionForceSolid(Nix, 1) - ForcesX
                      This%InteractionForceSolid(Niy, 1) = This%InteractionForceSolid(Niy, 1) - ForcesY
                    end if  
                    
                  end if
                end do
                  
              end if
   
            end do ! element loop
            
        end subroutine DetermineInteractionForces
  
        
        subroutine DetermineTwoLayerStatus(This)
        !**********************************************************************
        !
        !    Function: Determines phase status   (PhaseStatusSOLID or PhaseStatusLIQUID)
        !
        !**********************************************************************
        
        implicit none
        
          class(TwoLayerType), intent(inout) :: This
          
          ! local variables
          integer(INTEGER_TYPE) :: I, J, ParticleIndex
          real(REAL_TYPE) :: EVol, MassGrainsInContact, MassGrainsSeparated
          logical :: SolidGrainsInContact
          
        
          if (.not.(NFORMULATION==2)) RETURN
          if (.not.(CalParams%NumberOfPhases==2)) RETURN
                  
          do I = 1, Counters%NEl ! loop over all elements
              
            MassGrainsInContact = 0.0
            MassGrainsSeparated = 0.0
            SolidGrainsInContact = .false.
            TwoLayerData%Elements(I)%IsSolidMPwithPhaseLiquidInElement = .false.
          
            do J = 1, NPartEle(I) ! loop over all material points in element
              ParticleIndex = GetParticleIndex(J, I)
              
              !!---------------------------------- SOLID MATERIAL POINT ------------------------------!!            
              if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then ! SOLID Material Point
                Particles(ParticleIndex)%LiquidFreeSurface = 0.0

                  EVol =  (GetEpsStepI(Particles(ParticleIndex), 1) +  GetEpsStepI(Particles(ParticleIndex), 2) +  &
                                   GetEpsStepI(Particles(ParticleIndex), 3))
     
                  if ((Particles(ParticleIndex)%EffPorosity<CalParams%LimitPorosity).or. &
                      ((Particles(ParticleIndex)%EffPorosity==CalParams%LimitPorosity).and.(EVol<0.0))) then   
                    Particles(ParticleIndex)%PhaseStatus = PhaseStatusSOLID
                    MassGrainsInContact = MassGrainsInContact + MassArray(ParticleIndex)
                  else
                    Particles(ParticleIndex)%PhaseStatus = PhaseStatusLIQUID
                    MassGrainsSeparated = MassGrainsSeparated + MassArray(ParticleIndex)
                    TwoLayerData%Elements(I)%IsSolidMPwithPhaseLiquidInElement = .true.
                  end if
               end if ! solid points
                
              end do ! material point loop
          
     
                SolidGrainsInContact = ((MassGrainsInContact>0.0).or.(MassGrainsSeparated>0.0)).and. &
                             (MassGrainsInContact>MassGrainsSeparated)
          
              !!---------------------------------- LIQUID MATERIAL POINT ------------------------------!!  
              do J = 1, NPartEle(I) ! loop over all material points in element
                ParticleIndex = GetParticleIndex(J, I)
                if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then ! LIQUID Material Point
                  if (This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
                    Particles(ParticleIndex)%PhaseStatus = PhaseStatusLIQUID
                  else if(This%Elements(I)%ContainedMaterialTypes==ContainedMaterialTypeSOLIDLIQUID) then
     
                      if(SolidGrainsInContact) then ! Solid In Contact 
                        Particles(ParticleIndex)%PhaseStatus = PhaseStatusSOLID
                      else if(.not.SolidGrainsInContact) then ! Solid Separated                         
                        Particles(ParticleIndex)%PhaseStatus = PhaseStatusLIQUID
                      end if
                  end if
                end if
              end do    ! loop over all material points in element
              
            
          end do ! element loop
                     
        end subroutine DetermineTwoLayerStatus
 
        
        subroutine ComputeTwoLayerVolumetricStrainLiquid(DULiquid, DUSolid)
        !**********************************************************************
        !
        !    Function: compute volumetric strain of the liquid for double-point formulation 
        !
        !**********************************************************************        
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity), intent(in) :: DULiquid, DUSolid
          
          !! local variables
          integer(INTEGER_TYPE) :: IAElement, IElement, NElemPart, IParticle, ParticleIndex, I, J, K, NN, NAdjacentElements, IEl, II
          integer(INTEGER_TYPE) :: Ix, Nx, idx
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs, Ni
          real(REAL_TYPE) :: DetJac, WeiGP, DU_LiqSolTerm, PorosityValue
          real(REAL_TYPE) :: absValDistance
          real(REAL_TYPE), dimension(NTENSOR) :: EpsV, EpsW
          real(REAL_TYPE), dimension(NVECTOR) :: PosGP, DU_Liq, DU_LiqSolGP, GradN, LocPosition, GradientDu
          real(REAL_TYPE), dimension(NVECTOR) :: DU_Sol, DU_LiqGP, DU_SolGP, Distance
          real(REAL_TYPE), dimension(NVECTOR, NVECTOR) :: RJac, InvRJac
          real(REAL_TYPE), dimension(ELEMENTNODES) :: HS
          real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: dHS, CartD, ElementNodalSolidIncrDispl
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: CartDtemp
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE), dimension(Counters%N, Counters%NEntity) :: InterpDUSolid, TrialDUSolid
          real(REAL_TYPE), dimension(Counters%NodTot) :: NodalVolumeS, SumInvNormDistance
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR) :: NodalGradDU
          real(REAL_TYPE), dimension(Counters%NEl, NVECTOR) :: GradDU
          logical, dimension(Counters%NodTot) :: IsNodeAlreadyConsidered, IsNodeToGetSolidVelocity
          
          if (.not.(NFORMULATION==2)) RETURN
          
          TrialDUSolid = DUSolid
          InterpDUSolid = 0.0
          NodalGradDU = 0.0
          NodalVolumeS = 0.0
          SumInvNormDistance = 0.0
          GradDU = 0.0
          IsNodeAlreadyConsidered = .false.
          IsNodeToGetSolidVelocity = .false.
          
          ! Determine gradient of incremental displacement in element of boundary solid domain
          do IAElement = 1, Counters%NAEl ! loop over all active elements
            IElement = ActiveElement(IAElement)
            if(TwoLayerData%Elements(IElement)%IsBoundaryOfSolidElements) then
               NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)

               do I = 1, NVECTOR
                 Ni = ReducedDof(NodeIDs) + I ! global storage coordinate of i-val
                 ElementNodalSolidIncrDispl(:, i) = DUSolid(Ni, 1)
               end do
               
               call GaussPointLocalCoordinates(1, WeiGP, PosGP)
               call ShapeFunctionData(PosGP, ELEMENTNODES, HS, DHS) !shape function derivatives dHS for LocPos
               call DetJacob(PosGP, Counters%NEl, Counters%NodTot, NVECTOR,  &
                              IElement, ElementConnectivities, NodalCoordinates, RJac, InvRJac, DetJac)  ! Calculate Jacobian matrix RJac and the inverse of the Jacobian matrix InvRJac
               WeiGP = WeiGP * DetJac 
                
               ! Calculate the cartesian derivatives and transpose
               CartDtemp = 0.0
               do J = 1, NVECTOR
                do II = 1, ELEMENTNODES
                  do K = 1, NVECTOR
                    CartDtemp(K, II) = CartDtemp(K,II) + InvRJac(K, J) * DHS(II, J)
                  end do
                end do
               end do
               CartD = TRANSPOSE(CartDtemp)
               
               GradientDU = 0.0
               do J = 1, ELEMENTNODES
                  do i = 1, NVECTOR
                    GradientDu(i) = GradientDu(i) + CartD(J, i) * ElementNodalSolidIncrDispl(J, i)
                  end do
               end do

               do i = 1, NVECTOR
                GradDU(IElement, i) = GradientDu(i)
               end do
            end if
          end do
          

          ! Determine gradient of incremental displacement at node of boundary solid domain
          do I= 1, Counters%NodTot
                if(TwoLayerData%Nodes(I)%IsNodeOfBoundSolidDomain) then
                  NAdjacentElements = GetNElmOfNode(I)
                  do J = 1, NAdjacentElements
                    IEl = GetElmIOfNode(I, J)
                    if (IsActiveElement(IEl)) then
                      if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes/=ContainedMaterialTypeLIQUID) then
                         do idx = 1, NVECTOR
                           NodalGradDu(i, idx) = NodalGradDU(I, idx) + GradDU(IEl, idx) * ElementSpace(IEl)
                         end do
                      end if
                    end if
                  end do
                end if
          end do
          
          do I= 1, Counters%NodTot
                if(TwoLayerData%Nodes(I)%IsNodeOfBoundSolidDomain) then
                  NAdjacentElements = GetNElmOfNode(I)
                  do J = 1, NAdjacentElements
                    IEl = GetElmIOfNode(I, J)
                    if (IsActiveElement(IEl)) then
                      if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes/=ContainedMaterialTypeLIQUID) then
                         NodalVolumeS(I) =  NodalVolumeS(I) + ElementSpace(IEl)
                      end if
                    end if
                  end do
                end if
          end do   
          
          do I= 1, Counters%NodTot
                if(TwoLayerData%Nodes(I)%IsNodeOfBoundSolidDomain) then
                  do idx = 1, NVECTOR
                    NodalGradDu(i, idx) = NodalGradDU(I, idx) / NodalVolumeS(I)
                  end do
                end if
          end do  
          
          ! Determine IsNodeToGetSolidVelocity 
          do I= 1, Counters%NodTot
                if(TwoLayerData%Nodes(I)%IsNodeOfBoundSolidDomain) then
                  NAdjacentElements = GetNElmOfNode(I)
                  do J = 1, NAdjacentElements
                    IEl = GetElmIOfNode(I, J)
                    if (IsActiveElement(IEl)) then
                        if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
                            do K = 1,ELEMENTNODES
                                NN = ElementConnectivities(K,IEl)
                                if(.not.TwoLayerData%Nodes(NN)%IsNodeOfBoundSolidDomain) then
                                    IsNodeToGetSolidVelocity(NN) = .true.
                                end if
                            enddo
                        end if
                    end if
                  end do
                end if
          end do   
          
          ! Determine Increment solid displacement at node IsNodeToGetSolidVelocity
          do I= 1, Counters%NodTot
                if(IsNodeToGetSolidVelocity(I)) then
                  Ix = ReducedDof(I)
                  NAdjacentElements = GetNElmOfNode(I)
                  do J = 1, NAdjacentElements
                    IEl = GetElmIOfNode(I, J)
                    if (IsActiveElement(IEl)) then
                        if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) then
                            do K = 1,ELEMENTNODES
                                NN = ElementConnectivities(K,IEl)
                                if(TwoLayerData%Nodes(NN)%IsNodeOfBoundSolidDomain.and. &
                                  (.not.IsNodeAlreadyConsidered(NN))) then
                                  IsNodeAlreadyConsidered(NN) = .true.
                                  Distance = NodalCoordinates(I,:) - NodalCoordinates(NN,:)
                                  absValDistance = DotProduct(Distance, Distance, 3) ! generalise dimension
                                  Nx = ReducedDof(NN)

                                  do idx = 1, NVECTOR
                                    InterpDUSolid(Ix + idx, 1) = InterpDUSolid(Ix + idx, 1) + (DUSolid(Nx + idx, 1) + NodalGradDU(NN, idx)*Distance(idx)) / absValDistance
                                  end do
                                  
                                  
                                  SumInvNormDistance(I) = SumInvNormDistance(I) + 1.0 / absValDistance
                                end if
                             enddo
                             do K = 1,ELEMENTNODES 
                                 NN = ElementConnectivities(K,IEl)
                                 IsNodeAlreadyConsidered(NN) = .false. !set to false for a new point I
                             end do
                        end if
                    end if
                  end do
                end if
      end do       
      
      do I= 1, Counters%NodTot
             if(IsNodeToGetSolidVelocity(I)) then
                Ix = ReducedDof(I)
                do idx = 1, NVECTOR
                  InterpDUSolid(Ix + idx, 1) = InterpDUSolid(Ix + idx,1) / SumInvNormDistance(I)
                end do  
             end if
      end do
      
          do IAElement = 1, Counters%NAEl ! loop over all active elements
            IElement = ActiveElement(IAElement)
            
            call GaussPointLocalCoordinates(1, WeiGP, PosGP)
   
            
            NElemPart = NPartEle(IElement)
            if ( ISAXISYMMETRIC .and. .not.IsParticleIntegration(IElement) ) then
                NElemPart = ELEMENTGAUSSPOINTS ! Number of Gauss points per element
            end if
            
            do IParticle = 1, NElemPart ! loop over all material points in element
              
                ParticleIndex = GetParticleIndex(IParticle, IElement)
                
                if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then ! only LIQUID material point
                    
                    LocPosition = PosGP
                    PorosityValue = (1.0 - TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL)
                    
                  call ShapeFunctionData(LocPosition, ELEMENTNODES, HS, DHS) !shape function derivatives dHS for LocPos
                  
                  call BMatrix(LocPosition, ELEMENTNODES, Counters%NEl, Counters%NodTot, NDIM, &
                             IElement, ElementConnectivities, NodalCoordinatesUpd, BMatrixDeformed, DetJac)
     
                  call Get_Strain(IElement, IParticle, ElementConnectivities, &
                                BMatrixDeformed, DULiquid(1:Counters%N, 1), ReducedDof, EpsW)
     
                  !------------------------------------------------------
                  DU_Liq = 0.0
                  DU_Sol = 0.0
                  DU_LiqGP = 0.0
                  DU_SolGP = 0.0
                  DU_LiqSolGP = 0.0
                  DU_LiqSolTerm = 0.0
                  
                  do J= 1, ELEMENTNODES ! loop over element nodes
                    NN = ElementConnectivities(J,IElement)
                    do idx = 1, NVECTOR
                      DU_Liq(idx) = DULiquid(ReducedDof(NN) + idx,1)
                    end do

                    DU_LiqGP = DU_LiqGP + DU_Liq * HS(J) ! velocity at GP location

                    do idx = 1, NVECTOR
                      DU_Sol(idx) = DUSolid(ReducedDof(NN) + idx,1)
                    end do
                    
                    if ((IsNodeToGetSolidVelocity(NN)).and. &
                       ((TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID).and. &
                       (TwoLayerData%Elements(IElement)%ConcentrationRatioSolidL>0.0))) then
                        
                       do idx = 1, NVECTOR
                        if(PBoundary(ReducedDof(NN) + idx) /= 0.0) then
                            DU_Sol(idx) = InterpDUSolid(ReducedDof(NN) + idx,1)
                            TrialDUSolid(ReducedDof(NN) + idx, 1) = InterpDUSolid(ReducedDof(NN) + idx, 1)
                        end if
                       end do
                    end if
                    
                      DU_SolGP = DU_SolGP + DU_Sol * HS(J) ! velocity at GP location
                  end do
        
                  DU_LiqSolGP = DU_LiqGP - DU_SolGP
                  GradN = TwoLayerData%Elements(IElement)%GradientPorosityL

                  DU_LiqSolTerm = DotProduct(DU_LiqSolGP, GradN, NVECTOR) ! generalise dimension
                  
                  !------------------------------------------------------
                  
                  if((TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID).and. &
                      (TwoLayerData%Elements(IElement)%ConcentrationRatioSolidl>0.0)) then
                       call Get_Strain(IElement, IParticle, ElementConnectivities, &
                                BMatrixDeformed, TrialDUSolid(1:Counters%N, 1), ReducedDof, EpsV)
                  else 
                      call Get_Strain(IElement, IParticle, ElementConnectivities, &
                                BMatrixDeformed, DUSolid(1:Counters%N, 1), ReducedDof, EpsV)
                  end if
     
     
                    Particles(ParticleIndex)%WaterVolumetricStrain =  &
                      DU_LiqSolTerm  + &
                      (1.0-PorosityValue) * (EpsV(1) + EpsV(2) + EpsV(3)) + &
                      PorosityValue * (EpsW(1) + EpsW(2) + EpsW(3))
      
                  Particles(ParticleIndex)%WaterVolumetricStrain =  &
                      Particles(ParticleIndex)%WaterVolumetricStrain /  &
                      PorosityValue

                end if ! only LIQUID Material Point
              end do ! material point loop
              


          end  do ! element loop
         
        end subroutine ComputeTwoLayerVolumetricStrainLiquid

        
        subroutine DetermineDensityField(NodalDensity, MaterialPointType, IsLiquidModel, OnlyElementsWLiquidMPs)
        !**********************************************************************
        !
        !    Function: Determines density field
        !
        !**********************************************************************        
        implicit none
        
          real(REAL_TYPE), dimension(:), intent(inout) :: NodalDensity
          integer(INTEGER_TYPE), intent(in) :: MaterialPointType
          logical, intent(in) :: IsLiquidModel, OnlyElementsWLiquidMPs
          
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, I, J, ParticleIndex, NAdjacentElements
          logical, dimension(Counters%NodTot) :: LocalActiveNode
          real(REAL_TYPE) :: ConsideredMass, ConstDensityLiquid
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
          
          ! Determine nodal mass from material points of activated elements connected to respective nodes
          NodalDensity = 0.0
          LocalActiveNode = .false.
          
          !--- Calculate  density  ----
          ConstDensityLiquid = 0.0
                
          do J = 1, Counters%NLayers ! loop over all material points in element
              if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                  ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
              else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                  ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)  
              end if
          end do 
          ! ---------------------------------------------------
          
          do IAEl = 1, Counters%NAEl ! loop over all active elements
            IEl = ActiveElement(IAEl)
            
            if(OnlyElementsWLiquidMPs) then
              if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) CYCLE
            else
              if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID) CYCLE
            end if
            
            NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl)
            LocalActiveNode(NodeIDs) = .true.
          
            do I = 1, NPartEle(IEl) ! loop over all material points in element
              ParticleIndex = GetParticleIndex(I, IEl)
              
              ConsideredMass = 0.0
              
              if (IsLiquidModel) then ! Temporary solution
                if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                  ConsideredMass = MassArray(ParticleIndex)
                else
                  ConsideredMass = ConstDensityLiquid * Particles(ParticleIndex)%IntegrationWeight ! Replace solid density by liquid density
                end if
              else
                if (MaterialPointTypeArray(ParticleIndex)==MaterialPointType) then
                  ConsideredMass = MassArray(ParticleIndex)
                end if
              end if
              
              NodalDensity(NodeIDs) = NodalDensity(NodeIDs) + ConsideredMass * ShapeValuesArray(ParticleIndex,:)
              
            end do ! material point loop
            
          end do ! element loop
           
          ActiveNodeElementVolume = 0.0
          do I = 1, Counters%NodTot
            if (LocalActiveNode(I)) then
              NAdjacentElements = GetNElmOfNode(I)
              do J = 1, NAdjacentElements
                IEl = GetElmIOfNode(I, J)
                if (IsActiveElement(IEl)) then
                    if(OnlyElementsWLiquidMPs) then
                      if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes/=ContainedMaterialTypeSOLID) then                        
                        ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) + ElementSpace(IEl)
                      end if
                    else
                      if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes/=ContainedMaterialTypeLIQUID) then
                        ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) + ElementSpace(IEl)
                      end if
                    end if
                end if
              end do
            end if
          end do
          
          do I = 1, Counters%NodTot
            if (LocalActiveNode(I)) then
              ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) / dble(ELEMENTNODES)
            end if
          end do
          
          do I = 1, Counters%NodTot
            if (LocalActiveNode(I)) then
                if(ActiveNodeElementVolume(I)>0.0) then
                    NodalDensity(I) = NodalDensity(I) / ActiveNodeElementVolume(I)
                end if
            end if
          end do
                      
      end subroutine DetermineDensityField
      
      
      subroutine DetermineFillingRatioField()
        !**********************************************************************
        !
        !    Function: Determines filling ratio
        !
        !**********************************************************************        
        implicit none
        
          ! local variables
          integer(INTEGER_TYPE) :: IAEl, IEl, I, J, ParticleIndex, NAdjacentElements
          logical, dimension(Counters%NodTot) :: LocalActiveNode
          real(REAL_TYPE) :: ConsideredVolume, ConstDensityLiquid, ParticleConstVolume
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
          
          ! Determine nodal mass from material points of activated elements connected to respective nodes
          TwoLayerData%Nodes%FillingRatio = 0.0
          LocalActiveNode = .false.
          
          !--- Calculate  density  ----
          ConstDensityLiquid = 0.0
                
          do J = 1, Counters%NLayers ! loop over all material points in element
              if(trim(MatParams(J)%MATERIALTYPE)=='2-phase'.or.MatParams(j)%MaterialPhases=='2-phase') then
                  ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)
              else if(trim(MatParams(J)%MATERIALTYPE)=='1-phase-liquid'.or.MatParams(j)%MaterialPhases=='1-phase-liquid') then
                  ConstDensityLiquid = (MatParams(J)%DensityLiquid/1000)  
              end if
          end do 
          ! ---------------------------------------------------
          
          do IAEl = 1, Counters%NAEl ! loop over all active elements
            IEl = ActiveElement(IAEl)
            
            if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes==ContainedMaterialTypeSOLID) CYCLE
            
            NodeIDs = ElementConnectivities(1:ELEMENTNODES, IEl)
            LocalActiveNode(NodeIDs) = .true.
            
            do I = 1, NPartEle(IEl) ! loop over all material points in element
              ParticleIndex = GetParticleIndex(I, IEl)
              
              ParticleConstVolume = 0.0
              
              if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                  ParticleConstVolume = MassArray(ParticleIndex) / ConstDensityLiquid
              end if
              
              TwoLayerData%Nodes(NodeIDs)%FillingRatio = TwoLayerData%Nodes(NodeIDs)%FillingRatio +  &
                       ParticleConstVolume * ShapeValuesArray(ParticleIndex,:)
              
            end do ! material point loop
          
            ConsideredVolume = (TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL) * ElementSpace(IEl)/ dble(ELEMENTNODES)
          
            TwoLayerData%Nodes(NodeIDs)%FillingRatio = TwoLayerData%Nodes(NodeIDs)%FillingRatio +  &
                          ConsideredVolume     
            
          end do ! element loop
           
          
          ActiveNodeElementVolume = 0.0
          do I = 1, Counters%NodTot
            if (LocalActiveNode(I)) then
              NAdjacentElements = GetNElmOfNode(I)
              do J = 1, NAdjacentElements
                IEl = GetElmIOfNode(I, J)
                if (IsActiveElement(IEl)) then
                  if(TwoLayerData%Elements(IEl)%ContainedMaterialTypes/=ContainedMaterialTypeSOLID) then                        
                      ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) + ElementSpace(IEl)
                  end if
                end if
              end do
            end if
          end do
          
          do I = 1, Counters%NodTot
            if (LocalActiveNode(I)) then
              ActiveNodeElementVolume(I) = ActiveNodeElementVolume(I) / dble(ELEMENTNODES)
            end if
          end do
          
          do I = 1, Counters%NodTot
            if (LocalActiveNode(I)) then
                if(ActiveNodeElementVolume(I)>0.0) then
                    TwoLayerData%Nodes(I)%FillingRatio = TwoLayerData%Nodes(I)%FillingRatio / ActiveNodeElementVolume(I)
                end if
            end if
          end do
                      
      end subroutine DetermineFillingRatioField
      
         subroutine DetermineBoundaryElementLiquidDomain2LayerForm()
        !**********************************************************************
        !
        !  Function:  Determines :
        !             1) Nodes which are at the boundary between elements which contain 
        !                SOLID+LIQUID and LIQUID MPs, and elements which contain only SOLID MPs or are empty, 
        !                -> then these nodes are stored in a flag IsNodeOfBoundLiquidDomain
        !             2) If element contains SOLID+LIQUID and LIQUID MPs, 
        !                and has one of the nodes stored in IsNodeOfBoundLiquidDomain, 
        !                -> then the element is stored in a flag IsBoundaryOfLiquidElements
        !
        !**********************************************************************
        
        implicit none 
        
          ! Local variables
          integer(INTEGER_TYPE) :: IElement, I, IDNode
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
          integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NodeOfLiquidElm, NodeOfNonLiquidElm, FreeSurfNodes 
          
          ! Initialize
          NodeOfLiquidElm = 0
          NodeOfNonLiquidElm = 0
          FreeSurfNodes = 0
          
          do IElement = 1, Counters%NEl 
              
            TwoLayerData%Elements(IElement)%IsBoundaryOfLiquidElements = .false. ! Initialize
            
            if((IsActiveElement(IElement)).and. &
               (.not.(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeSOLID))) then ! Loop over all initially active elements with liquid material points
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              NodeOfLiquidElm(NodeIDs) = 1
            else !Inactive element or element without liquid MPs
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              NodeOfNonLiquidElm(NodeIDs) = 1  ! set to 1 each node of inactive elements
            end if

          end do ! elements
          
          !Determine Free Surface Nodes as the ones which have value 2
          FreeSurfNodes = NodeOfLiquidElm + NodeOfNonLiquidElm
          
          
          ! Determine Boundary Elements
          do IElement = 1, Counters%NEl 
            if((IsActiveElement(IElement)).and. &
               (.not.(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeSOLID))) then ! Loop over all initially active elements with Liquid Material Points
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              do I = 1, ELEMENTNODES
                  IDNode = NodeIDs(I)
                  if(FreeSurfNodes(IDNode)==2) then
                      TwoLayerData%Elements(IElement)%IsBoundaryOfLiquidElements = .true.
                  end if
              end do
            end if

          end do ! elements
             
        
        end subroutine DetermineBoundaryElementLiquidDomain2LayerForm
      
      
      subroutine DetermineBoundaryElementSolidDomain2LayerForm()
        !**********************************************************************
        !
        !  Function:  Determines :
        !             1) Nodes which are at the boundary of Elements which contain 
        !                SOLID+LIQUID and SOLID MPs, and elements which contain only LIQUID MPs)-> 
        !                -> then these nodes are stored in a flag IsNodeOfBoundSolidDomain
        !             2) If element contains SOLID+LIQUID and SOLID MPs, 
        !                and has one of the nodes stored in IsNodeOfBoundSolidDomain, 
        !                -> then the element is stored in a flag IsBoundaryOfSolidElements
        !
        !**********************************************************************
        
        implicit none 
        
          ! local variables
          integer(INTEGER_TYPE) :: IElement, I, IDNode
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: NodeIDs
          integer(INTEGER_TYPE), dimension(Counters%NodTot) :: NodeOfSolidElm, NodeOfNonSolidElm, FreeSurfNodes 
          
          ! Initialize
          NodeOfSolidElm = 0
          NodeOfNonSolidElm = 0
          FreeSurfNodes = 0
          
          do IElement = 1, Counters%NEl 
              
            TwoLayerData%Elements(IElement)%IsBoundaryOfSolidElements = .false. ! Initialize
            
            if((IsActiveElement(IElement)).and. &
               (.not.(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID))) then ! Loop over all initially active elements with Liquid Material Points
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              NodeOfSolidElm(NodeIDs) = 1
            else !Inactive element or element without liquid MPs
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              NodeOfNonSolidElm(NodeIDs) = 1  ! set to 1 each node of inactive elements
            end if
            TwoLayerData%Nodes(NodeIDs)%IsNodeOfBoundSolidDomain = .false. ! initialize as false

          end do ! elements
          
          ! Determine Free Surface Nodes as the ones which have value 2
          FreeSurfNodes = NodeOfSolidElm + NodeOfNonSolidElm
          
          
          ! Determine Boundary Elements and nodes at the boundary
          do IElement = 1, Counters%NEl 
            if((IsActiveElement(IElement)).and. &
               (.not.(TwoLayerData%Elements(IElement)%ContainedMaterialTypes==ContainedMaterialTypeLIQUID))) then ! Loop over all initially active elements with Liquid Material Points
              NodeIDs = ElementConnectivities(1:ELEMENTNODES, IElement)
              do I = 1, ELEMENTNODES
                  IDNode = NodeIDs(I)
                  if(FreeSurfNodes(IDNode)==2) then
                      TwoLayerData%Elements(IElement)%IsBoundaryOfSolidElements = .true.
                      TwoLayerData%Nodes(IDNode)%IsNodeOfBoundSolidDomain = .true.
                  end if
              end do
            end if

          end do ! elements
             
        
        end subroutine DetermineBoundaryElementSolidDomain2LayerForm

        subroutine AdjustMassLiquidMP()
        !**********************************************************************
        !
        !    Function: Modifies the mass of liquid material points
        !
        !**********************************************************************
        
        
        ! local
        integer(INTEGER_TYPE) :: IAEl, IEl,  IParticle, ParticleIndex

        
        
        if(NFORMULATION==1) RETURN
        
        do IAEl = 1, Counters%NAEl ! loop over all active elements
          IEl = ActiveElement(IAEl)
              
          if(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0) then
           
            do IParticle = 1, NPartEle(IEl) ! loop over total number of material points 
              ParticleIndex = GetParticleIndex(IParticle, IEl)
        
              ! Modify mass liquid material point
              if (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then !

                 MassArray(ParticleIndex) = (1.0 - Particles(ParticleIndex)%ConcentrationRatioSolidL) * &
                           Particles(ParticleIndex)%ConstDensity *   Particles(ParticleIndex)%IntegrationWeight
     
                 Particles(ParticleIndex)%Density = MassArray(ParticleIndex) / Particles(ParticleIndex)%IntegrationWeight
                 
                 Particles(ParticleIndex)%FBody = MassArray(ParticleIndex) *  &
                           CalParams%GravityData%GAccel * CalParams%GravityData%GravityVector
                 
              end if
            end do
          end if
        end do
            
                
                
        end subroutine AdjustMassLiquidMP
      
      end module ModTwoLayerFormulation
