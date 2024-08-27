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
	  
	  
	  module ModStrainSmoothing
      !**********************************************************************
      !
      !    Function:  This module contains all routines which are used for both quasi-static 
      !               and dynamic MPM with low-order elements for smoothing of strains to
      !               overcome locking.
      !
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               strain smoothing.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************

      use ModCounters
      use ModMeshAdjacencies
      use ModMPMData
      use ModMeshInfo
      use ModTwoLayerFormulation
        
      implicit none

      contains ! Routines of this module
      
            
    
        subroutine SmoothenStrains(DUTot, IMatSet)
        !**********************************************************************
        !
        !    Function:  Smoothens Gauss point strains for low-order elements.
        !
        !     DUTot : Incremental nodal displacements
        !     IMatSet : The no. of the material set
        !
        !**********************************************************************
       
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: DUTot
          integer(INTEGER_TYPE), intent (in) :: IMatSet

          call GetEnhancedVolumetricStrain(DUTot, IMatSet)
     
          call SmoothenVolumetricStrain(IMatSet)
          
        end subroutine SmoothenStrains
          
        subroutine GetEnhancedVolumetricStrain(DUTot, IMatSet)
        !**********************************************************************
        !
        !    Function: Returns enhanced volumetric strains
        !
        !     DUTotAll : Incremental nodal displacements of all entities
        !     IMatSet : The no. of the material set
        !
        !**********************************************************************
       
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: DUTot
          integer(INTEGER_TYPE), intent(in) :: IMatSet
          
          ! Local variables  
          real(REAL_TYPE), dimension(NTENSOR) :: Eps
          real(REAL_TYPE),  dimension(NVECTOR, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac, EV, EVD
          real(REAL_TYPE), dimension(NVECTOR) :: LI
          real(REAL_TYPE), dimension(ELEMENTNODES) :: volumeElNodes ! added for axisymmetric analysis
          integer(INTEGER_TYPE) :: IEl, IAEl, I, IParticle
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
          logical :: DoSmoothing2Layer, DoSmoothing

          LI = 0.25D0
          StrainSmooth = 0.0
    
          do IAEl = 1, Counters%NAEl
          IEl = ActiveElement(IAEl)
          
            DoSmoothing2Layer = .false. !initialize
                
            if(NFORMULATION==2) then
                if(.not.(TwoLayerData%Elements(IEl)%IsSolidMPwithPhaseLiquidInElement)) then
                  DoSmoothing2Layer = .true.
                end if
            end if

            DoSmoothing = (MaterialElements(IMatSet,IEl)==1).and. &
                          ((NFORMULATION==1).or.(DoSmoothing2Layer))
     
              if (DoSmoothing) then
                call BMatrix(LI, &
                             ELEMENTNODES, Counters%NEl,  &
                             Counters%NodTot, NVECTOR, &
                             IEl, ElementConnectivities, &
                             NodalCoordinatesUpd,  &
                             BMatrixDeformed, DetJac)
				
                if ( ISAXISYMMETRIC ) then
                  call getNodalVolumesOfAnElement(IEl, volumeElNodes)
			    end if                

                IParticle = 1
                call Get_Strain(IEl, IParticle, ElementConnectivities, & ! GetStrain.FOR
                                BMatrixDeformed, DUTot, ReducedDof, &
                                Eps)

              do I = 1, NTENSOR
                ElementStrain(IEl, I) = Eps(I)
              end do

              LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
        
              EV = Eps(1) + Eps(2) + Eps(3)
              EVD = EV * DetJac
              if ( ISAXISYMMETRIC ) then
                StrainSmooth(LJ(:), 1) = StrainSmooth(LJ(:), 1) + EV * volumeElNodes
                StrainSmooth(LJ(:), 2) = StrainSmooth(LJ(:), 2) + volumeElNodes
              else
                StrainSmooth(LJ(:), 1) = StrainSmooth(LJ(:), 1) + EVD
                StrainSmooth(LJ(:), 2) = StrainSmooth(LJ(:), 2) + DetJac
              end if
              end if
          end do

          do I = 1, Counters%NodTot
            if (StrainSmooth(I, 2)>0.0) then
              StrainSmooth(I, 1) = StrainSmooth(I, 1) / StrainSmooth(I, 2) ! Only consider active nodes
            end if
          end do
        end subroutine GetEnhancedVolumetricStrain

        subroutine SmoothenVolumetricStrain(IMatSet)
        !**********************************************************************
        !
        !    Function: Smoothens the volumetric strain
        !
        !**********************************************************************
       
        implicit none
        
          integer(INTEGER_TYPE), intent (in) :: IMatSet
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, IAEl, NElemPart, IParticle, ParticleIndex, I, MaterialID
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
          real(REAL_TYPE) :: EnhancedVolumetricStrain ! average enhanced volumetric strain in element
          real(REAL_TYPE) :: ev
          real(REAL_TYPE), dimension (NTENSOR) :: Eps
          logical :: DoSmoothing2Layer, DoSmoothing
  
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            
            DoSmoothing2Layer = .false. !initialize
                
            if(NFORMULATION==2) then
                if(.not.(TwoLayerData%Elements(IEl)%IsSolidMPwithPhaseLiquidInElement)) then
                  DoSmoothing2Layer = .true.
                end if
            end if
                
            DoSmoothing = (MaterialElements(IMatSet,IEl)==1).and. &
                    ((NFORMULATION==1).or.(DoSmoothing2Layer))
     
            if (DoSmoothing) then     !element belongs to the material IMatSet
              LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
              
              EnhancedVolumetricStrain = 0.0
              do I = 1, ELEMENTNODES
                EnhancedVolumetricStrain = EnhancedVolumetricStrain + StrainSmooth(LJ(I), 1) / dble(ELEMENTNODES)
              end do    
                  
              do I = 1, NTENSOR
                Eps(I) = ElementStrain(IEl, I)
              end do

              if((NDIM==3).or.(ISAXISYMMETRIC)) then
                EV = (Eps(1) + Eps(2) + Eps(3) - EnhancedVolumetricStrain ) / 3.D0 ! Correction for volumetric strain ! valid for isaxisymmetric and 3D
                Eps(1:NPRINCIPAL) =  Eps(1:NPRINCIPAL) - EV
              else if ((NDIM==2).and.(.not.ISAXISYMMETRIC)) then
                EV = (Eps(1) + Eps(2) - EnhancedVolumetricStrain ) / 2.D0 ! Correction for volumetric strain ! valid for 2D plane strain
                Eps(1:NDIM) =  Eps(1:NDIM) - EV
              else
                call GiveError('Dimensions not correctly defined for strain smoothening')
              end if
              
              ! Update particle strains

              NElemPart = NPartEle(IEl)
              do IParticle = 1, NElemPart ! Loop over all particles of the element
                ParticleIndex = GetParticleIndex(IParticle, IEl)
                if ((NFORMULATION==1).or. &
                      (MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) ) then !If NumbOfLayers = 1 or Solid MatPoint
                    MaterialID = MaterialIDArray(ParticleIndex)
                    if (MaterialID==IMatSet)then ! this particle belongs to material IMatSet. Assign its strains
                      call SetEpsStep(Particles(ParticleIndex), Eps)
                      if (.not.CalParams%ApplyImplicitQuasiStatic) then
                        call IncreaseEps(Particles(ParticleIndex), Eps)
                      end if
                    end if
                end if
              end do
            end if
        end  do
       
        end subroutine SmoothenVolumetricStrain
        
          subroutine SmoothenStrainsLiquid(DULiquid, IMatSet)
        !**********************************************************************
        !
        !    Function:  Smoothens Gauss point strains for low-order elements.
        !
        !     DUTot : Incremental nodal displacements
        !     IMatSet : The no. of the material set
        !
        !**********************************************************************
       
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: DULiquid
          integer(INTEGER_TYPE), intent (in) :: IMatSet

          call GetEnhancedVolumetricStrainLiquid(DULiquid, IMatSet)

          call SmoothenVolumetricStrainLiquid(IMatSet)

        end subroutine SmoothenStrainsLiquid

        subroutine GetEnhancedVolumetricStrainLiquid(DULiquid, IMatSet)
        !**********************************************************************
        !
        !    Function: Returns enhanced volumetric strain of the liquid phase
        !
        !     DUTotAll : Incremental nodal displacements of all entities
        !     IMatSet : The no. of the material set
        !
        !**********************************************************************
       
        implicit none
        
          real(REAL_TYPE), dimension(Counters%N), intent(in) :: DULiquid
          integer(INTEGER_TYPE), intent(in) :: IMatSet
          
          ! Local variables  
          real(REAL_TYPE), dimension(NTENSOR) :: Eps
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: BMatrixDeformed
          real(REAL_TYPE) :: DetJac, EV, EVD
          real(REAL_TYPE), dimension(NVECTOR) :: LI
          real(REAL_TYPE), dimension(ELEMENTNODES) :: volumeElNodes ! added for axisymmetric analysis
          integer(INTEGER_TYPE) :: IEl, IAEl, I, IParticle
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) ::  LJ
          logical :: FreeWater2Layer, DoSmoothing
        
          LI = 0.25D0
          StrainSmooth = 0.0
        

          do IAEl = 1, Counters%NAEl
          IEl = ActiveElement(IAEl)
              if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
                  
              FreeWater2Layer = .false. !initialize
                
              if(NFORMULATION==2) then
                  if(.not.(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0)) then
                      FreeWater2Layer = .true.
                  end if
              end if
                
              DoSmoothing = (NFORMULATION==1).or.(FreeWater2Layer.and.(.not.IsElemWithLiquidFreeSurfMP(IEl)))
     
                if(DoSmoothing) then
                  call BMatrix(LI, &
                             ELEMENTNODES, Counters%NEl,  &
                             Counters%NodTot, NVECTOR, &
                             IEl, ElementConnectivities, &
                             NodalCoordinatesUpd,  &
                             BMatrixDeformed, DetJac)
                  
                  if ( ISAXISYMMETRIC ) then
                    call getNodalVolumesOfAnElement(IEl, volumeElNodes)
			      end if                   

                  IParticle = 1
                  call Get_Strain(IEl, IParticle, ElementConnectivities, & ! GetStrain.FOR
                                BMatrixDeformed, DULiquid, ReducedDof, &
                                Eps)

                  do I = 1, NTENSOR
                   ElementStrain(IEl, I) = Eps(I)
                  end do

                  LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
        
                  EV = Eps(1) + Eps(2) + Eps(3)
                  EVD = EV * DetJac

                if ( ISAXISYMMETRIC ) then
                  StrainSmooth(LJ(:), 1) = StrainSmooth(LJ(:), 1) + EV * volumeElNodes
                  StrainSmooth(LJ(:), 2) = StrainSmooth(LJ(:), 2) + volumeElNodes
                else
                  StrainSmooth(LJ(:), 1) = StrainSmooth(LJ(:), 1) + EVD
                  StrainSmooth(LJ(:), 2) = StrainSmooth(LJ(:), 2) + DetJac
                end if
                 end if
              end if
          end do

          do I = 1, Counters%NodTot
            if(StrainSmooth(I, 2)>0.0) then
              StrainSmooth(I, 1) = StrainSmooth(I, 1) / StrainSmooth(I, 2) ! Only consider active nodes
            end if
          end do
        end subroutine GetEnhancedVolumetricStrainLiquid

        subroutine SmoothenVolumetricStrainLiquid(IMatSet)
        !**********************************************************************
        !
        !    Function:  Smoothen Volumetric Strain of Liquid
        !
        !**********************************************************************
       
        implicit none
        
          integer(INTEGER_TYPE), intent (in) :: IMatSet
          ! Local variables
          integer(INTEGER_TYPE) :: IEl, IAEl, NElemPart, IParticle, ParticleIndex, I, MaterialID
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
          real(REAL_TYPE) :: EnhancedVolumetricStrain ! Average enhanced volumetric strain in element
          real(REAL_TYPE) :: ev
          real(REAL_TYPE), dimension(NTENSOR) :: Eps
          logical :: FreeWater2Layer, DoSmoothing
        
          do IAEl = 1, Counters%NAEl
            IEl = ActiveElement(IAEl)
            if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
                
              FreeWater2Layer = .false. !initialize
                
              if(NFORMULATION==2) then
                  if(.not.(TwoLayerData%Elements(IEl)%ConcentrationRatioSolidL>0.0)) then
                      FreeWater2Layer = .true.
                  end if
              end if
                
              DoSmoothing = (NFORMULATION==1).or.(FreeWater2Layer.and.(.not.IsElemWithLiquidFreeSurfMP(IEl)))
                
              if(DoSmoothing) then
                LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
                
                EnhancedVolumetricStrain = 0.0
                do I = 1, ELEMENTNODES
                  EnhancedVolumetricStrain = EnhancedVolumetricStrain + StrainSmooth(LJ(I), 1) / dble(ELEMENTNODES)
                end do    

                do I = 1, NTENSOR
                  Eps(I) = ElementStrain(IEl, I)
                end do
                
              if((NDIM==3).or.(ISAXISYMMETRIC)) then
                EV = (Eps(1) + Eps(2) + Eps(3) - EnhancedVolumetricStrain ) / 3.D0 ! Correction for volumetric strain ! valid for isaxisymmetric and 3D
                Eps(1:NPRINCIPAL) =  Eps(1:NPRINCIPAL) - EV
              else if ((NDIM==2).and.(.not.ISAXISYMMETRIC)) then
                EV = (Eps(1) + Eps(2) - EnhancedVolumetricStrain ) / 2.D0 ! Correction for volumetric strain ! valid for 2D plane strain
                Eps(1:NDIM) =  Eps(1:NDIM) - EV
              else
                call GiveError('Dimensions not correctly defined for strain smoothening')
              end if
             
                
                ! Update particle strains

                NElemPart = NPartEle(IEl)
                do IParticle = 1, NElemPart ! Loop over all particles of the element
                  ParticleIndex = GetParticleIndex(IParticle, IEl)
                  if ((NFORMULATION==2).and.(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid)) then
                    ! 2 Constituents
                    MaterialID = MaterialIDArray(ParticleIndex)
                    if (MaterialID==IMatSet)then ! this particle belongs to material IMatSet. Assign its strains
                      Particles(ParticleIndex)%WaterVolumetricStrain = Eps(1) + Eps(2) + Eps(3) ! valid for 2D and 3D
                      call SetEpsStep(Particles(ParticleIndex), Eps)
                      call IncreaseEps(Particles(ParticleIndex), Eps)
                    end if
                  end if
                end do
              end if
            end if
        end  do
       
        end subroutine SmoothenVolumetricStrainLiquid
        
        
         subroutine SmoothenStrainsWater(IMatSet)
        !**********************************************************************
        !
        !    Function: To smoothen the volumetric strain of water phase
        !
        !     IMatSet : The no. of the material set
        !
        !**********************************************************************
       
        implicit none
 
          integer(INTEGER_TYPE), intent(in) :: IMatSet
          ! Local variables  
          real(REAL_TYPE), dimension(Counters%NEl) :: ElementsVolumetricStrain
          real(REAL_TYPE) :: EV, EVD, EnhancedVolumetricStrain
          integer(INTEGER_TYPE) :: IEl, IAEl, I, NElemPart, MaterialID, IParticle, ParticleIndex
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
        
          StrainSmooth = 0.0
          ElementsVolumetricStrain = 0.0

          do IAEl = 1, Counters%NAEl  ! active element
            IEl = ActiveElement(IAEl)
            if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
              NElemPart = NPartEle(IEl)
              do IParticle = 1, NElemPart
                ParticleIndex = GetParticleIndex(IParticle, IEl)
                MaterialID = MaterialIDArray(ParticleIndex)
                if (MaterialID==IMatSet)then ! This particle belongs to the considered entity (Material)
                  ! The volumetric strain in this element for the considered entity (Material)
                  EV = Particles(ParticleIndex)%WaterVolumetricStrain
                  Exit
                end if
              end do ! particles

              ElementsVolumetricStrain(IEl) = EV
              LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
              EVD = EV * ElementSpace(IEl)

              StrainSmooth(LJ(:), 1) = StrainSmooth(LJ(:), 1) + EVD
              StrainSmooth(LJ(:), 2) = StrainSmooth(LJ(:), 2) + ElementSpace(IEl)

              end if ! entity
          end do ! active element


          do I = 1, Counters%NodTot
            if(StrainSmooth(I, 2)>0.0) then
              StrainSmooth(I, 1) = StrainSmooth(I, 1) / StrainSmooth(I, 2) ! Only consider active nodes
            end if
          end do

            do IAEl = 1, Counters%NAEl  ! active element 
              IEl = ActiveElement(IAEl)
                if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
                  LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
                  
                  EnhancedVolumetricStrain = 0.0
                  do I = 1, ELEMENTNODES
                    EnhancedVolumetricStrain = EnhancedVolumetricStrain + StrainSmooth(LJ(I), 1) / dble(ELEMENTNODES)
                  end do    

!              EV = ElementsVolumetricStrain (IEl) - EnhancedVolumetricStrain
     
              EV = EnhancedVolumetricStrain
        
              ! Update particle strains

                NElemPart = NPartEle(IEl)
                do IParticle = 1, NElemPart ! Loop over all particles of the element
                  ParticleIndex = GetParticleIndex(IParticle, IEl)
                    MaterialID = MaterialIDArray(ParticleIndex)
                      if (MaterialID==IMatSet)then ! this particle belongs to material IMatSet. Assign its strains
                        Particles(ParticleIndex)%WaterVolumetricStrain = EV
                      end if
                end do ! particles
              end if ! entity
           end  do ! active eleemnt
          
        end subroutine SmoothenStrainsWater
     
         subroutine SmoothenStrainsGas(IMatSet)
        !**********************************************************************
        !
        !    Function: To smoothen the volumetric strain of gas phase
        !
        !     IMatSet : The no. of the material set
        !
        !**********************************************************************
       
        implicit none
 
          integer(INTEGER_TYPE), intent(in) :: IMatSet
          ! Local variables  
           real(REAL_TYPE), dimension( Counters%NEl) :: ElementsVolumetricStrain
          real(REAL_TYPE) :: EV, EVD, EnhancedVolumetricStrain
          integer(INTEGER_TYPE) :: IEl, IAEl, I, NElemPart, MaterialID, IParticle, ParticleIndex
          integer(INTEGER_TYPE), dimension(ELEMENTNODES) :: LJ
        
          StrainSmooth = 0.0
          ElementsVolumetricStrain = 0.0

          do IAEl = 1, Counters%NAEl  ! active element
            IEl = ActiveElement(IAEl)
            if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
              NElemPart = NPartEle(IEl)
              do IParticle = 1, NElemPart
                ParticleIndex = GetParticleIndex(IParticle, IEl)
                MaterialID = MaterialIDArray(ParticleIndex)
                if (MaterialID==IMatSet)then ! This particle belongs to the considered entity (Material)
                  ! The volumetric strain in this element for the considered entity (Material)
                  EV = Particles(ParticleIndex)%GasVolumetricStrain
                  Exit
                end if
              end do ! particles

              ElementsVolumetricStrain(IEl) = EV
              LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
              EVD = EV * ElementSpace(IEl)

              StrainSmooth(LJ(:), 1) = StrainSmooth(LJ(:), 1) + EVD
              StrainSmooth(LJ(:), 2) = StrainSmooth(LJ(:), 2) + ElementSpace(IEl)

            end if ! entity
          end do ! active element


          do I = 1, Counters%NodTot
            if(StrainSmooth(I, 2)>0.0) then
              StrainSmooth(I, 1) = StrainSmooth(I, 1) / StrainSmooth(I, 2) ! Only consider active nodes
            end if
          end do
          

            do IAEl = 1, Counters%NAEl  ! active element 
              IEl = ActiveElement(IAEl)
                if (MaterialElements(IMatSet,IEl)==1) then     !element belongs to the material IMatSet
                  LJ(:) = ElementConnectivities(1:ELEMENTNODES, IEl)
                  
                  EnhancedVolumetricStrain = 0.0
                  do I = 1, ELEMENTNODES
                    EnhancedVolumetricStrain = EnhancedVolumetricStrain + StrainSmooth(LJ(I), 1) / dble(ELEMENTNODES)
                  end do    

!              EV = ElementsVolumetricStrain (IEl) - EnhancedVolumetricStrain                      
     
              EV = EnhancedVolumetricStrain
        
              ! Update particle strains

                NElemPart = NPartEle(IEl)
                do IParticle = 1, NElemPart ! Loop over all particles of the element
                  ParticleIndex = GetParticleIndex(IParticle, IEl)
                    MaterialID = MaterialIDArray(ParticleIndex)
                      if (MaterialID==IMatSet)then ! this particle belongs to material IMatSet. Assign its strains
                        Particles(ParticleIndex)%GasVolumetricStrain = EV
                      end if
                end do ! particles
              end if ! entity
           end  do ! active element
          
        end subroutine SmoothenStrainsGas
          
        subroutine getNodalVolumesOfAnElement(IElement, volumeElNodes)
        !**********************************************************************
        !
        !    Function:  Determines the volumes per nodes for all active elements.
        !
        !    IElement : Number of Element
        !    volumeElNodes : Volume of each Elementnode
        ! 
        !    Implemented in Anura3D during Sprint #4.
        !
        !**********************************************************************
        implicit none
    
        ! arguments
        integer, intent(in) :: IElement
        real(REAL_TYPE), dimension(ELEMENTNODES), intent(out) ::volumeElNodes
    
        ! Local variables
        integer :: IGaussPoint, INode
        real(REAL_TYPE), dimension(NVECTOR) :: PosGP
        real(REAL_TYPE) :: WeiGP
        real(REAL_TYPE) :: DetJac, weight
        real(REAL_TYPE), dimension(NDIM, NDIM) :: RJac
        real(REAL_TYPE), dimension(NDIM, NDIM) :: RJacInv
    
        if ( .not. ISAXISYMMETRIC ) RETURN
      
        volumeElNodes = 0.0
        do IGaussPoint = 1, ELEMENTGAUSSPOINTS ! loop gauss points of each element
          ! determine the location and integration weight of IGaussPoint
          call GaussPointLocalCoordinates(IGaussPoint, WeiGP, PosGP)
          ! Calculate the determinante of the Jacobian matrix
          call DetJacob(PosGP, Counters%NEl, Counters%NodTot, NVECTOR, IElement, ElementConnectivities, NodalCoordinates, RJac, RJacInv, DetJac)
          if ( ISAXISYMMETRIC ) then
            weight = WeiGP * GPGlobalPositionElement(1, IGaussPoint, IElement)
          else
            weight = WeiGP
          end if
          do INode = 1, ELEMENTNODES
            volumeElNodes(INode) = volumeElNodes(INode) + ShapeValuesArray(IGaussPoint, INode) * weight * DetJac
          end do
        end do
      
      end subroutine getNodalVolumesOfAnElement
        
      end module ModStrainSmoothing