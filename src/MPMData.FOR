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
	
	
	  module ModMPMData
      !**********************************************************************
      !
      !    Function:  Contains the housekeeping data structure required for the MPM project.
      !
      !               In order to keep the size of this source file reasonably small,
      !               this module only contains routines that are directly related to
      !               the basic manipulation of contained data (maintenance of data structures).
      !
      !               First, routines are provided to access the material point data, followed by routines
      !               for generating the material point house-keeping data structure.
      !
      !               Routines for initialising or updating material point information or 
      !               routines for output of material point data are located in separate modules that access the
      !               data stored in this module.
      !
      !     $Revision: 9808 $
      !     $Date: 2022-10-13 16:48:47 +0200 (do, 13 okt 2022) $
      !
      !**********************************************************************

      use ModGlobalConstants
      use ModCounters
      use ModReadCalculationData
      use ModElementEvaluation
      use ModMeshAdjacencies
      use ModParticle
      use ModMeshInfo
      
      implicit none

        !> Decimal factor used for defining particle - element connectivities with the array EleParticles
        integer(kind = 8), private, save  :: PartFactor = 10
        
        ! implicit variables
        !> NPartEle(IElement) = Number of particles in element IElement, size = NEl
        integer(INTEGER_TYPE), dimension(:), allocatable :: NPartEle_N ! -> implicit bookmark
        integer(INTEGER_TYPE), dimension(:), allocatable :: NPartEle_NPlus1 ! -> implicit bookmark

        !> NSolidEle(IElement) = Number of solid material points in element IElement, size = NEl
        integer(INTEGER_TYPE), dimension(:), allocatable :: NSolidEle_N !-> implicit bookmark
        integer(INTEGER_TYPE), dimension(:), allocatable :: NSolidEle_NPlus1 !-> implicit bookmark

        !> NLiquidEle(IElement) = Number of liquid material points in element IElement, size = NEl
        integer(INTEGER_TYPE), dimension(:), allocatable :: NLiquidEle_N !-> not needed at the moment
              
        integer(INTEGER_TYPE), dimension(:), allocatable :: NLiquidEle_NPlus1 !-> not needed at the moment

        
        !> List which contains the particle connectivities (particle - element assigments), size = NParticles
        integer(kind = 8), dimension(:), allocatable :: EleParticles_N ! --> implicit bookmark
        !> List which contains the particle connectivities (particle - element assigments), size = NParticles
        integer(kind = 8), dimension(:), allocatable :: EleParticles_NPlus1 ! --> implicit bookmark
        
        !> Helper list for EleParticles, index of first particle of element, size = NEl
        integer(INTEGER_TYPE), dimension(:), allocatable :: EleParticlesHelp_N !-> implicit bookmark
        !> Helper list for EleParticles, index of first particle of element, size = NEl
        integer(INTEGER_TYPE), dimension(:), allocatable :: EleParticlesHelp_NPlus1 !-> implicit bookmark
        
        !> True, if particle based integration is used inside element IElement, size = NEl
        logical, dimension(:), allocatable :: IsParticleIntegration
        
        !> True, if liquid free surface Material Point is inside element IElement, size = NEl
        logical, dimension(:), allocatable :: IsElemWithLiquidFreeSurfMP
        
        type(ParticleType(:, :, :)), dimension(:), allocatable :: Particles ! List of particle objects
        
        !> The following arrays correspond to the stripped members of the Particle datatype to improve efficency (spatial locality)  
        integer(INTEGER_TYPE), allocatable          :: IDArray(:)                   !> Unique ID of a particle. This variable should only be set once through the routine InitParticle.
        
        ! implicit variable
        integer(INTEGER_TYPE), allocatable          :: ElementIDArray_N(:)            !> Unique ID of the element to whom the particle belongs. A value of -1 indicates no assignment.
        integer(INTEGER_TYPE), allocatable          :: ElementIDArray_NPlus1(:)            !> Unique ID of the element to whom the particle belongs. A value of -1 indicates no assignment.
        
        integer(INTEGER_TYPE), allocatable          :: EntityIDArray(:)             !> Unique ID indicating to which entity the particle belongs
        
        ! implicit varibale 
        integer(INTEGER_TYPE), allocatable          :: MaterialIDArray_N(:)           !> Material ID of the particle referring to the material data read from input files
        integer(INTEGER_TYPE), allocatable          :: MaterialIDArray_NPlus1(:)           !> Material ID of the particle referring to the material data read from input files
        
        integer(INTEGER_TYPE), allocatable          :: MaterialPointTypeArray(:)    !> Material properties (properties related to the matter represented by particles)
        real(REAL_TYPE), allocatable :: MassArray(:)                 !> Particle mass [kg] (corresponds to the dry density of soil) (fixed during simulation - ensures conservation of mass)
        real(REAL_TYPE), allocatable :: MassWaterArray(:)            !> Particle mass [kg] (corresponds to the water density)(fixed during simulation - ensures conservation of mass)
        
        ! Displacements (solid)
        real(REAL_TYPE), allocatable :: UArray(:,:)                  !> Total displacements in x, y, z direction of the particle
        real(REAL_TYPE), allocatable :: UPhaseArray(:,:)             !> Displacements in x, y, z direction of the particle during a phase
        real(REAL_TYPE), allocatable :: UStepArray(:,:)              !> Displacements in x, y, z direction of the particle during a load step
       
        ! implicit scheme variables
        real(REAL_TYPE), allocatable :: ShapeValuesArray_N(:,:) ! -> implicit bookmark       
                                        !> Shape function values for the particle, these values are updated when new local coordinates are calculated.
        real(REAL_TYPE), allocatable :: ShapeValuesArray_NPlus1(:,:) !-> implicit bookmark       
                                        !> Shape function values for the particle, these values are updated when new local coordinates are calculated.
        
        ! implicit scheme variables
        real(REAL_TYPE), allocatable :: GlobPosArray_N(:,:) ! -> implicit bookmark           
                                        !> Global particle position x = GlobPos(1), y = GlobPos(2), z = GlobPos(3)
        real(REAL_TYPE), allocatable :: GlobPosArray_NPlus1(:,:)
        
        ! implicit scheme variables
        real(REAL_TYPE), allocatable :: VelocityArray_N(:,:)           ! Velocity of the particle (soil)
        real(REAL_TYPE), allocatable :: VelocityArray_NPlus1(:,:)           ! Velocity of the particle (soil)

        
        real(REAL_TYPE), allocatable :: VelocityWaterArray(:,:)      ! Velocity of the particle (water)
        real(REAL_TYPE), allocatable :: VelocityGasArray(:,:)        ! Velocity of the particle (gas)
        
        ! Stresses
        real(REAL_TYPE), allocatable :: SigmaEffArray_N(:,:) ! -> implicit bookmark          
        real(REAL_TYPE), allocatable :: SigmaEffArray_NPlus1(:,:) ! -> implicit bookmark          
        
        !> Effective stresses (Sxx, Syy, Szz, Sxy, Syz, Sxz) of the particle 
        real(REAL_TYPE), allocatable :: SigmaEff0Array(:,:) ! -> implicit bookmark --> this might not be needed         
        !> Initial Effective stresses (Sxx, Syy, Szz, Sxy, Syz, Sxz) of the particle from the last load step
        
        ! implicit scheme variables
        real(REAL_TYPE), allocatable :: DShapeValuesArray_N(:,:,:)     !> Shape function derivatives for the particle, these values are updated when new local coordinates are calculated.
        real(REAL_TYPE), allocatable :: DShapeValuesArray_NPlus1(:,:,:)     !> Shape function derivatives for the particle, these values are updated when new local coordinates are calculated.
        
        real(REAL_TYPE), allocatable :: AccelerationArray(:,:)       !> Acceleration of the particle (soil)
        
        ! User defined model
        !real(REAL_TYPE), allocatable :: ESMstatevArray(:,:)         !> array of state variables of the user defined material model
        real(REAL_TYPE), allocatable :: ESMstatevArray_N(:,:)         !> array of state variables of the user defined material model
        real(REAL_TYPE), allocatable :: ESMstatevArray_NPlus1(:,:)         !> array of state variables of the user defined material model
        real(REAL_TYPE), allocatable :: ESMstatevArray_Iteration(:,:)         !> array of state variables of the user defined material model
        
        !real(REAL_TYPE), allocatable :: ESMpropsArray(:,:)         !> array of material model properties of the user defined material model (Particle,property)

      contains ! Routines of this module


        logical function IsMaterialParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  returns .true. if the material point with ID MaterialPointIndex is
        !               a material point
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsMaterialParticle =  Particles(MaterialPointIndex)%Kind==MATERIALPARTICLE
          
        end function IsMaterialParticle


        logical function IsVirtualParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the material point with ID MaterialPointIndex is
        !               a virtual material point
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsVirtualParticle = Particles(MaterialPointIndex)%Kind==VIRTUALPARTICLE
          
        end function IsVirtualParticle


        logical function IsRemovedParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the material point with ID MaterialPointIndex is
        !               a material point to be removed
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsRemovedParticle = Particles(MaterialPointIndex)%Kind==REMOVEDPARTICLE
          
        end function IsRemovedParticle


        logical function IsNewVirtualParticle(MaterialPointIndex)
        !**********************************************************************
        !
        !    Function:  Returns .true. if the material point with ID MaterialPointIndex is
        !               a newly initialised virtual material point
        !
        !    MaterialPointIndex : ID of the considered material point
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialPointIndex
          
          IsNewVirtualParticle = Particles(MaterialPointIndex)%Kind==NEWVIRTUALPARTICLE
          
        end function IsNewVirtualParticle


        integer(INTEGER_TYPE) function NMaterialParticlesEle(IElement, NPartEle)
        !**********************************************************************
        !
        !    Function:  Returns the number of material points in IElement
        !
        !    IElement : Considered element
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, MaterialPointIndex
          integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle

          NMaterialParticlesEle = 0
          do I = 1, NPartEle(IElement)
            MaterialPointIndex = GetParticleIndex_N(I, IElement)
            
            if (IsMaterialParticle(MaterialPointIndex)) then
              NMaterialParticlesEle = NMaterialParticlesEle + 1
            end if

          end do
        
        end function NMaterialParticlesEle


        integer(INTEGER_TYPE) function NVirtualParticlesEle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns the number of virtual material point in IElement
        !
        !    IElement : Considered element
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, ParticleIndex

          NVirtualParticlesEle = 0
          do I = 1, NPartEle_N(IElement)
            ParticleIndex = GetParticleIndexFunction(I, IElement)!, &
                                                     !EleParticlesHelp, EleParticles, NPartEle)
            
            if (IsVirtualParticle(ParticleIndex) ) then
              NVirtualParticlesEle = NVirtualParticlesEle + 1
            end if
          end do

        end function NVirtualParticlesEle

        integer(INTEGER_TYPE) function NRemovedParticlesEle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns the number of material points to be removed in IElement
        !
        !    IElement : Considered element
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: I, ParticleIndex

          NRemovedParticlesEle = 0
          do I = 1, NPartEle_N(IElement)
            ParticleIndex = GetParticleIndex_N(I, IElement)
            
            if (IsRemovedParticle(ParticleIndex) ) then
              NRemovedParticlesEle = NRemovedParticlesEle + 1
            end if
          end do

        end function NRemovedParticlesEle


        logical function HasBoundaryParticle(IElement)
        !**********************************************************************
        !
        !    Function:  Returns .true. if IElement contains a material point marked as
        !               boundary material point
        !
        !     IElement : Considered element
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex
        
          HasBoundaryParticle = .false.
          do IParticle = 1, NPartEle_N(IElement) ! Loop over particles in IElement
            ParticleIndex = GetParticleIndexFunction(IParticle, IElement)!, &
                                                     !EleParticlesHelp, EleParticles, NPartEle))
      
            if (Particles(ParticleIndex)%IsBoundaryParticle) then
              HasBoundaryParticle = .true.
              EXIT
            end if
          end do
        
        end function HasBoundaryParticle


        integer(INTEGER_TYPE) function GetParticleIndexFunction(IParticle, IElement) result(res) !EleParticlesHelp, EleParticles, NPartEle
        !**********************************************************************
        !
        !    Function:  Returns the index of the material point with the local number
        !               IParticle inside element IElement. The returned index
        !               is the location of the material point in the array 'Particles'.
        !               Note that the index of the material point might not be identical
        !               to the particle ID (if material points might be deleted ... ).
        !               If less than IParticle material points are inside IElement, -1 is
        !               returned.
        !
        !     IParticle : Local number of the material point inside IElement
        !     IElement : ID of the element, that the considered material point is located in
        !
        ! O   GetParticleIndex : Index of the considered material point in the array 'Particles'
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          integer(INTEGER_TYPE), intent(in) :: IElement
          ! Local variables
          integer(INTEGER_TYPE) :: EleParticleIndex
          integer(kind = 8) :: LocalIElement
          ! Input variables 
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: EleParticlesHelp
          !integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(in) :: EleParticles
          !integer(INTEGER_TYPE), dimension(Counters%NEl), intent(in) :: NPartEle
          
          if ( (allocated(EleParticlesHelp_N) ).and.(allocated(EleParticles_N) ).and. &
               (IParticle>=1).and.(IParticle<=NPartEle_N(IElement))            ) then
            EleParticleIndex = EleParticlesHelp_N(IElement) + IParticle - 1
            LocalIElement = IElement
            res = EleParticles_N(EleParticleIndex) - LocalIElement * PartFactor
          else
            res = -1
          end if

        end function GetParticleIndexFunction

        


        integer(INTEGER_TYPE) function GetElementIDFromList(IParticle, EleParticles)
        !**********************************************************************
        !
        !    Function:  Returns the Element ID of the particle stored at index
        !               IParticle in array EleParticles.
        !               
        !     IParticle : Index of list EleParticles
        !
        ! O   GetElementIDFromList : ID of the element that particle at index IParticle is located in
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: IParticle
          integer(kind = 8), dimension(Counters%NParticles), intent(in) :: EleParticles

          GetElementIDFromList = int(EleParticles(IParticle) /  &
                                     PartFactor)
        
        end function GetElementIDFromList

        integer(INTEGER_TYPE) function GetParticleIndexFromList(IParticle)
        !**********************************************************************
        !
        !    Function:  Returns the index of the particle for the location of the 
        !               particle inside the list EleParticles. The returned index
        !               is the location of the particle in the array 'Particles'.
        !               Note: This routine is used for updating the particles, where
        !                     a loop over the items of EleParticles is performed.
        !               
        !     IParticle : Index of list EleParticles
        !
        ! O   GetParticleIndexFromList : Index of the considered particle in the array 'Particles'
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          
          if ( (allocated(EleParticles_N) ).and. &
               (IParticle>=1) ) then
            GetParticleIndexFromList = mod(EleParticles_N(IParticle), &
                                           PartFactor)
          else
            GetParticleIndexFromList = -1
          end if
        
        end function GetParticleIndexFromList


        integer(kind = 8) function GetPartFactor()
        !**********************************************************************
        !
        !    Function:  Returns PartFactor, which is declared private as it should only
        !               be modified inside this module.
        !
        ! O   GetPartFactor : Decimal factor used internally for storing particle-element connectivities
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
          
          GetPartFactor = PartFactor
        
        end function GetPartFactor


        subroutine SetPartFactor(NewPartFactor)
        !**********************************************************************
        !
        !    Function:  Set PartFactor, which is declared private as it should only
        !               be modified inside this module. PartFactor stores the decimal factor 
        !               used internally for storing particle-element connectivities.
        !
        !     NewPartFactor : Decimal factor read from the BRF file of the previous load phase
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(kind = 8), intent(in) :: NewPartFactor
        
          PartFactor = NewPartFactor
        
        end subroutine SetPartFactor


        integer(INTEGER_TYPE) function NumberOfIntegrationPoints(IElement)
        !**********************************************************************
        !
        !    Function: Returns the number of particles in element IElement if the element is 
        !              partially filled and the number of Gauss points per element if the 
        !              element is fully filled.
        !
        !     IElement : ID of the considered element
        !
        ! O   NumberOfIntegrationPoints : Number of Gauss points or particles
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
     
        implicit none
        
          integer(INTEGER_TYPE) :: IElement
        
          if (IsParticleIntegration(IElement) ) then ! Partially filled element
            NumberOfIntegrationPoints = NPartEle_N(IElement) ! Number of particles in IElement
          else ! Fully filled element
            NumberOfIntegrationPoints = ELEMENTGAUSSPOINTS ! Number of Gauss points in element
          end if
        
        end function NumberOfIntegrationPoints


        subroutine SetEleParticles(IElement, IParticle, ParticleIndex, &
                                   EleParticles, EleParticlesHelp)
        !**********************************************************************
        !
        !    Function:  Sets EleParticles(ParticleIndex) and EleParticlesHelp(IElement)
        !               if ParticleIndex is larger than 0. Else EleParticlesHelp(IElement)
        !               is set to -1 to indicate that IElement is inactive and therefor
        !               contains no particles. The particle-element connectivity stored
        !               in EleParticles consists of the element ID multiplied by PartFactor
        !               and ParticleIndex. EleParticlesHelp contains the ParticleIndex of
        !               the first particle located in IElement. The following NPartEle(IElement)
        !               particles of EleParticles also belong to IElement.
        !
        !     IElement : ID of the considered element
        !     IParticle : Local ID of the considered particle inside IElement
        !     ParticleIndex : Index of the considered particle in the Particles array
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
          implicit none
          
            integer(INTEGER_TYPE), intent(in) :: IElement
            integer(INTEGER_TYPE), intent(in) :: IParticle
            integer(INTEGER_TYPE), intent(in) :: ParticleIndex
            ! Local variables
            integer(kind = 8) :: LocalIElement
            integer(kind = 8) :: LocalParticleIndex
            
            ! Input variables
            integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(out) :: EleParticles !-> implicit bookmark
            integer(INTEGER_TYPE), dimension(Counters%NParticles), intent(out) :: EleParticlesHelp !-> implicit bookmark

            
            if (ParticleIndex>0) then
              LocalIElement = IElement
              LocalParticleIndex = ParticleIndex

              EleParticles(ParticleIndex) = LocalIElement * PartFactor + &
                                            LocalParticleIndex
     
          
              if (IParticle==1) then ! First particle of element
                EleParticlesHelp(IElement) = ParticleIndex
              end if
            else ! Inactive element
              EleParticlesHelp(IElement) = -1
            end if
        
        end subroutine SetEleParticles


        subroutine SetParticleElementID(IParticle, ParticleIndex, NewElementID, ElementIDArray, EleParticles)
        !**********************************************************************
        !
        !    Function:  Sets EleParticles(IParticle) and Particles(ParticleIndex)%ElementID.
        !               NOTE: Only this routine should be used to update these two values in order
        !               to guarantee that they are identical!
        !
        !     IParticle : Index of the considered particle in the EleParticles array
        !     ParticleIndex : Index of the considered particle in the Particles array
        !     NewElementID : ID of the new element of the considered particle
        !
        !**********************************************************************

        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IParticle
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex
          integer(INTEGER_TYPE), intent(in) :: NewElementID
          
          integer(INTEGER_TYPE), intent(out), dimension(Counters%NParticles) :: ElementIDArray
          integer(kind = 8), intent(out), dimension(Counters%NParticles) :: EleParticles!Counters%NEl) :: 
          
          
          ! Local variables
          integer(kind = 8) :: LocalNewElementID
          integer(kind = 8) :: LocalParticleIndex
        
          LocalNewElementID = NewElementID
          LocalParticleIndex = ParticleIndex

          ElementIDArray(ParticleIndex)= NewElementID
          EleParticles(IParticle) = LocalNewElementID * PartFactor +  &
                                    LocalParticleIndex
        end subroutine SetParticleElementID


        subroutine InitialiseHouseKeepingArrays(&
            EleParticlesHelp_N, EleParticles_N, &
            EleParticlesHelp_NPlus1, EleParticles_NPlus1, &
            ElementIDArray_N, MaterialIDArray_N, &
            ElementIDArray_NPlus1, MaterialIDArray_NPlus1, &
            ShapeValuesArray_N, DShapeValuesArray_N, &
            ShapeValuesArray_NPlus1, DShapeValuesArray_NPlus1, &
            GlobPosArray_N, GlobPosArray_NPlus1, &
            SigmaEffArray_N, SigmaEffArray_NPlus1)
        !**********************************************************************
        !
        !    Function:  Allocates the arrays storing house-keeping data.
        !               Sets their initial values.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
          
          ! local variables
          integer(INTEGER_TYPE) :: IError
          integer(INTEGER_TYPE) :: I
          
          ! input 
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: EleParticlesHelp_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: EleParticlesHelp_NPlus1
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: EleParticles_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: EleParticles_NPlus1
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: ElementIDArray_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: ElementIDArray_NPlus1
          
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: MaterialIDArray_N
          integer(INTEGER_TYPE), allocatable, dimension(:), intent(inout) :: MaterialIDArray_NPlus1
          
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: ShapeValuesArray_N
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: ShapeValuesArray_NPlus1
          
          real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: DShapeValuesArray_N
          real(REAL_TYPE), allocatable, dimension(:,:,:), intent(inout) :: DShapeValuesArray_NPlus1
          
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: GlobPosArray_N
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: GlobPosArray_NPlus1
          
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: SigmaEffArray_N
          real(REAL_TYPE), allocatable, dimension(:,:), intent(inout) :: SigmaEffArray_NPlus1
          
          
          
          ! input - implicit bookmarks

          allocate(EleParticlesHelp_N(Counters%NEl), stat = IError)
          EleParticlesHelp_N = -1
          
          allocate(EleParticles_N(Counters%NParticles), stat = IError)
          EleParticles_N = -1
          
          ! implicit scheme variables
          allocate(EleParticlesHelp_NPlus1(Counters%NEl), stat = IError)
          EleParticlesHelp_NPlus1 = -1
          
          allocate(EleParticles_NPlus1(Counters%NParticles), stat = IError)
          EleParticles_NPlus1 = -1
          
          allocate(ParticleType(NTENSOR, NVECTOR, MAX_LOAD_SYSTEMS)::Particles(Counters%NParticles), stat = IError)
          do i = 1, Counters%NParticles 
              call SetParticleStructureDefault(Particles(i))
          end do
          
          allocate(ElementIDArray_N(Counters%NParticles), stat = IError)
          allocate(ElementIDArray_NPlus1(Counters%NParticles), stat = IError)
          
          allocate(IDArray(Counters%NParticles), stat = IError)
          allocate(EntityIDArray(Counters%NParticles), stat = IError)
          
          allocate(MaterialIDArray_N(Counters%NParticles), stat = IError)
          allocate(MaterialIDArray_NPlus1(Counters%NParticles), stat = IError)
          
          allocate(UArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(MaterialPointTypeArray(Counters%NParticles), stat = IError)
          allocate(UStepArray(Counters%NParticles,NVECTOR), stat = IError)
          
          ! implicit scheme variable
          allocate(ShapeValuesArray_N(Counters%NParticles,ELEMENTNODES), stat = IError) !-> implicit bookmark
          allocate(ShapeValuesArray_NPlus1(Counters%NParticles,ELEMENTNODES), stat = IError) ! -> implicit bookmark
          
          allocate(UPhaseArray(Counters%NParticles,NVECTOR), stat = IError)
          
          ! implicit scheme variable 
          allocate(GlobPosArray_N(Counters%NParticles,NVECTOR), stat = IError)
          allocate(GlobPosArray_NPlus1(Counters%NParticles,NVECTOR), stat = IError)

          ! implicit scheme variable
          allocate(VelocityArray_N(Counters%NParticles,NVECTOR), stat = IError)
          allocate(VelocityArray_NPlus1(Counters%NParticles,NVECTOR), stat = IError)

          ! implicit scheme variable
          allocate(SigmaEffArray_N(Counters%NParticles,NTENSOR), stat = IError)
          allocate(SigmaEff0Array(Counters%NParticles,NTENSOR), stat = IError)
          
          allocate(SigmaEffArray_NPlus1(Counters%NParticles,NTENSOR), stat = IError)
          !allocate(SigmaEff0Array(Counters%NParticles,NTENSOR), stat = IError)
          
          
          ! implicit scheme variable
          allocate(DShapeValuesArray_N(Counters%NParticles,ELEMENTNODES,NVECTOR), stat = IError) !-> implicit bookmark
          allocate(DShapeValuesArray_NPlus1(Counters%NParticles,ELEMENTNODES,NVECTOR), stat = IError) !-> implicit bookmark
          
          allocate(VelocityWaterArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(VelocityGasArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(AccelerationArray(Counters%NParticles,NVECTOR), stat = IError)
          allocate(MassArray(Counters%NParticles), stat = IError)
          
          ! implicit scheme variables 
          allocate(ESMstatevArray_N(Counters%NParticles, NSTATEVAR), stat = IError)
          allocate(ESMstatevArray_NPlus1(Counters%NParticles, NSTATEVAR), stat = IError)
          allocate(ESMstatevArray_Iteration(Counters%NParticles, NSTATEVAR), stat = IError)
          !allocate(ESMpropsArray(Counters%NParticles, NPROPERTIES), stat = IError)
          allocate(MassWaterArray(Counters%NParticles), stat = IError)
          allocate(IsParticleIntegration(Counters%NEl), stat = IError)
          
          ! Initialize
          
          ! implicit scheme variable 
          ElementIDArray_N = -1 ! -> implicit bookmark
          ElementIDArray_NPlus1 = -1 ! -> implicit bookmark
          
          IDArray = - 1
          EntityIDArray = -1
          
          MaterialIDArray_N = -1 ! -> implicit bookmark
          MaterialIDArray_NPlus1 = -1 ! -> implicit bookmark
          
          MaterialPointTypeArray = MaterialPointTypeUndefined
          UStepArray = 0.0
          UArray = 0.0
          
          ShapeValuesArray_N = 0.0 ! -> implicit bookmark
          ShapeValuesArray_NPlus1 = 0.0 ! -> implicit bookmark
          
          UPhaseArray = 0.0
          
          GlobPosArray_N =  0.0! -> implicit bookmark
          GlobPosArray_NPlus1 =  0.0! -> implicit bookmark
          
          VelocityArray_N = 0.0
          VelocityArray_NPlus1 = 0.0
          
          SigmaEffArray_N = 0.0
          SigmaEffArray_NPlus1 = 0.0
          
          SigmaEff0Array = 0.0
          
          DShapeValuesArray_N = 0.0
          DShapeValuesArray_NPlus1 = 0.0
          
          VelocityWaterArray = 0.0
          VelocityGasArray = 0.0
          AccelerationArray = 0.0
          MassArray = -1.0
          MassWaterArray = -1.0
          
          if (IsMPMWithMPIntegration()) then
            ! Initially (and always) all elements are partially filled, always use mass point integration
            IsParticleIntegration = .true.
          else
            ! At least initially all elements are fully filled
            IsParticleIntegration = .false.
          end if

        end subroutine InitialiseHouseKeepingArrays


        subroutine DestroyHouseKeeping(&
            NPartEle_N, &
            NPartEle_NPlus1, &
            EleParticlesHelp_N, &
            EleParticlesHelp_NPlus1, &
            EleParticles_N, &
            EleParticles_NPlus1, &
            ElementIDArray_N, &
            ElementIDArray_NPlus1, &
            MaterialIDArray_N, &
            MaterialIDArray_NPlus1, &
            ShapeValuesArray_N, &
            ShapeValuesArray_NPlus1, &
            DShapeValuesArray_N, &
            DShapeValuesArray_NPlus1, &
            GlobPosArray_N, &
            GlobPosArray_NPlus1, &
            SigmaEffArray_N, &
            SigmaEffArray_NPlus1)
        !**********************************************************************
        !
        !    Function:  Frees the arrays belonging to this module.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IError
          
          ! Inout variables
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: NPartEle_N
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: NPartEle_NPlus1
          
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: EleParticlesHelp_N
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: EleParticlesHelp_NPlus1
          
          integer(kind = 8), dimension(:), allocatable, intent(inout) :: EleParticles_N      
          integer(kind = 8), dimension(:), allocatable, intent(inout) :: EleParticles_NPlus1
          
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: ElementIDArray_N
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: ElementIDArray_NPlus1
          
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: MaterialIDArray_N
          integer(INTEGER_TYPE), dimension(:), allocatable, intent(inout) :: MaterialIDArray_NPlus1
          
          real(REAL_TYPE), dimension(:,:), allocatable, intent(inout) :: ShapeValuesArray_N
          real(REAL_TYPE), dimension(:,:), allocatable, intent(inout) :: ShapeValuesArray_NPlus1
          
          real(REAL_TYPE), dimension(:,:,:), allocatable, intent(inout) :: DShapeValuesArray_N
          real(REAL_TYPE), dimension(:,:,:), allocatable, intent(inout) :: DShapeValuesArray_NPlus1
          
          real(REAL_TYPE), dimension(:,:), allocatable, intent(inout) :: GlobPosArray_N
          real(REAL_TYPE), dimension(:,:), allocatable, intent(inout) :: GlobPosArray_NPlus1
          
          real(REAL_TYPE), dimension(:,:), allocatable, intent(inout) :: SigmaEffArray_N
          real(REAL_TYPE), dimension(:,:), allocatable, intent(inout) :: SigmaEffArray_NPlus1
          
          
          
          ! implicit variables 
          if (allocated(NPartEle_N) ) then !--> implicit bookmark
            deallocate(NPartEle_N, stat = IError)
          end if
          if (allocated(NPartEle_NPlus1) ) then !--> implicit bookmark
            deallocate(NPartEle_NPlus1, stat = IError)
          end if
        
          ! implicit variables 
          if (allocated(EleParticlesHelp_N) ) then
            deallocate(EleParticlesHelp_N, stat = IError)
          end if
          if (allocated(EleParticlesHelp_NPlus1) ) then
            deallocate(EleParticlesHelp_NPlus1, stat = IError)
          end if
          
          ! implicit variables
          if (allocated(EleParticles_N) ) then
            deallocate(EleParticles_N, stat = IError)
          end if
          if (allocated(EleParticles_NPlus1) ) then
            deallocate(EleParticles_NPlus1, stat = IError)
          end if
          
          if (allocated(Particles) ) then !-> this is a huge structure 
            deallocate(Particles, stat = IError)
          end if
         
         !Fileds stripped from Particle datatype
         if (allocated(IDArray) ) then
            deallocate(IDArray, stat = IError)
         end if
                
         ! implicit variable 
          if (allocated(ElementIDArray_N) ) then !--> implicit bookmark 
            deallocate(ElementIDArray_N, stat = IError)
          end if
          if (allocated(ElementIDArray_NPlus1) ) then !--> implicit bookmark
            deallocate(ElementIDArray_NPlus1, stat = IError)
          end if

         if (allocated(EntityIDArray) ) then
            deallocate(EntityIDArray, stat = IError)
         end if
                 
         ! implicit variable
         if (allocated(MaterialIDArray_N) ) then !--> implicit bookmark
            deallocate(MaterialIDArray_N, stat = IError)
         end if  
         if (allocated(MaterialIDArray_NPlus1) ) then !--> implicit bookmark
            deallocate(MaterialIDArray_NPlus1, stat = IError)
         end if  
         
         if (allocated(UArray) ) then
            deallocate(UArray, stat = IError)
         end if  
         
        if (allocated(MaterialPointTypeArray) ) then
            deallocate(MaterialPointTypeArray, stat = IError)
        end if
        
        if (allocated(UStepArray) ) then
            deallocate(UStepArray, stat = IError)
        end if 
        
        ! implicit variable
        if (allocated(ShapeValuesArray_N) ) then !--> implicit bookmark
            deallocate(ShapeValuesArray_N, stat = IError)
        end if
        if (allocated(ShapeValuesArray_NPlus1) ) then !--> implicit bookmark
            deallocate(ShapeValuesArray_NPlus1, stat = IError)
        end if
        
         if (allocated(UPhaseArray) ) then
            deallocate(UPhaseArray, stat = IError)
         end if
                 
         ! implicit variable
         if (allocated(GlobPosArray_N) ) then !-> implicit bookmark
            deallocate(GlobPosArray_N, stat = IError)
         end if
         if (allocated(GlobPosArray_NPlus1) ) then !-> implicit bookmark
            deallocate(GlobPosArray_NPlus1, stat = IError)
         end if
         
         ! implicit scheme variables 
         if (allocated(VelocityArray_N) ) then
            deallocate(VelocityArray_N, stat = IError)
         end if
         
         if (allocated(VelocityArray_NPlus1) ) then
            deallocate(VelocityArray_NPlus1, stat = IError)
         end if
         
         ! implicit variable
         if (allocated(SigmaEffArray_N) ) then !--> implicit bookmark
            deallocate(SigmaEffArray_N, stat = IError)
         end if
         if (allocated(SigmaEffArray_NPlus1) ) then
            deallocate(SigmaEffArray_NPlus1, stat = IError)
         end if
         

         if (allocated(SigmaEff0Array) ) then
            deallocate(SigmaEff0Array, stat = IError) !--> I am not sure if we are needing this...
         end if

         ! implicit variable
         if (allocated(DShapeValuesArray_N) ) then
            deallocate(DShapeValuesArray_N, stat = IError)
         end if
         if (allocated(DShapeValuesArray_NPlus1) ) then
            deallocate(DShapeValuesArray_NPlus1, stat = IError)
         end if
         
         if (allocated(VelocityWaterArray) ) then
            deallocate(VelocityWaterArray, stat = IError)
         end if
         
         if (allocated(VelocityGasArray) ) then
            deallocate(VelocityGasArray, stat = IError)
         end if
         
         if (allocated(AccelerationArray) ) then
            deallocate(AccelerationArray, stat = IError)
         end if
         
         if (allocated(MassArray)) then
            deallocate(MassArray, stat = IError)
         end if
         
         ! implicit scheme variables 
         !if (allocated(ESMstatevArray)) then
         !   deallocate(ESMstatevArray, stat = IError)
         !end if
         if (allocated(ESMstatevArray_N)) then
            deallocate(ESMstatevArray_N, stat = IError)
         end if
         if (allocated(ESMstatevArray_NPlus1)) then
            deallocate(ESMstatevArray_NPlus1, stat = IError)
         end if
         if (allocated(ESMstatevArray_Iteration)) then
            deallocate(ESMstatevArray_Iteration, stat = IError)
         end if
         
         !if (allocated(ESMpropsArray)) then
         !   deallocate(ESMpropsArray, stat = IError)
         !end if
         
         if (allocated(MassWaterArray)) then
            deallocate(MassWaterArray, stat = IError)
         end if
         
          if (allocated(IsParticleIntegration) ) then
            deallocate(IsParticleIntegration, stat = IError)
          end if
          
        end subroutine DestroyHouseKeeping


        subroutine DetermineDecimalFactor()
        !**********************************************************************
        !
        !    Function:  Determines the decimal factor PartFactor which is used
        !               to set up the particle house-keeping arrays. Requires
        !               Counters%NParticles. Writes PartFactor.
        !               Determines the next highest decimal number.
        !               For example: Value = 345 return 1000 etc
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          PartFactor = DecimalFactor(Counters%NParticles)
        
        end subroutine DetermineDecimalFactor


        integer(kind = 8) function DecimalFactor(Number)
        !**********************************************************************
        !
        !    Function:  Determines the next highest decimal number of Number
        !               For example: Value = 345 return 1000 etc
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          integer(INTEGER_TYPE), intent(in) :: Number

          DecimalFactor = 10
          do while (DecimalFactor<=Number)
            DecimalFactor = 10 * DecimalFactor
          end do
        
        end function DecimalFactor


        subroutine CoordLocalToGlobal(IElement, LNodalCoordinates)
        !**********************************************************************
        !
        !    Function:  Determines the global coordinates associated with the local
        !               coordinates of particles inside a considered element.
        !               Updates the global coordinates of particles inside IElement.
        !
        !     IElement : ID of the considered element
        !     LNodalCoordinates : Nodal coordinates used for updating the integration point locations
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        use ModCounters

        implicit none 

          integer(INTEGER_TYPE), intent(in) :: IElement
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR) :: LNodalCoordinates
          ! Local variables
          integer(INTEGER_TYPE) :: NodeID, ParticleIndex
          integer(INTEGER_TYPE) :: I, IParticle, IDim

          do IParticle = 1, NPartEle_N(IElement) ! Loop over particles of IElement
            ParticleIndex = GetParticleIndexFunction(IParticle, IElement)

            GlobPosArray_N(ParticleIndex,:) = 0.0 ! Reset coordinates to zero
            do I = 1, ELEMENTNODES ! Loop over nodes of IElement
              NodeID = iabs(ElementConnectivities(I, IElement) )
              do IDim = 1, NVECTOR ! Loop over dimensions of IElement
                GlobPosArray_N(ParticleIndex,IDim) = GlobPosArray_N(ParticleIndex,IDim) + LNodalCoordinates(NodeID, IDim) * ShapeValuesArray_N(ParticleIndex,I)
              end do
            end do
            
          end do

        end subroutine CoordLocalToGlobal

        subroutine SetParticleShapeFunctionData(Particle, ParticleIndex, ShapeValuesArray, DShapeValuesArray)
        !**********************************************************************
        !
        !    Function:  Calculate the shape function values and derivatives
        !               for the new local position of Particle.
        !
        !     Particle : Considered particle
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

          type(ParticleType(NTENSOR, NVECTOR,MAX_LOAD_SYSTEMS)) :: Particle
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex
          ! Local variables
          real(REAL_TYPE), dimension(ELEMENTNODES) :: ShapeValues !, intent(out)
          real(REAL_TYPE), dimension(ELEMENTNODES, NVECTOR) :: DShapeValues !, intent(out)
          
          		
          ! implicit scheme - shape values array inout 
         real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES), intent(inout) :: ShapeValuesArray 
         real(REAL_TYPE), dimension(Counters%NParticles, ELEMENTNODES, NVECTOR), intent(inout) :: DShapeValuesArray
		 

          call ShapeFunctionData(Particle%LocPos, ELEMENTNODES, ShapeValues, DShapeValues)
          ShapeValuesArray(ParticleIndex, : ) = ShapeValues
          DShapeValuesArray(ParticleIndex, :, : ) = DShapeValues

        end subroutine SetParticleShapeFunctionData


        subroutine DetermineAdjacentParticles(ISide, NParticles, ParticleStatus)
        !**********************************************************************
        !
        !    Function:  Determines which particles from 1 to NParticles
        !               of the element lie next to ISide and how many.
        !               Note: ParticleStatus is not initialised to .false. in order
        !                     to allow for a more flexible usage! (Which particles of
        !                     an element lie adjacent to any side of the element? )
        !
        !     ISide : Local number of considered element side
        !     NParticles : Number of particles inside the considered element
        !
        ! O   ParticleStatus : Set to .true., if the particle lies next to ISide
        !
        !**********************************************************************
        implicit none

          integer(INTEGER_TYPE), intent(in) :: ISide, NParticles
          logical, dimension(NParticles), intent(inout) :: ParticleStatus

            select case(ELEMENTTYPE)
                case(TRI3)
                    select case (NParticles)  
                        case(1) ! 1 material point per elmeent / FEM calculation
                          call DetermineAdjacentParticlesTRI_MP1(NParticles, ParticleStatus)
                        case(3) ! 3 material points per element
                          call DetermineAdjacentParticlesTRI_MP3(ISide, NParticles, ParticleStatus)
                        case(6) ! 6 material points per element
                            call DetermineAdjacentParticlesTRI_MP6(ISide, NParticles, ParticleStatus)
                        case(12) ! 12 material points per element
                            call DetermineAdjacentParticlesTRI_MP12(ISide, NParticles, ParticleStatus)
                        case(16) ! 16 material points per element
                            call DetermineAdjacentParticlesTRI_MP16(ISide, NParticles, ParticleStatus)
                        case(25) ! 25 material points per element
                            call DetermineAdjacentParticlesTRI_MP25(ISide, NParticles, ParticleStatus)
                        case(46) ! 46 material points per element
                            call DetermineAdjacentParticlesTRI_MP46(ISide, NParticles, ParticleStatus)
                        case(88) ! 88 material points per element
                            call DetermineAdjacentParticlesTRI_MP88(ISide, NParticles, ParticleStatus)
                        case default
                          call GiveError("Undefined number of material points per element in [subroutine DetermineAdjacentParticles()], 2D case.")
                    end select
                case(TETRAOLD)
                    select case (NParticles)
                        case(1) ! 1 material point per element / FEM calculation
                            call DetermineAdjacentParticlesTETRA_MP1(NParticles, ParticleStatus)
                        case(4) ! 4 material points per element
                            call DetermineAdjacentParticlesTETRA_MP4(ISide, NParticles, ParticleStatus)
                        case(7) ! 7 material points per element
                            call DetermineAdjacentParticlesTETRA_MP7(ISide, NParticles, ParticleStatus)
                        case(8) ! 8 material points per element
                            call DetermineAdjacentParticlesTETRA_MP8(ISide, NParticles, ParticleStatus)
                        case(10) ! 10 material points per element
                            call DetermineAdjacentParticlesTETRA_MP10(ISide, NParticles, ParticleStatus)
                        case(13) ! 13 material points per element
                            call DetermineAdjacentParticlesTETRA_MP13(ISide, NParticles, ParticleStatus)
                        case(20) ! 20 material points per element
                            call DetermineAdjacentParticlesTETRA_MP20(ISide, NParticles, ParticleStatus)
                        case default
                          call GiveError("Undefined number of material points per element in [subroutine DetermineAdjacentParticles()], 3D case.")
                    end select
                end select

        end subroutine DetermineAdjacentParticles

        
        integer(INTEGER_TYPE) function DetermineNAdjacentParticles(ParticleStatus, NParticles)
        !**********************************************************************
        !
        !    Function:  Determines how many of the particles from 1 to Flags%NMaterialPoints
        !               of the element lie next to a side by evaluating ParticleStatus.
        !
        !     ParticleStatus : Array storing for each particle whether it was found
        !                      to lie next to a checked element side
        !
        ! O   DetermineNAdjacentParticles : Number of particles adjacent to side(s) of 
        !                                   the element
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
         integer(INTEGER_TYPE), intent (in):: NParticles
          logical, dimension(NParticles), intent(in) :: ParticleStatus
          ! Local variables
          integer(INTEGER_TYPE) :: I

          DetermineNAdjacentParticles = 0
          
          do I = 1, NParticles
            if (ParticleStatus(I) ) then
              DetermineNAdjacentParticles = DetermineNAdjacentParticles + 1
            end if
          end do

        end function DetermineNAdjacentParticles

        subroutine DefineVariableData()
        !**********************************************************************
        !
        !    Function:  Sets Particles(I)%VariableData(1 .. 3) fields. This fields can be
        !               used for debugging purposes to display some particle data
        !               with the Output program.
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          ! Local variables
          integer(INTEGER_TYPE) :: IParticle, ParticleIndex, IElement
        
          do IElement = 1, Counters%NEl
        
            do IParticle = 1, NPartEle_N(IElement)

              ParticleIndex = GetParticleIndex_N(IParticle, IElement)

              ! Indicate partially and fully filled elements -> VariableData(1)
              if (IsParticleIntegration(IElement) ) then
                call SetVariableDataI(Particles(ParticleIndex), 1, 1.d0)
              else
                call SetVariableDataI(Particles(ParticleIndex), 1, 0.d0)
              end if

              call SetVariableDataI(Particles(ParticleIndex), 2, Particles(ParticleIndex)%IntegrationWeight)

              call SetVariableDataI(Particles(ParticleIndex), 3, (Particles(ParticleIndex)%Density_N))

            end do
          
          end do
        
        end subroutine DefineVariableData        

      end module ModMPMData
