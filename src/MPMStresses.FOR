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
	  
	  
	  module ModMPMStresses
      !**********************************************************************
      !
      !    Function: Contains the routines related to calculating the stresses 
      !              All routines here are shared by both the quasi-static and dynamic MPM.  
      !
      ! Implemented in the frame of the MPM project.
      !
      !     $Revision: 9928 $
      !     $Date: 2023-04-05 16:14:43 +0200 (wo, 05 apr 2023) $
      !
      !**********************************************************************
      use ModCounters
      use ModElementEvaluation
      use ModMPMData
      use ModMeshInfo
      use ModGeometryMath
      use ModReadMaterialData
      use ModGlobalConstants
      use ModReadGeometryData
      
      contains ! Routines of this module

        subroutine Hill(IEl, IElTyp, Disp, NDoFEx, ICon, B, Sig0, SigE, DEpsV)
      !**********************************************************************
      !
      !  Function:  apply objective stresses
      !             
      !
      !**********************************************************************
        implicit none
          
          integer(INTEGER_TYPE), intent(in) :: IEl, IElTyp
          real(REAL_TYPE), dimension(:), intent(in) :: Disp
          integer(INTEGER_TYPE), dimension(:), intent(in) :: NDoFEx
          integer(INTEGER_TYPE), dimension(:, :), intent(in) :: ICon
          real(REAL_TYPE), dimension(:,:), intent(in) :: B
          real(REAL_TYPE), dimension(:), intent(in) :: Sig0
          real(REAL_TYPE), intent(in) :: DEpsV
          real(REAL_TYPE), dimension(:), intent(inout) :: SigE
          ! Local variables
          real(REAL_TYPE) :: EpsV, Om12, Om23, Om31, Ux, Uy, Uz
          integer(INTEGER_TYPE) :: J, NN
          
          Om12 = 0.0
          Om23 = 0.0
          Om31 = 0.0
          
		  select case(NTENSOR)
            case(4)
                do J = 1, IElTyp
                    NN = ICon(J, IEl)
                    Ux = Disp(NDofEx(NN) + 1)
                    Uy = Disp(NDofEx(NN) + 2)
                    Om12 = Om12 + B(1, J) * Uy - B(2, J) * Ux  
                end do
                
                    Om12 = Om12 / 2.0
                    EpsV = DEpsV

                    SigE(1) = SigE(1) - EpsV * Sig0(1) - 2 * Sig0(4) * Om12
                    SigE(2) = SigE(2) - EpsV * Sig0(2) + 2 * Sig0(4) * Om12
                    SigE(3) = SigE(3) - EpsV * Sig0(3)
                    SigE(4) = SigE(4) - EpsV * Sig0(4) + (Sig0(1) - Sig0(2) ) * Om12
                
            case(6)
				do J = 1, IElTyp
            		NN = ICon(J, IEl)
            		Ux = Disp(NDofEx(NN) + 1)
            		Uy = Disp(NDofEx(NN) + 2)
            		Uz = Disp(NDofEx(NN) + 3)
            		Om12 = Om12 + B(1, J) * Uy - B(2, J) * Ux
            		Om23 = Om23 + B(2, J) * Uz - B(3, J) * Uy
            		Om31 = Om31 + B(3, J) * Ux - B(1, J) * Uz
          		end do
    
          			Om12 = Om12 / 2.0
          			Om23 = Om23 / 2.0
          			Om31 = Om31 / 2.0
          			EpsV = DEpsV

          			SigE(1) = SigE(1) - EpsV * Sig0(1) - 2 * Sig0(4) * Om12 + 2 * Sig0(6) * Om31
          			SigE(2) = SigE(2) - EpsV * Sig0(2) + 2 * Sig0(4) * Om12 - 2 * Sig0(5) * Om23
          			SigE(3) = SigE(3) - EpsV * Sig0(3) - 2 * Sig0(6) * Om31 + 2 * Sig0(5) * Om23
          			SigE(4) = SigE(4) - EpsV * Sig0(4) + (Sig0(1) - Sig0(2) ) * Om12 - Sig0(6) * Om23 + Sig0(5) * Om31
          			SigE(5) = SigE(5) - EpsV * Sig0(5) + Sig0(6) * Om12 + (Sig0(2) - Sig0(3) ) * Om23 - Sig0(4) * Om31
          			SigE(6) = SigE(6) - EpsV * Sig0(6) - Sig0(5) * Om12 + Sig0(4) * Om23 + (Sig0(3) - Sig0(1) ) * Om31
			end select
        end subroutine Hill

        subroutine DetermineCoordRelativeToSurface(SoilPoint, Depth)
        !**********************************************************************
        !
        !    Function: SoilPoint is an arbitrary point of the discretisation
        !              For the simplified computation of initial stresses the depth of the
        !              point below the soil surface is required and, in case the soil surface
        !              has the shape of a wave, also the distance on the soil surface from some
        !              point specified as an origin for the wave function.
        !
        !              While trivial for the case that the gravity direction is aligned with the y-direction,
        !              here also the case is taken into consideration that the gravity direction is rotated.
        !              Here, it is assumed that the surface is rotated around the z-axis, so that the direction
        !              vector of the surface lies in the x-y-plane. Furthermore, it is assumed that the soil surface
        !              is oriented perpendicular to the gravity direction.
        !
        !   SoilPoint : Arbitrary point (material point, node)
        !
        ! O Depth : Depth of SoilPoint relative to the surface specified as described above
        ! 
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(2), intent(in) :: SoilPoint
          real(REAL_TYPE), intent(out) :: Depth
          ! Local variables
          real(REAL_TYPE) :: Determinante, SurfaceCoordinate
          real(REAL_TYPE), dimension(2) :: SurfaceNormVector, DepthVector, SurfaceVector, SurfacePoint
          integer(INTEGER_TYPE) :: I
          logical :: IsZeroGravity
          
          IsZeroGravity = .true.
          do I = 1, NVECTOR
            IsZeroGravity = IsZeroGravity .and. (CalParams%GravityData%GravityVector(I) == 0.0)
          end do
          if (IsZeroGravity) RETURN          
          
          Depth = 0.0
          
          ! Determine vector parallel to surface (it might not be horizontal) assuming it lies in the x-y-plane
          ! and is perpendicular to the gravity vector
          SurfaceNormVector(1) = -1.0 * CalParams%GravityData%GravityVector(2)
          SurfaceNormVector(2) = CalParams%GravityData%GravityVector(1)
          
          ! Compute distance of material column above material point from specified origin
          SurfaceVector = SoilPoint - CalParams%SoilSurfacePoint
          SurfaceCoordinate = SurfaceVector(1) * SurfaceNormVector(1) + SurfaceVector(2) * SurfaceNormVector(2)
          ! Compute height of material column above material point
          SurfacePoint = CalParams%SoilSurfacePoint + SurfaceCoordinate * SurfaceNormVector
          DepthVector = SurfacePoint - SoilPoint
          Depth = sqrt(DepthVector(1) * DepthVector(1) + DepthVector(2) * DepthVector(2))
          
          ! Check whether point lies above soil surface (make its depth negative in this case)
          ! TODO: This does not work yet properly for wave-shaped surfaces (rather exceptional case)
          Determinante = SurfaceNormVector(1) * (SoilPoint(2) - CalParams%SoilSurfacePoint(2)) -  &
                         SurfaceNormVector(2) * (SoilPoint(1) - CalParams%SoilSurfacePoint(1))
          if (Determinante<=0.0) then
            Depth = -1.0 * Depth
          end if
          
          if ((Depth>=0.0)) then
            Depth = 0.0
          end if

          end subroutine DetermineCoordRelativeToSurface
          
          subroutine DetermineCoordRelativeToDefinedSurface(SoilPoint,CoordMatrix, Distance, NewSurfaceReference)
        !**********************************************************************
        !
        !    Function: SoilPoint is an arbitrary point of the discretisation
        !              For the simplified computation of initial stresses the depth of the
        !              point below the soil surface is required.
        !
        !   SoilPoint : Arbitrary point (material point, node)
        !   CoordMatrix: nodes coordinates belonging to Soil Surface or Phreatic surface, defined in the GOM file
        !   NewSurfaceReference: reference point on the soil/phreatic surface, necessary to compute the distance from the  SoilPoint;
        !                        it is obtained through linear interpolation
        !
        ! Distance : Depth of SoilPoint relative to the surface point, specified as described above
        ! 
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(2), intent(in) :: SoilPoint
          real(REAL_TYPE),dimension(:,:), intent(in) :: CoordMatrix
              real(REAL_TYPE), intent(out) :: Distance
              real(REAL_TYPE), dimension(2), intent(out):: NewSurfaceReference
          ! Local variables
          real(REAL_TYPE) ::  SurfaceCoordinate 
          
          real(REAL_TYPE) :: Elevation
          real(REAL_TYPE), dimension(2) :: SurfaceNormVector, DepthVector, SurfaceVector, SurfacePoint
          integer(INTEGER_TYPE) :: I, SurfaceNodesCounter
          logical :: IsZeroGravity
          
          IsZeroGravity = .true.
          do I = 1, NVECTOR
            IsZeroGravity = IsZeroGravity .and. (CalParams%GravityData%GravityVector(I) == 0.0)
          end do
          if (IsZeroGravity) RETURN          
          
          Distance = 0.0
          Elevation = 0.0
          NewSurfaceReference = 0.0       
         
          ! Determine vector parallel to surface (it might not be horizontal) assuming it lies in the x-y-plane
          ! and is perpendicular to the gravity vector
          SurfaceNormVector(1) = -1.0 * CalParams%GravityData%GravityVector(2)
          SurfaceNormVector(2) = CalParams%GravityData%GravityVector(1)
                        
          SurfaceNodesCounter = size(CoordMatrix,1)
          do I = 1,  SurfaceNodesCounter
              
                  if (NDIM == 2) then
                  if (CoordMatrix(I,1) > SoilPoint(1)) then
                      ! interpolated delta-y value (2D) (INTERPOLATION considering a line)
                  if (I==1) then
                  Elevation = CoordMatrix(I,2) 
                                    else
                                        Elevation = CoordMatrix(I-1,2) + (CoordMatrix(I,2)-CoordMatrix(I-1,2)) / (CoordMatrix(I,1)- &
                          CoordMatrix(I-1,1)) * (SoilPoint(1) - CoordMatrix(I-1,1))
                          end if
                      exit
                      else
                    Elevation = CoordMatrix(SurfaceNodesCounter,2)   
                 end if
                   else if (NDIM == 3) then
                   if (CoordMatrix(I,1) > SoilPoint(1)) then
                      ! interpolated delta-z value (3D) (INTERPOLATION considering a plane)
                    if (I==1) then
                  Elevation =  CoordMatrix(I,3) 
                                    else
                      Elevation = CoordMatrix(i-1,3) + & 
                      ((SoilPoint(2)-CoordMatrix(i-1,2)) * ( (CoordMatrix(i,1)-CoordMatrix(i-1,1)) * (CoordMatrix(i+1,3)-CoordMatrix(i-1,3))- &
                      (CoordMatrix(i+1,1)-CoordMatrix(i-1,1)) * (CoordMatrix(i,3)-CoordMatrix(i-1,3)))- &
                      (SoilPoint(1)-CoordMatrix(i-1,1)) *  ( (CoordMatrix(i,2)-CoordMatrix(i-1,2)) * (CoordMatrix(i+1,3)-CoordMatrix(i-1,3))- &
                      (CoordMatrix(i,3)-CoordMatrix(i-1,3)) *  (CoordMatrix(i+1,2)-CoordMatrix(i-1,2))))/ &
                      ((CoordMatrix(i,1)-CoordMatrix(i-1,1))*(CoordMatrix(i+1,2)-CoordMatrix(i-1,2))- (CoordMatrix(i,2)-CoordMatrix(i-1,2))*(CoordMatrix(i+1,1)-CoordMatrix(i-1,1)))
                      end if
                      exit
                      else
                    Elevation = CoordMatrix(SurfaceNodesCounter,3)   
                  end if
              end if
          end do
          
          NewSurfaceReference(1) = 0.0
          NewSurfaceReference(2) = Elevation
          
          ! Compute distance of material column above material point from specified origin
          SurfaceVector = SoilPoint - NewSurfaceReference
          SurfaceCoordinate = SurfaceVector(1) * SurfaceNormVector(1) + SurfaceVector(2) * SurfaceNormVector(2)
          ! Compute height of material column above material point
          SurfacePoint = NewSurfaceReference + SurfaceCoordinate * SurfaceNormVector
          DepthVector = SurfacePoint - SoilPoint
          Distance = sqrt(DepthVector(1) * DepthVector(1) + DepthVector(2) * DepthVector(2))
          
               
        end subroutine DetermineCoordRelativeToDefinedSurface
          
        subroutine DetermineCoordRelativeToSurfaceDefinedInExternalFile(MaterialIndex, SoilPoint, Depth, LiquidSurfacePoint)
        !**********************************************************************
        !
        !    Function: SoilPoint is an arbitrary point of the discretisation
        !              For the simplified computation of initial stresses the depth of the
        !              point below the soil surface is required.
        !              
        !              Surface coordinates are read from a text file whose name and location are specified in GOM file
        !
        !   SoilPoint : Arbitrary point (material point, node)
        !
        ! O Depth : Depth of SoilPoint relative to the surface specified as described above
        ! 
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialIndex
          real(REAL_TYPE), dimension(2), intent(in) :: SoilPoint
          real(REAL_TYPE), intent(out) :: Depth
          real(REAL_TYPE), dimension(2), intent(out) ::LiquidSurfacePoint
          ! Local variables
          real(REAL_TYPE) :: SurfaceCoordinate
          real(REAL_TYPE), dimension(2) :: SurfaceNormVector, DepthVector, SurfaceVector, SurfacePoint
          integer(INTEGER_TYPE) :: I
          logical :: IsZeroGravity
          
          character(len = 255) :: FileName
          character(len = 1) :: FileNumber
          integer(INTEGER_TYPE) :: Npoint
          real(REAL_TYPE), allocatable, dimension(:,:) :: DumI
          
          IsZeroGravity = .true.
          do I = 1, NVECTOR
            IsZeroGravity = IsZeroGravity .and. (CalParams%GravityData%GravityVector(I) == 0.0)
          end do
          if (IsZeroGravity) RETURN          
          
          Depth = 0.0
          LiquidSurfacePoint = 0.0
          
          ! Determine vector parallel to surface (it might not be horizontal) assuming it lies in the x-y-plane
          ! and is perpendicular to the gravity vector
          SurfaceNormVector(1) = -1.0 * CalParams%GravityData%GravityVector(2)
          SurfaceNormVector(2) = CalParams%GravityData%GravityVector(1)
          
          ! Compute distance of material column above material point from surface defined in external file
          do I = 1, GeoParams%NumberWaterSurfaceMaterials
               if (MaterialIndex == GeoParams%WaterSurfaceMaterialID(I))then  
                   FileNumber = GeoParams%WaterSurfaceFileNumber(I)
               end if
          end do

          FileName=trim(CalParams%FileNames%ProjectName)//PHREATIC_SURFACE_FILE_EXTENSION//trim(FileNumber)
            
            open(PSFUnit, FILE = FileName)
            read(PSFUnit,*) Npoint
            allocate(DumI(2,Npoint))
            do I = 1, Npoint                
                read(PSFUnit,*) DumI(1,I), DumI(2,I)
                if (I > 1)then
                if ( (SoilPoint(1) <= DumI(1,I)) .and. SoilPoint(1) >= DumI(1,I-1) )then
                    LiquidSurfacePoint(1) = 0.0d0 
                    LiquidSurfacePoint(2) = (( DumI(2,I) - DumI(2,I-1) )/( DumI(1,I) - DumI(1,I-1) )) *(SoilPoint(1) - DumI(1,I-1)) + DumI(2,I-1)
                end if 
                else
                if ( (SoilPoint(1) >= DumI(1,I)) )then
                    LiquidSurfacePoint(1) = 0.0d0 
                    LiquidSurfacePoint(2) = DumI(2,I)
                end if
                end if
            end do 
            deallocate(DumI)
            close(PSFUnit) 
          
          SurfaceVector = SoilPoint - LiquidSurfacePoint
          SurfaceCoordinate = SurfaceVector(1) * SurfaceNormVector(1) + SurfaceVector(2) * SurfaceNormVector(2)  
          
          ! Compute height of water column above material point
          SurfacePoint = LiquidSurfacePoint + SurfaceCoordinate * SurfaceNormVector
          DepthVector = SurfacePoint - SoilPoint
          Depth = sqrt(DepthVector(1) * DepthVector(1) + DepthVector(2) * DepthVector(2))
          
         end subroutine DetermineCoordRelativeToSurfaceDefinedInExternalFile

          subroutine DetermineCoordRelativeToSurface2Layers(SoilPoint, Depth, ParticleIndex)
        !**********************************************************************
        !
        !    Function: SoilPoint is an arbitrary point of the discretisation
        !              For the simplified computation of initial stresses the depth of the
        !              point below the soil surface is required and, in case the soil surface
        !              has the shape of a wave, also the distance on the soil surface from some
        !              point specified as an origin for the wave function.
        !
        !              While trivial for the case that the gravity direction is aligned with the y-direction,
        !              here also the case is taken into consideration that the gravity direction is rotated.
        !              Here, it is assumed that the surface is rotated around the z-axis, so that the direction
        !              vector of the surface lies in the x-y-plane. Furthermore, it is assumed that the soil surface
        !              is oriented perpendicular to the gravity direction.
        !
        !   SoilPoint : Arbitrary point (material point, node)
        !
        ! O Depth : Depth of SoilPoint relative to the surface specified as described above
        ! 
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(2), intent(in) :: SoilPoint
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex
          real(REAL_TYPE), intent(out) :: Depth
          ! Local variables
          real(REAL_TYPE) :: Determinante, SurfaceCoordinate
          real(REAL_TYPE), dimension(2) :: SurfaceNormVector, DepthVector, SurfaceVector, SurfacePoint
          integer(INTEGER_TYPE) :: I
          logical :: IsZeroGravity
          
          IsZeroGravity = .true.
          do I = 1, NVECTOR
            IsZeroGravity = IsZeroGravity .and. (CalParams%GravityData%GravityVector(I) == 0.0)
          end do
          if (IsZeroGravity) RETURN 
          
          Depth = 0.0
          
          ! Determine vector parallel to surface (it might not be horizontal) assuming it lies in the x-y-plane
          ! and is perpendicular to the gravity vector
          SurfaceNormVector(1) = -1.0 * CalParams%GravityData%GravityVector(2)
          SurfaceNormVector(2) = CalParams%GravityData%GravityVector(1)


              if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
                ! Compute distance of material column above material point from specified origin
                SurfaceVector = SoilPoint - CalParams%LiquidSurfacePoint
                SurfaceCoordinate = SurfaceVector(1) * SurfaceNormVector(1) + SurfaceVector(2) * SurfaceNormVector(2)
                ! Compute height of material column above material point
                SurfacePoint = CalParams%LiquidSurfacePoint + SurfaceCoordinate * SurfaceNormVector
                DepthVector = SurfacePoint - SoilPoint
                Depth = sqrt(DepthVector(1) * DepthVector(1) + DepthVector(2) * DepthVector(2))
                Determinante = SurfaceNormVector(1) * (SoilPoint(2) - CalParams%LiquidSurfacePoint(2)) -  &
                         SurfaceNormVector(2) * (SoilPoint(1) - CalParams%LiquidSurfacePoint(1))   
              end if
                
              if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then
                ! Compute distance of material column above material point from specified origin
                SurfaceVector = SoilPoint - CalParams%SoilSurfacePoint
                SurfaceCoordinate = SurfaceVector(1) * SurfaceNormVector(1) + SurfaceVector(2) * SurfaceNormVector(2)
                ! Compute height of material column above material point
                SurfacePoint = CalParams%SoilSurfacePoint + SurfaceCoordinate * SurfaceNormVector
                DepthVector = SurfacePoint - SoilPoint
                Depth = sqrt(DepthVector(1) * DepthVector(1) + DepthVector(2) * DepthVector(2))
                Determinante = SurfaceNormVector(1) * (SoilPoint(2) - CalParams%SoilSurfacePoint(2)) -  &
                         SurfaceNormVector(2) * (SoilPoint(1) - CalParams%SoilSurfacePoint(1))   
              end if
                
              if (Determinante<=0.0) then
                Depth = -1.0 * Depth
              end if

          end subroutine DetermineCoordRelativeToSurface2Layers


          subroutine DetermineCoordRelativeToSurface2LayersSolid(SoilPoint, Depth)
        !**********************************************************************
        !
        !    Function: SoilPoint is an arbitrary point of the discretisation
        !              For the simplified computation of initial stresses the depth of the
        !              point below the soil surface is required and, in case the soil surface
        !              has the shape of a wave, also the distance on the soil surface from some
        !              point specified as an origin for the wave function.
        !
        !              While trivial for the case that the gravity direction is aligned with the y-direction,
        !              here also the case is taken into consideration that the gravity direction is rotated.
        !              Here, it is assumed that the surface is rotated around the z-axis, so that the direction
        !              vector of the surface lies in the x-y-plane. Furthermore, it is assumed that the soil surface
        !              is oriented perpendicular to the gravity direction.
        !
        !   SoilPoint : Arbitrary point (material point, node)
        !
        ! O Depth : Depth of SoilPoint relative to the surface specified as described above
        ! 
        !**********************************************************************
        
        implicit none
        
          real(REAL_TYPE), dimension(2), intent(in) :: SoilPoint
          real(REAL_TYPE), intent(out) :: Depth
          ! Local variables
          real(REAL_TYPE) :: Determinante, SurfaceCoordinate
          real(REAL_TYPE), dimension(2) :: SurfaceNormVector, DepthVector, SurfaceVector, SurfacePoint
          integer(INTEGER_TYPE) :: I
          logical :: IsZeroGravity
          
          IsZeroGravity = .true.
          do I = 1, NVECTOR
            IsZeroGravity = IsZeroGravity .and. (CalParams%GravityData%GravityVector(I) == 0.0)
          end do
          if (IsZeroGravity) RETURN 
          
          Depth = 0.0
                 
              ! Determine vector parallel to surface (it might not be horizontal) assuming it lies in the x-y-plane
              ! and is perpendicular to the gravity vector
              SurfaceNormVector(1) = -1.0 * CalParams%GravityData%GravityVector(2)
              SurfaceNormVector(2) = CalParams%GravityData%GravityVector(1)
          
              ! Compute distance of material column above material point from specified origin
              SurfaceVector = SoilPoint - CalParams%SoilSurfacePoint
              SurfaceCoordinate = SurfaceVector(1) * SurfaceNormVector(1) + SurfaceVector(2) * SurfaceNormVector(2)
              ! Compute height of material column above material point
              SurfacePoint = CalParams%SoilSurfacePoint + SurfaceCoordinate * SurfaceNormVector
              DepthVector = SurfacePoint - SoilPoint
              Depth = sqrt(DepthVector(1) * DepthVector(1) + DepthVector(2) * DepthVector(2))
              Determinante = SurfaceNormVector(1) * (SoilPoint(2) - CalParams%SoilSurfacePoint(2)) -  &
                         SurfaceNormVector(2) * (SoilPoint(1) - CalParams%SoilSurfacePoint(1))   
              
                
              if (Determinante<=0.0) then
                Depth = -1.0 * Depth
              end if
  
          end subroutine DetermineCoordRelativeToSurface2LayersSolid
          
          subroutine GenerateBoreholeInfo()
		!**********************************************************************
        !
        !    Function: It fills a system of matrices that will contain compacted information
        !              about layer thickness, densities, and stresses required for a correct
        !              K0 initialization.
        !**********************************************************************
		  Implicit none
		  !Variable initialization
		  integer(INTEGER_TYPE) :: I, MatIndex, Iel, IPart, NodeID, NElemPart, ParticleIndex
		  real(REAL_TYPE) :: TopCoord, BotCoord, GammaMixture, GammaLiquid
		  real(REAL_TYPE), Dimension(NDIM) :: PartCoord
		  logical:: IsUndrEffectiveStress, MatFound
		  		  
		   do I = 1 ,CalParams%NumberSoilLayers
			   CalParams%BoreHoleData(I,1)=CalParams%ThicknessSoilLayer(I)
			   ! Top Coordinate
			   if (I==1) then
			    TopCoord=CalParams%SoilSurfacePoint(2)
			   else
				TopCoord=CalParams%BoreHoleData(I-1,2)-CalParams%BoreHoleData(I-1,1)
			   end if
			   CalParams%BoreHoleData(I,2)=TopCoord
			   ! Bottom coordinate
			   BotCoord=CalParams%BoreHoleData(I,2)-CalParams%BoreHoleData(I,1)
			   CalParams%BoreHoleData(I,3)=BotCoord
			   ! Get the Gamma mixture
			   
			   !get an element inside the zone of interest
			   Matfound=.false.
			   do Iel= 1 , Counters%NEl! loop trough elements
				   NElemPart = NPartEle(Iel)
				   do IPart= 1 ,  NElemPart ! loop over material points of the element
					  ParticleIndex = GetParticleIndex(IPart, Iel)
					  PartCoord=GlobPosArray(ParticleIndex,1:NVECTOR)
					  if ( PartCoord(2)<CalParams%BoreHoleData(I,2) .and. &
						  PartCoord(2)>CalParams%BoreHoleData(I,3) ) then !Node is inside layer
						  if (EntityIDArray(ParticleIndex) /= CalParams%RigidBody%RigidEntity) then
						  MatIndex=ElementMaterialID(Iel)
						  Matfound=.true.
						  EXIT
						  end if
					  end if					  
				   end do
				   if (Matfound) EXIT
			   end do
			   
			if(NFORMULATION==1) then
							  IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MatIndex)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MatIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
							  
				if (CalParams%ApplySubmergedCalculation) then
				GammaMixture = MatParams(MatIndex)%WeightSubmerged
				GammaLiquid = 0.0
				else if (MatParams(MatIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MatIndex)%MaterialPhases=='1-phase-liquid') then
				GammaMixture = MatParams(MatIndex)%WeightLiquid 
				!GammaLiquid = MatParams(MaterialIndex)%WeightLiquid
				GammaLiquid = 0.0 !fluid is the only one material and therefore thought as MIXTURE with only 1 material
				else
				GammaLiquid = MatParams(MatIndex)%WeightLiquid
				if ((CalParams%NumberOfPhases==1).and.(.not.IsUndrEffectiveStress)) then
				 GammaMixture = MatParams(MatIndex)%DryWeight
				else
					GammaMixture = MatParams(MatIndex)%WeightMixture
				end if
				end if
			else ! 2 Constituents
				if(MatParams(MatIndex)%MaterialType=='2-phase'.or.MatParams(MatIndex)%MaterialPhases=='2-phase') then
					GammaMixture = MatParams(MatIndex)%DryWeight  + MatParams(MatIndex)%WeightLiquid *  &
							   MatParams(MatIndex)%InitialPorosity   ! Soil Density * SoilConcentrationRatio
					GammaLiquid = MatParams(MatIndex)%WeightLiquid
				end if
				if(MatParams(MatIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MatIndex)%MaterialPhases=='1-phase-liquid') then
				 GammaMixture = MatParams(MatIndex)%WeightLiquid
				 GammaLiquid = MatParams(MatIndex)%WeightLiquid
				end if
			end if
			CalParams%BoreHoleData(I,4)=GammaMixture
			CalParams%BoreHoleData(I,5)=GammaLiquid
			! compute total effective stress at the interfaces
			CalParams%BoreHoleData(I,6)=-CalParams%BoreHoleData(I,1) * &
				GammaMixture * CalParams%Multipliers%GravityRealised
			CalParams%BoreHoleData(I,7)=-CalParams%BoreHoleData(I,1) * &
				GammaLiquid * CalParams%Multipliers%GravityRealised
		   end do
		   
		   do I= 1 , CalParams%NumberSoilLayers
			   ! Total stresses
			   if (I==1) then
				   CalParams%BoreHoleData(I,6)=CalParams%BoreHoleData(I,6) + &
					   CalParams%InitialVerticalLoadK0 * CalParams%Multipliers%SolidARealised(1)
			   else
				   CalParams%BoreHoleData(I,6)=CalParams%BoreHoleData(I,6)+&
					   CalParams%BoreHoleData(I-1,6)
			   end if
				! Pore pressures
			   if (I>1) then
			   CalParams%BoreHoleData(I,7)=CalParams%BoreHoleData(I,7)+&
				   CalParams%BoreHoleData(I-1,7)
			   end if
			   
			   ! Effective stresses
			   CalParams%BoreHoleData(I,8)=CalParams%BoreHoleData(I,6)-&
				   CalParams%BoreHoleData(I,7)
		   end do		  

  end subroutine GenerateBoreholeInfo
              
        
        subroutine ComputeVerticalStressForK0(MaterialIndex, SoilPoint, UpdatedFluidDensity, SigY, SigWP, SoilSurfCoord, PhreatSurfCoord)
        !**********************************************************************
        !
        !    Function: compute vertical stresses for K0 procedure
        ! 
        !**********************************************************************        
       
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: MaterialIndex ! Material Index 
          real(REAL_TYPE), dimension(2), intent(in) :: SoilPoint ! Global position of material point / node
          real(REAL_TYPE), intent(out) :: UpdatedFluidDensity, SigY, SigWP
          real(REAL_TYPE), dimension(:,:), intent(in) :: SoilSurfCoord, PhreatSurfCoord
          ! Local variables
          
          integer(INTEGER_TYPE) :: I, MatInLayer
          real(REAL_TYPE) :: GammaMixture, GammaLiquid
          real(REAL_TYPE) :: Depth, SoilSurfDepth, PhreaticSurfDepth,  Elevation, SDepth
          real(REAL_TYPE), dimension(2):: NewReferenceSoil, NewReferencePhreatic
          !real(REAL_TYPE) :: Depth
          !integer(INTEGER_TYPE) :: SoilSurfaceNumberofSides, PhreaticSurfaceNumberofSides
          logical :: IsUndrEffectiveStress, IsUndrTotalStress
          
          
          SoilSurfDepth = 0.0
          PhreaticSurfDepth = 0.0
          SigY = 0.0
          SigWP = 0.0
          
          IsUndrEffectiveStress = &
              !code version 2016 and previous
              ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(MaterialIndex)%MaterialType)=='2-phase')) .or. &
              !code version 2017.1 and following
              (trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
          
          IsUndrTotalStress = trim(MatParams(MaterialIndex)%MaterialType)==SATURATED_SOIL_UNDRAINED_TOTAL
    
          ! Determine depth of point below
          call DetermineCoordRelativeToSurface(SoilPoint, Depth)
          
          if (CalParams%NumberSoilLayers>1) then !using multilayer K0 for CPT
            if (Depth/=0.0) then		  
		      !Determine soil layer
		      Elevation=CalParams%SoilSurfacePoint(2)+Depth
		      do I= 1 , CalParams%NumberSoilLayers
			      if (Elevation <= CalParams%BoreHoleData(I,2) .and. &
				      Elevation > CalParams%BoreHoleData(I,3)) then
				      MatInLayer=I
				      EXIT
			      end if				  
		      end do
		  
		      SDepth=Elevation-CalParams%BoreHoleData(MatInLayer,2)
              ! Retrieve possibly depth-dependent material parameters
              if(NFORMULATION==1) then
			      GammaMixture=CalParams%BoreHoleData(MatInLayer,4)
			      GammaLiquid=CalParams%BoreHoleData(MatInLayer,5)  
              end if
          
              ! Compute vertical total stress for gravity loading
              ! Zero stresses for points above the specified material surface (positive depth is point above surface)
              if(MatParams(MaterialIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') then
                  UpdatedFluidDensity  = GammaMixture / CalParams%GravityData%GAccel * &
                            exp((CalParams%Multipliers%GravityRealised * CalParams%GravityData%GAccel) /  &
                             MatParams(MaterialIndex)%BulkModulusLiquid * (- Depth))
     
                 SigY = - MatParams(MaterialIndex)%BulkModulusLiquid * ( UpdatedFluidDensity -  &
                            GammaMixture / CalParams%GravityData%GAccel)
		      else
			      If (MatInLayer>1) then
                 SigY = CalParams%BoreHoleData(MatInLayer-1,6)+&
				     SDepth * GammaMixture * CalParams%Multipliers%GravityRealised
			      else
				     SigY = SDepth * GammaMixture * CalParams%Multipliers%GravityRealised 
			      end if
			  
		      end if
     
              ! Compute vertical total stress for constant load
              If (MatInLayer==1) then
		      SigY = SigY + CalParams%InitialVerticalLoadK0 * CalParams%Multipliers%SolidARealised(1)
		      End if
              ! Calculate water pressure
		      If (MatInLayer>1) then
               SigWP = CalParams%BoreHoleData(MatInLayer-1,7)+&
			      SDepth * GammaLiquid * CalParams%Multipliers%GravityRealised
		      else
		       SigWP = SDepth * GammaLiquid * CalParams%Multipliers%GravityRealised
		      end if
              if (CalParams%ConsiderPTInLoadK0) then
                SigY = SigY - MatParams(MaterialIndex)%HypoPt
                if ((MatParams(MaterialIndex)%MaterialModel)==ESM_EXTERNAL_SOIL_MODEL) then
                SigY = SigY - MatParams(MaterialIndex)%ESM_Solid(2)
                end if
		      end if
		    end if            
        endif
        
        !Otherwise calculate using UNSAT definitions
    
          ! Determine depth of point below a defined (not necessarily perpendicular to gravity vector) soil surface
          if (Counters%SoilSurfaceNumberofSides > 0) then
              call DetermineCoordRelativeToDefinedSurface(SoilPoint, SoilSurfCoord, SoilSurfDepth,NewReferenceSoil)
              call DetermineCoordRelativeToDefinedSurface(SoilPoint, PhreatSurfCoord, PhreaticSurfDepth,NewReferencePhreatic)
          end if          
         
        
          ! Determine depth of point below a phreatic surface defined in an external file     
          if (GeoParams%ApplyInitialWaterSurfaceFromFile == .true.)then
               do I = 1, GeoParams%NumberWaterSurfaceMaterials
               if (MaterialIndex == GeoParams%WaterSurfaceMaterialID(I))then               
			   call DetermineCoordRelativeToSurfaceDefinedInExternalFile(MaterialIndex, SoilPoint, PhreaticSurfDepth,NewReferencePhreatic)
               end if
               end do
		  end if
          
          ! Determine depth of point below a defined phreatic surface
          
          ! point above phreatic surface (dry/partially sat conditions)
          if (SoilPoint(2) > NewReferencePhreatic(2)) then
          PhreaticSurfDepth = - PhreaticSurfDepth
          end if
          
          ! Retrieve possibly depth-dependent material parameters
          if(NFORMULATION==1) then
              if (CalParams%ApplySubmergedCalculation) then
                  GammaMixture = MatParams(MaterialIndex)%WeightSubmerged
                  GammaLiquid = 0.0
              else if (MatParams(MaterialIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') then
                  GammaMixture = MatParams(MaterialIndex)%WeightLiquid
                  GammaLiquid = 0.0 !fluid is the only one material and therefore thought as MIXTURE with only 1 material
              else
                  GammaLiquid = MatParams(MaterialIndex)%WeightLiquid
                  if ((CalParams%NumberOfPhases==1).and.(.not.IsUndrEffectiveStress)) then
                      GammaMixture = MatParams(MaterialIndex)%DryWeight
                  else
                      GammaMixture = MatParams(MaterialIndex)%WeightMixture
                  end if
              end if
          end if
          
          ! Compute vertical total stress for gravity loading
          ! Zero stresses for points above the specified material surface (positive depth is point above surface)
          if(MatParams(MaterialIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') then
              UpdatedFluidDensity  = GammaMixture / CalParams%GravityData%GAccel * &
                  exp((CalParams%Multipliers%GravityRealised * CalParams%GravityData%GAccel) /  &
                  MatParams(MaterialIndex)%BulkModulusLiquid * (- Depth))

              SigY = - MatParams(MaterialIndex)%BulkModulusLiquid * ( UpdatedFluidDensity -  &
                  GammaMixture / CalParams%GravityData%GAccel)
         end if
         
         if (Counters%SoilSurfaceNumberofSides > 0) then 
             ! case 1: soil surface overlaps phreatic surface
             if (SoilSurfDepth == PhreaticSurfDepth) then
                 SigY = - SoilSurfDepth * GammaMixture * CalParams%Multipliers%GravityRealised
                 ! case 2: phreatic surface is above soil surface ("ponding" or river)
             else if (SoilSurfDepth < PhreaticSurfDepth) then
                 SigY = - (SoilSurfDepth * GammaMixture + (PhreaticSurfDepth - SoilSurfDepth)*GammaLiquid) * CalParams%Multipliers%GravityRealised
                 ! case 3a:   phreatic surface below soil surface and point below phreatic
             else if ((SoilSurfDepth > PhreaticSurfDepth) .and. (PhreaticSurfDepth > 0)) then
                 SigY = - (PhreaticSurfDepth * GammaMixture + (SoilSurfDepth - PhreaticSurfDepth ) * MatParams(MaterialIndex)%DryWeight) &
                     * CalParams%Multipliers%GravityRealised
                 ! case 3b: phreatic surface below soil surface and point above phreatic
             else if (PhreaticSurfDepth < 0) then
                 SigY = - SoilSurfDepth* MatParams(MaterialIndex)%DryWeight * CalParams%Multipliers%GravityRealised
             end if
         else
             SigY = Depth * GammaMixture * CalParams%Multipliers%GravityRealised
         end if
              
          ! Compute vertical total stress for constant load
          SigY = SigY + CalParams%InitialVerticalLoadK0 * CalParams%Multipliers%SolidARealised(1)

          !Calculate water pressure
          if ((Counters%SoilSurfaceNumberofSides > 0) .or. (GeoParams%ApplyInitialWaterSurfaceFromFile == .true.)) then
              SigWP = - PhreaticSurfDepth * GammaLiquid * CalParams%Multipliers%GravityRealised
          else
              SigWP = Depth * GammaLiquid * CalParams%Multipliers%GravityRealised
          end if
         !add initial water pressure
          SigWP = SigWP + CalParams%InitialWaterPressure*CalParams%Multipliers%WaterARealised(1)

        if (SigWP> CalParams%K0MaxSurfaceSuction) then
          SigWP = CalParams%K0MaxSurfaceSuction
        end if
        
        if (MatParams(MaterialIndex)%MaterialType==DRY_SOIL) then
            SigWP = 0.0
        end if

        end subroutine ComputeVerticalStressForK0
        
        
        subroutine ComputeVerticalStressForK02Layers(MatPoint, UpdatedFluidDensity, SigY, SigWP, ParticleIndex, MaterialSetLayer)
        !**********************************************************************
        !
        !    Function: compute vertical stresses for K0 initialisation for double-point formulation (2 layers of material points)
        ! 
        !**********************************************************************        
        implicit none
        
          real(REAL_TYPE), dimension(2), intent(in) :: MatPoint ! Global position of material point / node
          integer(INTEGER_TYPE), dimension(CalParams%NumberSoilLayers, 1), intent(in) :: MaterialSetLayer
          real(REAL_TYPE), intent(out) :: UpdatedFluidDensity, SigY, SigWP
          integer(INTEGER_TYPE), intent(in) :: ParticleIndex ! Material Point Index 
          ! Local variables
          integer(INTEGER_TYPE) :: MaterialIndex, ISoilLayer
          real(REAL_TYPE) :: GammaMixture, GammaLiquid
          real(REAL_TYPE) :: EVol
          real(REAL_TYPE) :: Depth, WaterThick, BottomLayerDepth, TopLayerDepth, ThicknessLayer
               
          SigY = 0.0
          SigWP = 0.0
          Depth = 0.0
          WaterThick = 0.0
          BottomLayerDepth = 0.0
          TopLayerDepth = 0.0
          ThicknessLayer = 0.0

          ! Determine depth of point below
          call DetermineCoordRelativeToSurface2Layers(MatPoint, Depth, ParticleIndex)

        if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeLiquid) then
           MaterialIndex = MaterialIDArray(ParticleIndex)
           GammaMixture = MatParams(MaterialIndex)%WeightLiquid
           GammaLiquid = MatParams(MaterialIndex)%WeightLiquid

           ! Compute vertical total stress for gravity loading
           ! Zero stresses for points above the specified material surface (positive depth is point above surface)
           SigY = (Depth * GammaMixture - WaterThick * GammaLiquid) * CalParams%Multipliers%GravityRealised
          
           ! Compute vertical total stress for constant load
           !SigY = SigY + CalParams%InitialVerticalLoadK0 * CalParams%Multipliers%SolidARealised

           ! Calculate water pressure
           SigWP = (Depth - WaterThick)  * GammaLiquid
           SigWP = SigWP + CalParams%InitialWaterPressure*CalParams%Multipliers%WaterARealised(1)
           
        end if ! End Liquid MP
        
        if(MaterialPointTypeArray(ParticleIndex)==MaterialPointTypeSolid) then
            
          do ISoilLayer = 1, CalParams%NumberSoilLayers
              
            ThicknessLayer =  CalParams%ThicknessSoilLayer(ISoilLayer)
 
            ! Assign new bottom layer depth
            BottomLayerDepth = BottomLayerDepth - ThicknessLayer
            
            if(Depth<TopLayerDepth) then
            
              MaterialIndex = MaterialSetLayer(ISoilLayer,1) !Set material intex
              
              if(MatParams(MaterialIndex)%MaterialType=='2-phase'.or.MatParams(MaterialIndex)%MaterialPhases=='2-phase') then
                GammaMixture = MatParams(MaterialIndex)%DryWeight  + MatParams(MaterialIndex)%WeightLiquid *  &
                             MatParams(MaterialIndex)%InitialPorosity   
                GammaLiquid = MatParams(MaterialIndex)%WeightLiquid
              end if
              if(MatParams(MaterialIndex)%MaterialType=='1-phase-solid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-solid') then
                GammaMixture = MatParams(MaterialIndex)%DryWeight
                GammaLiquid = 0.0
              end if          
          if(MatParams(MaterialIndex)%MaterialType=='1-phase-solid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-solid') then
                GammaMixture = MatParams(MaterialIndex)%DryWeight
                GammaLiquid = 0.0
          end if
          
              ! Compute vertical total stress for gravity loading
              ! Zero stresses for points above the specified material surface (positive depth is point above surface)
              if(Depth<BottomLayerDepth) then
                 SigY = SigY - ThicknessLayer * GammaMixture
              else if (Depth>=BottomLayerDepth) then
                 SigY = SigY + (Depth - TopLayerDepth) * GammaMixture
              end if

            end if
              
            ! Assign new top layer depth
            TopLayerDepth = TopLayerDepth - ThicknessLayer             
          
          end do
          
          
          call DetermineCoordRelativeToSurface2Layers(CalParams%LiquidSurfacePoint, WaterThick, ParticleIndex)
          ! add surface load
          SigY = SigY + CalParams%InitialVerticalLoadK0 * CalParams%Multipliers%SolidARealised(1)
          ! Zero stresses for points above the specified material surface (positive depth is point above surface)
          SigY = (SigY - WaterThick * GammaLiquid) * CalParams%Multipliers%GravityRealised
              
          ! Calculate water pressure
             SigWP = (Depth - WaterThick)  * GammaLiquid
             SigWP = SigWP + CalParams%InitialWaterPressure*CalParams%Multipliers%WaterARealised(1)
             
          ! Update density in case of fluid
          if (MatParams(MaterialIndex)%MaterialType=='1-phase-liquid'.or.MatParams(MaterialIndex)%MaterialPhases=='1-phase-liquid') then
            EVol = SigY / MatParams(MaterialIndex)%BulkModulusLiquid
            if (CalParams%GravityData%GAccel>0.0) then
              UpdatedFluidDensity = GammaMixture / CalParams%GravityData%GAccel * (1.0 - EVol / (1.0 + EVol))
            else
               call GiveError("No gravity acceleration specified for fluid stress initialisation")
            end if
          end if

          ! for hypoplastic model reduce total vertical stress by pt
          SigY = SigY - MatParams(MaterialIndex)%HypoPt

        end if ! End Solid MP

        end subroutine ComputeVerticalStressForK02Layers

        
        subroutine GetMaterialData(IntGlo, ISet, &
                                   XNu, BFac, SPhi, &
                                   SPsi, GG, Cohec, Tens)
        !**********************************************************************
        !
        !    Function: Determines material data for the considered integration point
        !
        !              Note: In case of fully filled elements, the material data of the
        !                    first particle inside the element is used. (There must be
        !                    at least one particle because the element is active.)
        !                    Of course, this only works when assuming one material and thereby
        !                    forms a temporary solution. If more than one material is permitted
        !                    inside an element, the solution might be a weighted averaging,
        !                    determining material values from particle data at Gauss points.
        !
        !     IEl : ID of the element, in which the integration point is located
        !     IntGlo : Global ID of the considered integration point
        ! 
        ! O   ISet : ID of the material set of the integration point
        ! O   XNu : Poisson's ratio
        ! O   BFac : ...
        ! O   SPhi : Sinus Phi
        ! O   SPsi : Sinus Psi
        ! O   GG : Shear modulus
        ! O   Cohec : Cohesion
        ! O   Tens : Tension-cut-off value
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          integer(INTEGER_TYPE), intent(in) :: IntGlo
          integer(INTEGER_TYPE), intent(out) :: ISet
          real(REAL_TYPE), intent(out) :: XNu, BFac, &
                                           SPhi, SPsi, GG, Cohec, &
                                           Tens
                                
          ! Local variables
          real(REAL_TYPE) :: VNu_Un
          character(len=64) :: SoilModel ! name of the constitutive model

            ISet = MaterialIDArray(IntGlo) ! is the material number stored in $$MATERIAL_INDEX in the GOM-file
            SoilModel = MatParams(ISet)%MaterialModel ! name of constitutive model as specified in GOM-file

            if (SoilModel == ESM_HYPOPLASTICITY_SAND) then
              ! the rest of this subroutine is not needed in case of HP model
              RETURN
            endif

            XNu = MatParams(ISet)%PoissonRatio
            VNu_Un = MatParams(ISet)%UndrainedPoissonRatio
            BFac = (1.0 + VNu_Un) / (1.0 - 2.0 * VNu_Un) - (1.0 + XNu) / (1.0 - 2.0 * XNu)
            BFac = 2.0D0 / 3.0D0 * BFac
            SPsi = SIN(MatParams(ISet)%DilatancyAngle*(Pi/180.0))
            SPhi = SIN(MatParams(ISet)%FrictionAngle*(Pi/180.0))
            GG = Particles(IntGlo)%ShearModulus
            Cohec = Particles(IntGlo)%CohesionCosPhi

            Tens = MatParams(ISet)%TensileStrength
            if (Tens>Particles(IntGlo)%SFail) then
              Tens = Particles(IntGlo)%SFail
            end if
            
 
            


        end subroutine GetMaterialData


        subroutine AssignStressStrainToLocalArray(IntGlo, IDim, Sig0, SigC, DEpsPrevious)
        !**********************************************************************
        !
        !    Function: Assigns the stresses of the considered integration point (Gauss point
        !              or particle) to a local copy.
        !
        !     IEl : ID of the element, in which the integration point is located
        !     IntGlo : Global ID of the considered integration point
        !     ID : Number of stress / strain components (6 in 3D, 4 in 2D)
        ! 
        ! O   Sig0 : Local copy of the initial stresses of the integration point of the considered load step
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IntGlo, IDim
          real(REAL_TYPE), dimension(IDim), intent(out) :: Sig0, SigC
          real(REAL_TYPE), dimension(IDim), intent(out) :: DEpsPrevious
          ! Local variables
          integer(INTEGER_TYPE) :: I
       
          do I = 1, IDim
              Sig0(I) = SigmaEff0Array(IntGlo, I)
              SigC(I) = SigmaEffArray(IntGlo, I)
              DEpsPrevious(I) = GetEpsStepPreviousI(Particles(IntGlo), I)
          end do
          
        end subroutine AssignStressStrainToLocalArray 
          
        subroutine GetBMatrix(IEl, LocalID, Co, B, WtN)
        !**********************************************************************
        !
        !    Function: Returns the B matrix and the global weight of the considered
        !              integration point (Gauss point for fully filled elements, particle
        !              for partially filled elements).
        !
        !     IEl : ID of the element, in which the integration point is located
        !     LocalID : Local ID of the Gauss point or particle inside the element
        !     Co : Nodal coordinates
        ! 
        ! O   B : B matrix at the considered integration point
        ! O   WtN : Global weight at the considered integration point
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IEl, LocalID
          real(REAL_TYPE), dimension(Counters%NodTot, NVECTOR), intent(in) :: Co
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES), intent(out) :: B
          real(REAL_TYPE), intent(out) :: WtN
          ! Local variables
          integer(INTEGER_TYPE) :: IParticleIndex
          real(REAL_TYPE) :: Det
                                            
          if (IsParticleIntegration(IEl) ) then ! Particle
            IParticleIndex = GetParticleIndex(LocalID, IEl)   
            call BMatrix(Particles(IParticleIndex)%LocPos, & ! Get the B matrix for the particle
                         ELEMENTNODES, Counters%NEl,  Counters%NodTot, NVECTOR, IEl, ElementConnectivities, Co, B, Det)
            WtN = Particles(IParticleIndex)%IntegrationWeight 
          else ! Gauss point
            call FormB3(LocalID, IEl, ElementConnectivities, Co, B, Det, WtN)
          end if
            
        end subroutine GetBMatrix
           
        integer(INTEGER_TYPE) function GetIPL(IntGlo)
        !**********************************************************************
        !
        !    Function: Returns the plasticity state of the integration point
        !              identified by IntGlo 
        !
        !     IntGlo : ID of the integration point
        !     IEl : ID of the element, in which the integration point is located
        ! 
        ! O   GetIPL : Plasticity state of the integration point
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IntGlo
          
            GetIPL = Particles(IntGlo)%IPL
        
        end function GetIPL
        
        subroutine SetIPL(IntGlo, IEl, PlasticityState)
        !**********************************************************************
        !
        !    Function: Updated the plasticity state of the integration point
        !              identified by IntGlo (Gauss point or particle).
        !
        !     IntGlo : ID of the integration point (Gauss point or particle)
        !     IEl : ID of the element, in which the integration point is located
        !     Plasticity state : New value for the plasticity state of IntGlo
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IntGlo, IEl
          integer(INTEGER_TYPE), intent(in) :: PlasticityState
          !local variable
          integer(INTEGER_TYPE):: NElemPart, IPart, ParticleIndex
          
          if (IsParticleIntegration(IEl) ) then ! Partially filled element
            Particles(IntGlo)%IPL = PlasticityState
          else ! full filled element
            NElemPart = NPartEle(IEl)  ! Number of particles in element
            do IPart = 1, NElemPart
              ParticleIndex = GetParticleIndex(IPart, IEl)
              Particles(ParticleIndex)%IPL = PlasticityState
            end do
          end if
        
        end subroutine SetIPL
        

          


        subroutine UpdateParticleSoftParam(c,phi,psi,IntGlo)
        !**********************************************************************
        !
        !    Function: update the softening parameters of the particle
        !
        !     IntGlo: particle ID
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

        !In variables
        integer(INTEGER_TYPE), intent(in):: IntGlo
        real(REAL_TYPE), intent(in) :: c,phi,psi

            Particles(IntGlo)%CohesionStSoft = c
            Particles(IntGlo)%PhiStSoft = phi
            Particles(IntGlo)%PsiStSoft = psi

        end subroutine UpdateParticleSoftParam

        subroutine UpdateParticleMCCParam(PreconsolidationPressure, IntGlo)
        !**********************************************************************
        !
        !    Function: update the preconsolidation pressure of the particle
        !                 used for Modified Cam Clay model
        !
        !     IntGlo: particle ID
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        
        implicit none

        !In variables
        integer(INTEGER_TYPE), intent(in):: IntGlo
        real(REAL_TYPE) :: PreconsolidationPressure
          
            Particles(IntGlo)%pp = PreconsolidationPressure
                     
        end subroutine UpdateParticleMCCParam
        
        subroutine AssignStressStrainToGlobalArray(IntGlo, IDim, DSigWP, DSigGP, DSig, SigPrin, DEps)
        !**********************************************************************
        !
        !  Function: Updates the global stress and strain data of the integration
        !            point identified by IntGlo (Gauss point or material point) from
        !            the local stress and strain arrays.
        !
        !  IEl : ID of the element, in which the integration point is located
        !  IntGlo : ID of the Gauss point or material point inside the element
        !  ID : Number of stress / strain components (6 in 3D, 4 in 2D)
        !  DSigWP : Change of water pressures
        !  DSigGP : Change of gas pressures
        !  DSig : Change of stresses
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IntGlo, IDim
          real(REAL_TYPE), intent(in) :: DSigWP
          real(REAL_TYPE), intent(in) :: DSigGP
          real(REAL_TYPE), dimension(IDim), intent(in) :: DEps
          real(REAL_TYPE), dimension(IDim), intent(in) :: DSig
          real(REAL_TYPE), dimension(IDim), intent(in) :: SigPrin
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: Value
     
          do I = 1, IDim
              Value = SigmaEff0Array(IntGlo, I) + DSig(I)
              SigmaEffArray(IntGlo, I) = Value
              if (CalParams%ApplyImplicitQuasiStatic) then
                call SetEpsStepPreviousI(Particles(IntGlo), I, DEps(I))
              end if
          end do
          
          Particles(IntGlo)%WaterPressure = Particles(IntGlo)%WaterPressure0 + DSigWP
          Particles(IntGlo)%GasPressure = Particles(IntGlo)%GasPressure0 + DSigGP
          
          ! Ensure that the gas pressure is = 0 when the material point is saturated
          if (Particles(IntGlo)%DegreeSaturation>=1) then
              Particles(IntGlo)%GasPressure=0.0
          end if 

          ! Important assumption: gas can not be under tension. This is realistic in most of the cases.
          if (Particles(IntGlo)%GasPressure>0) then
              Particles(IntGlo)%GasPressure=0.0
          end if 
          
          call SetSigmaPrin(Particles(IntGlo),SigPrin)
          
        end subroutine AssignStressStrainToGlobalArray
          
        subroutine AssignStressStrainToGlobalArrayESM(IntGlo, IDim, DSig, SigPrin, DEps)
        !**********************************************************************
        !
        !  Function: Updates the global stress and strain data of the integration
        !            point identified by IntGlo (Gauss point or material point) from
        !            the local stress and strain arrays.
        !
        !  IEl : ID of the element, in which the integration point is located
        !  IntGlo : ID of the Gauss point or material point inside the element
        !  ID : Number of stress / strain components (6 in 3D, 4 in 2D)
        !  DSig : Change of stresses
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IntGlo, IDim
          real(REAL_TYPE), dimension(IDim), intent(in) :: DEps
          real(REAL_TYPE), dimension(IDim), intent(in) :: DSig
          real(REAL_TYPE), dimension(IDim), intent(in) :: SigPrin
          ! Local variables
          integer(INTEGER_TYPE) :: I
          real(REAL_TYPE) :: Value
     
          do I = 1, IDim
              Value = SigmaEff0Array(IntGlo, I) + DSig(I)
              SigmaEffArray(IntGlo, I) = Value
              if (CalParams%ApplyImplicitQuasiStatic) then
                call SetEpsStepPreviousI(Particles(IntGlo), I, DEps(I))
              end if
          end do
          
          call SetSigmaPrin(Particles(IntGlo),SigPrin)
          
        end subroutine AssignStressStrainToGlobalArrayESM
        
        subroutine AssignWatandgasPressureToGlobalArray(IntGlo, DSigWP, DSigGP)
        !**********************************************************************
        !
        !  Function: Updates the global water and gas pressure of the integration
        !            point identified by IntGlo (Gauss point or material point) from
        !            the local stress and strain arrays.
        !
        !  IEl : ID of the element, in which the integration point is located
        !  IntGlo : ID of the Gauss point or material point inside the element
        !  ID : Number of stress / strain components (6 in 3D, 4 in 2D)
        !  DSigWP : Change of water pressures
        !  DSigGP : Change of gas pressures
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: IntGlo
          real(REAL_TYPE), intent(in) :: DSigWP
          real(REAL_TYPE), intent(in) :: DSigGP
          ! Local variables
          
          Particles(IntGlo)%WaterPressure = Particles(IntGlo)%WaterPressure0 + DSigWP
          Particles(IntGlo)%GasPressure = Particles(IntGlo)%GasPressure0 + DSigGP
          
          ! Ensure that the gas pressure is = 0 when the material point is saturated
          if (Particles(IntGlo)%DegreeSaturation>=1) then
              Particles(IntGlo)%GasPressure=0.0
          end if 

          ! Important assumption: gas can not be under tension. This is realistic in most of the cases.
          if (Particles(IntGlo)%GasPressure>0) then
              Particles(IntGlo)%GasPressure=0.0
          end if 
          

          
        end subroutine AssignWatandgasPressureToGlobalArray
        
        subroutine SolveBalanceEquations(IntGlo,DEpsVol,dPw,dPg,dT)  
        !**********************************************************************
        !
        !    Function: Solves the sytem of equations (M*x = b) formed by the balances of mass and energy 
        !              
        !               SYSTEM OF EQUATIONS:
        !               A1*dPw + A2*dPg + A3*dT = A  ---> Air mass balance eq.
        !               W1*dPw + W2*dPg + W3*dT = W  ---> Water mass balance eq.
        !               H1*dPw + H2*dPg + H3*dT = H  ---> Energy balance eq.
        !                
        !               Terms of the air balance equation ----> { A1, A2, A3 }
        !               Terms of the water balance equation --> { W1, W2, W3 } ---> M
        !               Terms of the energy balance equation -> { H1, H2, H3 }
        ! 
        !               COLUMN VECTOR ------------------------> b = { A, W, H }
        !
        !               SOLUTION OF THE SYSTEM ---------------> x = { dPw, dPg, dT }
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          integer(INTEGER_TYPE), parameter :: IDim = 3 ! vectorsize of A, W, H
          integer(INTEGER_TYPE), intent(in) :: IntGlo
          real(REAL_TYPE), intent(inout) :: DEpsVol
          real(REAL_TYPE), intent(inout) :: dPw
          real(REAL_TYPE), intent(inout) :: dPg
          real(REAL_TYPE), intent(inout) :: dT
        ! Local variables
          real(REAL_TYPE), dimension(:,:), allocatable :: M !Matrix of the system
          real(REAL_TYPE), dimension(:), allocatable :: b   !Coumn vector
          real(REAL_TYPE), dimension(:), allocatable :: x   !Solution vector
          integer(INTEGER_TYPE) :: NEq ! vectorsize of Air, Water, Heat
          integer(INTEGER_TYPE) :: IError
          
            if (Particles(IntGlo)%DegreeSaturation>=1.0) then
                NEq = 1
            else
                NEq = 2
            end if

            allocate(M(NEq,NEq), stat = IError)
            allocate(b(NEq), stat = IError)
            allocate(x(NEq), stat = IError)
            
            !there are 3 unknowns (dPw, dPg, dT)
            x = 0.0d0
            M = 0.0d0
            b = 0.0d0
            
            call CalculateTermsMassBalanceEquations(IntGlo,DEpsVol,M,b,NEq)   !This subroutine is in the MPMStresses
            call SolveSystem_Gauss(NEq,M,b,x)                         !This subroutine is in the MPMStresses

            if (NEq==1) then            
                dPw = x(1)
                dPg = 0.0d0
                dT  = 0.0d0
            else if (NEq==2) then
                dPw = x(1)
                dPg = x(2)
                dT  = 0.0d0
            else if (NEq==3) then
                dPw = x(1)
                dPg = x(2)
                dT  = x(3)
            end if

            deallocate(M, stat = IError)
            deallocate(b, stat = IError)
            deallocate(x, stat = IError)
            
        end subroutine SolveBalanceEquations

        subroutine CalculateTermsMassBalanceEquations(IntGlo,DEpsVol,M,b,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the terms of the mass balance equations for a particle
        !              
        !               MASS BALANCE EQUATIONS
        !               W1*dPw + W2*dPg + W3*dT = W  ---> Water mass balance eq.
        !               A1*dPw + A2*dPg + A3*dT = A  ---> Air mass balance eq.
        !                
        !               Terms of the water balance equation --> W1, W2, W3, W
        !               Terms of the air balance equation ----> A1, A2, A3, A 
        !
        ! 
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        integer(INTEGER_TYPE), intent(in) :: IntGlo, NEq
        real(REAL_TYPE) :: DEpsVol
        real(REAL_TYPE), dimension(NEq,NEq) :: M
        real(REAL_TYPE), dimension(NEq) :: b
    
        ! Local variables
        real(REAL_TYPE) ::  g, N, WD, GD, Sr, Sg, wAW, wVG
        real(REAL_TYPE) ::  WaterAdvectiveFlux, AirAdvectiveFlux
        real(REAL_TYPE), dimension(NEq) :: dSr, dWD, dGD, dwAW, dwVG
        real(REAL_TYPE) ::  q1, q2, q3, t1, t2, t3
        real(REAL_TYPE) ::  xVG, xAW, xW, xG

        real(REAL_TYPE):: DivergenceNonAdvectiveFluxAirInWater = 0.0
        real(REAL_TYPE):: DivergenceNonAdvectiveVapourInGas = 0.0
        
          g = CalParams%GravityData%GAccel   !gravity
            
          N = Particles(IntGlo)%Porosity            ! Porosity of the particle
          WD = Particles(IntGlo)%WaterWeight/g      ! Water density of the particle
          GD = Particles(IntGlo)%GasWeight/g        ! Gas density of the particle
          Sr = Particles(IntGlo)%DegreeSaturation   ! Degree of Saturation (liquid) of the particle
          Sg = 1.0d0 - Sr                           ! Degree of Saturation (gas) of the particle
          wAW = Particles(IntGlo)%AirInWaterMassFraction    ! Mass Fraction of Air in Water of the particle
          wVG = Particles(IntGlo)%VapourInGasMassFraction   ! Mass Fraction of Vapour in Gas of the particle
          WaterAdvectiveFlux = Particles(IntGlo)%WaterAdvectiveFlux
          AirAdvectiveFlux = Particles(IntGlo)%AirAdvectiveFlux
          
          dSr = 0.0
          dWD = 0.0
          dGD = 0.0
          dwAW =0.0
          dwVG =0.0

        ! In order to calculate all terms, to calculate different derivatives respect to Pw, Pg, T are needed
            ! Derivatives of Degree of Saturation                  ---> dSr 
                call CalculateDerivDegreeSaturation(IntGlo,dSr,NEq)
            ! Derivatives of Liquid density                        ---> dWD
                call CalculateDerivLiquidDensity(IntGlo,dWD,NEq)
            ! Derivatives of Gas density                           ---> dGD
                call CalculateDerivGasDensity(IntGlo,dGD,NEq)
            ! Derivatives of Mass fraction of the Air in the Water         ---> dwAW   
                    call CalculateDerivAirInWaterMassFraction(dwAW, NEq)
            ! Derivatives of Mass fraction of the Vapour (water) in the Gas    ---> dwVG    
                    call CalculateDerivVapourInGasMassFraction(IntGlo,dwVG, NEq)

        ! The Divergence of the Non Advective fluxes is needed (AirInWater and VapourInGas)    
                call GetDivergenceNonAdvectiveFluxes(IntGlo, DivergenceNonAdvectiveFluxAirInWater, DivergenceNonAdvectiveVapourInGas ) 
                
            t1 = GD*Sg
            t2 = wVG*Sg
            t3 = (WD-wVG*GD)
          
            q1 = WD*Sr
            q2 = wAW*Sr
            q3 = (wAW*WD-GD)
          
            xVG = wVG*GD*Sg
            xAW = wAW*WD*Sr
            xW = WD*Sr
            xG = GD*Sg
          
            !Determine the terms
            if (NEq==1) then
                M(1,1) = N * (t1*dwVG(1) + t2*dGD(1) + Sr*dWD(1) + t3*dSr(1))
                b(1) = -(xVG + xW)*DEpsVol - WaterAdvectiveFlux - DivergenceNonAdvectiveVapourInGas 
                
            else if (NEq==2) then
                M(1,1) = N * (t1*dwVG(1) + t2*dGD(1) + Sr*dWD(1) + t3*dSr(1))
                M(1,2) = N * (t1*dwVG(2) + t2*dGD(2) + Sr*dWD(2) + t3*dSr(2))
                b(1) = -(xVG + xW)*DEpsVol - WaterAdvectiveFlux - DivergenceNonAdvectiveVapourInGas 

                M(2,1) = N * (q1*dwAW(1) + q2*dWD(1) + Sg*dGD(1) + q3*dSr(1))
                M(2,2) = N * (q1*dwAW(2) + q2*dWD(2) + Sg*dGD(2) + q3*dSr(2))
                b(2) = -(xAW + xG)*DEpsVol - AirAdvectiveFlux - DivergenceNonAdvectiveFluxAirInWater
            
            else
                M(1,1) = N * (t1*dwVG(1) + t2*dGD(1) + Sr*dWD(1) + t3*dSr(1))
                M(1,2) = N * (t1*dwVG(2) + t2*dGD(2) + Sr*dWD(2) + t3*dSr(2))
                M(1,3) = N * (t1*dwVG(3) + t2*dGD(3) + Sr*dWD(3) + t3*dSr(3))
                b(1) = -(xVG + xW)*DEpsVol - WaterAdvectiveFlux - DivergenceNonAdvectiveVapourInGas 

                M(2,1) = N * (q1*dwAW(1) + q2*dWD(1) + Sg*dGD(1) + q3*dSr(1))
                M(2,2) = N * (q1*dwAW(2) + q2*dWD(2) + Sg*dGD(2) + q3*dSr(2))
                M(2,3) = N * (q1*dwAW(3) + q2*dWD(3) + Sg*dGD(3) + q3*dSr(3))
                b(2) = -(xAW + xG)*DEpsVol - AirAdvectiveFlux - DivergenceNonAdvectiveFluxAirInWater
            end if
            
        end subroutine CalculateTermsMassBalanceEquations
        
        subroutine CalculateDerivDegreeSaturation(IntGlo,dSr,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the Degree of Saturation (Sr)
        !
        !**********************************************************************
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo, NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dSr
        integer(INTEGER_TYPE) :: ISet
        
        ISet = MaterialIDArray(IntGlo)

        if (MatParams(ISet)%RetentionCurve == SWRC_VANGENUCHTEN) then
            call CalculateDerivDegreeSaturationVanGenuchten(IntGlo,ISet,dSr,NEq)
        else if (MatParams(ISet)%RetentionCurve == SWRC_LINEAR) then
            call CalculateDerivDegreeSaturationLinear(IntGlo,ISet,dSr,NEq)
        end if
        
        end subroutine CalculateDerivDegreeSaturation

        
        subroutine CalculateDerivDegreeSaturationVanGenuchten(IntGlo,ISet,dSr,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the Degree of Saturation (Sr)
        !               with respect de Pw, Pg, and T
        !               The expression of the Degree of Saturation is the RETENTION CURVE, based
        !               on the VAN GENUCHTEN MODEL (1980)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo, ISet, NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dSr
        ! Local variables
        real(REAL_TYPE) ::  Pw, Pg, T, N
        real(REAL_TYPE) ::  Suc, Smin, Smax
        real(REAL_TYPE) ::  L, P, P0, dPdT
        real(REAL_TYPE) ::  n00, n1, n2, n3, n4, n5, n6, n7, n8, n9


            Pw = Particles(IntGlo)%WaterPressure
            Pg = Particles(IntGlo)%GasPressure
            T = Particles(IntGlo)%Temperature 
            N = Particles(IntGlo)%Porosity
            
            n00 = 1.0d0
            n1 = 0.625d0            
            n2 = 0.2358d0
            n3 = 1.256d0
            n4 = 0.0019406d0
            n5 = 0.05d0
            n6 = 360.0d0
            n7 = 374.15d0
            n8 = 647.3d0
            n9 = 0.2961648d0
           
            Suc = Pw-Pg     !Suction (compressive pressure <0)
            Smin = MatParams(ISet)%Smin_SWRC   !If Smin=0.0 and Smax=1.0, the Effective Degree of Saturation (Se) is equivalent to the "liquid phase" Degree of Saturation (Sr)   
            Smax = MatParams(ISet)%Smax_SWRC   !If Smin=0.0 and Smax=1.0, the Effective Degree of Saturation (Se) is equivalent to the "liquid phase" Degree of Saturation (Sr)   
            P0 = MatParams(ISet)%P0_SWRC       !P0=15000 kPa (mudstone model) / P0=18000 kPa (bentonite model) / P0=15 kPa (sand) 
            L = MatParams(ISet)%Lambda_SWRC         !L=0.36 (mudstone model) / P0=0.38 (bentonite model)

            P = P0
            dPdT = 0.0d0
            
            if (Suc<=0.0) then !Saturated
                dSr = 0.0d0 
                return
            end if
            
           
            
            select case(NEq)
            case(1)
                !Derivative with respect to water pressure
                dSr(1) = (Smax-Smin)*(-L/(1.-L))*((1.+(Suc/P)**(1./(1.-L)))**(-L-1.))* ((Suc/P)**(L/(1.-L)))*(1./P)
            case(2)
                !Derivative with respect to water pressure
                dSr(1) = (Smax-Smin)*(-L/(1.-L))*((1.+(Suc/P)**(1./(1.-L)))**(-L-1.))* ((Suc/P)**(L/(1.-L)))*(1./P)
                !Derivative with respect to gas pressure
                dSr(2) = -dSr(1)
            end select
            
        end subroutine CalculateDerivDegreeSaturationVanGenuchten

        
        subroutine CalculateDerivDegreeSaturationLinear(IntGlo,ISet,dSr,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the Degree of Saturation (Sr)
        !               with respect Pw, Pg, and T
        !               The expression of the Degree of Saturation is a LINEAR Model
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo, ISet, NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dSr
        ! Local variables
        real(REAL_TYPE) :: av, Pw, Pg, Suc

            Pw = Particles(IntGlo)%WaterPressure
            Pg = Particles(IntGlo)%GasPressure
            Suc = Pw-Pg     !Suction (compressive pressure <0)
            
            if (Suc<=0.0) then !Saturated
                dSr = 0.0d0 
                return
            end if
            
            av = MatParams(ISet)%av_SWRC   ! m2/KN
            
            dSr = 0.0
            select case(NEq)
            case(1)
              dSr(1)=-av
            case(2)
              dSr(1)=-av
              dSr(2)=av
            case(3)
              dSr(1)=-av
              dSr(2)=av
            end select

            
        end subroutine CalculateDerivDegreeSaturationLinear
                    

        subroutine CalculateDerivLiquidDensity(IntGlo,dWD,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the liquid density (WD)
        !               with respect de Pw, Pg, and T
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo, NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dWD
        ! Local variables
        real(REAL_TYPE) ::  Pw, T, Kw
        real(REAL_TYPE) ::  Alpha, Beta, Pw0, WD0, T0
        real(REAL_TYPE) ::  c1

            Pw = Particles(IntGlo)%WaterPressure
            T = Particles(IntGlo)%Temperature   !Temperature (ºC)
            Kw = Particles(IntGlo)%BulkWater
            
            Alpha = -0.00034d0  !Default=-0.00034, volumetric thermal expansion coefficient for water (1/ºC)
            Beta = -1/Kw         !Default=0.00045, water compressibility (1/MPa) ---> Inverse of the bulk modulus
            Pw0 =  0.0d0        !Default=0.1, Reference water pressure (MPa) 
            WD0 = 1.0026d0      !Default=1002.6, Reference water density (kg/m3)
            T0 = 20.0d0         !Default=20, Reference temperature (ºC)
        	c1=exp(Beta*(Pw-Pw0) + Alpha*(T-T0))

        select case(NEq)
            case(1)
                !Derivative with respect to water pressure
                dWD(1)= WD0*Beta*c1
            case(2)
                !Derivative with respect to water pressure
                dWD(1)= WD0*Beta*c1
                !Derivative with respect to gas pressure
                dWD(2)= 0.0d0
            case(3)
                !Derivative with respect to water pressure
                dWD(1)= WD0*Beta*c1
                !Derivative with respect to gas pressure
                dWD(2)= 0.0d0
                !Derivative with respect to temperature
                dWD(3)= WD0*Alpha*c1
            end select

        end subroutine CalculateDerivLiquidDensity


        subroutine CalculateDerivGasDensity(IntGlo,dGD,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the gas density (GD)
        !               with respect de Pw, Pg, and T
        !               The expression of the Gas Density is based
        !               on the IDEAL GASES LAW (PV=nRT)
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo, NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dGD
        ! Local variables
        real(REAL_TYPE) ::  Pw, Pg, T, WD
        real(REAL_TYPE) ::  g, Suc, Ma, Mw, R, Pv0, F, Pv
        real(REAL_TYPE) ::  n1, n2

            g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
            Pw = Particles(IntGlo)%WaterPressure
            Pg = Particles(IntGlo)%GasPressure 
            T = Particles(IntGlo)%Temperature           !Temperature (ºC)
            WD = Particles(IntGlo)%WaterWeight/g        !Water density of the particle
      
            Suc = Pw-Pg  !Suction (compressive pressure < 0)
             
            if (Suc<=0.0) then  !Saturated
                dGD = 0.0d0
                return
            end if

            T = T+273.15d0              !Temperature (ºK)
            Ma = 0.028d0                !Default=0.02895, Molecular gas of dry air (kg/mol)
            Mw = 0.018d0                !Default=0.01801528, Molecular mass of water (=vapour) (kg/mol)
            R = 0.008314d0              !Default=8.3144621, Gas constant value (J/mol*K)


            !Calculation of the vapour Pressure (Pg = Pv + Pa) using the psychrometric law 
            n1 = 136.075d0
            n2 = -5239.7d0
            Pv0 = n1*exp(n2/T)
            F = exp(-(Mw*Suc)/(R*T*WD))     !Psychrometric Law
            Pv = Pv0*F                      !Vapour Pressure

        select case(NEq)
            case(1)
                !Derivative with respect to water pressure
                dGD(1)= 0.0d0
            case(2)
                !Derivative with respect to water pressure
                dGD(1)= 0.0d0
                !Derivative with respect to gas pressure
                dGD(2)= -Ma /(R*T)  
            case(3)
                !Derivative with respect to water pressure
                dGD(1)= 0.0d0
                !Derivative with respect to gas pressure
                dGD(2)= -Ma /(R*T)  
                !Derivative with respect to temperature
                dGD(3)= (Pv*(Mw-Ma)+Pg*Ma)/R
            end select

        end subroutine CalculateDerivGasDensity


        subroutine CalculateDerivAirInWaterMassFraction(dwAW,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the Mass Fraction of Air in Water (wAW)
        !               with respect de Pw, Pg, and T
        !               The expression of the Air in Water mass fraction is the HENRY's LAW 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dwAW
        ! Local variables
        real(REAL_TYPE) ::  Ma, Mw, H
!        real(REAL_TYPE) ::  T, wAW, Pa, H0, T0, C
        
            Ma = 0.028d0               !Default=0.02895, Molecular gas of dry air (kg/mol)
            Mw = 0.018d0               !Default=0.01801528 Molecular mass of water (kg/mol)
            
            !Taking Henry's coefficient as a constant. Default=10000000 (KPa)
            H=10000000
            
        select case(NEq)
            case(1)
                !Derivative with respect to water pressure
                dwAW(1)= 0.0d0
            case(2)
                !Derivative with respect to water pressure
                dwAW(1)= 0.0d0
                !Derivative with respect to gas pressure
                dwAW(2)= Ma/(Mw*H)
            case(3)
                !Derivative with respect to water pressure
                dwAW(1)= 0.0d0
                !Derivative with respect to gas pressure
                dwAW(2)= Ma/(Mw*H)
                !Derivative with respect to temperature
                dwAW(3)= 0.0d0
            end select


        end subroutine CalculateDerivAirInWaterMassFraction


        subroutine CalculateDerivVapourInGasMassFraction(IntGlo,dwVG,NEq)
        !**********************************************************************
        !
        !    Function: Calculates the derivatives of the Mass Fraction of Vapour in the Gass (wVG)
        !               with respect de Pw, Pg, and T
        !               The expression of the Vapour in Gas mass fraction is the PSYCHROMETRIC LAW 
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
        integer(INTEGER_TYPE), intent(in) :: IntGlo, NEq
        real(REAL_TYPE), dimension(NEq), intent(inout) :: dwVG
        ! Local variables
        real(REAL_TYPE) ::  Pw, Pg, T, wVG, WD
        real(REAL_TYPE) ::  g, Mw, R, Suc

            g = CalParams%GravityData%GAccel  !Gravity (m/s2)
            
            Pw = Particles(IntGlo)%WaterPressure
            Pg = Particles(IntGlo)%GasPressure
            T = Particles(IntGlo)%Temperature           !Temperature (ºC)
            wVG = Particles(IntGlo)%VapourInGasMassFraction
            WD = Particles(IntGlo)%WaterWeight/g        !Water density of the particle
      
            Suc = Pw-Pg     !Suction

            if (Suc<=0.0) then  !Saturated
                dwVG = 0.0d0
                return
            end if

            T = T+273.15d0              !Temperature (ºK)
            Mw = 0.018d0                !Default=0.01801528, Molecular mass of water (kg/mol)
            R = 0.008314d0              !Default= 8.3144621, Gas constant value (J/mol*K)
                                                
        select case(NEq)
            case(1)                                              
                !Derivative respect the water pressure
                dwVG(1)= (Mw/(R*T*WD))*exp(-(Mw*Suc)/(R*T*WD))   
            case(2)                                              
                !Derivative respect the water pressure
                dwVG(1)= (Mw/(R*T*WD))*exp(-(Mw*Suc)/(R*T*WD))
                !Derivative respect the gas pressure
                dwVG(2)= -dwVG(1) 
            case(3)                                              
                !Derivative respect the water pressure
                dwVG(1)= (Mw/(R*T*WD))*exp(-(Mw*Suc)/(R*T*WD))
                !Derivative respect the gas pressure
                dwVG(2)= -dwVG(1)
                !Derivative respect the temperature
                dwVG(3)= Mw*Suc/(R*T*T*WD)
        end select
            

        end subroutine CalculateDerivVapourInGasMassFraction
        
 
       subroutine GetDivergenceNonAdvectiveFluxes(IP,DivergenceNAAirInWater,DivergenceNAVapourInGas)
        !**********************************************************************
        !
        !    Function:  Calculate Divergence of the Non Advective Fluxes in the particles
        !                   (Fluxes are in the nodes!!)
        !
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none

          ! Local variables
          integer(INTEGER_TYPE) :: IP
          real(REAL_TYPE), dimension(NVECTOR, ELEMENTNODES) :: B
          real(REAL_TYPE) :: WtN, Det
          integer(INTEGER_TYPE) :: ElementID, iEntity, INode
          integer(INTEGER_TYPE) :: nn, I, J
          integer(INTEGER_TYPE), dimension(NVECTOR) :: IDof 
          real(REAL_TYPE), intent(inout)  :: DivergenceNAAirInWater 
          real(REAL_TYPE), intent(inout)  :: DivergenceNAVapourInGas 
            
            ElementID = ElementIDArray(IP)
          
            call FormB3(1, ElementID, ElementConnectivities, NodalCoordinatesUpd, B, Det, WTN) ! get the B-matrix once per element

            !get particle entity ID
            if (CalParams%ApplyContactAlgorithm) then
              iEntity = EntityIDArray(IP) 
            else
              iEntity = 1
            end if 
              
            DivergenceNAAirInWater = 0.0d0
            DivergenceNAVapourInGas = 0.0d0
                               
            do INode = 1, ELEMENTNODES  ! Loop over all nodes in element

                nn=ElementConnectivities(iNode,ElementID) ! get global node number
                do I = 1, NVECTOR
                  IDof(I) = ReducedDof(nn) + I  
                end do  

                do J = 1, NVECTOR
                  DivergenceNAAirInWater = DivergenceNAAirInWater + B(J,INode) * NonAdvectiveFluxAirInWater(IDof(J),iEntity) 
                  DivergenceNAVapourInGas = DivergenceNAVapourInGas + B(J,INode) * NonAdvectiveFluxVapourInGas(IDof(J),iEntity) 
                end do

            end do ! Loop over all nodes in element
     
      end subroutine GetDivergenceNonAdvectiveFluxes

      

      subroutine SolveSystem_Gauss(n,M,b,x)
        !**********************************************************************
        !
        !    Function: Calculates the solution of a linear system of equations using Gauss Method
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
    
        implicit none
          integer(INTEGER_TYPE), intent(in) :: n 
          real(REAL_TYPE), dimension(n,n), intent(inout) :: M
          real(REAL_TYPE), dimension(n), intent(inout) :: b
          real(REAL_TYPE), dimension(n), intent(inout) :: x

            !Gauss method
            call elimgauss(n,M,b)
      
            !Solution
            call solve_tu(n,M,b,x)
           
      end subroutine SolveSystem_Gauss
  

      subroutine elimgauss (n,M,b)
        !**********************************************************************
        !
        !    Function: Calculates The tiangular Gauss Matrix
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
      
        implicit none
          integer(INTEGER_TYPE), intent(in) :: n 
          real(REAL_TYPE), dimension(n,n), intent(inout) :: M
          real(REAL_TYPE), dimension(n), intent(inout) :: b
          ! Local variables
          integer(INTEGER_TYPE) :: k, i, j
          real(REAL_TYPE) :: factorm
          
            do  k=1,n-1
                do  i=k+1,n
                    factorm =  M(i,k)/M(k,k)
                    b(i) = b(i)-factorm*b(k)
                    M(k,i) = 0.0d0
                    do j=k+1,n
                        M(i,j) = M(i,j)-factorm* M(k,j)
                    end do
                end do  
            end do
      
        end subroutine elimgauss


        subroutine solve_tu(n, tumat, b, x)
        !**********************************************************************
        !
        !    Function: Solution The Triangular System of equations
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************

        implicit none
          integer(INTEGER_TYPE), intent(in) :: n 
          real(REAL_TYPE), dimension(n,n), intent(in) :: tumat
          real(REAL_TYPE), dimension(n), intent(in) :: b
          real(REAL_TYPE), dimension(n), intent(inout) :: x
          ! Local variables
          integer(INTEGER_TYPE) :: i, j
        
            x=0.0d0
            x(1) = b(1) / tumat(1,1)
            
            if (tumat(1,1)==0.0) then
                x(1) = 0.0
            end if
            
            do i=2,n
                x(i) = b(i)
                do j=1,i-1
                    x(i) = x(i) - tumat(i,j)*x(j)
                end do
                x(i) = x(i) / tumat(i,i)
                if (tumat(i,i)==0.0d0) then
                    x(i) = 0.0d0
                end if
            end do

        end subroutine solve_tu  

        subroutine InitialiseStVarMCC()
        !**********************************************************************
        !
        !    Function: Initialize the preconsolidation pressure for Modified Cam Clay model 
        !              
        !
        !     xM : slope of the CSL in the p-q plane
        !     Sig0 : initial stress coming from K0 procedure
        !     OCR : overconsolidation ratio
        !     xN1, xN2, xN3: direction of principal stresses
        !     S1, S2, S3: principal stresses
        !     P: mean effective stress
        !     Q: deviatoric stress
        !     pp: preconsolidation pressure
        !
        ! Implemented in the frame of the MPM project.
        !
        !**********************************************************************
        implicit none
        
        real(REAL_TYPE) :: xM, S1, S2, S3, P, Q, OCR, K0nc
        real(REAL_TYPE) :: xK, xG, xnu, xkappa, xlambda, xe0, nu_u, r
        real(REAL_TYPE), dimension(NVECTOR) :: xN1, xN2, xN3
        real(REAL_TYPE), dimension(NTENSOR) :: Sig0, SigYield
        integer(INTEGER_TYPE) :: I, iOpt, ISet 
        character(len=64) :: SoilModel ! name of the constitutive model
        logical :: IsUndrEffectiveStress
         
              if (.not.(Calparams%IStep==1)) Return
              
              do I=1,Counters%NParticles
                ISet = MaterialIDArray(I)
                SoilModel = MatParams(ISet)%MaterialModel ! name of constitutive model as specified in GOM-file
                
                IsUndrEffectiveStress = &
                !code version 2016 and previous
                ((CalParams%ApplyEffectiveStressAnalysis.and.(trim(MatParams(ISet)%MaterialType)=='2-phase')) .or. &
                !code version 2017.1 and following
                (trim(MatParams(ISet)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE))
          
                if (SoilModel==ESM_MODIFIED_CAM_CLAY) then
                  xM = MatParams(ISet)%MCC_M
                  OCR = MatParams(ISet)%OCR
                  Sig0 = SigmaEff0Array(I,:)
                  SigYield = 0.0
                  if (OCR>1.0) then !overconsolidated
                    K0nc = 1. - 3.*xM / (6.+xM) !compute K0 for normal consolidated clay, not equal to input K0
                    SigYield(1) = K0nc * OCR * Sig0(2) !compute horizontal preconsolidation stress
                    SigYield(2) = OCR * Sig0(2) !compute vertical preconsolidation stress
                    SigYield(3) = SigYield(1)
                  else !normally consolidated
                    SigYield = Sig0
                  end if
                
                                 iOpt = 0
                  call PrnSig(iOpt, SigYield, xN1, xN2, xN3, S1, S2, S3, P, Q)
                  P = MAX(-P,1d0)
                  Particles(I)%PP = Q**2 / (xM**2 * P) + P  ! normal consolidated
         
                  
                  if (IsUndrEffectiveStress) then
                      xNu = MatParams(ISet)%PoissonRatio
                      xkappa = MatParams(ISet)%MCC_Kappa
                      xlambda = MatParams(ISet)%MCC_Lambda
                      xe0 = MatParams(ISet)%InitialVoidRatio
                      nu_u = MatParams(ISet)%UndrainedPoissonRatio

                      xK = (xe0 + 1) / xkappa * P  ! bulk modulus at start of time step, assumed constant
                      r = 3. * ( 1. - 2. * xNu) / ( 2. * (1. + xNu))
                      xG = r * xK

                      Particles(I)%BulkWater = 2. * xG / 3. * &
                          (((1. + nu_u)/(1. + xnu))/(1. - 2. * nu_u) - 1./(1. - 2. * xnu) )
                  end if
      	        end if
               end do
        end subroutine InitialiseStVarMCC

      end module ModMPMStresses