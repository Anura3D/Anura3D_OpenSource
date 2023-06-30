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
    Module WriteOutPut_GiD
    !**********************************************************************
    !
    ! Function: Contains routines for writing output data in Binary and ASCII format
    !           to be read by GiD Posprocess
    !
    !     $Revision: 9408 $
    !     $Date: 2022-02-25 03:53:12 -0500 (Fri, 25 Feb 2022) $
    !
    !**********************************************************************
    !***************Mesh: .mesh archive extension**************************
   
    use ModReadCalculationData
    use ModReadGeometryData
    use ModReadMaterialData
    use ModGlobalConstants
    use ModCounters
    use ModParticle
    use ModMPMData
    use ModWriteVTKASCII
    use ModWriteVTKBinary
    use ModFileIO
	use ModMPMInit
    use ModMPMStresses
      
    #pragma comment(lib, "gidpost.lib")

    use gidpost
    contains 
       
    !************************* Subroutines *******************************
    !*********************************************************************
    
	Subroutine WriteGiDMeshASCII
    !**********************************************************************
    !
    !    This subroutine print the mesh of the model in a unique time step
    !   
    !   No output of this subroutine
    !**********************************************************************

	implicit none 
	
	!Local variables
    integer(INTEGER_TYPE) :: I,J,NumberMaterialPoints       !ID material points
	integer(INTEGER_TYPE) :: TimeStep, IDim, NVECTOR, NumberElements, NNodes, NoMPs
    integer(INTEGER_TYPE), dimension(3) :: NMATElem
    integer(INTEGER_TYPE), dimension(4) :: MatType = 0.0
    real(REAL_TYPE), dimension(3) :: MPCo = 0.0 ! dimensions 3D and 2D
   
    NumberMaterialPoints = Counters%NParticles      ! Total number of material points
    TimeStep = CalParams%IStep                      ! Time step
    NumberElements = Counters%NEl                   ! Total number of elements
    NNodes = Counters%NodTot                        ! Total number of nodes
    NVECTOR = NDIM                                  ! Dimension
    NoMPs = Counters%NParticles                     ! Number of MPs

    CALL GID_OPENPOSTMESHFILE(trim(CalParams%FileNames%ProjectName)//'.POST.MSH',GiD_PostAscii)
    
     CALL GiD_BeginMesh('Material Points', GiD_2D,GiD_Point,1)
        CALL GiD_BeginMeshColor('trimesh',GiD_2D,GiD_Point,1,0.7d0,0.7d0,0.4d0)
        CALL GID_MESHUNIT('m')
        
        CALL GID_BEGINCOORDINATES  
        do I=1, NumberMaterialPoints 
	        do IDim = 1, NVECTOR 
                MPCo(IDim) = GlobPosArray(I,IDim)    
            end do
            CALL GID_WRITECOORDINATES(I,  MPCo(1), MPCo(2), MPCo(3))
        end do
      CALL GID_ENDCOORDINATES
      
      CALL GID_BEGINELEMENTS
      do J = 1, NumberMaterialPoints 
          NMATElem(1) = J
          CALL GID_WRITEELEMENT(J, NMATElem(1))
      enddo
      CALL GID_ENDELEMENTS

    CALL GiD_BeginMesh('Fixed Mesh', GiD_2D,GiD_Triangle,3)
    CALL GiD_BeginMeshColor('trimesh',GiD_2D,GiD_Triangle,3,0.7d0,0.7d0,0.4d0)
    CALL GID_MESHUNIT('m')

    !   **** printing coordinates of the mesh ****
    CALL GID_BEGINCOORDINATES  
      do I=1, NNodes
	    do IDim = 1, NVECTOR 
             MPCo(IDim) = NodalCoordinates(I, IDim)     
        end do
        CALL GID_WRITECOORDINATES(I+NoMPs,  MPCo(1), MPCo(2), MPCo(3))
      end do
      CALL GID_ENDCOORDINATES
      
      ! **** printing elements ****
      CALL GID_BEGINELEMENTS
      
      do J = 1, NumberElements ! Elements empty and filled
          MatType(1) = ElementConnectivities(1, J)+NoMPs  
          MatType(2) = ElementConnectivities(2, J)+NoMPs   
          MatType(3) = ElementConnectivities(3, J)+NoMPs  
          MatType(4) = ElementMaterialID(J)
          if (MatType(4)<=0.0) then 
              MatType(4)=0
          endif
          
          CALL GID_WRITEELEMENTMAT(J, MatType)
      enddo 
      CALL GID_ENDELEMENTS
      
      CALL GID_ENDMESH 
      
     CALL GID_CLOSEPOSTMESHFILE()
    End Subroutine WriteGiDMeshASCII
      
    
    Subroutine WriteGiDResultsASCII
    !**********************************************************************
    !
    !    This subroutine write the results using GiD posprocess module
    !    using ASCII format
    !
    !    Scalars, vectors, and tensor
    !**********************************************************************   
	implicit none 
	
	!Local variables
    integer(INTEGER_TYPE) :: I,J,K,NumberMaterialPoints,TimeStepInt    
	real(REAL_TYPE) :: TimeStep
    real(REAL_TYPE) :: EpsD, SigD, EpsI, MLiquid, MSolid, PWP, Wvolstrain, Liqw, Solidw
    real(REAL_TYPE) :: EpsD_Vol, Eps(3), mean_stress, Strain(6), Stress(6)
    Real(REAL_TYPE), dimension(3) :: Acc = 0.0               
    Real(REAL_TYPE), dimension(3) :: BodyForce = 0.0, &
                                     BodyForceliq = 0.0, &
                                     BodyForcemix = 0.0, &
                                     BodyForcegas = 0.0, &
                                     Uliq = 0.0, &
                                     Usolid = 0.0, & 
                                     Ugas = 0.0, &
                                     Extforce = 0.0, & 
                                     Extforceliq = 0.0, & 
                                     Extforcegas = 0.0, &
                                     Globalpos = 0.0, &
                                     Localpos = 0.0, &
                                     Vliquid = 0.0, &
                                     Vsolid = 0.0, &
                                     Vgas = 0.0

	Character (len=200) :: ARCH_POST_RES, ActualStep 
    Logical :: Hasvalue, file_exists
    
    NumberMaterialPoints = Counters%NParticles   
    TimeStep = CalParams%IStep                   
	TimeStepInt = TimeStep
    
    ARCH_POST_RES = trim(CalParams%FileNames%ProjectName)//'.POST.RES'
    write(ActualStep, '(i3)') TimeStepInt
	
     INQUIRE(FILE= ARCH_POST_RES, EXIST=file_exists)
     

     CALL GID_OPENPOSTRESULTFILE(trim(CalParams%FileNames%ProjectName)//'.POST.RES',GiD_PostAscii)


    If ((CalParams%IStep==1).and.(CalParams%TimeStep==1)) then   
       TimeStep = 0.0
    end if
    
    ! **********************************************************************
    ! ************************ SCALAR RESULTS ******************************
    ! **********************************************************************

    ! ***** Degree_saturation_Liquid *****
    CALL GID_FLUSHPOSTFILE
    CALL GID_BEGINRESULTHEADER('Degree_Saturation','Scalar Results',TimeStep,GiD_Scalar,GiD_onNodes,'Material_Points_Mesh')
    CALL GID_BEGINSCALARRESULT('Degree_Saturation','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL)
        DO I=1,NumberMaterialPoints                    
	     CALL GID_WRITESCALAR(I,Particles(I)%DegreeSaturation)         
        END DO
    CALL GID_ENDRESULT
    ! ***** Density_Liquid *****
    CALL GID_BEGINSCALARRESULT('Density_Liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL)
        DO I=1,NumberMaterialPoints                    
	     CALL GID_WRITESCALAR(I,Particles(I)%Density)         
        END DO
    CALL GID_ENDRESULT
    
    ! ***** Deviatoric_strain_solid *****
    CALL GID_BEGINSCALARRESULT('Dev_strain_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          EpsD = getEpsq(GetEps(Particles(I)))              
          CALL GID_WRITESCALAR(I,EpsD)
      END DO
    CALL GID_ENDRESULT

    ! ***** Deviatoric_stress *****
    CALL GID_BEGINSCALARRESULT('Dev_stress','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          SigD = getq(GetSigmaPrin(Particles(I)))               
          CALL GID_WRITESCALAR(I,SigD)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Incremental_deviatoric_strain_solid *****
    CALL GID_BEGINSCALARRESULT('Incre_dev_strain_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          EpsD = getEpsq(GetEpsStep(Particles(I)))               
          CALL GID_WRITESCALAR(I,EpsD)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Incremental_volumetric_strain_solid *****
    CALL GID_BEGINSCALARRESULT('Incr_vol_strain_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          Eps(1) = GetEpsStepI(Particles(I),1)
          Eps(2) = GetEpsStepI(Particles(I),2)
          Eps(3) = GetEpsStepI(Particles(I),3)
          EpsD_Vol = Eps(1)+Eps(2)+Eps(3)              
          CALL GID_WRITESCALAR(I,EpsD_Vol)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Initial_porosity *****
    CALL GID_BEGINSCALARRESULT('Initial_porosity','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                                           
          CALL GID_WRITESCALAR(I,Particles(I)%InitialPorosity)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Mass_gas *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                                           
          CALL GID_WRITESCALAR(I,Particles(I)%MassGas)
    END DO
    CALL GID_ENDRESULT

    ! ***** Mass_Liquid *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints  
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              MLiquid = 0.0d0
          else
              MLiquid = Particles(I)%MaterialWeight
          end if
          CALL GID_WRITESCALAR(I,MLiquid)
    END DO
    CALL GID_ENDRESULT
   
    ! ***** Mass_mixture *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_mixture','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                                           
          CALL GID_WRITESCALAR(I,Particles(I)%MassMixed)
    END DO
    CALL GID_ENDRESULT

    ! ***** Mass_solid *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints  
            if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              MSolid = MassArray(I)
            else
              MSolid = MassWaterArray(I)
            end if
          CALL GID_WRITESCALAR(I,MSolid)
    END DO
    CALL GID_ENDRESULT
   
    ! ***** 'Mean_eff_stress' *****
    CALL GID_BEGINSCALARRESULT('Mean_eff_stress','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          mean_stress = (SigmaEffArray(I,1) + SigmaEffArray(I,2) + SigmaEffArray(I,3)) / 3.0
          CALL GID_WRITESCALAR(I,mean_stress)
    END DO
    CALL GID_ENDRESULT
            
    ! ***** Porosity *****
    CALL GID_BEGINSCALARRESULT('Porosity','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%Porosity)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Pressure_gas *****
    CALL GID_BEGINSCALARRESULT('Pressure//Gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%GasPressure)
    END DO
    CALL GID_ENDRESULT

    ! ***** Pressure_liquid *****
    CALL GID_BEGINSCALARRESULT('Pressure//Liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
          PWP = (SigmaEffArray(I,1) + SigmaEffArray(I,2) + SigmaEffArray(I,3) ) / 3.0 
        else
          PWP = Particles(I)%WaterPressure
        end if
          CALL GID_WRITESCALAR(I,PWP)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Volumetric_strain_gas *****
    CALL GID_BEGINSCALARRESULT('Vol_strain//gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%GasVolumetricStrain)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Volumetric_strain_liquid *****
    CALL GID_BEGINSCALARRESULT('Vol_strain//liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
            if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              Wvolstrain = GetEpsI(Particles(I),1) + GetEpsI(Particles(I),2) + GetEpsI(Particles(I),3) 
            else
               Wvolstrain = Particles(I)%watervolumetricstrain
            end if          
          CALL GID_WRITESCALAR(I,Wvolstrain)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Volumetric_strain_solid *****    
    CALL GID_BEGINSCALARRESULT('Vol_strain//solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          EpsI = GetEpsI(Particles(I),1) + GetEpsI(Particles(I),2) + GetEpsI(Particles(I),3)
          CALL GID_WRITESCALAR(I,EpsI)
    END DO
    CALL GID_ENDRESULT

    ! ***** Weight_gas *****
    CALL GID_BEGINSCALARRESULT('Weight//gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%GasWeight)
    END DO
    CALL GID_ENDRESULT
    

    ! ***** Weight_liquid *****
    CALL GID_BEGINSCALARRESULT('Weight//liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
            if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              LiqW = Particles(I)%MaterialWeight
            else
              LiqW = Particles(I)%WaterWeight
            end if          
          CALL GID_WRITESCALAR(I,LiqW)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Weight_mixture *****
    CALL GID_BEGINSCALARRESULT('Weight//mixture','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%MixedWeight)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Weight_solid *****    
    CALL GID_BEGINSCALARRESULT('Weight//solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
               SolidW = 0.0d0
            else
               SolidW = Particles(I)%MaterialWeight
            end if
            
          CALL GID_WRITESCALAR(I,SolidW)
    END DO
    CALL GID_ENDRESULT

    ! **********************************************************************
    ! ************************ VECTOR RESULTS ******************************
    ! **********************************************************************
    
    ! ***** Acceleration solid *****
    CALL GID_BEGINRESULTHEADER('Acceleration solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('Accel X','Accel Y','Accel Z','Magnitude')
    CALL GID_RESULTUNIT('m/s2')    
    CALL GiD_BeginVectorResult('Acceleration solid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Accel_x','Accel_y','Accel_z')
      DO I=1,NumberMaterialPoints
          Acc(1) = AccelerationArray(I, 1)
          Acc(2) = AccelerationArray(I, 2)
          if (NDIM==3) Acc(3) = AccelerationArray(I, 3) 
          CALL GID_WRITEVECTOR(I,Acc)
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces *****
    CALL GID_BEGINRESULTHEADER('Body Forces//total'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body force','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BF_x','BF_y','BF_z')
      DO I=1,NumberMaterialPoints
          BodyForce(1) = Particles(I)%FBody(1)
          BodyForce(2) = Particles(I)%FBody(2)
          if (NDIM==3) BodyForce(3) = Particles(I)%FBody(3) 
          CALL GID_WRITEVECTOR(I,BodyForce)
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces liquid *****
    CALL GID_BEGINRESULTHEADER('Body Forces//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body Forces Liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BFW_x','BFW_y','BFW_z')
      DO I=1,NumberMaterialPoints
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
           BodyForceliq(1) = Particles(I)%FBody(1)
           BodyForceliq(2) = Particles(I)%FBody(2)
           if (NDIM==3) BodyForceliq(3) = Particles(I)%FBody(3) 
           CALL GID_WRITEVECTOR(I,BodyForceliq) 
          else
           BodyForceliq(1) = Particles(I)%FBodyWater(1)
           BodyForceliq(2) = Particles(I)%FBodyWater(2)
           if (NDIM==3) BodyForceliq(3) = Particles(I)%FBodyWater(3) 
           CALL GID_WRITEVECTOR(I,BodyForceliq) 
          end if	    
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces mixture ***** 
    CALL GID_BEGINRESULTHEADER('Body Forces//mixture'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body force mixture','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BFM_x','BFM_y','BFM_z')
      DO I=1,NumberMaterialPoints
         BodyForcemix(1) = Particles(I)%FBodyMixed(1)
         BodyForcemix(2) = Particles(I)%FBodyMixed(2)
         if (NDIM==3) BodyForcemix(3) = Particles(I)%FBodyMixed(3) 
          CALL GID_WRITEVECTOR(I,BodyForcemix)
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces gas *****      
    CALL GID_BEGINRESULTHEADER('Body Forces//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body force gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BFG_x','BFG_y','BFG_z')
      DO I=1,NumberMaterialPoints
         BodyForcegas(1) = Particles(I)%FBodyGas(1)
         BodyForcegas(2) = Particles(I)%FBodyGas(2)
         if (NDIM==3) BodyForcegas(3) = Particles(I)%FBodyGas(3) 
          CALL GID_WRITEVECTOR(I,BodyForcegas)
      ENDDO
     CALL GID_ENDRESULT      
 
    ! ***** Displacement liquid ***** 
    CALL GID_BEGINRESULTHEADER('Displacement//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('WX','WY','WZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Displacement liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Displ_x','Displ_y','Displ_z')
      DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then           
            Uliq(1) = UArray(I,1)
            Uliq(2) = UArray(I,2)
            if (NDIM==3) Uliq(3) = UArray(I,3)
          CALL GID_WRITEVECTOR(I,Uliq)
         else
            Uliq(1) = Particles(I)%UW(1)
            Uliq(2) = Particles(I)%UW(2)
            if (NDIM==3) Uliq(3) = Particles(I)%UW(3)   
          CALL GID_WRITEVECTOR(I,Uliq)
        end if     
      ENDDO
     CALL GID_ENDRESULT        

     ! ***** Displacement solid *****
    CALL GID_BEGINRESULTHEADER('Displacement//solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('UX','UY','UZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Displacement solid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Displ_x','Displ_y','Displ_z')
      DO I=1,NumberMaterialPoints
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              Usolid = 0.0
               CALL GID_WRITEVECTOR(I,Usolid)
          else
            Usolid(1) = UArray(I,1)
            Usolid(2) = UArray(I,2)
            if (NDIM==3) Usolid(3) = UArray(I,3)
                CALL GID_WRITEVECTOR(I,Usolid)
          end if
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Displacement gas *****
    CALL GID_BEGINRESULTHEADER('Displacement//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('UgX','UgY','UgZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Displacement gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Displ_x','Displ_y','Displ_z')
      DO I=1,NumberMaterialPoints
         Ugas(1) = Particles(I)%UG(1)
         Ugas(2) = Particles(I)%UG(2)
         if (NDIM==3) Ugas(3) = Particles(I)%UG(3)
         CALL GID_WRITEVECTOR(I,Ugas)   
      ENDDO
    CALL GID_ENDRESULT
    
    ! ***** External force *****    
    CALL GID_BEGINRESULTHEADER('External Forces//mixture'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('ExtFX','ExtFY','ExtFZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('External force','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','ExtF_x','ExtF_y','ExtF_z')
      DO I=1,NumberMaterialPoints
        HasValue = .false.
        do j = 1, NVECTOR
          do k=1,MAX_LOAD_SYSTEMS  
            HasValue = HasValue .or. (Particles(I)%FExt(j,k) /= 0.0)     
          end do  
        end do
        
        if (HasValue) then
            Extforce = 0.0
          do k=1,MAX_LOAD_SYSTEMS
            Extforce(1) = Extforce(1) + Particles(I)%FExt(1,k)
            Extforce(2) = Extforce(2) + Particles(I)%FExt(2,k)
            if (NDIM==3) Extforce(3) = Extforce(3) + Particles(I)%FExt(3,k)
          end do
            CALL GID_WRITEVECTOR(I,Extforce)
        else
           Extforce = 0.0
           CALL GID_WRITEVECTOR(I,Extforce)
        end if
      ENDDO
    CALL GID_ENDRESULT
 
    ! ***** External force liquid *****  
    CALL GID_BEGINRESULTHEADER('External Forces//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('ExtFWX','ExtFWY','ExtFWZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Extforce liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','ExtF_x','ExtF_y','ExtF_z')
      DO I=1,NumberMaterialPoints
        HasValue = .false.
        do j = 1, NVECTOR
          do k=1, MAX_LOAD_SYSTEMS
            HasValue = HasValue .or. (Particles(I)%FExtWater(j,k) /= 0.0)
          END DO
        end do
        if (HasValue) then
          Extforceliq = 0.0
          do k=1, MAX_LOAD_SYSTEMS
            Extforceliq(1) = Extforceliq(1)+Particles(I)%FExtWater(1,k)
            Extforceliq(2) = Extforceliq(2)+Particles(I)%FExtWater(2,k)
            if (NDIM==3) Extforceliq(3) = Extforceliq(3)+Particles(I)%FExtWater(3,k)
          end do
            CALL GID_WRITEVECTOR(I,Extforceliq)
        else
            Extforceliq = 0.0
            CALL GID_WRITEVECTOR(I,Extforceliq)
        end if
      ENDDO
    CALL GID_ENDRESULT

    ! ***** External force gas *****     
    CALL GID_BEGINRESULTHEADER('External Forces//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('ExtFGX','ExtFGY','ExtFGZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Extforce gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','ExtF_x','ExtF_y','ExtF_z')
      DO I=1,NumberMaterialPoints
          HasValue = .false.
          do j = 1, NVECTOR
             do k=1,MAX_LOAD_SYSTEMS 
             HasValue = HasValue .or. (Particles(I)%FExtGas(j,k) /= 0.0)
             END DO
          end do
          if (HasValue) then
              Extforcegas = 0.0
             do k=1,MAX_LOAD_SYSTEMS 
             Extforcegas(1) = Extforcegas(1)+Particles(I)%FExtGas(1,k)
             Extforcegas(2) = Extforcegas(2)+Particles(I)%FExtGas(2,k)
             if (NDIM==3) Extforcegas(3) = Extforcegas(3)+Particles(I)%FExtGas(3,k)
             end do
             CALL GID_WRITEVECTOR(I,Extforcegas)
          else
             Extforcegas = 0.0
             CALL GID_WRITEVECTOR(I,Extforcegas) 
          end if
      ENDDO
    CALL GID_ENDRESULT   

    ! ***** Global position *****         
    CALL GID_BEGINRESULTHEADER('Position//Global'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('GPosX','GPosY','GPosZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Global position','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','PosGl_x','PosGl_y','PosGl_z')
      DO I=1,NumberMaterialPoints
          Globalpos(1) = GlobPosArray(I,1)
          Globalpos(2) = GlobPosArray(I,2)
          if (NDIM==3) Globalpos(3) = GlobPosArray(I,3)  
          CALL GID_WRITEVECTOR(I,Globalpos) 
      ENDDO
    CALL GID_ENDRESULT 
   
    ! ***** Local position *****
    CALL GID_BEGINRESULTHEADER('Position//Local'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('LPosX','LPosY','LPosZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Local position','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','LocPos_x','LocPos_y','LocPos_z')
      DO I=1,NumberMaterialPoints
          Localpos(1) = Particles(I)%LocPos(1)
          Localpos(2) = Particles(I)%LocPos(2)
          if (NDIM==3) Localpos(3) = Particles(I)%LocPos(3)
          CALL GID_WRITEVECTOR(I,Localpos) 
      ENDDO
    CALL GID_ENDRESULT

    ! ***** Velocity liquid *****
    CALL GID_BEGINRESULTHEADER('Velocity//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('WX','WY','WZ','Magnitude')
    CALL GID_RESULTUNIT('m/s')    
    CALL GiD_BeginVectorResult('Velocity liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Veloc_x','Veloc_y','Veloc_z')
      DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid)))then
            Vliquid(1) = VelocityArray(I,1)
            Vliquid(2) = VelocityArray(I,2)
            if (NDIM==3) Vliquid(3) = VelocityArray(I,3)
            CALL GID_WRITEVECTOR(I,Vliquid)
         else
            Vliquid(1) = VelocityWaterArray(I,1)
            Vliquid(2) = VelocityWaterArray(I,2)
            if (NDIM==3) Vliquid(3) = VelocityWaterArray(I,3)
            CALL GID_WRITEVECTOR(I,Vliquid)
        end if   
      ENDDO
    CALL GID_ENDRESULT  

    ! ***** Velocity solid *****
    CALL GID_BEGINRESULTHEADER('Velocity//solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('VX','VY','VZ','Magnitude')
    CALL GID_RESULTUNIT('m/s')    
    CALL GiD_BeginVectorResult('Velocity solid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Veloc_x','Veloc_y','Veloc_z')
      DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
            Vsolid = 0.0
            CALL GID_WRITEVECTOR(I,Vsolid)
         else
            Vsolid(1) = VelocityArray(I,1)
            Vsolid(2) = VelocityArray(I,2)
            if (NDIM==3) Vsolid(3) = VelocityArray(I,3)
            CALL GID_WRITEVECTOR(I,Vsolid)
         end if      
      ENDDO
    CALL GID_ENDRESULT  

    ! ***** Velocity gas *****
    CALL GID_BEGINRESULTHEADER('Velocity//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('VX','VY','VZ','Magnitude')
    CALL GID_RESULTUNIT('m/s')    
    CALL GiD_BeginVectorResult('Velocity gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Veloc_x','Veloc_y','Veloc_z')
      DO I=1,NumberMaterialPoints
        Vgas(1) = VelocityGasArray(I,1)
        Vgas(2) = VelocityGasArray(I,2)
        if (NDIM==3) Vgas(3) = VelocityGasArray(I,3)
        CALL GID_WRITEVECTOR(I,Vgas)
      ENDDO
    CALL GID_ENDRESULT  
    
    ! **********************************************************************
    ! ************************ TENSOR RESULTS ******************************
    ! **********************************************************************
    
    ! ***** Strains ***** 
    CALL GID_BEGINRESULTHEADER('Strains'//char(0),'Tensor Results'//char(0), TimeStep, GiD_Matrix, GiD_onNodes, GiD_NULL)
    CALL GiD_3DMatrixComp('EpsXX','EpsYY','EpsZZ','EpsXY','EpsYZ','EpsXZ')
    CALL GID_RESULTUNIT('-') 
    CALL GiD_Begin3DMatResult('Strains','Tensor Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,'Comp1','Comp2','Comp3','Comp4','Comp5','Comp6')
      DO I=1,NumberMaterialPoints
          if (NDIM == 3) then 
          Strain(1) = Particles(I)%Eps(1)
          Strain(2) = Particles(I)%Eps(2)
          Strain(3) = Particles(I)%Eps(3)
          Strain(4) = Particles(I)%Eps(4)
          Strain(5) = Particles(I)%Eps(5)
          Strain(6) = Particles(I)%Eps(6)
          call GiD_Write3DMatrix(I,Strain(1),Strain(2),Strain(3),Strain(4),Strain(5),Strain(6))
          else if (NDIM == 2) then 
          Strain(1) = Particles(I)%Eps(1)
          Strain(2) = Particles(I)%Eps(2)
          Strain(3) = Particles(I)%Eps(3)
          Strain(4) = Particles(I)%Eps(4)
          Strain(5) = 0.0
          Strain(6) = 0.0
          call GiD_Write3DMatrix(I,Strain(1),Strain(2),Strain(3),Strain(4),Strain(5),Strain(6))
          endif
      ENDDO
    CALL GID_ENDRESULT 

    ! ***** Effective stress solid ***** 
    CALL GID_BEGINRESULTHEADER('Eff_stress_solid'//char(0),'Tensor Results'//char(0), TimeStep, GiD_Matrix, GiD_onNodes, GiD_NULL)
    CALL GiD_3DMatrixComp('SXX','SYY','SZZ','SXY','SYZ','SXZ')
    CALL GID_RESULTUNIT('kPa') 
    CALL GiD_Begin3DMatResult('Eff_stress_solid','Tensor Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,'Comp1','Comp2','Comp3','Comp4','Comp5','Comp6')
      DO I=1,NumberMaterialPoints
          if ((MaterialPointTypeArray(I)==MaterialPointTypeSolid).or. &
            (MaterialPointTypeArray(I)==MaterialPointTypeMixture)) then
          if (NDIM == 3) then 
            Stress(1) = SigmaEffArray(I,1)
            Stress(2) = SigmaEffArray(I,2)
            Stress(3) = SigmaEffArray(I,3)
            Stress(4) = SigmaEffArray(I,4)
            Stress(5) = SigmaEffArray(I,5)
            Stress(6) = SigmaEffArray(I,6)
            call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))    
          else if (NDIM == 2) then 
            Stress(1) = SigmaEffArray(I,1)
            Stress(2) = SigmaEffArray(I,2)
            Stress(3) = SigmaEffArray(I,3)
            Stress(4) = SigmaEffArray(I,4)
            Stress(5) = 0.0
            Stress(6) = 0.0
            call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))              
          end if
        end if
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
            Stress(1) = 0.0
            Stress(2) = 0.0
            Stress(3) = 0.0
            Stress(4) = 0.0
            Stress(5) = 0.0
            Stress(6) = 0.0
            call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))  
          end if
      ENDDO
    CALL GID_ENDRESULT 
    
    ! ***** Stress liquid *****
    CALL GID_BEGINRESULTHEADER('Stress_liquid'//char(0),'Tensor Results'//char(0), TimeStep, GiD_Matrix, GiD_onNodes, GiD_NULL)
    CALL GiD_3DMatrixComp('SXX','SYY','SZZ','SXY','SYZ','SXZ')
    CALL GID_RESULTUNIT('kPa') 
    CALL GiD_Begin3DMatResult('Stress_liquid','Tensor Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,'Comp1','Comp2','Comp3','Comp4','Comp5','Comp6')
      DO I=1,NumberMaterialPoints
         if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
            if (NDIM == 3) then 
                Stress(1) = SigmaEffArray(I,1)
                Stress(2) = SigmaEffArray(I,2)
                Stress(3) = SigmaEffArray(I,3)
                Stress(4) = SigmaEffArray(I,4)
                Stress(5) = SigmaEffArray(I,5)
                Stress(6) = SigmaEffArray(I,6)                
                call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))             
            elseif (NDIM == 2) then
                Stress(1) = SigmaEffArray(I,1)
                Stress(2) = SigmaEffArray(I,2)
                Stress(3) = SigmaEffArray(I,3)
                Stress(4) = SigmaEffArray(I,4)
                Stress(5) = 0.0
                Stress(6) = 0.0
                call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))  
            end if
        end if
                
        if ((MaterialPointTypeArray(I)==MaterialPointTypeSolid).or. &
                (MaterialPointTypeArray(I)==MaterialPointTypeMixture)) then
                Stress(1) = 0.0
                Stress(2) = 0.0
                Stress(3) = 0.0
                Stress(4) = 0.0
                Stress(5) = 0.0
                Stress(6) = 0.0
                call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6)) 
        end if      
      ENDDO
    CALL GID_ENDRESULT 

    if (TimeStep==CalParams%NLoadSteps) then
    CALL GID_CLOSEPOSTRESULTFILE() 
    endif
    
    End Subroutine WriteGiDResultsASCII
    

    
    Subroutine WriteGiDResultsBIN
    !**********************************************************************
    !
    !    This subroutine write the results using GiD posprocess module
    !    using Binary format
    !
    !    Scalars, vectors, and tensor
    !**********************************************************************   
 
	implicit none 
	
	!Local variables
    integer(INTEGER_TYPE) :: I,J, K, NumberMaterialPoints       
    integer(INTEGER_TYPE) :: IDim, NVECTOR, NumberElements, NNodes, NoMPs
	real(REAL_TYPE) :: TimeStep, TimeStepEnd
    real(REAL_TYPE) :: EpsD, SigD, EpsI, MLiquid, MSolid, PWP, Wvolstrain, Liqw, Solidw
    real(REAL_TYPE) :: EpsD_Vol, Eps(3), mean_stress, Strain(6), Stress(6)
    integer(INTEGER_TYPE), dimension(3) :: NMATElem
    integer(INTEGER_TYPE), dimension(4) :: MatType = 0.0 
    integer(INTEGER_TYPE), dimension(5) :: MatType3D = 0.0 
    Real(REAL_TYPE), dimension(3) :: Acc = 0.0, MPCo = 0.0           
    Real(REAL_TYPE), dimension(3) :: BodyForce = 0.0, &
                                     BodyForceliq = 0.0, &
                                     BodyForcemix = 0.0, &
                                     BodyForcegas = 0.0, &
                                     Uliq = 0.0, &
                                     Usolid = 0.0, & 
                                     Ugas = 0.0, &
                                     Extforce = 0.0, & 
                                     Extforceliq = 0.0, & 
                                     Extforcegas = 0.0, &
                                     Globalpos = 0.0, &
                                     Localpos = 0.0, &
                                     Vliquid = 0.0, &
                                     Vsolid = 0.0, &
                                     Vgas = 0.0
    Logical :: Hasvalue
    Character (len=200) :: ARCH_POST_RES
    
    NumberMaterialPoints = Counters%NParticles   
    TimeStep = CalParams%IStep                   
    NumberElements = Counters%NEl                   
    NNodes = Counters%NodTot                        
    NVECTOR = NDIM                                 
    NoMPs = Counters%NParticles                     
    
    ARCH_POST_RES = trim(CalParams%FileNames%ProjectName)//'.POST.lst'


    CALL GID_OPENPOSTRESULTFILE(trim(CalParams%FileNames%ProjectName)//trim(CalParams%FileNames%LoadStepExt)//'.POST.bin',GiD_PostBinary)    
    
    ! ***** Opening list for print results ***** 
    If ((CalParams%IStep==1).and.(CalParams%TimeStep==1)) then
    OPEN (10,FILE=ARCH_POST_RES, STATUS ='UNKNOWN')
    WRITE(10,12)
    WRITE(10,14)
    WRITE(10,'(A)') trim(CalParams%FileNames%ProjectName)//trim(CalParams%FileNames%LoadStepExt)//'.POST.bin'
    CLOSE (10) 
    else

        OPEN (10,FILE=ARCH_POST_RES, STATUS ='OLD', POSITION = 'APPEND')   ! Open the exixting post-proccesing file
        WRITE(10,'(A)') trim(CalParams%FileNames%ProjectName)//trim(CalParams%FileNames%LoadStepExt)//'.POST.bin'
        CLOSE (10) ! close post-proccesing file

    endif

    If ((CalParams%IStep==1).and.(CalParams%TimeStep==1)) then   
       TimeStep = 0.0
    end if
    TimeStepEnd = 0.0 
    
    if(CalParams%PreviouslyRealisedLoadStep/=0) then
        k=CalParams%PreviouslyRealisedLoadStep+1.0
    else
        k=0
    endif
    
    ! ********** PRINTING THE MESH **********
    ! ***************************************
    if (NDIM==2) then
     CALL GiD_BeginMesh('Material Points', GiD_2D,GiD_Point,1)
        CALL GiD_BeginMeshColor('Material Points',GiD_2D,GiD_Point,1,0.7d0,0.7d0,0.4d0)
        CALL GID_MESHUNIT('m')
        
        CALL GID_BEGINCOORDINATES  
        do I=1, NumberMaterialPoints 
	        do IDim = 1, NVECTOR 
                MPCo(IDim) = GlobPosArray(I,IDim)    
            end do
            CALL GID_WRITECOORDINATES(I,  MPCo(1), MPCo(2), MPCo(3))
        end do
      CALL GID_ENDCOORDINATES
      
      CALL GID_BEGINELEMENTS
      do J = 1, NumberMaterialPoints !
          NMATElem(1) = J
          CALL GID_WRITEELEMENT(J, NMATElem(1))
      enddo
      CALL GID_ENDELEMENTS

    CALL GiD_BeginMesh('Fixed Mesh', GiD_2D,GiD_Triangle,3)
    CALL GiD_BeginMeshColor('Fixed Mesh',GiD_2D,GiD_Triangle,3,0.7d0,0.7d0,0.4d0)
    CALL GID_MESHUNIT('m')
    
    ! ***** printing coordinates of the mesh *****
    CALL GID_BEGINCOORDINATES  
      do I=1, NNodes 
	    do IDim = 1, NVECTOR ! Loop over dimensions of IElement
             MPCo(IDim) = NodalCoordinates(I, IDim) 
        end do
        CALL GID_WRITECOORDINATES(I+NoMPs,  MPCo(1), MPCo(2), MPCo(3))
      end do
      CALL GID_ENDCOORDINATES
      
      ! ***** printing elements *****
      CALL GID_BEGINELEMENTS
      
      do J = 1, NumberElements 
          MatType(1) = ElementConnectivities(1, J)+NoMPs  
          MatType(2) = ElementConnectivities(2, J)+NoMPs   
          MatType(3) = ElementConnectivities(3, J)+NoMPs  
          MatType(4) = ElementMaterialID(J)
          if (MatType(4)<=0.0) then 
              MatType(4)=0
          endif
          
          CALL GID_WRITEELEMENTMAT(J, MatType)
      enddo 
      CALL GID_ENDELEMENTS
      
      CALL GID_ENDMESH
    
    else ! 3D
    
     CALL GiD_BeginMesh('Material Points', GiD_3D,GiD_Point,1)
        CALL GiD_BeginMeshColor('Material Points',GiD_3D,GiD_Point,1,0.7d0,0.7d0,0.4d0)
        CALL GID_MESHUNIT('m')
        
        CALL GID_BEGINCOORDINATES  
        do I=1, NumberMaterialPoints 
	        do IDim = 1, NVECTOR 
                MPCo(IDim) = GlobPosArray(I,IDim)    
            end do
            CALL GID_WRITECOORDINATES(I,  MPCo(1), MPCo(2), MPCo(3))
        end do
      CALL GID_ENDCOORDINATES
      
      CALL GID_BEGINELEMENTS
      do J = 1, NumberMaterialPoints 
          NMATElem(1) = J
          CALL GID_WRITEELEMENT(J, NMATElem(1))
      enddo
      CALL GID_ENDELEMENTS
    
    CALL GiD_BeginMesh('Fixed Mesh', GiD_3D,GiD_Tetrahedra,4)
    CALL GiD_BeginMeshColor('Fixed Mesh',GiD_3D,GiD_Tetrahedra,10,0.7d0,0.7d0,0.4d0)
    CALL GID_MESHUNIT('m')
    
    ! ***** printing coordinates of the mesh *****
    CALL GID_BEGINCOORDINATES  
      do I=1, NNodes
	    do IDim = 1, NVECTOR 
             MPCo(IDim) = NodalCoordinates(I, IDim)  
        end do
        CALL GID_WRITECOORDINATES(I+NoMPs,  MPCo(1), MPCo(2), MPCo(3))
      end do
      CALL GID_ENDCOORDINATES
      
    ! ***** printing elements *****
      CALL GID_BEGINELEMENTS
      
      do J = 1, NumberElements 
          MatType3D(1) = ElementConnectivities(1, J)+NoMPs  
          MatType3D(2) = ElementConnectivities(2, J)+NoMPs   
          MatType3D(3) = ElementConnectivities(3, J)+NoMPs  
          MatType3D(4) = ElementConnectivities(4, J)+NoMPs 
          MatType3D(5) = ElementMaterialID(J)
          
          CALL GID_WRITEELEMENT(J, MatType3D)
      enddo 
      CALL GID_ENDELEMENTS
      
      CALL GID_ENDMESH
     endif
    
    ! **********************************************************************
    ! ************************ SCALAR RESULTS ******************************
    ! **********************************************************************
    
    ! ***** Degree_Saturation_Liquid *****
    CALL GID_BEGINRESULTHEADER('Degree_Saturation','Scalar Results',TimeStep,GiD_Scalar,GiD_onNodes,'Material_Points_Mesh')
    CALL GID_BEGINSCALARRESULT('Degree_Saturation','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL)
        DO I=1,NumberMaterialPoints                    
	     CALL GID_WRITESCALAR(I,Particles(I)%DegreeSaturation)         
        END DO
    CALL GID_ENDRESULT
    
    ! ***** Density_Liquid *****
    CALL GID_BEGINSCALARRESULT('Density_Liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL)
        DO I=1,NumberMaterialPoints                    ! loop over material points
	     CALL GID_WRITESCALAR(I,Particles(I)%Density)         
        END DO
    CALL GID_ENDRESULT
    
    ! ***** Deviatoric_strain_solid *****
    CALL GID_BEGINSCALARRESULT('Dev_strain_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          EpsD = getEpsq(GetEps(Particles(I)))              
          CALL GID_WRITESCALAR(I,EpsD)
      END DO
    CALL GID_ENDRESULT

    ! ***** Deviatoric_stress *****
    CALL GID_BEGINSCALARRESULT('Dev_stress','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          SigD = getq(GetSigmaPrin(Particles(I)))               
          CALL GID_WRITESCALAR(I,SigD)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Incremental_deviatoric_strain_solid *****
    CALL GID_BEGINSCALARRESULT('Incre_dev_strain_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          EpsD = getEpsq(GetEpsStep(Particles(I)))               
          CALL GID_WRITESCALAR(I,EpsD)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Incremental_volumetric_strain_solid *****
    CALL GID_BEGINSCALARRESULT('Incr_vol_strain_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                             
          Eps(1) = GetEpsStepI(Particles(I),1)
          Eps(2) = GetEpsStepI(Particles(I),2)
          Eps(3) = GetEpsStepI(Particles(I),3)
          EpsD_Vol = Eps(1)+Eps(2)+Eps(3)              
          CALL GID_WRITESCALAR(I,EpsD_Vol)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Initial_porosity *****
    CALL GID_BEGINSCALARRESULT('Initial_porosity','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                                           
          CALL GID_WRITESCALAR(I,Particles(I)%InitialPorosity)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Mass_gas *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                                           
          CALL GID_WRITESCALAR(I,Particles(I)%MassGas)
    END DO
    CALL GID_ENDRESULT

    ! ***** Mass_Liquid *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints  
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              MLiquid = 0.0d0
          else
              MLiquid = Particles(I)%MaterialWeight
          end if
          CALL GID_WRITESCALAR(I,MLiquid)
    END DO
    CALL GID_ENDRESULT
   
    ! ***** Mass_mixture *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_mixture','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints                                           
          CALL GID_WRITESCALAR(I,Particles(I)%MassMixed)
    END DO
    CALL GID_ENDRESULT

    ! ***** Mass_solid *****
    CALL GID_BEGINSCALARRESULT('Mass//Mass_solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints  
            if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              MSolid = MassArray(I)
            else
              MSolid = MassWaterArray(I)
            end if
          CALL GID_WRITESCALAR(I,MSolid)
    END DO
    CALL GID_ENDRESULT
   
    ! ***** 'Mean_eff_stress' *****
    CALL GID_BEGINSCALARRESULT('Mean_eff_stress','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          mean_stress = (SigmaEffArray(I,1) + SigmaEffArray(I,2) + SigmaEffArray(I,3)) / 3.0
          CALL GID_WRITESCALAR(I,mean_stress)
    END DO
    CALL GID_ENDRESULT
            
    ! ***** Porosity *****
    CALL GID_BEGINSCALARRESULT('Porosity','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%Porosity)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Pressure_gas *****
    CALL GID_BEGINSCALARRESULT('Pressure//Gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%GasPressure)
    END DO
    CALL GID_ENDRESULT

    ! ***** Pressure_liquid *****
    CALL GID_BEGINSCALARRESULT('Pressure//Liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
          PWP = (SigmaEffArray(I,1) + SigmaEffArray(I,2) + SigmaEffArray(I,3) ) / 3.0 
        else
          PWP = Particles(I)%WaterPressure
        end if
          CALL GID_WRITESCALAR(I,PWP)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Volumetric_strain_gas *****
    CALL GID_BEGINSCALARRESULT('Volumetric_strain//gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%GasVolumetricStrain)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Volumetric_strain_liquid *****
    CALL GID_BEGINSCALARRESULT('Volumetric_strain//liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
            if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              Wvolstrain = GetEpsI(Particles(I),1) + GetEpsI(Particles(I),2) + GetEpsI(Particles(I),3) 
            else
               Wvolstrain = Particles(I)%watervolumetricstrain
            end if          
          CALL GID_WRITESCALAR(I,Wvolstrain)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Volumetric_strain_solid *****   
    CALL GID_BEGINSCALARRESULT('Volumetric_strain//solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          EpsI = GetEpsI(Particles(I),1) + GetEpsI(Particles(I),2) + GetEpsI(Particles(I),3)
          CALL GID_WRITESCALAR(I,EpsI)
    END DO
    CALL GID_ENDRESULT

    ! ***** Weight_gas *****
    CALL GID_BEGINSCALARRESULT('Weight//gas','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%GasWeight)
    END DO
    CALL GID_ENDRESULT
    

    ! ***** Weight_liquid *****
    CALL GID_BEGINSCALARRESULT('Weight//liquid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
            if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              LiqW = Particles(I)%MaterialWeight
            else
              LiqW = Particles(I)%WaterWeight
            end if          
          CALL GID_WRITESCALAR(I,LiqW)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Weight_mixture *****
    CALL GID_BEGINSCALARRESULT('Weight//mixture','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,Particles(I)%MixedWeight)
    END DO
    CALL GID_ENDRESULT
    
    ! ***** Weight_solid *****   
    CALL GID_BEGINSCALARRESULT('Weight//solid','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
               SolidW = 0.0d0
            else
               SolidW = Particles(I)%MaterialWeight
            end if
            
          CALL GID_WRITESCALAR(I,SolidW)
    END DO
    CALL GID_ENDRESULT

    ! ***** State Variables *****    
    ! * 1 *
    CALL GID_BEGINSCALARRESULT('State Variables//1','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,1))
      END DO
    CALL GID_ENDRESULT
    ! * 2 *
    CALL GID_BEGINSCALARRESULT('State Variables//2','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,2))
      END DO
    CALL GID_ENDRESULT
    ! * 3 *
    CALL GID_BEGINSCALARRESULT('State Variables//3','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,3))
      END DO
    CALL GID_ENDRESULT
    ! * 4 *
    CALL GID_BEGINSCALARRESULT('State Variables//4','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,4))
      END DO
    CALL GID_ENDRESULT
    ! * 5 *
    CALL GID_BEGINSCALARRESULT('State Variables//5','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,5))
      END DO
    CALL GID_ENDRESULT
    ! * 6 *
    CALL GID_BEGINSCALARRESULT('State Variables//6','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,6))
      END DO
    CALL GID_ENDRESULT
    ! * 7 *
    CALL GID_BEGINSCALARRESULT('State Variables//7','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,7))
      END DO
    CALL GID_ENDRESULT
    ! * 8 *
    CALL GID_BEGINSCALARRESULT('State Variables//8','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,8))
      END DO
    CALL GID_ENDRESULT
    ! * 9 *
    CALL GID_BEGINSCALARRESULT('State Variables//9','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,9))
      END DO
    CALL GID_ENDRESULT
    ! * 10 *
    CALL GID_BEGINSCALARRESULT('State Variables//10','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,10))
      END DO
    CALL GID_ENDRESULT
    ! * 11 *
    CALL GID_BEGINSCALARRESULT('State Variables//11','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,11))
      END DO
    CALL GID_ENDRESULT
    ! * 12 *
    CALL GID_BEGINSCALARRESULT('State Variables//12','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,12))
      END DO
    CALL GID_ENDRESULT
    ! * 13 *
    CALL GID_BEGINSCALARRESULT('State Variables//13','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,13))
      END DO
    CALL GID_ENDRESULT
    ! * 14 *
    CALL GID_BEGINSCALARRESULT('State Variables//14','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,14))
      END DO
    CALL GID_ENDRESULT
    ! * 15 *
    CALL GID_BEGINSCALARRESULT('State Variables//15','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,15))
      END DO
    CALL GID_ENDRESULT
    ! * 16 *
    CALL GID_BEGINSCALARRESULT('State Variables//16','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,16))
      END DO
    CALL GID_ENDRESULT
    ! * 17 *
    CALL GID_BEGINSCALARRESULT('State Variables//17','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,17))
      END DO
    CALL GID_ENDRESULT
    ! * 18 *
    CALL GID_BEGINSCALARRESULT('State Variables//18','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,18))
      END DO
    CALL GID_ENDRESULT
    ! * 19 *
    CALL GID_BEGINSCALARRESULT('State Variables//19','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,19))
      END DO
    CALL GID_ENDRESULT
    ! * 20 *
    CALL GID_BEGINSCALARRESULT('State Variables//20','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,20))
      END DO
    CALL GID_ENDRESULT
    ! * 21 *
    CALL GID_BEGINSCALARRESULT('State Variables//21','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,21))
      END DO
    CALL GID_ENDRESULT
    ! * 22 *
    CALL GID_BEGINSCALARRESULT('State Variables//22','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,22))
      END DO
    CALL GID_ENDRESULT
    ! * 23 *
    CALL GID_BEGINSCALARRESULT('State Variables//23','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,23))
      END DO
    CALL GID_ENDRESULT
    ! * 24 *
    CALL GID_BEGINSCALARRESULT('State Variables//24','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,24))
      END DO
    CALL GID_ENDRESULT
    ! * 25 *
    CALL GID_BEGINSCALARRESULT('State Variables//25','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,25))
      END DO
    CALL GID_ENDRESULT
    ! * 26 *
    CALL GID_BEGINSCALARRESULT('State Variables//26','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,26))
      END DO
    CALL GID_ENDRESULT
    ! * 27 *
    CALL GID_BEGINSCALARRESULT('State Variables//27','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,27))
      END DO
    CALL GID_ENDRESULT
    ! * 28 *
    CALL GID_BEGINSCALARRESULT('State Variables//28','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,28))
      END DO
    CALL GID_ENDRESULT
    ! * 29 *
    CALL GID_BEGINSCALARRESULT('State Variables//29','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,29))
      END DO
    CALL GID_ENDRESULT
    ! * 30 *
    CALL GID_BEGINSCALARRESULT('State Variables//30','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,30))
      END DO
    CALL GID_ENDRESULT
    ! * 31 *
    CALL GID_BEGINSCALARRESULT('State Variables//31','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,31))
      END DO
    CALL GID_ENDRESULT
    ! * 32 *
    CALL GID_BEGINSCALARRESULT('State Variables//32','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,32))
      END DO
    CALL GID_ENDRESULT
    ! * 33 *
    CALL GID_BEGINSCALARRESULT('State Variables//33','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,33))
      END DO
    CALL GID_ENDRESULT
    ! * 34 *
    CALL GID_BEGINSCALARRESULT('State Variables//34','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,34))
      END DO
    CALL GID_ENDRESULT
    ! * 35 *
    CALL GID_BEGINSCALARRESULT('State Variables//35','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,35))
      END DO
    CALL GID_ENDRESULT
    ! * 36 *
    CALL GID_BEGINSCALARRESULT('State Variables//36','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,36))
      END DO
    CALL GID_ENDRESULT
    ! * 37 *
    CALL GID_BEGINSCALARRESULT('State Variables//37','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,37))
      END DO
    CALL GID_ENDRESULT
    ! * 38 *
    CALL GID_BEGINSCALARRESULT('State Variables//38','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,38))
      END DO
    CALL GID_ENDRESULT
    ! * 39 *
    CALL GID_BEGINSCALARRESULT('State Variables//39','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,39))
      END DO
    CALL GID_ENDRESULT
    ! * 40 *
    CALL GID_BEGINSCALARRESULT('State Variables//40','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,40))
      END DO
    CALL GID_ENDRESULT
    ! * 41 *
    CALL GID_BEGINSCALARRESULT('State Variables//41','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,41))
      END DO
    CALL GID_ENDRESULT
    ! * 42 *
    CALL GID_BEGINSCALARRESULT('State Variables//42','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,42))
      END DO
    CALL GID_ENDRESULT
    ! * 43 *
    CALL GID_BEGINSCALARRESULT('State Variables//43','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,43))
      END DO
    CALL GID_ENDRESULT
    ! * 44 *
    CALL GID_BEGINSCALARRESULT('State Variables//44','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,44))
      END DO
    CALL GID_ENDRESULT
    ! * 45 *
    CALL GID_BEGINSCALARRESULT('State Variables//45','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,45))
      END DO
    CALL GID_ENDRESULT
    ! * 46 *
    CALL GID_BEGINSCALARRESULT('State Variables//46','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,46))
      END DO
    CALL GID_ENDRESULT
    ! * 47 *
    CALL GID_BEGINSCALARRESULT('State Variables//47','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,47))
      END DO
    CALL GID_ENDRESULT
    ! * 48 *
    CALL GID_BEGINSCALARRESULT('State Variables//48','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,48))
      END DO
    CALL GID_ENDRESULT
    ! * 49 *
    CALL GID_BEGINSCALARRESULT('State Variables//49','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,49))
      END DO
    CALL GID_ENDRESULT
    ! * 50 *
    CALL GID_BEGINSCALARRESULT('State Variables//50','Scalar Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
	  DO I=1,NumberMaterialPoints
          CALL GID_WRITESCALAR(I,ESMstatevArray(I,50))
      END DO
    CALL GID_ENDRESULT    
    ! **********************************************************************
    ! ************************ VECTOR RESULTS ******************************
    ! **********************************************************************
    
    ! ***** Acceleration solid *****
    CALL GID_BEGINRESULTHEADER('Acceleration solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('Accel X','Accel Y','Accel Z','Magnitude')
    CALL GID_RESULTUNIT('m/s2')    
    CALL GiD_BeginVectorResult('Acceleration solid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Accel_x','Accel_y','Accel_z')
      DO I=1,NumberMaterialPoints
          Acc(1) = AccelerationArray(I, 1)
          Acc(2) = AccelerationArray(I, 2)
          if (NDIM==3) Acc(3) = AccelerationArray(I, 3) 
          CALL GID_WRITEVECTOR(I,Acc)
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces *****
    CALL GID_BEGINRESULTHEADER('Body Forces//solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body Force//total','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BF_x','BF_y','BF_z')
      DO I=1,NumberMaterialPoints
          BodyForce(1) = Particles(I)%FBody(1)
          BodyForce(2) = Particles(I)%FBody(2)
          if (NDIM==3) BodyForce(3) = Particles(I)%FBody(3) 
          CALL GID_WRITEVECTOR(I,BodyForce)
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces liquid *****
    CALL GID_BEGINRESULTHEADER('Body Forces//Liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body Forces Liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BFW_x','BFW_y','BFW_z')
      DO I=1,NumberMaterialPoints
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
           BodyForceliq(1) = Particles(I)%FBody(1)
           BodyForceliq(2) = Particles(I)%FBody(2)
           if (NDIM==3) BodyForceliq(3) = Particles(I)%FBody(3) 
           CALL GID_WRITEVECTOR(I,BodyForceliq) 
          else
           BodyForceliq(1) = Particles(I)%FBodyWater(1)
           BodyForceliq(2) = Particles(I)%FBodyWater(2)
           if (NDIM==3) BodyForceliq(3) = Particles(I)%FBodyWater(3) 
           CALL GID_WRITEVECTOR(I,BodyForceliq) 
          end if	    
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces mixture ***** 
    CALL GID_BEGINRESULTHEADER('Body Forces//mixture'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body force mixture','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BFM_x','BFM_y','BFM_z')
      DO I=1,NumberMaterialPoints
         BodyForcemix(1) = Particles(I)%FBodyMixed(1)
         BodyForcemix(2) = Particles(I)%FBodyMixed(2)
         if (NDIM==3) BodyForcemix(3) = Particles(I)%FBodyMixed(3) 
          CALL GID_WRITEVECTOR(I,BodyForcemix)
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Body forces gas *****      
    CALL GID_BEGINRESULTHEADER('Body Forces//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('FX','FY','FZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Body force gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','BFG_x','BFG_y','BFG_z')
      DO I=1,NumberMaterialPoints
         BodyForcegas(1) = Particles(I)%FBodyGas(1)
         BodyForcegas(2) = Particles(I)%FBodyGas(2)
         if (NDIM==3) BodyForcegas(3) = Particles(I)%FBodyGas(3) 
          CALL GID_WRITEVECTOR(I,BodyForcegas)
      ENDDO
     CALL GID_ENDRESULT      
 
    ! ***** Displacement liquid ***** 
    CALL GID_BEGINRESULTHEADER('Displacement//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('WX','WY','WZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Displacement liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Disp_x','Disp_y','Disp_z')
      DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then           
            Uliq(1) = UArray(I,1)
            Uliq(2) = UArray(I,2)
            if (NDIM==3) Uliq(3) = UArray(I,3)
          CALL GID_WRITEVECTOR(I,Uliq)
         else
            Uliq(1) = Particles(I)%UW(1)
            Uliq(2) = Particles(I)%UW(2)
            if (NDIM==3) Uliq(3) = Particles(I)%UW(3)   
          CALL GID_WRITEVECTOR(I,Uliq)
        end if     
      ENDDO
     CALL GID_ENDRESULT        

     ! ***** Displacement solid *****
    CALL GID_BEGINRESULTHEADER('Displacement//solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('UX','UY','UZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Displacement solid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Disp_x','Disp_y','Disp_z')
      DO I=1,NumberMaterialPoints
          if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
              Usolid = 0.0
               CALL GID_WRITEVECTOR(I,Usolid)
          else
            Usolid(1) = UArray(I,1)
            Usolid(2) = UArray(I,2)
            if (NDIM==3) Usolid(3) = UArray(I,3)
                CALL GID_WRITEVECTOR(I,Usolid)
          end if
      ENDDO
     CALL GID_ENDRESULT 

    ! ***** Displacement gas *****
    CALL GID_BEGINRESULTHEADER('Displacement//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('UgX','UgY','UgZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Displacement gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Disp_x','Disp_y','Disp_z')
      DO I=1,NumberMaterialPoints
         Ugas(1) = Particles(I)%UG(1)
         Ugas(2) = Particles(I)%UG(2)
         if (NDIM==3) Ugas(3) = Particles(I)%UG(3)
         CALL GID_WRITEVECTOR(I,Ugas)   
      ENDDO
    CALL GID_ENDRESULT
    
    ! ***** External forces solid *****    
    CALL GID_BEGINRESULTHEADER('External Forces//solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('ExtFX','ExtFY','ExtFZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('External Forces//Total','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','ExtF_x','ExtF_y','ExtF_z')
      DO I=1,NumberMaterialPoints
        HasValue = .false.
        do j = 1, NVECTOR
          do k=1,MAX_LOAD_SYSTEMS  
            HasValue = HasValue .or. (Particles(I)%FExt(j,k) /= 0.0)     
          end do         
        end do
        
        if (HasValue) then
            Extforce = 0.0
          do k=1,MAX_LOAD_SYSTEMS
            Extforce(1) = Extforce(1) + Particles(I)%FExt(1,k)
            Extforce(2) = Extforce(2) + Particles(I)%FExt(2,k)
            if (NDIM==3) Extforce(3) = Extforce(3) + Particles(I)%FExt(3,k)
          end do
            CALL GID_WRITEVECTOR(I,Extforce)
        else
           Extforce = 0.0
           CALL GID_WRITEVECTOR(I,Extforce)
        end if
      ENDDO
    CALL GID_ENDRESULT
 
    CALL GID_FLUSHPOSTFILE
    ! ***** External forces liquid *****  
    CALL GID_BEGINRESULTHEADER('External Forces//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('ExtFWX','ExtFWY','ExtFWZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Extforce liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','ExtF_x','ExtF_y','ExtF_z')
      DO I=1,NumberMaterialPoints
        HasValue = .false.
        do j = 1, NVECTOR
          do k=1, MAX_LOAD_SYSTEMS  
            HasValue = HasValue .or. (Particles(I)%FExtWater(j,k) /= 0.0)
          end do
        end do
        if (HasValue) then
          Extforceliq = 0.0
          do k=1, MAX_LOAD_SYSTEMS
            Extforceliq(1) = Extforceliq(1)+ Particles(I)%FExtWater(1,K)
            Extforceliq(2) = Extforceliq(1) + Particles(I)%FExtWater(2,K)
            if (NDIM==3) Extforceliq(3) =  Extforceliq(3) + Particles(I)%FExtWater(3,K)
          END DO  
            CALL GID_WRITEVECTOR(I,Extforceliq)
        else
            Extforceliq = 0.0
            CALL GID_WRITEVECTOR(I,Extforceliq)
        end if
      ENDDO
    CALL GID_ENDRESULT

    ! ***** External forces gas *****     
    CALL GID_BEGINRESULTHEADER('External Forces//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('ExtFGX','ExtFGY','ExtFGZ','Magnitude')
    CALL GID_RESULTUNIT('kN')    
    CALL GiD_BeginVectorResult('Extforce gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','ExtF_x','ExtF_y','ExtF_z')
      DO I=1,NumberMaterialPoints
          HasValue = .false.
          do j = 1, NVECTOR
            do k=1, MAX_LOAD_SYSTEMS  
             HasValue = HasValue .or. (Particles(I)%FExtGas(j,k) /= 0.0)
            end do
          end do
          if (HasValue) then
           Extforcegas = 0.0
           do k=1,MAX_LOAD_SYSTEMS
             Extforcegas(1) = Extforcegas(1) + Particles(I)%FExtGas(1,k)
             Extforcegas(2) = Extforcegas(2) + Particles(I)%FExtGas(2,k)
             if (NDIM==3) Extforcegas(3) = Extforcegas(3) + Particles(I)%FExtGas(3,k)
           end do
           CALL GID_WRITEVECTOR(I,Extforcegas)
          else
             Extforcegas = 0.0
             CALL GID_WRITEVECTOR(I,Extforcegas) 
          end if
      ENDDO
    CALL GID_ENDRESULT   

    ! ***** Global position *****         
    CALL GID_BEGINRESULTHEADER('Position//Global'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('GPosX','GPosY','GPosZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Global position','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','PosGl_x','PosGl_y','PosGl_z')
      DO I=1,NumberMaterialPoints
          Globalpos(1) = GlobPosArray(I,1)
          Globalpos(2) = GlobPosArray(I,2)
          if (NDIM==3) Globalpos(3) = GlobPosArray(I,3)  
          CALL GID_WRITEVECTOR(I,Globalpos) 
      ENDDO
    CALL GID_ENDRESULT 
   
    ! ***** Local position *****
    CALL GID_BEGINRESULTHEADER('Position//Local'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('LPosX','LPosY','LPosZ','Magnitude')
    CALL GID_RESULTUNIT('m')    
    CALL GiD_BeginVectorResult('Local position','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','LocPos_x','LocPos_y','LocPos_z')
      DO I=1,NumberMaterialPoints
          Localpos(1) = Particles(I)%LocPos(1)
          Localpos(2) = Particles(I)%LocPos(2)
          if (NDIM==3) Localpos(3) = Particles(I)%LocPos(3)
          CALL GID_WRITEVECTOR(I,Localpos) 
      ENDDO
    CALL GID_ENDRESULT

    ! ***** Velocity liquid *****
    CALL GID_BEGINRESULTHEADER('Velocity//liquid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('WX','WY','WZ','Magnitude')
    CALL GID_RESULTUNIT('m/s')    
    CALL GiD_BeginVectorResult('Velocity liquid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Veloc_x','Veloc_y','Veloc_z')
      DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid)))then
            Vliquid(1) = VelocityArray(I,1)
            Vliquid(2) = VelocityArray(I,2)
            if (NDIM==3) Vliquid(3) = VelocityArray(I,3)
            CALL GID_WRITEVECTOR(I,Vliquid)
         else
            Vliquid(1) = VelocityWaterArray(I,1)
            Vliquid(2) = VelocityWaterArray(I,2)
            if (NDIM==3) Vliquid(3) = VelocityWaterArray(I,3)
            CALL GID_WRITEVECTOR(I,Vliquid)
        end if   
      ENDDO
    CALL GID_ENDRESULT  

    ! ***** Velocity solid *****
    CALL GID_BEGINRESULTHEADER('Velocity//solid'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('VX','VY','VZ','Magnitude')
    CALL GID_RESULTUNIT('m/s')    
    CALL GiD_BeginVectorResult('Velocity solid','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Veloc_x','Veloc_y','Veloc_z')
      DO I=1,NumberMaterialPoints
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((NFORMULATION==2).and.(MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
            Vsolid = 0.0
            CALL GID_WRITEVECTOR(I,Vsolid)
         else
            Vsolid(1) = VelocityArray(I,1)
            Vsolid(2) = VelocityArray(I,2)
            if (NDIM==3) Vsolid(3) = VelocityArray(I,3)
            CALL GID_WRITEVECTOR(I,Vsolid)
         end if      
      ENDDO
    CALL GID_ENDRESULT  

    ! ***** Velocity gas *****
    CALL GID_BEGINRESULTHEADER('Velocity//gas'//char(0),'Vector Results'//char(0), TimeStep, GiD_Vector, GiD_onNodes, GiD_NULL)
    CALL GID_VECTORCOMP('VX','VY','VZ','Magnitude')
    CALL GID_RESULTUNIT('m/s')    
    CALL GiD_BeginVectorResult('Velocity gas','Vector Results',TimeStep,GiD_onNodes, GiD_NULL, GiD_NULL,'MPs','Veloc_x','Veloc_y','Veloc_z')
      DO I=1,NumberMaterialPoints
        Vgas(1) = VelocityGasArray(I,1)
        Vgas(2) = VelocityGasArray(I,2)
        if (NDIM==3) Vgas(3) = VelocityGasArray(I,3)
        CALL GID_WRITEVECTOR(I,Vgas)
      ENDDO
    CALL GID_ENDRESULT  
    
    ! **********************************************************************
    ! ************************ TENSOR RESULTS ******************************
    ! **********************************************************************
    
    ! ***** Strains ***** 
    
    CALL GID_BEGINRESULTHEADER('Strains'//char(0),'Tensor Results'//char(0), TimeStep, GiD_Matrix, GiD_onNodes, GiD_NULL)
    CALL GiD_3DMatrixComp('EpsXX','EpsYY','EpsZZ','EpsXY','EpsYZ','EpsXZ')
    CALL GID_RESULTUNIT('-') 
    CALL GiD_Begin3DMatResult('Strains','Tensor Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,'Comp1','Comp2','Comp3','Comp4','Comp5','Comp6')
      DO I=1,NumberMaterialPoints
        if (NDIM == 3) then ! 3D case
        Strain(1) = Particles(I)%Eps(1)
        Strain(2) = Particles(I)%Eps(2)
        Strain(3) = Particles(I)%Eps(3)
        Strain(4) = Particles(I)%Eps(4)
        Strain(5) = Particles(I)%Eps(5)
        Strain(6) = Particles(I)%Eps(6)
        call GiD_Write3DMatrix(I,Strain(1),Strain(2),Strain(3),Strain(4),Strain(5),Strain(6))
        else  !2D case
        Strain(1) = Particles(I)%Eps(1)
        Strain(2) = Particles(I)%Eps(2)
        Strain(3) = Particles(I)%Eps(3)
        Strain(4) = Particles(I)%Eps(4)
        Strain(5) = 0.0
        Strain(6) = 0.0
        call GiD_Write3DMatrix(I,Strain(1),Strain(2),Strain(3),Strain(4),Strain(5),Strain(6))
        endif
      ENDDO
    CALL GID_ENDRESULT 
    CALL GID_FLUSHPOSTFILE
    ! ***** Effective stress solid *****  
    CALL GID_BEGINRESULTHEADER('Eff_stress_solid'//char(0),'Tensor Results'//char(0), TimeStep, GiD_Matrix, GiD_onNodes, GiD_NULL)
    CALL GiD_3DMatrixComp('SXX','SYY','SZZ','SXY','SYZ','SXZ')
    CALL GID_RESULTUNIT('kPa') 
    CALL GiD_Begin3DMatResult('Eff_stress_solid','Tensor Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,'Comp1','Comp2','Comp3','Comp4','Comp5','Comp6')
      DO I=1,NumberMaterialPoints
        if ((MaterialPointTypeArray(I)==MaterialPointTypeSolid).or. &
            (MaterialPointTypeArray(I)==MaterialPointTypeMixture)) then
        if (NDIM == 3) then ! 3D case
        Stress(1) = SigmaEffArray(I,1)
        Stress(2) = SigmaEffArray(I,2)
        Stress(3) = SigmaEffArray(I,3)
        Stress(4) = SigmaEffArray(I,4)
        Stress(5) = SigmaEffArray(I,5)
        Stress(6) = SigmaEffArray(I,6)
        call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))    
        else !2D case
        Stress(1) = SigmaEffArray(I,1)
        Stress(2) = SigmaEffArray(I,2)
        Stress(3) = SigmaEffArray(I,3)
        Stress(4) = SigmaEffArray(I,4)
        Stress(5) = 0.0
        Stress(6) = 0.0
        call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))              
        end if
        end if
        if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
            ((MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
        Stress(1) = 0.0
        Stress(2) = 0.0
        Stress(3) = 0.0
        Stress(4) = 0.0
        Stress(5) = 0.0
        Stress(6) = 0.0
        call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))  
        end if
      ENDDO
    CALL GID_ENDRESULT 
    
    ! ***** Stress liquid *****
    CALL GID_BEGINRESULTHEADER('Stress_liquid'//char(0),'Tensor Results'//char(0), TimeStep, GiD_Matrix, GiD_onNodes, GiD_NULL)
    CALL GiD_3DMatrixComp('SXX','SYY','SZZ','SXY','SYZ','SXZ')
    CALL GID_RESULTUNIT('kPa') 
    CALL GiD_Begin3DMatResult('Stress_liquid','Tensor Results',TimeStep,GiD_onNodes,GiD_NULL,GiD_NULL,'Comp1','Comp2','Comp3','Comp4','Comp5','Comp6')
      DO I=1,NumberMaterialPoints
         if ((MatParams(MaterialIDArray(I))%MaterialType=='1-phase-liquid').or.MatParams(MaterialIDArray(I))%MaterialPhases=='1-phase-liquid'.or. &
                ((MaterialPointTypeArray(I)==MaterialPointTypeLiquid))) then
            if (NDIM == 3) then ! 3D case
                Stress(1) = SigmaEffArray(I,1)
                Stress(2) = SigmaEffArray(I,2)
                Stress(3) = SigmaEffArray(I,3)
                Stress(4) = SigmaEffArray(I,4)
                Stress(5) = SigmaEffArray(I,5)
                Stress(6) = SigmaEffArray(I,6)                
                call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))             
            else ! 2D Case
                Stress(1) = SigmaEffArray(I,1)
                Stress(2) = SigmaEffArray(I,2)
                Stress(3) = SigmaEffArray(I,3)
                Stress(4) = SigmaEffArray(I,4)
                Stress(5) = 0.0
                Stress(6) = 0.0
                call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6))  
            end if
        end if
                
        if ((MaterialPointTypeArray(I)==MaterialPointTypeSolid).or. &
                (MaterialPointTypeArray(I)==MaterialPointTypeMixture)) then
          Stress(1) = 0.0
          Stress(2) = 0.0
          Stress(3) = 0.0
          Stress(4) = 0.0
          Stress(5) = 0.0
          Stress(6) = 0.0
          call GiD_Write3DMatrix(I,Stress(1),Stress(2),Stress(3),Stress(4),Stress(5),Stress(6)) 
        end if      
      ENDDO
    CALL GID_ENDRESULT 
 

    if (TimeStep==CalParams%NLoadSteps) then

    endif
     12    FORMAT('Multiple')
     14    FORMAT('# postprocess files')
     CALL GID_CLOSEPOSTRESULTFILE() ! Close the file
    End Subroutine WriteGiDResultsBIN
    
    
  End Module WriteOutPut_GiD