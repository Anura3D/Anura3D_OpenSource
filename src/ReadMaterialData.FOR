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
    !	Copyright (C) 2024  Members of the Anura3D MPM Research Community 
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
	  
	  
	  module ModReadMaterialData
      !**********************************************************************
      !
      ! Function : Contains routines for intialisation and reading material data from GOM file
      !
      !     $Revision: 10175 $
      !     $Date: 2024-06-26 15:36:42 +0200 (Wed, 26 Jun 2024) $
      !
      !**********************************************************************

      
      use ModGlobalConstants
      use ModReadCalculationData
      use ModString
      use kernel32
      use ModFileIO
      use ModLinearElasticity
      use ModMohrCoulomb
      use ModBingham
      implicit none




      
        type MaterialParameterType
          integer(INTEGER_TYPE) :: &
                 MaterialIndex = -1 ! unique identifier of the material
          character(len=64) :: &
                 MaterialType = "UNDEFINED", & ! material type (number of phases): 1-phase-solid, 1-phase-liquid, 2-phase, 3-phase (v2016)
                 MaterialPhases = "UNDEFINED", & ! number of phasess: 1-phase-solid, 1-phase-liquid, 2-phase, 3-phase (v2017)
                 MaterialName = "UNDEFINED", & ! user-given name of the material
                 MaterialModel = "UNDEFINED", & ! name of constitutive model as defined in GlobalConstants.FOR
                 RetentionCurve = "UNDEFINED", & ! van_genuchten, linear
                 HydraulicConductivityCurve = "UNDEFINED" ! costant, Hillel, Mualem. Hydraulic conductivity k=k(Sr)
          character(len=64) :: &
                 SoilModelDLL = "UNDEFINED"   ! name of the external DLL containing the soil consititutive model in UMAT/VUMAT format
          integer(handle) :: SoilModelDLLHandle ! This is handle to dll
          
          real(REAL_TYPE) :: &
                 InitialPorosity = 0.0, & ! initial porosity n_0
                 DensitySolid = 0.0, & ! density of the solid, rho_s
                 K0Value = 0.0, & ! K0-value of the solid, K_0
                 DensityLiquid = 0.0, & ! density of the liquid, rho_l
                 BulkModulusLiquid = 0.0,  & ! bulk modulus of the liquid, K_l
                 ViscosityLiquid = 0.0, & ! dynamic viscosity of the liquid, nu_l
                 IntrinsicPermeabilityLiquid = 0.0, & ! intrinsic permeability of the liquid, kappa_l
                 DensityGas = 0.0, & ! density of the gas, rho_g
                 BulkModulusGas = 0.0, & ! bulk modulus of the gas, K_g
                 ViscosityGas = 0.0, & ! dynamic viscosity of the gas, nu_g
                 IntrinsicPermeabilityGas = 0.0, & ! intrinsic permeability of the gas, kappa_g
                 YoungModulus = 0.0, & ! Linear Elasticity: Young modulus, E
                 PoissonRatio = 0.0, & ! Linear Elasticity: Poisson ratio, nu
                 FrictionAngle = 0.0, & ! Mohr-Coulomb: friction angle, phi
                 Cohesion = 0.0, & ! Mohr-Coulomb: cohesion, c
                 DilatancyAngle = 0.0, & ! Mohr-Coulomb: dilatancy angle, psi
                 TensileStrength = 0.0, & ! Mohr-Coulomb: tensile strenght for tension cut-off
                 PeakCohesion = 0.0, & ! Mohr-Coulomb Strain-Softening: peak cohesion, c_p
                 ResidualCohesion = 0.0, & ! Mohr-Coulomb Strain-Softening: residual cohesion, c_res
                 PeakFrictionAngle = 0.0, & ! Mohr-Coulomb Strain-Softening: peak friction angle, phi_p
                 ResidualFrictionAngle = 0.0, & ! Mohr-Coulomb Strain-Softening: residual friction angle, phi_res
                 PeakDilatancyAngle = 0.0, & ! Mohr-Coulomb Strain-Softening: peak dilatancy angle, psi_p
                 ResidualDilatancyAngle = 0.0, & ! Mohr-Coulomb Strain-Softening: residual dilatancy angle, psi_res
                 ShapeFactor = 0.0, & ! Mohr-Coulomb Strain-Softening: shape factor
                 SwellingIndexSuc, & ! Swelling Index with respect to suction 
                 InitialDegreeOfSaturation = 1.0, & ! initial degree of saturation, S_l0
                 DensityMixture = 0.0, & ! density of the mixture, rho_mix
                 DensitySubmerged = 0.0, & ! density of the submerged mixture, rho_sub
                 ShearModulus = 0.0, & ! shear modulus, G
                 UndrainedPoissonRatio = 0.0, & ! undrained Poisson ratio, nu_undr
                 WeightMixture = 0.0, & ! weight of the mixture, gamma_mix
                 DryWeight = 0.0, & ! dry weight, gamma_dry
                 WeightGas = 0.0, & ! weight of the gas, gamma_gas
                 WeightLiquid = 0.0, & ! weight of the liquid, gamma_liquid
                 WeightSubmerged = 0.0, & ! weight of the submerged mixture, gamma_sub
                 HydraulicConductivityLiquid = 0.0, & ! hydraulic conductivity (Darcy) of the liquid, k_l
                 HydraulicConductivityGas = 0.0, & ! hydraulic conductivity (Darcy) of the gas, k_g
                 BishopsAlpha = 1.0, & ! Bishops parameter, alpha
                 FluidThresholdDensity = 0.0, & ! temporary fixed value of 1.0
                 MCC_Lambda = 0.0,& 
                 MCC_Kappa = 0.0, &
                 MCC_M = 0.0, &
                 MCC_N = 0.0, &
                 MCC_PC = 0.0, &
                 InitialVoidRatio = 0.0, &
                 HypoPt = 0.0, &
                 IGSmr = 0.0, &
                 IGSmt = 0.0, &
                 BinghamYieldStress , & ! Bingham fluid: yield shear stress                   
                 OCR= 0., & ! overconsolidation ratio
                 Smin_SWRC= 0., &               ! minimum degree of saturation (Retention curve van Genuchten Model)
                 Smax_SWRC= 1., &               ! maximum degree of saturation (Retention curve van Genuchten Model)
                 P0_SWRC= 0., &                 ! Reference pressure (Retention curve van Genuchten Model); P0=15000 kPa (mudstone model) / P0=18000 kPa (bentonite model) / P0=15 kPa (sand)
                 Lambda_SWRC= 0., &                  ! Lambda (Retention curve van Genuchten Model); L=0.36 (mudstone model) / P0=0.38 (bentonite model)
                 av_SWRC= 0., &                   ! Linear factor (Retention curve Linear Model)
                 rexponentHillel_HCC= 1.0   ! exponent of Hillel equation (Hydraulic conductivity curve)
         real(REAL_TYPE), dimension(NPROPERTIES) :: &
                 ESM_Solid = 0.0         ! array of NPROPERTIES material parameters solid for the external soil model in UMAT/VUMAT format
         real(REAL_TYPE), dimension(NSTATEVAR) :: &
                 ESM_Statvar_in = 0.0    ! array of NSTATEVAR initial value of the State Variables for the external soil model in UMAT/VUMAT format
         integer(INTEGER_TYPE) :: UMATDimension = 6 !User can define if the provided UMAT is written for 3D (6 componet stress tensor, default value) or 2D plane strain (4 component stress tensor)        
         procedure(DUMMYESM), pointer, nopass :: ESM_POINTER
        end type MaterialParameterType
        
        type(MaterialParameterType), dimension(:), allocatable, public, save :: MatParams ! stores material parameters

    contains ! subroutines of this module
  

     Subroutine DUMMYESM(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER,&
        DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)


        implicit double precision (a-h, o-z)
        INTEGER NPT, NOEL, IDSET, NSTATEV, NADDVAR, NPROPS, NUMBEROFPHASES, NTENS
        CHARACTER*80 CMNAME     
        DIMENSION STRESS(NTENS), DSTRAN(NTENS),STATEV(NSTATEV), &
        ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
        real(REAL_TYPE) :: Eunloading, PlasticMultiplier
        
      end subroutine DUMMYESM
        

    
        subroutine ReadMaterialParameters()
        !**********************************************************************
        !
        ! Function: Determines GOM file version and calls respective GOM reader
        !
        !**********************************************************************
        implicit none
        
          ! local variables       
          character(len=MAX_FILENAME_LENGTH) :: FileName, FileVersion
          integer(INTEGER_TYPE) :: FileUnit
          character(len=255) :: BName
          integer(INTEGER_TYPE) :: ios

          FileName = trim(CalParams%FileNames%ProjectName)//GOM_FILE_EXTENSION
          FileUnit = TMP_UNIT
          
          ! check if GOM file exists in project folder, otherwise give error and stop execution
          if ( FExist(trim(FileName)) ) then
            call GiveMessage('Reading GOM file (Materials): ' // trim(FileName) )  
          else
            call GiveError('GOM file does not exist!' // NEW_LINE('A') // 'required GOM file: ' // trim(FileName) )
          end if
          
          ! open GOM file
          call FileOpen(FileUnit, trim(FileName))
          
          ! determine current version of GOM file 
          read(FileUnit, '(A)', iostat=ios) BName ! NB: if no version is specified in the header of the GOM file, the default case will be chosen
          call Assert( ios == 0, 'GOM file: Can''t read flag from GOM file.' )
          FileVersion = trim(BName)
          
          ! read GOM data 
          select case (FileVersion) ! read GOM data depending on file version
            case (Anura3D_v2024)               
              call ReadMaterial_v2021(FileUnit,FileVersion)
            case (Anura3D_v2023)               
              call ReadMaterial_v2021(FileUnit,FileVersion)
            case (Anura3D_v2022)               
              call ReadMaterial_v2021(FileUnit,FileVersion)  
            case (Anura3D_v2021)               
              call ReadMaterial_v2021(FileUnit,FileVersion)         
            case (Anura3D_v2019_2)               
              call ReadMaterial_v2021(FileUnit,FileVersion)   
                             
           
            case default  
              call GiveError('Wrong version of GOM file! ' // NEW_LINE('A') // 'Supported GOM versions: ' &
                // trim(Anura3D_v2019_2) // ', ' &
                // trim(Anura3D_v2021) // ', ' &
                // trim(Anura3D_v2022) // ', ' &
                // trim(Anura3D_v2023) // ', ' &
                // trim(Anura3D_v2024)     // '.' )
end select


          ! close GOM file    
          close(FileUnit)
      
        end subroutine ReadMaterialParameters      
      
      
      
        subroutine ReadMaterial_v2021(FileUnit,FileVersion)
        !**********************************************************************
        !
        ! Function : Reads material data from the GOM file and fills MatParams(I)%...
        ! 
        !            GOM version 2021
        !            GOM version 2019.2
        !
        !**********************************************************************
        
        implicit none

          character(len=MAX_FILENAME_LENGTH), intent(in) :: FileVersion
          integer(INTEGER_TYPE), intent(in) :: FileUnit

          ! local variables        
          integer(INTEGER_TYPE) :: ios ! used for error control
          character(len=21) :: messageIOS = 'GOM file: Can''t read '
          character(len=255) :: TName, BName, SelectMaterial, SelectConstitutiveModel
          character(len=255) :: SelectHydraulicModel, SelectHydraulicConductivityModel
          integer(INTEGER_TYPE) :: I, J
          integer(INTEGER_TYPE) :: DumI
          character(len=255) :: DumS
          real(REAL_TYPE) :: DumR
          integer(INTEGER_TYPE), dimension(:), allocatable :: NumberofPhases
          real(REAL_TYPE) :: n, nu_u, nu_eff, E
          integer :: umat
          integer(handle) :: lib_handle
          pointer (p, umat)
          integer(INTEGER_TYPE) :: K = 1          
          
          ! set GOM version number          
            select case (FileVersion)
            case (Anura3D_v2024) 
                call GiveMessage('Reading... ' // Anura3D_v2024)
                CalParams%GOMversion = Anura3D_v2021
            case (Anura3D_v2023) 
                call GiveMessage('Reading... ' // Anura3D_v2023)
                CalParams%GOMversion = Anura3D_v2021
            case (Anura3D_v2022) 
                call GiveMessage('Reading... ' // Anura3D_v2022)
                CalParams%GOMversion = Anura3D_v2021
            case (Anura3D_v2021) 
                call GiveMessage('Reading... ' // Anura3D_v2021)
                CalParams%GOMversion = Anura3D_v2021 
            case (Anura3D_v2019_2) 
                call GiveMessage('Reading... ' // Anura3D_v2021)
                CalParams%GOMversion = Anura3D_v2021 
            end select          

          do ! read GOM-file
            read(FileUnit,'(A)') TName
            BName = TName
            if (trim(BName)=='$$NUMBER_OF_MATERIALS') then ! read number of materials
              read(FileUnit, *, iostat=ios) DumI
              call Assert( ios == 0, messageIOS//trim(BName) )
              call Assert( DumI > 0, 'GOM file: ' //trim(BName)// ' must be larger than 0.' )
              CalParams%NumberOfMaterials = DumI 
              
              call InitialiseMaterialParameters() ! to initialise MatParams(I)%...
              allocate ( NumberOfPhases(CalParams%NumberOfMaterials) )
              NumberOfPhases = 0

              do I = 1, CalParams%NumberOfMaterials ! loop over number of materials

                !------ BASIC MATERIAL PROPERTIES ------ 
                  
                call ReadInteger(FileUnit, BName, DumI) ! read material index $$MATERIAL_INDEX
                call Assert( DumI > 0, 'GOM file: ' //trim(BName)// ' must be larger than 0.' )
                MatParams(I)%MaterialIndex = DumI
                
                call ReadString(FileUnit, BName, DumS) ! read material name $$MATERIAL_NAME
                MatParams(I)%MaterialName = DumS
                
                call ReadString(FileUnit, BName, DumS) ! read material type $$MATERIAL_TYPE
                MatParams(I)%MaterialType = DumS
                
                SelectMaterial = trim(MatParams(I)%MaterialType)
                select case (SelectMaterial) ! read basic properties for solid and liquid phase depending on material type
                    
                case (DRY_SOIL)
                    
                  call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR
                  
                  if (NFORMULATION==2) then
                    ! overwrite intrinsic permeability liquid, kappa_l, $$INTRINSIC_PERMEABILITY_LIQUID with value from CPS file
                    MatParams(I)%IntrinsicPermeabilityLiquid = CalParams%IntrinsicPermeability
                  
                    !Decide material index for the given particles diameters in the CPS 
                    if ( K == 1 ) then 
                        CalParams%FirstSolidMaterialIndex = I
                        K = K + 1
                    else 
                        CalParams%SecondSolidMaterialIndex = I
                    end if 
                  end if

                  ! calculate density of the mixture
                  MatParams(I)%DensityMixture = ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid
                    
                  ! set number of phases
                  NumberOfPhases(I) = 1
                  MatParams(I)%MaterialPhases = '1-phase-solid'

                case (SATURATED_SOIL_DRAINED)
                    
                  call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR

                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR
                  
                  ! calculate density of the mixture
                  MatParams(I)%DensityMixture = ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid +  MatParams(I)%InitialPorosity * MatParams(I)%DensityLiquid 
                  
                  ! set number of phases
                  NumberOfPhases(I) = 1  
                  MatParams(I)%MaterialPhases = '1-phase-solid'
                  
                  ! set calculation type to 'submerged calculation'
                  CalParams%ApplySubmergedCalculation = .true.
                  
                case (SATURATED_SOIL_UNDRAINED_EFFECTIVE)

                  call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR

                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR
                  
                  ! calculate density of the mixture
                  MatParams(I)%DensityMixture = ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid + MatParams(I)%InitialPorosity * MatParams(I)%DensityLiquid 
                    
                  ! set number of phases
                  NumberOfPhases(I) = 1  
                  MatParams(I)%MaterialPhases = '1-phase-solid'
                                    
                  ! set calculation type to 'effective stress analysis'
                  CalParams%ApplyEffectiveStressAnalysis = .true.

                case (SATURATED_SOIL_UNDRAINED_TOTAL)

                  call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR

                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR
                  
                  ! calculate density of the mixture
                  MatParams(I)%DensityMixture = ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid + MatParams(I)%InitialPorosity * MatParams(I)%DensityLiquid 
                    
                  ! set number of phases
                  NumberOfPhases(I) = 1
                  MatParams(I)%MaterialPhases = '1-phase-solid'

                case (SATURATED_SOIL_COUPLED)
                    
                  call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read intrinsic permeability liquid, kappa_l, $$INTRINSIC_PERMEABILITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 1.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%IntrinsicPermeabilityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read bulk modulus liquid, K_l, $$BULK_MODULUS_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 2.1e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR > 2.2e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%BulkModulusLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read dynamic viscosity liquid, mu_l, $$DYNAMIC_VISCOSITY_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%ViscosityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR
                  
                  ! calculate density of the mixture
                  MatParams(I)%DensityMixture = ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid + MatParams(I)%InitialPorosity * MatParams(I)%DensityLiquid 
                  
                  ! calculate hydraulic (Darcy) conductivities
                  MatParams(I)%HydraulicConductivityLiquid = MatParams(I)%DensityLiquid * CalParams%GravityData%GAccel * &
                                                             MatParams(I)%IntrinsicPermeabilityLiquid / &
                                                             (1000 * MatParams(I)%ViscosityLiquid)

                  call GiveMessage('Darcy permeability is ' //trim(String(MatParams(I)%HydraulicConductivityLiquid))// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert(.not.(MatParams(I)%BulkModulusLiquid < 0), 'bulk modulus of liquid must be positive')
                  call AssertWarning(MatParams(I)%BulkModulusLiquid < BULK_MODULUS_WATER, 'bulk modulus of liquid seems to be high, K:' // trim(String(MatParams(I)%BulkModulusLiquid)))
                
                  ! calculate liquid threshold density (temporary solution)
                  MatParams(I)%FluidThresholdDensity = MatParams(I)%DensityLiquid * &
                           ( 1 - CalParams%LiquidPressureCavitationThreshold / MatParams(I)%BulkModulusLiquid ) ! temporary if p=0 (at 1atm) 

                  if ((IsMPMComputation()).and.(MatParams(I)%FluidThresholdDensity<=0.0)) then
                    call GiveError('FluidThresholdDensity is zero or negative!'           // NEW_LINE('A') // &
                                   '$$LIQUID_PRESSURE_CAVITATION_THRESHOLD is too large,' // &
                                   ' or $$BULK_MODULUS_LIQUID is too small.'              // NEW_LINE('A') // &
                                   'Note the requirement: Cavitation pressure < Bulk modulus, and change accordingly.')
                  end if
                  CalParams%FirstSolidMaterialIndex = 1
                  if (NFORMULATION==2) then
                    ! overwrite intrinsic permeability liquid, kappa_l, $$INTRINSIC_PERMEABILITY_LIQUID with value from CPS file
                    MatParams(I)%IntrinsicPermeabilityLiquid = CalParams%IntrinsicPermeability
                  
                    !Decide material index for the given particles diameters in the CPS 
                    if ( K == 1 ) then 
                        CalParams%FirstSolidMaterialIndex = I
                        K = K + 1
                    else 
                        CalParams%SecondSolidMaterialIndex = I
                    end if 
                  end if
                  
                  ! set number of phases
                  NumberOfPhases(I) = 2  
                  MatParams(I)%MaterialPhases = '2-phase'
                  
                CASE (UNSATURATED_SOIL_TWOPHASE)
                   call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read intrinsic permeability liquid, kappa_l, $$INTRINSIC_PERMEABILITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 1.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%IntrinsicPermeabilityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read bulk modulus liquid, K_l, $$BULK_MODULUS_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 2.1e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR > 2.2e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%BulkModulusLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read dynamic viscosity liquid, mu_l, $$DYNAMIC_VISCOSITY_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%ViscosityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR
                  
                                   
                  ! calculate density of the SATURATED mixture
                  MatParams(I)%DensityMixture = ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid + MatParams(I)%InitialPorosity * MatParams(I)%DensityLiquid 
                  
                  ! calculate hydraulic (Darcy) conductivities SATURATED
                  MatParams(I)%HydraulicConductivityLiquid = MatParams(I)%DensityLiquid * CalParams%GravityData%GAccel * &
                                                             MatParams(I)%IntrinsicPermeabilityLiquid / &
                                                             (1000 * MatParams(I)%ViscosityLiquid)

                  call GiveMessage('Darcy permeability is ' //trim(String(MatParams(I)%HydraulicConductivityLiquid))// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert(.not.(MatParams(I)%BulkModulusLiquid < 0), 'bulk modulus of liquid must be positive')
                  call AssertWarning(MatParams(I)%BulkModulusLiquid < BULK_MODULUS_WATER, 'bulk modulus of liquid seems to be high, K:' // trim(String(MatParams(I)%BulkModulusLiquid)))

                  CalParams%FirstSolidMaterialIndex = 1
                  if (NFORMULATION==2) then
                    call GiveWarning('unsaturated material not allowed in double-point formulation')
                  end if
                  
                  ! set number of phases
                  NumberOfPhases(I) = 2  
                  MatParams(I)%MaterialPhases = '2-phase'
                  Calparams%ApplyPartialSaturation = .true.       
          
                case (LIQUID)
                    
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read bulk modulus liquid, K_l, $$BULK_MODULUS_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 2.1e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR > 2.2e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%BulkModulusLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read dynamic viscosity liquid, mu_l, $$DYNAMIC_VISCOSITY_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%ViscosityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read liquid cavitation, $$LIQUID_CAVITATION
                  call AssertWarning( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' should be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  CalParams%LiquidPressureCavitationThreshold = DumR

                  call ReadInteger(FileUnit, BName, DumI) ! read apply detect free surface, $$APPLY_DETECT_LIQUID_SURFACE
                  if ( DumI == 1 ) CalParams%ApplyDetectLiquidFreeSurface = .true.
                  
                  ! set some necessary standard parameters
                  MatParams(I)%PoissonRatio = 0.45 ! set Poisson ratio, relevant for calculation of shear modulus
                  MatParams(I)%ShearModulus = 3.0 * MatParams(I)%BulkModulusLiquid * (1.0 - 2.0 * MatParams(I)%PoissonRatio) / &
                                                 (2 * (1.0 + MatParams(I)%PoissonRatio)) ! shear modulus computed from Poisson ratio and specified bulk modulus
                  MatParams(I)%K0Value = 1.0

                  ! calculate liquid threshold density (temporary solution)
                  MatParams(I)%FluidThresholdDensity = MatParams(I)%DensityLiquid * &
                           ( 1 - CalParams%LiquidPressureCavitationThreshold / MatParams(I)%BulkModulusLiquid ) ! temporary if p=0 (at 1atm) 

                  if ((IsMPMComputation()).and.(MatParams(I)%FluidThresholdDensity<=0.0)) then
                    call GiveError('FluidThresholdDensity is zero or negative!'           // NEW_LINE('A') // &
                                   '$$LIQUID_PRESSURE_CAVITATION_THRESHOLD is too large,' // &
                                   ' or $$BULK_MODULUS_LIQUID is too small.'              // NEW_LINE('A') // &
                                   'Note the requirement: Cavitation pressure < Bulk modulus, and change accordingly.')
                  end if

                  ! calculate density of the mixture (temporary solution)
                  MatParams(I)%DensitySolid = MatParams(I)%DensityLiquid
                  
                  if(NFORMULATION==1) then ! one set of material points (1-point formulation)
                    MatParams(I)%InitialPorosity = 0.0
                    MatParams(I)%DensityMixture = (1-MatParams(I)%InitialPorosity) * MatParams(I)%DensitySolid
                  else if (NFORMULATION==2) then ! two sets of material points (2-point formulation)
                    MatParams(I)%InitialPorosity = 1.0
                    MatParams(I)%DensityMixture = 0.0
                  else
                    call GiveError('CPS-ERROR: Allowed values for $$NUMBER_OF_LAYERS: 1 or 2.')
                  end if

                  ! set number of phases
                  NumberOfPhases(I) = 1
                  MatParams(I)%MaterialPhases = '1-phase-liquid'

                case (UNSATURATED_SOIL_THREEPHASE) 
                  call ReadReal(FileUnit, BName, DumR) ! read initial porosity solid, n_0, $$POROSITY_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%InitialPorosity = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density solid, rho_s, $$DENSITY_SOLID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensitySolid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read density liquid, rho_l, $$DENSITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityLiquid = DumR
                    
                  call ReadReal(FileUnit, BName, DumR) ! read density gas, rho_g, $$DENSITY_GAS
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 10000.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%DensityGas = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read intrinsic permeability liquid, kappa_l, $$INTRINSIC_PERMEABILITY_LIQUID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 1.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%IntrinsicPermeabilityLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read intrinsic permeability liquid, kappa_l, $$INTRINSIC_PERMEABILITY_GAS
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 1.0, 'GOM file: Unrealistic value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%IntrinsicPermeabilityGas = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read bulk modulus liquid, K_l, $$BULK_MODULUS_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 2.1e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR > 2.2e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%BulkModulusLiquid = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read bulk modulus gas, K_g, $$BULK_MODULUS_GAS
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR < 2.1e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call AssertWarning( DumR > 2.2e6, 'GOM file: Non-default value for ' //trim(BName)// ' (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%BulkModulusGas = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read dynamic viscosity liquid, mu_l, $$DYNAMIC_VISCOSITY_LIQUID
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%ViscosityLiquid = DumR                                    
                                    
                  call ReadReal(FileUnit, BName, DumR) ! read dynamic viscosity gas, mu_g, $$DYNAMIC_VISCOSITY_GAS
                  call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%ViscosityGas = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! read swelling Index with respect to suction, Ksuc, $$SWELLING_INDEX
                  !call Assert( DumR > 0.0, 'GOM file: ' //trim(BName)// ' must be larger than 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%SwellingIndexSuc = DumR       
                  
                  call ReadReal(FileUnit, BName, DumR) ! read K0-value, K0, $$K0_VALUE_SOLID
                  call Assert( DumR >= 0.0, 'GOM file: ' //trim(BName)// ' must be larger than or equal to 0.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  call Assert( DumR <= 1.0, 'GOM file: ' //trim(BName)// ' must be smaller than or equal to 1.0 (MATERIAL: ' //trim(MatParams(I)%MaterialName)// ').' )
                  MatParams(I)%K0Value = DumR 
                  
                  MatParams(I)%InitialDegreeOfSaturation = 1.0 ! initially the soil is assumed to be fully saturated
                  
                  ! calculate density of the mixture
                  MatParams(I)%DensityMixture = (1-MatParams(I)%InitialPorosity) * MatParams(I)%DensitySolid &
                        + MatParams(I)%InitialPorosity * MatParams(I)%InitialDegreeOfSaturation * MatParams(I)%DensityLiquid &
                        + MatParams(I)%InitialPorosity * (1-MatParams(I)%InitialDegreeOfSaturation) * MatParams(I)%DensityGas
                  
                  ! calculate hydraulic (Darcy) conductivities
                  MatParams(I)%HydraulicConductivityLiquid = MatParams(I)%DensityLiquid * CalParams%GravityData%GAccel * &
                                                      MatParams(I)%IntrinsicPermeabilityLiquid /  &
                                                      (1000 * MatParams(I)%ViscosityLiquid)
                  MatParams(I)%HydraulicConductivityGas = MatParams(I)%DensityGas * CalParams%GravityData%GAccel * &
                                                      MatParams(I)%IntrinsicPermeabilityGas / MatParams(I)%ViscosityGas
                  
                  CalParams%FirstSolidMaterialIndex = 1
                  if (NFORMULATION==2) then
                    call GiveWarning('unsaturated material not allowed in double-point formulation')
                  end if
                  
                  ! set number of phases
                  NumberOfPhases(I) = 3  
                  MatParams(I)%MaterialPhases = '3-phase'
                   
				case default
                
                  call GiveError('GOM file: $$MATERIAL_TYPE is ' //trim(MatParams(I)%MaterialType)// ' and is therefore not properly specified.')
                
                end select ! material type
                
                !------ CONSTITUTIVE MODEL MATERIAL PROPERTIES ------ 

                call ReadString(FileUnit, BName, DumS) ! read constitutive model type $$MATERIAL_MODEL_SOLID or $$MATERIAL_MODEL_LIQUID
                MatParams(I)%MaterialModel = DumS

                SelectConstitutiveModel = trim(MatParams(I)%MaterialModel)
                
                call GiveMessage('Loading constitutive model... ' // trim(MatParams(I)%MaterialModel) )
                
                select case (SelectConstitutiveModel) ! read constitutive model properties for solid or liquid 
                    
				case (ESM_RIGID_BODY)!Rigid body
					CalParams%RigidBody%IsRigidBody=.true. !Activate Rigid body algorithm
					CalParams%RigidBody%RigidEntity=2 !Set 2 as the entity for rigid material so it works with the contact a.
					if (.not.CalParams%ApplyContactAlgorithm) then !check contact algorithm is enabled
						call GiveError('The contact algorithm is not activated!')
					end if
					call ReadInteger(FileUnit, BName, DumI) !Read x constrain
					CalParams%RigidBody%Constrains(1)=DumI
					call ReadInteger(FileUnit, BName, DumI) !Read y constrain
					CalParams%RigidBody%Constrains(2)=DumI		
					if (NVECTOR==3) then ! 3D case
					call ReadInteger(FileUnit, BName, DumI) !Read z constrain
					CalParams%RigidBody%Constrains(3)=DumI						
					end if					
					
				case (ESM_LINEAR_ELASTICITY)
                    
                  call ReadReal(FileUnit, BName, DumR) ! Young Modulus, E, $$YOUNG_MODULUS
                  MatParams(I)%YoungModulus = DumR
                  MatParams(I)%ESM_Solid(1) = DumR 
                  call ReadReal(FileUnit, BName, DumR) ! Poisson ratio, nu, $$POISSON_RATIO
                  MatParams(I)%PoissonRatio = DumR
                  MatParams(I)%ESM_Solid(2) = DumR
                  
                  if (trim(MatParams(I)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE) then
                    call ReadReal(FileUnit, BName, DumR) ! Poisson ratio undrained, nu_undr, $$POISSON_RATIO_UNDRAINED
                    MatParams(I)%UndrainedPoissonRatio = DumR
                    n = MatParams(I)%InitialPorosity
                    nu_u = MatParams(I)%UndrainedPoissonRatio
                    nu_eff = MatParams(I)%PoissonRatio
                    E = MatParams(I)%YoungModulus
                    MatParams(I)%BulkModulusLiquid = n*(nu_u-nu_eff)*E/   &
                        ((1-2.*nu_u)*(1+nu_eff)*(1-2.*nu_eff))
                  else
                    MatParams(I)%UndrainedPoissonRatio = 0.495   
                  end if                      
                
                  MatParams(I)%ESM_POINTER => ESM_LINEAR
                  
                  if ( NDIM == 2 ) then
                      MatParams(I)%UMATDimension = 4
                  end if 
                case (ESM_MOHR_COULOMB_STANDARD)
                    
                  call ReadReal(FileUnit, BName, DumR) ! Young Modulus, E, $$YOUNG_MODULUS
                  MatParams(I)%YoungModulus = DumR
                    
                  call ReadReal(FileUnit, BName, DumR) ! Poisson ratio, nu, $$POISSON_RATIO
                  MatParams(I)%PoissonRatio = DumR
                    
                  if (trim(MatParams(I)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE) then
                    call ReadReal(FileUnit, BName, DumR) ! Poisson ratio undrained, nu_undr, $$POISSON_RATIO_UNDRAINED
                    MatParams(I)%UndrainedPoissonRatio = DumR
                    n = MatParams(I)%InitialPorosity
                    nu_u = MatParams(I)%UndrainedPoissonRatio
                    nu_eff = MatParams(I)%PoissonRatio
                    E = MatParams(I)%YoungModulus
                    MatParams(I)%BulkModulusLiquid = n*(nu_u-nu_eff)*E/   &
                        ((1-2.*nu_u)*(1+nu_eff)*(1-2.*nu_eff))
                  else 
                    MatParams(I)%UndrainedPoissonRatio = 0.495 
                  end if                      

                  call ReadReal(FileUnit, BName, DumR) ! friction angle, phi, $$FRICTION_ANGLE
                  MatParams(I)%FrictionAngle = DumR
                    
                  call ReadReal(FileUnit, BName, DumR) ! cohesion, c, $$COHESION
                  MatParams(I)%Cohesion = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! dilatancy angle, psi, $$DILATANCY_ANGLE
                  MatParams(I)%DilatancyAngle = DumR
                    
                  call ReadReal(FileUnit, BName, DumR) ! tensile strength, T, $$TENSILE_STRENGTH
                  MatParams(I)%TensileStrength = DumR

                  ! initialise DLL
                  MatParams(I)%ESM_POINTER => ESM_MC              
                  
                  if ( NDIM == 2 ) then
                      MatParams(I)%UMATDimension = 4
                  end if 
                  
                case (ESM_MOHR_COULOMB_STRAIN_SOFTENING)
                    
                  call ReadReal(FileUnit, BName, DumR) ! Young Modulus, E, $$YOUNG_MODULUS
                  MatParams(I)%YoungModulus = DumR
                  MatParams(I)%ESM_Solid(1) = DumR  
                  call ReadReal(FileUnit, BName, DumR) ! Poisson ratio, nu, $$POISSON_RATIO
                  MatParams(I)%PoissonRatio = DumR
                  MatParams(I)%ESM_Solid(2) = DumR 
                  if (trim(MatParams(I)%MaterialType)==SATURATED_SOIL_UNDRAINED_EFFECTIVE) then
                    call ReadReal(FileUnit, BName, DumR) ! Poisson ratio undrained, nu_undr, $$POISSON_RATIO_UNDRAINED
                    MatParams(I)%UndrainedPoissonRatio = DumR
                    n = MatParams(I)%InitialPorosity
                    nu_u = MatParams(I)%UndrainedPoissonRatio
                    nu_eff = MatParams(I)%PoissonRatio
                    E = MatParams(I)%YoungModulus
                    MatParams(I)%BulkModulusLiquid = n*(nu_u-nu_eff)*E/   &
                        ((1-2.*nu_u)*(1+nu_eff)*(1-2.*nu_eff))
                  else
                    MatParams(I)%UndrainedPoissonRatio = 0.495 
                  end if                      

                  call ReadReal(FileUnit, BName, DumR) ! peak cohesion, c_p, $$PEAK_COHESION
                  MatParams(I)%PeakCohesion = DumR
                  MatParams(I)%ESM_Solid(3) = DumR 
                  call ReadReal(FileUnit, BName, DumR) ! residual cohesion, c_r, $$RESIDUAL_COHESION
                  MatParams(I)%ResidualCohesion = DumR
                  MatParams(I)%ESM_Solid(4) = DumR  
                  call ReadReal(FileUnit, BName, DumR) ! peak friction angle, phi_p, $$PEAK_FRICTION_ANGLE
                  MatParams(I)%PeakFrictionAngle = DumR
                  MatParams(I)%ESM_Solid(5) = DumR 
                  call ReadReal(FileUnit, BName, DumR) ! residual friction angle, phi_r, $$RESIDUAL_FRICTION_ANGLE
                  MatParams(I)%ResidualFrictionAngle = DumR
                  MatParams(I)%ESM_Solid(6) = DumR  
                  call ReadReal(FileUnit, BName, DumR) ! peak dilatancy angle, psi_p, $$PEAK_DILATANCY_ANGLE
                  MatParams(I)%PeakDilatancyAngle = DumR
                  MatParams(I)%ESM_Solid(7) = DumR  
                  call ReadReal(FileUnit, BName, DumR) ! residual dilatancy angle, psi_r, $$RESIDUAL_DILATANCY_ANGLE
                  MatParams(I)%ResidualDilatancyAngle = DumR
                  MatParams(I)%ESM_Solid(8) = DumR 
                  call ReadReal(FileUnit, BName, DumR) ! shape factor, $$SHAPE_FACTOR
                  MatParams(I)%ShapeFactor = DumR
                  MatParams(I)%ESM_Solid(9) = DumR
                  
				  !MatParams(I)%ESM_POINTER => ESM_MCSS                    
                   
                  
                case (ESM_EXTERNAL_SOIL_MODEL)
                    
                  !if (MatParams(I)%MaterialName=="MyModel")  MatParams(I)%ESM_POINTER => ESM_Mymodel  !Example of the pointer definition
                    
                  call ReadString(FileUnit, BName, DumS) ! soil model DLL, $$SOIL_MODEL_DLL
                  MatParams(I)%SoilModelDLL = DumS


                  
                  call ReadString(FileUnit, BName, DumS) ! Umat Dimension, $$UMAT_DIMENSION
                  if (DumS=='2D_plane_strain') then
                   MatParams(I)%UMATDimension = 4                  ! 3D Umat
                  else !full_3D_(default)
                   MatParams(I)%UMATDimension = 6                  ! 2D Umat
                  end if

                  if (MatParams(I)%UMATDimension < NTENSOR) call GiveError('UMAT-ERROR: UMAT Dimension is lower then Geometyry Dimension ') ! error                   
                  
                  do J = 1, NPROPERTIES ! reading material parameters (n=NPROPERTIES) solid the dimension of the array of properties is defined as a global constant
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%ESM_Solid(J) = DumR
                  end do
                  
                  do J = 1, NSTATEVAR ! reading initial value of the state variables (n=NSTATEVAR) solid the dimension of the array of properties is defined as a global constant
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%ESM_Statvar_in(J) = DumR
                  end do
                  
                case (ESM_NEWTONIAN_LIQUID) 
                    
                  ! no additional parameters to read  
                    
                case (ESM_BINGHAM_LIQUID)
                    
                  call ReadReal(FileUnit, BName, DumR) ! Bingham liquid: yield shear stress
                  MatParams(I)%ESM_Solid(3) = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! Bingham liquid: Young Modulus, E
                  MatParams(I)%ESM_Solid(1) = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! Bingham liquid: Poisson ratio, nu
                  MatParams(I)%ESM_Solid(2) = DumR
                  
                  
                  
                  MatParams(I)%ESM_POINTER => ESM_BINGHAM   
                  
                  if ( NDIM == 2 ) then
                      MatParams(I)%UMATDimension = 4
                  end if 
                  
                  
                  MatParams(I)%ESM_Solid(4)=MatParams(I)%ViscosityLiquid
                  MatParams(I)%ESM_Solid(5)=MatParams(I)%BulkModulusLiquid
                  MatParams(I)%ESM_Solid(6)=CalParams%LiquidPressureCavitationThreshold
                    
                case (ESM_FRICTIONAL_LIQUID)  
                    
                  call ReadReal(FileUnit, BName, DumR) ! frictional liquid: yield shear stress, $$LIQUID_FRICTION_ANGLE
                  MatParams(I)%FrictionAngle = DumR
                  !CalParams%FrictionalFluidFrictionAngle =DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! frictional liquid: Young Modulus, E, $$LIQUID_YOUNG_MODULUS
                  MatParams(I)%YoungModulus = DumR
                  !CalParams%BinghamYoungModulus = DumR
                  
                  call ReadReal(FileUnit, BName, DumR) ! frictional liquid: Poisson ratio, nu, $$LIQUID_POISSON_RATIO
                  MatParams(I)%PoissonRatio = DumR
                  !CalParams%BinghamPoissonRatio = DumR

                  ! set material model
                  !MatParams(I)%MaterialModel = ESM_MOHR_COULOMB_STANDARD ! material model: Mohr-Coulomb
                
                case default
                
                  call GiveError('GOM-ERROR: $$MATERIAL_MODEL is ' //trim(MatParams(I)%MaterialModel)// ' and is therefore not properly specified.')
                
                end select ! constitutive models for solid or liquid
                
          !---------------------------------------------------------------------------------------------
          ! READING UNSATURATED MATERIAL PROPERTIES
          ! NOT YED IN GID
          !------------------------------------------------------------------------------------

                If (MatParams(I)%MaterialType==UNSATURATED_SOIL_TWOPHASE.or.MatParams(I)%MaterialType==UNSATURATED_SOIL_THREEPHASE) then
                    
                  call ReadString(FileUnit,BName, DumS) 
                  MatParams(I)%RetentionCurve = DumS
                    
                  SelectHydraulicModel = trim(MatParams(I)%RetentionCurve)
                
                  call GiveMessage('Loading water retention curve model... ' // trim(MatParams(I)%RetentionCurve) )
                
                select case (SelectHydraulicModel)
                    
                case (SWRC_VANGENUCHTEN)

                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%Smin_SWRC = DumR ! Minimum degree of saturation   
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%Smax_SWRC = DumR ! Maximum degree of saturation                
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%P0_SWRC = DumR ! Reference Pressure
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%Lambda_SWRC = DumR ! Parameter 

                 case (SWRC_LINEAR)    
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%av_SWRC = DumR ! Linear factor 

                 case default
                
                  call GiveError('GOM-ERROR: $$WATER_RETENTION_CURVE is ' //trim(MatParams(I)%RetentionCurve)// ' and is therefore not properly specified.')
                
                 end select
                 
                 call ReadString(FileUnit,BName, DumS) 
                  MatParams(I)%HydraulicConductivityCurve = DumS
                    
                  SelectHydraulicConductivityModel = trim(MatParams(I)%HydraulicConductivityCurve)
                
                  call GiveMessage('Loading hydraulic conductivity curve model... ' // trim(MatParams(I)%HydraulicConductivityCurve) )
                
                select case (SelectHydraulicConductivityModel)
                    
                case (HCC_CONSTANT)
                     ! no additional parameters to read
                case (HCC_HILLEL) 
                  call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%rexponentHillel_HCC = DumR ! Linear factor 
                case(HCC_MUALEM)
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%Smin_SWRC = DumR ! Minimum degree of saturation   
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%Smax_SWRC = DumR ! Maximum degree of saturation                
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%P0_SWRC = DumR ! Reference Pressure
                    call ReadReal(FileUnit, BName, DumR)
                    MatParams(I)%Lambda_SWRC = DumR ! Parameter 
                case default
                
                  call GiveError('GOM-ERROR: $$HYDR_CONDUCTIVITY_CURVE is ' //trim(MatParams(I)%HydraulicConductivityCurve)// ' and is therefore not properly specified.')
                
                end select 
                MatParams(I)%InitialDegreeOfSaturation = 1.0 ! initially the soil is assumed to be fully saturated
                MatParams(I)%BishopsAlpha = 1.0 ! initially the soil is assumed to be fully saturated
            end if        
                    
          !---------------------------------------------------------------------------------------------
          ! END READING UNSATURATED MATERIAL PROPERTIES
          ! NOT YED IN GID
          !------------------------------------------------------------------------------------
                 
                 
                !------ DERIVED PARAMETERS ------ 
               
                ! calculate G = E/2(1+nu) of the solid (for LE, MC, ICH, SSMC; not applicable for HP)
                MatParams(I)%ShearModulus = MatParams(I)%YoungModulus / 2 / ( 1 + MatParams(I)%PoissonRatio )
                
                ! calculate dry weight, gamma_dry = DryDensity*g/1000
                MatParams(I)%DryWeight =  ( 1 - MatParams(I)%InitialPorosity ) * MatParams(I)%DensitySolid * CalParams%GravityData%GAccel / 1000.0 
                
                ! calculate weight of the mixture, gamma_mix = DensityMixture*g/1000
                MatParams(I)%WeightMixture = MatParams(I)%DensityMixture * CalParams%GravityData%GAccel / 1000.0 
                
                ! calculate weight of the liquid, gamma_liquid = DensityLiquid*g/1000
                MatParams(I)%WeightLiquid = MatParams(I)%DensityLiquid * CalParams%GravityData%GAccel / 1000.0 
                
                ! calculate weight of the gas, gamma_gas = DensityGas*g/1000
                MatParams(I)%WeightGas = MatParams(I)%DensityGas * CalParams%GravityData%GAccel / 1000.0 
                
                ! calculate density of the submerged mixture ( = density_mixture - density_liquid ) and weight
                if (CalParams%ApplySubmergedCalculation) then ! only if submerged calculation
                  MatParams(I)%DensitySubmerged = MatParams(I)%DensityMixture - MatParams(I)%DensityLiquid
                  MatParams(I)%WeightSubmerged = MatParams(I)%DensitySubmerged * CalParams%GravityData%GAccel / 1000.0
                  MatParams(I)%WeightMixture = MatParams(I)%WeightSubmerged
                end if

              end do ! loop over number of materials

            else if (trim(BName)=='$$FINISH') then
                
              EXIT
              
            end if 

          end do ! loop for reading GOM-file  
                   
          ! determine global number of phases
          CalParams%NumberOfPhases = maxval(NumberOfPhases)
          
          if (NFORMULATION==2) then
            CalParams%NumberOfPhases = 2
          end if
            
        end subroutine ReadMaterial_v2021        
        
        
        subroutine ReadInteger(unit, name, value)
        !**********************************************************************
        !
        ! Function : Read integer from GOM file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(len=255), intent(out) :: name
          integer(INTEGER_TYPE), intent(out) :: value
        
          integer(INTEGER_TYPE) :: ios
          character(len=21) :: messageIOS = 'GOM file: Can''t read '

          read(unit, *, iostat=ios) name 
          call Assert( ios == 0, messageIOS//trim(name) )
          read(unit, *, iostat=ios) value
          call Assert( ios == 0, messageIOS//trim(name) )

        end subroutine ReadInteger

        
        subroutine ReadString(unit, name, value)
        !**********************************************************************
        !
        ! Function : Read string from GOM file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(len=255), intent(out) :: name
          character(*), intent(out) :: value
        
          integer(INTEGER_TYPE) :: ios
          character(len=21) :: messageIOS = 'GOM file: Can''t read '

          read(unit, *, iostat=ios) name 
          call Assert( ios == 0, messageIOS//trim(name) )
          read(unit, *, iostat=ios) value
          call Assert( ios == 0, messageIOS//trim(name) )

        end subroutine ReadString

        
        subroutine ReadReal(unit, name, value)
        !**********************************************************************
        !
        ! Function : Read real from GOM file
        !
        !**********************************************************************
        
        implicit none
        
          integer(INTEGER_TYPE), intent(in) :: unit
          character(len=255), intent(out) :: name
          real(REAL_TYPE), intent(out) :: value
        
          integer(INTEGER_TYPE) :: ios
          character(len=21) :: messageIOS = 'GOM file: Can''t read '

          read(unit, *, iostat=ios) name 
          call Assert( ios == 0, messageIOS//trim(name) )
          read(unit, *, iostat=ios) value
          call Assert( ios == 0, messageIOS//trim(name) )

        end subroutine ReadReal

        
        subroutine InitialiseMaterialParameters()
        !**********************************************************************
        !
        ! Function : Allocates the array(s) belonging to this module storing material parameters 
        !
        !**********************************************************************
        
        implicit none
          
          ! local variables
          integer(INTEGER_TYPE) :: IError
        
          if (allocated(MatParams)) then
            deallocate(MatParams, stat = IError)
          end if
          allocate(MatParams(CalParams%NumberOfMaterials), stat = IError)

        end subroutine InitialiseMaterialParameters
 
       
        subroutine DestroyMaterialParameters()
        !**********************************************************************
        !
        ! Function : Deallocates the array(s) belonging to this module storing material parameters 
        !
        !**********************************************************************
        
        implicit none
          
          ! local variables
          integer(INTEGER_TYPE) :: IError
        
          if (allocated(MatParams)) then
            deallocate(MatParams, stat = IError)
          end if
          
        end subroutine DestroyMaterialParameters
        
        subroutine Error_NoDLL(name)
        !**********************************************************************
        !
        ! Function:  Error function
        !            If the required DLL "name" is not available an error message
        !            is printed on the screen and the program execution is stopped
        !
        !**********************************************************************
        
        implicit none

            character(*), intent(in) :: name

            call GiveError('UMAT-ERROR: DLL could not be loaded, DLL: ' // trim(name))

        end subroutine Error_NoDLL


        subroutine Error_NoUMAT(name)
        !**********************************************************************
        !
        ! Function:  Error function
        !            If the required UMAT is not available in the DLL an error message
        !            is printed on the screen and the program execution is stopped
        !
        !**********************************************************************
        
        implicit none
          
            character(*), intent(in) :: name

            call GiveError('UMAT-ERROR: UMAT not found in: ' // trim(name))

        end subroutine Error_NoUMAT

   
      end module ModReadMaterialData
