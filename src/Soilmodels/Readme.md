# How to integrate External Soil Models (ESM) into Anura3D

1. Start by uncommenting line 870 in `ReadMaterialData.for`. Use a name that best describes your material constitutive model. Example:

    ```fortran
    if (MatParams(I)%SoilModelName=="CamClay")  MatParams(I)%ESM_POINTER => ESM_CamClay
    ```
2. Load the constitutive model into the ```Soilmodels``` folder:            
    i. Copy the files in the Folder     
    ii. In Visual Studio right click on the ```MaterialModels``` folder         
    iii. Click on ```Add```⟶```Existing item```        
    iv. Search the UMAT file and click ```Add```

3. Add the module structure as noted below:
    ```fortran
    module ModNameOfYourModel
    contains

    !Your code here

    end module
    ```

    
4. Add the ```ESM``` subroutine interface to your model if it doesn't have it. You can copy it from an existing model file such as the ```A3DLinearElasticity.f``` file. *Note: The ```ESM``` should be named consistently with the name given in step 1. e.g., ```ESM_CamClay```*

5. Add your module as a dependency on the ```ExternalSoilModel.for``` and the ```ReadMaterialData.for``` modules as follows:
    ```fortran
    use ModLinearElasticity
    use ModMohrCoulomb
    use ModBingham
    use ModNameOfYourModel !<-- This is your model
    ```
6. Compile, debug, and verify

# How to create a custom GUI for your constitutive model

1. Locate the file ```Materials.xml``` inside the ```xml``` folder of the Problemtype.

2. Locate the ```Solid Material Model``` nodule and add the name of your model as shown:
    ```xml
    <value n="_material_model_solid_" pn="Solid material model" values="Rigid body,Linear Elasticity,Mohr-Coulomb,YourModelName,External Material Model" v="Rigid body" state="[hide_show_const_model_mat_type_n Liquid %W]" actualize_tree="1"/>
    ```

3. Add a nodule for each material property needed in the same order used to define the variable ```PROPS()``` in the ```UMAT```. Example:
    ```xml
    <value n="mat_property_1" pn="My property 1 label" v="0.0" state="[hide_show_const_model {YourModelName} %W]"/>

    <value n="mat_property_2" pn="My property 2 label" v="0.0" state="[hide_show_const_model {YourModelName} %W]"/>
    ```

4. Locate the file ```createGOM.tcl``` inside the folder ```scripts``` of the ProblemType.

5. Add the following condition to the type of materials so that selecting your model also writes ```external_soil_model```:     
    ```tcl
        } elseif {$typemodel == "External Material Model"} {
            set model "external_soil_model"
        } else {                            #<-- New line
            set model "external_soil_model" #<-- New line
        }
    ```
5. Search for ```"Mohr-Coulomb"``` in the code and find a place to put the new lines of code that will write the input for Anura as illustrated below:

    ```tcl
    GiD_WriteCalculationFile puts $type
    } elseif {$typemodel == "YourModelName"} {
        GiD_WriteCalculationFile puts {$$MATERIAL_MODEL_NAME} 
        GiD_WriteCalculationFile puts {YourModelName_as_in_Anura3D}
        GiD_WriteCalculationFile puts {$$UMAT_DIMENSION} 
        GiD_WriteCalculationFile puts {3D} 
        #for 3D or
        # GiD_WriteCalculationFile puts {2D_plane_strain} 

        # The rest of your code will follow here

        
    } elseif {$typemodel == "External Material Model"} 
    ```

6. Add each of your properties or state variables as shown below (i.e., the rest of your code):
    ```tcl
    # Material parameters: 
    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="mat_property_1"]}]]
    set type [$node getAttribute "v"]
    GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_01} 
    GiD_WriteCalculationFile puts $type
    set node [$gNode selectNodes [format_xpath {container[@n="_material_constitutive_model"]/value[@n="mat_property_2"]}]]
    set type [$node getAttribute "v"]
    GiD_WriteCalculationFile puts {$$MATERIAL_PARAMETER_SOLID_02} 
    GiD_WriteCalculationFile puts $type
    ```

7. Complete up to 50 mat properties with zero value adding the following lines below:
    ```tcl
    #complete the additional 50 parameters with 0.0
    set var "\$\$MATERIAL_PARAMETER_SOLID"  
    for { set j nvar+1 } { $j <= 50 } { incr j } {
        if { $j<= 9 } {
            set jj "0$j" 
        } else {
            set jj "$j"
        }  
        set nn ${var}_${jj}
        set vv 0.0
        GiD_WriteCalculationFile puts $nn
        GiD_WriteCalculationFile puts $vv                        
    }
    ```
    Replace ```nvar+1``` with the total number of variables you have plus one.

8. Add the following block of code for the state variables. A total of 50 state variables is required.
    ```tcl
    # Complete the 50 state variables with 0.0
    set var "\$\$INITIAL_STATE_VARIABLE_SOLID"
    for { set j 1 } { $j <= 50 } { incr j } {
        if { $j<= 9 } {
            set jj "0$j" 
        } else {
            set jj "$j"
        }  
        set nn ${var}_${jj}
        set vv 0.0
        GiD_WriteCalculationFile puts $nn
        GiD_WriteCalculationFile puts $vv                        
    }
    ```
9. Replace the ProblemType back into your GID ProblemType folder. Test and debug as needed.

