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
