      program Main_CreateCPS
      
      use ModReadCalculationData
      use ModGlobalConstants
      use ModFileIO
      
      implicit none
      
        ! local variables
        character(len=MAX_FILENAME_LENGTH) :: CPSTemplateExternal
        character(len=MAX_FILENAME_LENGTH) :: CPSTemplateInternal
        integer :: FileUnit
        
        !InitialiseCalculationParameters()
        
        CPSTemplateExternal = 'template_v2017.2.CPS_001'
        CPSTemplateInternal = 'template_v2017.2i.CPS_001'
      
        FileUnit = CPSUnit
        
        ! Note that only the templates for the latest version are updated. The previous CPS files should not be changed after release date.
        
        ! write CPS-template for external release version
        call FileOpen(FileUnit, CPSTemplateExternal) 
        CalParams%IsExternalCPS = .true. ! .true. for external release version 
        call WriteCPS_v2017_2(FileUnit)
        close(FileUnit)
       
        ! write CPS-template for internal release/debug version
        call FileOpen(FileUnit, CPSTemplateInternal) 
        CalParams%IsExternalCPS = .false. ! .false. for internal release/debug version
        call WriteCPS_v2017_2(FileUnit)
        close(FileUnit) 
      
      end program Main_CreateCPS