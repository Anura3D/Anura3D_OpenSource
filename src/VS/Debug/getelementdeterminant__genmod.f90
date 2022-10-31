        !COMPILER-GENERATED INTERFACE MODULE: Mon Oct 31 15:12:14 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETELEMENTDETERMINANT__genmod
          INTERFACE 
            FUNCTION GETELEMENTDETERMINANT(IEL,ICON,CO,DSHAPEVALUESARRAY&
     &) RESULT(DET)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: IEL
              INTEGER(KIND=4), INTENT(IN) :: ICON(:,:)
              REAL(KIND=8), INTENT(IN) :: CO(:,:)
              REAL(KIND=8), INTENT(IN) :: DSHAPEVALUESARRAY(COUNTERS%   &
     &NPARTICLES,ELEMENTNODES,NVECTOR)
              REAL(KIND=8) :: DET
            END FUNCTION GETELEMENTDETERMINANT
          END INTERFACE 
        END MODULE GETELEMENTDETERMINANT__genmod
