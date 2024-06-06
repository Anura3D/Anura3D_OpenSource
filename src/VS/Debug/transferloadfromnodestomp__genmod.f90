        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun  5 15:45:18 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANSFERLOADFROMNODESTOMP__genmod
          INTERFACE 
            SUBROUTINE TRANSFERLOADFROMNODESTOMP(RLOAD,NDOF,            &
     &ELEMENTUPONWHICHLOADISAPPLIED)
              USE MODREADCALCULATIONDATA
              REAL(KIND=8), INTENT(IN) :: RLOAD(COUNTERS%N)
              INTEGER(KIND=4), INTENT(IN) :: NDOF(COUNTERS%NODTOT)
              INTEGER(KIND=4), INTENT(IN) ::                            &
     &ELEMENTUPONWHICHLOADISAPPLIED
            END SUBROUTINE TRANSFERLOADFROMNODESTOMP
          END INTERFACE 
        END MODULE TRANSFERLOADFROMNODESTOMP__genmod
