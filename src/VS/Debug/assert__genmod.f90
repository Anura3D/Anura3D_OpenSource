        !COMPILER-GENERATED INTERFACE MODULE: Mon Oct 16 14:28:55 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ASSERT__genmod
          INTERFACE 
            SUBROUTINE ASSERT(CONDITION,MESSAGE)
              LOGICAL(KIND=4), INTENT(IN) :: CONDITION
              CHARACTER(*), INTENT(IN) :: MESSAGE
            END SUBROUTINE ASSERT
          END INTERFACE 
        END MODULE ASSERT__genmod
