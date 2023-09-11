        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 11 11:08:18 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALCULATEPRINCIPALSTRESSES__genmod
          INTERFACE 
            SUBROUTINE CALCULATEPRINCIPALSTRESSES(INTGLO,SIG,SIGPRIN)
              INTEGER(KIND=4), INTENT(IN) :: INTGLO
              REAL(KIND=8), INTENT(IN) :: SIG(6)
              REAL(KIND=8), INTENT(OUT) :: SIGPRIN(6)
            END SUBROUTINE CALCULATEPRINCIPALSTRESSES
          END INTERFACE 
        END MODULE CALCULATEPRINCIPALSTRESSES__genmod
