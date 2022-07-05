        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  5 17:04:45 2022
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
