        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 21 23:58:55 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BUILDDELASTICINVERSE__genmod
          INTERFACE 
            SUBROUTINE BUILDDELASTICINVERSE(G,XNU,DI)
              USE MODGLOBALCONSTANTS
              REAL(KIND=8), INTENT(IN) :: G
              REAL(KIND=8), INTENT(IN) :: XNU
              REAL(KIND=8) :: DI(NTENSOR,NTENSOR)
            END SUBROUTINE BUILDDELASTICINVERSE
          END INTERFACE 
        END MODULE BUILDDELASTICINVERSE__genmod
