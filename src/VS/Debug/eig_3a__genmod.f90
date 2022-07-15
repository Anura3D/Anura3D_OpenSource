        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 15 10:53:44 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EIG_3A__genmod
          INTERFACE 
            SUBROUTINE EIG_3A(IOPT,ST,S1,S2,S3,P,Q)
              USE MODGLOBALCONSTANTS
              INTEGER(KIND=4), INTENT(IN) :: IOPT
              REAL(KIND=8), INTENT(IN) :: ST(NTENSOR)
              REAL(KIND=8), INTENT(OUT) :: S1
              REAL(KIND=8), INTENT(OUT) :: S2
              REAL(KIND=8), INTENT(OUT) :: S3
              REAL(KIND=8), INTENT(OUT) :: P
              REAL(KIND=8), INTENT(OUT) :: Q
            END SUBROUTINE EIG_3A
          END INTERFACE 
        END MODULE EIG_3A__genmod
