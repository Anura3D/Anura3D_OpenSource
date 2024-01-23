        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 23 14:16:40 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RJACINV__genmod
          INTERFACE 
            SUBROUTINE RJACINV(IDIMJ,RJAC,RJAC1,DET,DET1)
              INTEGER(KIND=4), INTENT(IN) :: IDIMJ
              REAL(KIND=8), INTENT(IN) :: RJAC(IDIMJ,IDIMJ)
              REAL(KIND=8), INTENT(OUT) :: RJAC1(IDIMJ,IDIMJ)
              REAL(KIND=8), INTENT(OUT) :: DET
              REAL(KIND=8), INTENT(OUT) :: DET1
            END SUBROUTINE RJACINV
          END INTERFACE 
        END MODULE RJACINV__genmod
