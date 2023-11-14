        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 14 13:12:30 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOAD2D_QUAD__genmod
          INTERFACE 
            SUBROUTINE LOAD2D_QUAD(RLOAD,NDOF,COORD,NINT,NNOD,ILOADCON, &
     &ELEMENTUPONWHICHLOADISAPPLIED)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: NNOD
              INTEGER(KIND=4), INTENT(IN) :: NINT
              REAL(KIND=8), INTENT(INOUT) :: RLOAD(COUNTERS%NODTOT*NDOFL&
     &)
              INTEGER(KIND=4) :: NDOF(COUNTERS%NODTOT)
              REAL(KIND=8), INTENT(IN) :: COORD(COUNTERS%NODTOT,NVECTOR)
              INTEGER(KIND=4) :: ILOADCON(NNOD)
              INTEGER(KIND=4), INTENT(IN) ::                            &
     &ELEMENTUPONWHICHLOADISAPPLIED
            END SUBROUTINE LOAD2D_QUAD
          END INTERFACE 
        END MODULE LOAD2D_QUAD__genmod
