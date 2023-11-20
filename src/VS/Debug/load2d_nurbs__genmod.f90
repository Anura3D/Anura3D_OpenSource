        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 17 18:18:14 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOAD2D_NURBS__genmod
          INTERFACE 
            SUBROUTINE LOAD2D_NURBS(RLOAD,NDOF,COORD,NINT,NNOD,ILOADCON,&
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
            END SUBROUTINE LOAD2D_NURBS
          END INTERFACE 
        END MODULE LOAD2D_NURBS__genmod
