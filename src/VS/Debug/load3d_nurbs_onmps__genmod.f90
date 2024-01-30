        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 22 17:07:07 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOAD3D_NURBS_ONMPS__genmod
          INTERFACE 
            SUBROUTINE LOAD3D_NURBS_ONMPS(RLOAD,NDOF,COORD,NINT,NNOD,   &
     &ILOADCON,ELEMENTUPONWHICHLOADISAPPLIED)
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
            END SUBROUTINE LOAD3D_NURBS_ONMPS
          END INTERFACE 
        END MODULE LOAD3D_NURBS_ONMPS__genmod
