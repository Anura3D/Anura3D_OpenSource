        !COMPILER-GENERATED INTERFACE MODULE: Sun Feb  9 10:25:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOAD2D__genmod
          INTERFACE 
            SUBROUTINE LOAD2D(RLOAD,NDOF,COORD,NUMBEROFINTEGRATIONPOINTS&
     &,ILOADCON,LOADVALUE,LOADTYPE)
              USE MODELEMENTEVALUATION
              REAL(KIND=8), INTENT(INOUT) :: RLOAD(COUNTERS%SUM_NODTOT* &
     &NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: NDOF(COUNTERS%SUM_NODTOT)
              REAL(KIND=8), INTENT(IN) :: COORD(COUNTERS%SUM_NODTOT,    &
     &NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: NUMBEROFINTEGRATIONPOINTS
              INTEGER(KIND=4), INTENT(IN) :: ILOADCON(                  &
     &ELEMENTBOUNDARYNODES)
              REAL(KIND=8), INTENT(IN) :: LOADVALUE(ELEMENTBOUNDARYNODES&
     &,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: LOADTYPE
            END SUBROUTINE LOAD2D
          END INTERFACE 
        END MODULE LOAD2D__genmod
