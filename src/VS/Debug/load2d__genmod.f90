        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 12 09:22:37 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOAD2D__genmod
          INTERFACE 
            SUBROUTINE LOAD2D(RLOAD,NDOF,COORD,NUMBEROFINTEGRATIONPOINTS&
     &,ILOADCON,LOADVALUE,LOADTYPE)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: NUMBEROFINTEGRATIONPOINTS
              REAL(KIND=8), INTENT(INOUT) :: RLOAD(COUNTERS%NODTOT*     &
     &NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: NDOF(COUNTERS%NODTOT)
              REAL(KIND=8), INTENT(IN) :: COORD(COUNTERS%NODTOT,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: ILOADCON(                  &
     &ELEMENTBOUNDARYNODES)
              REAL(KIND=8), INTENT(IN) :: LOADVALUE(ELEMENTBOUNDARYNODES&
     &,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: LOADTYPE
            END SUBROUTINE LOAD2D
          END INTERFACE 
        END MODULE LOAD2D__genmod
