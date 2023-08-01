        !COMPILER-GENERATED INTERFACE MODULE: Tue Aug  1 10:37:06 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NORM__genmod
          INTERFACE 
            SUBROUTINE NORM(NNODES,NDOF,V1,NDOFEX,VNORM)
              INTEGER(KIND=4), INTENT(IN) :: NDOF
              INTEGER(KIND=4), INTENT(IN) :: NNODES
              REAL(KIND=8), INTENT(IN) :: V1(NNODES*NDOF)
              INTEGER(KIND=4), INTENT(IN) :: NDOFEX(NNODES)
              REAL(KIND=8), INTENT(OUT) :: VNORM
            END SUBROUTINE NORM
          END INTERFACE 
        END MODULE NORM__genmod
