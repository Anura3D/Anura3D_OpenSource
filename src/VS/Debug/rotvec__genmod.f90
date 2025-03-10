        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb 27 18:39:45 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROTVEC__genmod
          INTERFACE 
            SUBROUTINE ROTVEC(IROT,NROTNODESLOC,RMAT,NDOFEX,V_IN,V_OUT)
              INTEGER(KIND=4) :: IROT(COUNTERS%SUM_NODTOT)
              INTEGER(KIND=4) :: NROTNODESLOC
              REAL(KIND=8) :: RMAT(3,3,*)
              INTEGER(KIND=4) :: NDOFEX(COUNTERS%SUM_NODTOT)
              REAL(KIND=8) :: V_IN(COUNTERS%N)
              REAL(KIND=8) :: V_OUT(COUNTERS%N)
            END SUBROUTINE ROTVEC
          END INTERFACE 
        END MODULE ROTVEC__genmod
