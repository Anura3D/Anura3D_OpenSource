        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun  9 15:47:33 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROTVEC__genmod
          INTERFACE 
            SUBROUTINE ROTVEC(IROT,NROTNODESLOC,RMAT,NDOFEX,V_IN,V_OUT)
              INTEGER(KIND=4) :: IROT(COUNTERS%NODTOT)
              INTEGER(KIND=4) :: NROTNODESLOC
              REAL(KIND=8) :: RMAT(3,3,*)
              INTEGER(KIND=4) :: NDOFEX(COUNTERS%NODTOT)
              REAL(KIND=8) :: V_IN(COUNTERS%N)
              REAL(KIND=8) :: V_OUT(COUNTERS%N)
            END SUBROUTINE ROTVEC
          END INTERFACE 
        END MODULE ROTVEC__genmod
