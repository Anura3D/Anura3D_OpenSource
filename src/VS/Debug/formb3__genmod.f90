        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 13 00:35:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FORMB3__genmod
          INTERFACE 
            SUBROUTINE FORMB3(INT,IEL,ICON,CO,B,DET,WTN)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4) :: INT
              INTEGER(KIND=4) :: IEL
              INTEGER(KIND=4) :: ICON(ELEMENTNODES,COUNTERS%NEL)
              REAL(KIND=8) :: CO(COUNTERS%NODTOT,NDIM)
              REAL(KIND=8) :: B(NDIM,ELEMENTNODES)
              REAL(KIND=8) :: DET
              REAL(KIND=8) :: WTN
            END SUBROUTINE FORMB3
          END INTERFACE 
        END MODULE FORMB3__genmod
