        !COMPILER-GENERATED INTERFACE MODULE: Wed Nov 16 15:15:04 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATVEC__genmod
          INTERFACE 
            SUBROUTINE MATVEC(XMAT,IM,VEC,N,VECR)
              INTEGER(KIND=4) :: IM
              REAL(KIND=8) :: XMAT(IM,*)
              REAL(KIND=8) :: VEC(*)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: VECR(*)
            END SUBROUTINE MATVEC
          END INTERFACE 
        END MODULE MATVEC__genmod
