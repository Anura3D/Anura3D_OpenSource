        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 14 10:11:31 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATMAT__genmod
          INTERFACE 
            SUBROUTINE MATMAT(XMAT1,ID1,XMAT2,ID2,NR1,NC2,NC1,XMATR,IDR)
              INTEGER(KIND=4) :: IDR
              INTEGER(KIND=4) :: ID2
              INTEGER(KIND=4) :: ID1
              REAL(KIND=8) :: XMAT1(ID1,*)
              REAL(KIND=8) :: XMAT2(ID2,*)
              INTEGER(KIND=4) :: NR1
              INTEGER(KIND=4) :: NC2
              INTEGER(KIND=4) :: NC1
              REAL(KIND=8) :: XMATR(IDR,*)
            END SUBROUTINE MATMAT
          END INTERFACE 
        END MODULE MATMAT__genmod
