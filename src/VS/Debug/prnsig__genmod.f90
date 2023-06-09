        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun  9 17:35:24 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PRNSIG__genmod
          INTERFACE 
            SUBROUTINE PRNSIG(IOPT,S,XN1,XN2,XN3,S1,S2,S3,P,Q)
              USE MODGLOBALCONSTANTS
              INTEGER(KIND=4), INTENT(IN) :: IOPT
              REAL(KIND=8), INTENT(IN) :: S(NTENSOR)
              REAL(KIND=8), INTENT(OUT) :: XN1(NPRINCIPAL)
              REAL(KIND=8), INTENT(OUT) :: XN2(NPRINCIPAL)
              REAL(KIND=8), INTENT(OUT) :: XN3(NPRINCIPAL)
              REAL(KIND=8), INTENT(OUT) :: S1
              REAL(KIND=8), INTENT(OUT) :: S2
              REAL(KIND=8), INTENT(OUT) :: S3
              REAL(KIND=8), INTENT(OUT) :: P
              REAL(KIND=8), INTENT(OUT) :: Q
            END SUBROUTINE PRNSIG
          END INTERFACE 
        END MODULE PRNSIG__genmod
