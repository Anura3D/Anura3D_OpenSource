        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 11 23:33:55 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FORMB3__genmod
          INTERFACE 
            SUBROUTINE FORMB3(INT,IEL,ICON,CO,B,WTN,DSHAPEVALUESARRAY,  &
     &MAXPARTICLE,MAXEL,IPATCH)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: IPATCH
              INTEGER(KIND=4), INTENT(IN) :: INT
              INTEGER(KIND=4) :: IEL
              INTEGER(KIND=4) :: ICON(ELEMENTNODES,COUNTERS%NEL((IPATCH)&
     &))
              REAL(KIND=8) :: CO(NURBS%MAXIMUM_NCONTROLPOINTS,NDIM,     &
     &COUNTERS%NPATCHES)
              REAL(KIND=8) :: B(NDIM,ELEMENTNODES)
              REAL(KIND=8) :: WTN
              REAL(KIND=8), INTENT(IN) :: DSHAPEVALUESARRAY(ELEMENTNODES&
     &,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: MAXPARTICLE
              INTEGER(KIND=4), INTENT(IN) :: MAXEL
            END SUBROUTINE FORMB3
          END INTERFACE 
        END MODULE FORMB3__genmod
