        !COMPILER-GENERATED INTERFACE MODULE: Sun Aug  6 12:35:23 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FORMB3__genmod
          INTERFACE 
            SUBROUTINE FORMB3(INT,IEL,ICON,CO,B,DET,WTN,                &
     &DSHAPEVALUESARRAY,MAXPARTICLE,MAXEL,GETPARTICLEINDEX)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: MAXEL
              INTEGER(KIND=4), INTENT(IN) :: MAXPARTICLE
              INTEGER(KIND=4), INTENT(IN) :: INT
              INTEGER(KIND=4) :: IEL
              INTEGER(KIND=4) :: ICON(ELEMENTNODES,COUNTERS%NEL)
              REAL(KIND=8) :: CO(COUNTERS%NODTOT,NDIM)
              REAL(KIND=8) :: B(NDIM,ELEMENTNODES)
              REAL(KIND=8) :: DET
              REAL(KIND=8) :: WTN
              REAL(KIND=8), INTENT(IN) :: DSHAPEVALUESARRAY(COUNTERS%   &
     &NPARTICLES,ELEMENTNODES,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: GETPARTICLEINDEX(          &
     &MAXPARTICLE,MAXEL)
            END SUBROUTINE FORMB3
          END INTERFACE 
        END MODULE FORMB3__genmod
