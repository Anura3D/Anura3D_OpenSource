        !COMPILER-GENERATED INTERFACE MODULE: Fri May 24 09:30:13 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FORMB3__genmod
          INTERFACE 
            SUBROUTINE FORMB3(INT,IEL,ICON,CO,B,DET,WTN,                &
     &DSHAPEVALUESARRAY,MAXPARTICLE,MAXEL,GETPARTICLEINDEX,IPATCH)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: IPATCH
              INTEGER(KIND=4), INTENT(IN) :: MAXEL
              INTEGER(KIND=4), INTENT(IN) :: INT
              INTEGER(KIND=4) :: IEL
              INTEGER(KIND=4) :: ICON(ELEMENTNODES,NEL_NURBS((IPATCH)))
              REAL(KIND=8) :: CO(MAXIMUM_NCONTROLPOINTS,NDIM)
              REAL(KIND=8) :: B(NDIM,ELEMENTNODES)
              REAL(KIND=8) :: DET
              REAL(KIND=8) :: WTN
              REAL(KIND=8), INTENT(IN) :: DSHAPEVALUESARRAY(ELEMENTNODES&
     &,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: MAXPARTICLE
              INTEGER(KIND=4), INTENT(IN) :: GETPARTICLEINDEX(COUNTERS% &
     &NPARTICLES,MAXEL)
            END SUBROUTINE FORMB3
          END INTERFACE 
        END MODULE FORMB3__genmod
