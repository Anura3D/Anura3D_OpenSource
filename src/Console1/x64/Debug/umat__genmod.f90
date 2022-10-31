        !COMPILER-GENERATED INTERFACE MODULE: Fri Oct  7 11:21:49 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UMAT__genmod
          INTERFACE 
            SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,&
     &DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,    &
     &CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,    &
     &CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
              INTEGER(KIND=4) :: NPROPS
              INTEGER(KIND=4) :: NSTATEV
              INTEGER(KIND=4) :: NTENS
              REAL(KIND=8) :: STRESS(NTENS)
              REAL(KIND=8) :: STATEV(NSTATEV)
              REAL(KIND=8) :: DDSDDE(NTENS,NTENS)
              REAL(KIND=8) :: SSE
              REAL(KIND=8) :: SPD
              REAL(KIND=8) :: SCD
              REAL(KIND=8) :: RPL
              REAL(KIND=8) :: DDSDDT(NTENS)
              REAL(KIND=8) :: DRPLDE(NTENS)
              REAL(KIND=8) :: DRPLDT
              REAL(KIND=8) :: STRAN(NTENS)
              REAL(KIND=8) :: DSTRAN(NTENS)
              REAL(KIND=8) :: TIME(2)
              REAL(KIND=8) :: DTIME
              REAL(KIND=8) :: TEMP
              REAL(KIND=8) :: DTEMP
              REAL(KIND=8) :: PREDEF(1)
              REAL(KIND=8) :: DPRED(1)
              CHARACTER(LEN=80) :: CMNAME
              INTEGER(KIND=4) :: NDI
              INTEGER(KIND=4) :: NSHR
              REAL(KIND=8) :: PROPS(NPROPS)
              REAL(KIND=8) :: COORDS(3)
              REAL(KIND=8) :: DROT(3,3)
              REAL(KIND=8) :: PNEWDT
              REAL(KIND=8) :: CELENT
              REAL(KIND=8) :: DFGRD0(3,3)
              REAL(KIND=8) :: DFGRD1(3,3)
              INTEGER(KIND=4) :: NOEL
              INTEGER(KIND=4) :: NPT
              INTEGER(KIND=4) :: LAYER
              INTEGER(KIND=4) :: KSPT
              INTEGER(KIND=4) :: KSTEP
              INTEGER(KIND=4) :: KINC
            END SUBROUTINE UMAT
          END INTERFACE 
        END MODULE UMAT__genmod
