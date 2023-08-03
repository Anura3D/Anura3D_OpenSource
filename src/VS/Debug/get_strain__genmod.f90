        !COMPILER-GENERATED INTERFACE MODULE: Thu Aug  3 16:51:30 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_STRAIN__genmod
          INTERFACE 
            SUBROUTINE GET_STRAIN(IEL,IPOINT,ICON,B,DISP,NDOFEX,EPS)
              USE MODMPMDATA
              USE MODMESHINFO
              INTEGER(KIND=4), INTENT(IN) :: IEL
              INTEGER(KIND=4), INTENT(IN) :: IPOINT
              INTEGER(KIND=4), INTENT(IN) :: ICON(ELEMENTNODES,COUNTERS%&
     &NEL)
              REAL(KIND=8), INTENT(IN) :: B(NVECTOR,ELEMENTNODES)
              REAL(KIND=8), INTENT(IN) :: DISP(COUNTERS%N)
              INTEGER(KIND=4), INTENT(IN) :: NDOFEX(COUNTERS%NODTOT+1)
              REAL(KIND=8), INTENT(OUT) :: EPS(NTENSOR)
            END SUBROUTINE GET_STRAIN
          END INTERFACE 
        END MODULE GET_STRAIN__genmod
