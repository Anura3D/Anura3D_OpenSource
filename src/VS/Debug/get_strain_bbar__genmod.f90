        !COMPILER-GENERATED INTERFACE MODULE: Tue Feb 18 14:29:14 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_STRAIN_BBAR__genmod
          INTERFACE 
            SUBROUTINE GET_STRAIN_BBAR(IEL,IPOINT,ICON,B,BBAR,DISP,     &
     &NDOFEX,EPS,IPATCH)
              USE MODMPMDATA
              USE MODMESHINFO
              INTEGER(KIND=4), INTENT(IN) :: IPATCH
              INTEGER(KIND=4), INTENT(IN) :: IEL
              INTEGER(KIND=4), INTENT(IN) :: IPOINT
              INTEGER(KIND=4), INTENT(IN) :: ICON(ELEMENTNODES,COUNTERS%&
     &NEL((IPATCH)))
              REAL(KIND=8), INTENT(IN) :: B(NVECTOR,ELEMENTNODES)
              REAL(KIND=8), INTENT(IN) :: BBAR(NVECTOR,ELEMENTNODES)
              REAL(KIND=8), INTENT(IN) :: DISP(COUNTERS%N)
              INTEGER(KIND=4), INTENT(IN) :: NDOFEX(                    &
     &NUMBEROFGLOBALCONTROLPOINTSUNIQUEMULTIPATCH+1)
              REAL(KIND=8), INTENT(OUT) :: EPS(NTENSOR)
            END SUBROUTINE GET_STRAIN_BBAR
          END INTERFACE 
        END MODULE GET_STRAIN_BBAR__genmod
