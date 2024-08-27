*USER SUBROUTINES
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


        ! ARGUMENTS:
        !          I/O  TYPE
        !  PROPS    I   R()  : LIST WITH MODEL PARAMETERS
        !  DSTRAN   I   R()  : STRAIN INCREMENT
        !  DDSDDE   O   R(,) : MATERIAL STIFFNESS MATRIX
        !  STRESS  I/O  R()  : STRESSES
        !  STATEV  I/O  R()  : STATE VARIABLES
 
        !---  LOCAL VARIABLES
        DIMENSION DSIG(NTENS), SIG(NTENS), STRAINRATE(NTENS),
     1  DEVSTRAINRATE(NTENS), SIGMA(NTENS)
        
        LOGICAL COMPUTEPRESSURE
        double precision:: EEL, ENU, YIELDSTRESS, VISCOSITY, BULKMODULUSLIQUID, LIQUIDPRESSURECAVITATION, ONE, TWO, G
        double precision:: SHEARSTRESS, FAC, D1,D2
        double precision:: VOLSTRAINRATE, VOLSTRAINRATECOMPONENT, LIQUIDPRESSURE, LIQUIDPRESSUREINCREMENT


        EEL = PROPS(1)   
        ENU = PROPS(2)
        YIELDSTRESS=PROPS(3)
        VISCOSITY = PROPS(4)
        BULKMODULUSLIQUID=PROPS(5)
        LIQUIDPRESSURECAVITATION=PROPS(6)
        ONE = 1.0D0
        TWO = 2.0D0
        G = EEL/(TWO*(ONE + ENU))
      IF (NTENS == 6) THEN ! 3D  
       SHEARSTRESS =SQRT (( (STRESS(1) + 
     1                         -STRESS(2))**2 + 
     2                  (STRESS(1) + 
     3                      -STRESS(3))**2 + 
     4                  (STRESS(2)+ 
     5                      -STRESS(3))**2 ) / 3.0 + 
     6                ( STRESS(4)**2 + 
     7                  STRESS(5)**2 + 
     8                  STRESS(6)**2 ) * 2.0)
      ELSE IF (NTENS == 4) THEN ! 2D
          SHEARSTRESS =SQRT (( (STRESS(1) + 
     1           -STRESS(2))**2 + 
     2           (STRESS(1)+ 
     3           -STRESS(3))**2 + 
     4           (STRESS(2)+ 
     5           -STRESS(3))**2 ) / 3.0 + 
     6           ( STRESS(4)**2) * 2.0)
      END IF    
        
      IF   (SHEARSTRESS.LE.SQRT(2.0)*YIELDSTRESS) THEN
        ! CALCULATE ELASTIC STRESS INCREMENT (DSIGE = ELASTIC STIFFNESS D * STRAIN INCREMENT DEPS)
        FAC = TWO * G / ( ONE - TWO * ENU )
        D1 = FAC * ( ONE - ENU )
        D2 = FAC * ENU
        DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
        DSIG(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
        DSIG(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
        DSIG(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
        DSIG(4) = G * DSTRAN(4)
        IF (NTENS == 6) THEN
        DSIG(5) = G * DSTRAN(5)
        DSIG(6) = G * DSTRAN(6)
        END IF
        ! ELASTIC STRESS
        SIG = STRESS + DSIG
        
        ! STRESS STATE PARAMETERS UPDATE
        DO I = 1, NTENS
          STRESS(I) = SIG(I)
      END DO
        
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = D2
        DDSDDE(1,1) = D1
        DDSDDE(2,2) = D1
        DDSDDE(3,3) = D1
        DDSDDE(4,4) = G
        IF (NTENS == 6) THEN
          DDSDDE(5,5) = G
          DDSDDE(6,6) = G
        END IF
      
        RETURN
      ELSE
                  
       SIGMA = STRESS


                    
          DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
          STRAINRATE = DSTRAN / DTIME
          
          VOLSTRAINRATE = STRAINRATE(1) + STRAINRATE(2) + STRAINRATE(3)
          VOLSTRAINRATECOMPONENT = VOLSTRAINRATE / 3.0
          DEVSTRAINRATE(1:3) = STRAINRATE(1:3) - VOLSTRAINRATECOMPONENT
          DEVSTRAINRATE(4:NTENS) = STRAINRATE(4:NTENS)
          
          
          ! CONSIDER STRESS STATE OF PREVIOUS TIME STEP    
          LIQUIDPRESSURE = (STRESS(1)+STRESS(2)+STRESS(3))/3
          

     

          

          COMPUTEPRESSURE = (LIQUIDPRESSURE.LT.LIQUIDPRESSURECAVITATION)

          
          IF(COMPUTEPRESSURE) THEN
             
            LIQUIDPRESSUREINCREMENT = BULKMODULUSLIQUID  * DSTRANVOL

            LIQUIDPRESSURE = LIQUIDPRESSURE + LIQUIDPRESSUREINCREMENT

          
            IF(LIQUIDPRESSURE .GT. LIQUIDPRESSURECAVITATION) THEN
              LIQUIDPRESSURE=LIQUIDPRESSURECAVITATION
              STRESS(1:3) = LIQUIDPRESSURE
              STRESS(4:NTENS) = 0.0
              RETURN
            END IF  
            
            
     
            
          SIGMA(1:3) = LIQUIDPRESSURE
          SIGMA(4:NTENS) = 0.0
          

          SIGMA(1:3) = LIQUIDPRESSURE +  
     1                    TWO * VISCOSITY * DEVSTRAINRATE(1:3)
          SIGMA(4:NTENS) = VISCOSITY * DEVSTRAINRATE(4:NTENS)

          ELSE 
            LIQUIDPRESSURE=LIQUIDPRESSURECAVITATION
            STRESS(1:3) = LIQUIDPRESSURE
            STRESS(4:NTENS) = 0.0
            RETURN
          
          END IF
            
      STRESS = SIGMA
     

   
      END IF
      END SUBROUTINE UMAT            