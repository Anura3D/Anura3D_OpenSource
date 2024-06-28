    module ModBingham
    contains 
    
    Subroutine ESM_BINGHAM(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, &
      DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)


      implicit double precision (a-h, o-z) 
      integer :: NTENS, NSTATEV, NADDVAR, NPROPS, NPT, NOEL, IDSET, NUMBEROFPHASES
      double precision :: EUNLOADING, PLASTICMULTIPLIER
      CHARACTER*80 CMNAME   
      DIMENSION STRESS(NTENS), DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      

!---Local variables required in standard UMAT
        integer :: IStep, TimeStep
        double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
        double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
        double precision, dimension(:), allocatable :: stran
        double precision, dimension(:), allocatable :: time
        double precision, dimension(:), allocatable :: predef
        double precision, dimension(:), allocatable :: dpred    
        double precision, dimension(:), allocatable :: coords
        double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
        double precision, dimension(:,:), allocatable :: drot
        double precision, dimension(:,:), allocatable :: dfgrd0
        double precision, dimension(:,:), allocatable :: dfgrd1
        double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
        double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
        double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
        double precision :: pnewdt, dtime, temp, dtemp, celent
        double precision :: Value ! auxiliary variable holding any real valued number
        double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  
              
  
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

!---Local variables defned by the user
! e.g. integer :: var_local	  
!---User can define here additional variables     
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
!Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0
          
!Rename additional variables
        Porosity = AdditionalVar(1)
        WaterPressure = AdditionalVar(2)
        WaterPressure0 = AdditionalVar(3)
        GasPressure = AdditionalVar(4)
        GasPressure0 = AdditionalVar(5)
        DegreeSaturation = AdditionalVar(6)
        time(1) = AdditionalVar(7)   !TotalRealTime
        time(2) = AdditionalVar(8)   !OverallTotalTime
        dtime = AdditionalVar(9)     !TimeIncrement
        IStep = AdditionalVar(10)    
        TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1   
!Call the UMAT
          call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, & 
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = props(5)
!---Always define this value to run the simulation
    
! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

    end subroutine ESM_BINGHAM
    
     SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)


      implicit double precision (a-h, o-z) 

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV), &
      DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS), &
      STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1), &
      PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


        ! ARGUMENTS:
        !          I/O  TYPE
        !  PROPS    I   R()  : LIST WITH MODEL PARAMETERS
        !  DSTRAN   I   R()  : STRAIN INCREMENT
        !  DDSDDE   O   R(,) : MATERIAL STIFFNESS MATRIX
        !  STRESS  I/O  R()  : STRESSES
        !  STATEV  I/O  R()  : STATE VARIABLES
 
        !---  LOCAL VARIABLES
        DIMENSION DSIG(NTENS), SIG(NTENS), STRAINRATE(NTENS), &
       DEVSTRAINRATE(NTENS), SIGMA(NTENS)
        
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
       SHEARSTRESS =SQRT (( (STRESS(1) + &
                              -STRESS(2))**2 + &
                       (STRESS(1) + &
                           -STRESS(3))**2 + &
                       (STRESS(2)+ &
                           -STRESS(3))**2 ) / 3.0 + &
                     ( STRESS(4)**2 + &
                       STRESS(5)**2 + &
                       STRESS(6)**2 ) * 2.0)
      ELSE IF (NTENS == 4) THEN ! 2D
          SHEARSTRESS =SQRT (( (STRESS(1) + &
                -STRESS(2))**2 + &
                (STRESS(1)+ &
                -STRESS(3))**2 + &
                (STRESS(2)+ &
                -STRESS(3))**2 ) / 3.0 + &
                ( STRESS(4)**2) * 2.0)
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
          

          SIGMA(1:3) = LIQUIDPRESSURE +  &
                         TWO * VISCOSITY * DEVSTRAINRATE(1:3)
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
     
    end module ModBingham 