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


        ! Arguments:
        !          I/O  Type
        !  PROPS    I   R()  : List with model parameters
        !  DSTRAN   I   R()  : Strain increment
        !  DDSDDE   O   R(,) : Material stiffness matrix
        !  STRESS  I/O  R()  : stresses
        !  STATEV  I/O  R()  : state variables
 
        !---  Local variables
        dimension dSig(NTENS), Sig(NTENS)

        ! Contents of PROPS(2)
        !  1 : G       shear modulus
        !  2 : ENU     Poisson's ratio

        G = PROPS(1)   
        ENU = PROPS(2)
        one = 1.0d0
        two = 2.0d0 
        ! calculate elastic stress increment (DSigE = elastic stiffness D * strain increment DEps)
        FAC = two * G / ( one - two * ENU )
        D1 = FAC * ( one - ENU )
        D2 = FAC * ENU
        DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
        dSig(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
        dSig(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
        dSig(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
        dSig(4) = G * DSTRAN(4)
        if (NTENS == 6) then
        dSig(5) = G * DSTRAN(5)
        dSig(6) = G * DSTRAN(6)
        end if
        ! elastic stress
        Sig = STRESS + dSig
        
        ! stress state parameters update
        do i = 1, NTENS
          STRESS(i) = Sig(i)
      end do
        
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = D2
        DDSDDE(1,1) = D1
        DDSDDE(2,2) = D1
        DDSDDE(3,3) = D1
        DDSDDE(4,4) = G
        if (NTENS == 6) then
          DDSDDE(5,5) = G
          DDSDDE(6,6) = G
        end if
      
        return
      
      end subroutine umat            