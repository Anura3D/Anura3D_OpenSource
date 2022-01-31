!*****************************************************************************  
!
!   This file incorporates work covered by the following copyright and 
!   permission notice:
!  
!   Copyright © 2017-2021 Bentley Systems, Incorporated. All rights reserved.
!
!   Permission is hereby granted, free of charge, to any person obtaining a copy 
!   of this software and associated documentation files (the "Software"), to deal 
!   in the Software without restriction, including without limitation the rights 
!   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
!   copies of the Software, and to permit persons to whom the Software is furnished 
!   to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be included in all 
!   copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
!   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
!   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
!   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!*****************************************************************************  
      
C***********************************************************************
      Subroutine MC_Tens( iArea, G, xNu, sPhi, sPsi, cCosPhi, sTens,
     *                    Prs_E, Prs, ipl
     *                    ,fme,fte,fhe,fm,ft,fh ,xLamM,xLamT
     *                   )
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
!
! Routine to solve stresses according to Coulomb / tension criterion
! Caution : compression positive
!
!     Sig = Sig^E - xLamM * DdGmdS - xLamT * DdGtdS
!     f_m = f_m^e - xLamM * a_mm - xLamT * a_mt
!     f_t = f_t^e - xLamM * a_tm - xLamT * a_tt
!
!        I/O Type
! iArea   I    I    : 1 : Triax compression corner     Sig1 > Sig2 = Sig3
!                     2 : Regular yield surface        Sig1 > Sig2 > Sig3
!                     3 : Triax extension              Sig1 = Sig2 > Sig3
! G       I    R    : Shear modulus
! xNu     I    R    : Poisson's ratio
! sPhi    I    R    : Sine of friction angle
! sPsi    I    R    : Sine of dilation angle
! cCosPhi I    R    : Cohesion * Cos(phi)
! sTens   I    R    : Allowable tensile stress (negative value)
! Prs_E   I    R(3) : Elastic principal stresses   Sig1E >= Sig2E >= Sig3E
! Prs     O    R(3) : Resulting principal stresses
! ipl     O    I    : Plasticity indicator
!                     0 : elastic
!                     1 : Coulomb surface
!                     2 : Tension surface (also corner with Coulomb)
! Local:
!   fme,fte,fm,ft : Value of yield functions m=Coulomb, t=tension, e=elastic
!   xLamM,xLamT   : Plastic multipliers
!   DdG#dS(3)     : Elastic D matrix times derivative of plastic potential
!   PrsE(3)       : Copy of Prs_E with possible correction for corners
!
      Dimension Prs_E(3),Prs(3)
      Dimension DdGmdS(3),DdGtdS(3),dPrs(3),PrsE(3)
      Do i=1,3
        PrsE(i) = Prs_E(i)
        Prs (i) = Prs_E(i)
      End Do
      ipl = 0
      xLamM = 0
      xLamT = 0
      ft=-1
      fm=-1
      Select Case (IArea)
        Case (1) ! Triax compression ! Sig1 > Sig2 = Sig3
          !!!!!!!!!!!!!!!!!!! COMPRESSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! first correct sig2=sig3
          PrsE(2) = (Prs_E(2) + Prs_E(3) ) / 2
          PrsE(3) = PrsE(2)
          Prs (2) = PrsE(2)
          Prs (3) = PrsE(3)

          fme = ( 2*PrsE(1) - PrsE(2) - PrsE(3) ) /4
     *         -( 2*PrsE(1) + PrsE(2) + PrsE(3) ) /4 * sPhi - cCosPhi
          fte = sTens - ( PrsE(2) + PrsE(3) ) /2
          fmTol = 1d-10*abs(fme)+1d-10
          ftTol = 1d-10*abs(fte)+1d-10
          If (fme.Lt.fmTol .And. fte.Lt.ftTol ) Goto 999 ! Return Both elastic

          ! Try Coulomb only
          DdGmdS(1) = G/(1-2*xNu)/2 * ( 2*(1-2*xNu)  -(2      )*sPsi )
          DdGmdS(2) = G/(1-2*xNu)/2 * ( - (1-2*xNu)  -(1+2*xNu)*sPsi )
          DdGmdS(3) = G/(1-2*xNu)/2 * ( - (1-2*xNu)  -(1+2*xNu)*sPsi )
          amm =G/4 * (3  - sPsi - sPhi + sPhi*sPsi*(3+2*xNu)/(1-2*xNu) )
          If (fme.Gt.fmtol) Then
            xLamM = fme / amm  ! assume only coulomb surface
            xLamT = 0
            ! Calculate stresses Prs = Prs - xLamM*DdGmdS
            Call AddVec( PrsE, DdGmdS, 1d0, -xLamM, 3, Prs )
            ! check yield functions
            fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *          -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens - ( Prs(2) + Prs(3) ) /2
            ! when correct return
            ipl = 1 ! Coulomb
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
!         coming here : tension surface active
          ipl = 2 ! tens
          DdGtdS(1) =  -G/(1-2*xNu) * ( 2* xNu  )
          DdGtdS(2) =  -G/(1-2*xNu) * (  1      )
          DdGtdS(3) =  -G/(1-2*xNu) * (  1      )
          atm =G/2 * ( 1 + sPsi*(1+2*xNu)/(1-2*xNu) )
          amt =G/2 * ( 1 + sPhi*(1+2*xNu)/(1-2*xNu) )
          att = G /(1-2*xNu)
          If (fte.Gt.fttol) Then
            xLamT = fte/att
            xLamM = 0
            ! Prs = PrsE - xLamT * DdGtdS
            Call AddVec( PrsE, DdGtdS, 1d0, -xLamT, 3, Prs )
            fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *          -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens - ( Prs(2) + Prs(3) ) /2
            ft1= sTens - Prs(1)
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Then ! Goto 999 ! Return
              If (ft1.Gt.0) Prs(1) = sTens
              Goto 999 ! Return
            End If
          End If
          ! combined: coulomb + tension
          ipl = 2 ! tens (+compr)
          Det = amm*att-amt*atm
          xLamM = ( att*fme - amt*fte ) / det
          xLamT = (-atm*fme + amm*fte ) / det
          ! Prs = PrsE - xLamM * DdGmdS - xLamT * DdGtdS
          Call AddVec( DdGmdS, DdGtdS, xLamM, xLamT, 3, dPrs )
          Call AddVec( PrsE, dPrs, 1d0, -1d0, 3, Prs )
          fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *        -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
          ft = sTens - ( Prs(2) + Prs(3) ) /2
          ft1= sTens - Prs(1)
          if (ft1.gt.1d-10) write(1,*)'FT1 positive !!!'
          If (ft1.Gt.0) Then
            Prs(1) = sTens
          End If
          If ( fm .Lt. fmTol .And.
     *         ft .Lt. ftTol       ) Goto 999 ! Return
        Case (2) ! Regular
          !!!!!!!!!!!!!!!!!!! REGULAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          fme = ( PrsE(1) - PrsE(3) ) /2
     *         -( PrsE(1) + PrsE(3) ) /2 * sPhi - cCosPhi
          fte = sTens - PrsE(3)
          fmTol = 1d-10*abs(fme)+1d-10
          ftTol = 1d-10*abs(fte)+1d-10
          If (fme.Lt.fmtol .And. fte.Lt.fttol ) Goto 999 ! Return ! Both elastic
          ! Try Coulomb only
          DdGmdS(1) = G/(1-2*xNu) * ( 1-2*xNu -     sPsi )
          DdGmdS(2) = G/(1-2*xNu) * (        -2*xNu*sPsi )
          DdGmdS(3) = G/(1-2*xNu) * (-1+2*xNu -     sPsi )
          amm = G * ( 1 + sPhi*sPsi/(1-2*xNu) )
          If (fme.Gt.fmtol) Then
            xLamM = fme / amm
            xLamT = 0
            ! Prs = PrsE - xLamM * DdGmdS
            Call AddVec( PrsE, DdGmdS, 1d0, -xLamM, 3, Prs )
            fm = ( Prs(1) - Prs(3) ) /2
     *          -( Prs(1) + Prs(3) ) /2 * sPhi - cCosPhi
            ft = sTens - Prs(3)
            ipl = 1 ! Coulomb
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
          ! tension surface
          ipl = 2 ! tens
          DdGtdS(1) = 2*G/(1-2*xNu) * ( -  xNu   )
          DdGtdS(2) = 2*G/(1-2*xNu) * ( -  xNu   )
          DdGtdS(3) = 2*G/(1-2*xNu) * ( -(1-xNu) )
          atm =   G * ( 1 + sPsi/(1-2*xNu) )
          amt =   G * ( 1 + sPhi/(1-2*xNu) )
          att = 2*G * ( 1-xNu)/(1-2*xNu)
          If (fte.Gt.fttol) Then
            xLamT = fte/att
            xLamM = 0
            ! Prs = PrsE - xLamT * DdGtdS
            Call AddVec( PrsE, DdGtdS, 1d0, -xLamT, 3, Prs )
            fm = ( Prs(1) - Prs(3) ) /2
     *          -( Prs(1) + Prs(3) ) /2 * sPhi - cCosPhi
            ft = sTens - Prs(3)
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
          ! coulomb + tension
          Det = amm*att-amt*atm
          xLamM = ( att*fme - amt*fte ) / det
          xLamT = (-atm*fme + amm*fte ) / det
          ! Prs = PrsE - xLamM * DdGmdS - xLamT * DdGtdS
          Call AddVec( DdGmdS, DdGtdS, xLamM, xLamT, 3, dPrs )
          Call AddVec( PrsE, dPrs, 1d0, -1d0, 3, Prs )
          fm = ( Prs(1) - Prs(3) ) /2
     *        -( Prs(1) + Prs(3) ) /2 * sPhi - cCosPhi
          ft = sTens - Prs(3)
          ft1= sTens - Prs(1)
          if (ft1.gt.1d-10) write(1,*)'FT1 positive !!!'
          If ( fm .Lt. fmTol .And.
     *         ft .Lt. ftTol       ) Goto 999 ! Return
        Case (3) ! Triax extension
          !!!!!!!!!!!!!!!!!!! EXTENSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! first correct sig1=sig2
          PrsE(2) = (Prs_E(1) + Prs_E(2) ) / 2
          PrsE(1) = PrsE(2)
          Prs (1) = PrsE(1)
          Prs (2) = PrsE(2)

          fme = ( PrsE(1) + PrsE(2) - 2*PrsE(3) ) /4
     *         -( PrsE(1) + PrsE(2) + 2*PrsE(3) ) /4 * sPhi - cCosPhi
          fte = sTens - PrsE(3)
          fmTol = 1d-10*abs(fme)+1d-10
          ftTol = 1d-10*abs(fte)+1d-10
          If (fme.Lt.fmtol .And. fte.Lt.fttol ) Goto 999 ! Return ! Both elastic

          ! Try Coulomb only
          DdGmds(1) = G/(1-2*xNu)/2 * (   (1-2*xNu)  -(1+2*xNu)*sPsi )
          DdGmds(2) = G/(1-2*xNu)/2 * (   (1-2*xNu)  -(1+2*xNu)*sPsi )
          DdGmds(3) = G/(1-2*xNu)/2 * (-2*(1-2*xNu)  -(2      )*sPsi )
          amm =G/4 * (3 + sPsi + sPhi + sPhi*sPsi*(3+2*xNu)/(1-2*xNu) )
          If (fme.Gt.fmtol) Then
            xLamM = fme / amm
            xLamT = 0
            ! Prs = PrsE - xLamM * DdGmdS
            Call AddVec( PrsE, DdGmdS, 1d0, -xLamM, 3, Prs )
            fm = ( Prs(1) + Prs(2) - 2*Prs(3) ) /4
     *          -( Prs(1) + Prs(2) + 2*Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens -  Prs(3)
            ipl = 1 ! Coulomb
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
          ! tension surface
          ipl = 2 ! tens
          DdGtds(1) =  -2*G/(1-2*xNu) * (  xNu  )
          DdGtds(2) =  -2*G/(1-2*xNu) * (  xNu  )
          DdGtds(3) =  -2*G/(1-2*xNu) * ( 1-xNu )

          amt =G * ( 1 + sPhi /(1-2*xNu) )
          atm =G * ( 1 + sPsi /(1-2*xNu) )
          att = 2*G*(1-xNu) /(1-2*xNu)

          If (fte.Gt.fttol) Then
            xLamT = fte/att
            xLamM = 0
            ! Prs = PrsE - xLamT * DdGtdS
            Call AddVec( PrsE, DdGtdS, 1d0, -xLamT, 3, Prs )
            fm = ( Prs(1) + Prs(2) - 2*Prs(3) ) /4
     *          -( Prs(1) + Prs(2) + 2*Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens -  Prs(3)

            ft1= sTens - Prs(1)
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Then ! Goto 999 ! Return
              If (ft1.Gt.0) Then
                Prs(1) = sTens
                Prs(2) = sTens
              End If
              Goto 999 ! Return
            End If
          End If
          ! coulomb + tension
          Det = amm*att-amt*atm
          xLamM = ( att*fme - amt*fte ) / det
          xLamT = (-atm*fme + amm*fte ) / det
          ! Prs = PrsE - xLamM * DdGmdS - xLamT * DdGtdS
          Call AddVec( DdGmdS, DdGtdS, xLamM, xLamT, 3, dPrs )
          Call AddVec( PrsE, dPrs, 1d0, -1d0, 3, Prs )
          fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *        -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
          ft = sTens - ( Prs(2) + Prs(3) ) /2
          ft1= sTens - Prs(1)
          if (ft1.gt.1d-10) write(1,*)'FT1 positive !!!'
          If ( fm .Lt. fmTol .And.
     *         ft .Lt. ftTol       ) Goto 999 ! Return
        Case Default
          Stop ' incorrect iArea in MC_Tens '
      End Select
      ! Normally, no situation should come here
      write(1,*)'Area: ',iArea
      write(1,*)'fme : ',fme
      write(1,*)'fte : ',fte
      write(1,*)'fm  : ',fm
      write(1,*)'ft  : ',ft
!      Call Error ( ' how did we get here ?' )
!      Stop ' how did we get here ?'
  999 Continue
      If (fm.Gt.fmTol .Or. ft.Gt.ftTol) Then
        Write(1,902) iarea,fm,ft
  902   Format('fm>0 or ft>0, iarea:',i5,2f10.5)
      End If
      Return

      End ! MC_Tens
