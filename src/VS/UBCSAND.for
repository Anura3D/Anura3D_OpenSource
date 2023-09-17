      !*****************************************************************************
	  !                                       ____  _____  
      !           /\                         |___ \|  __ \ 
      !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
      !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
      !        / ____ \| | | | |_| | | | (_| |___) | |__| |
      !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
      !
	  !
	  !	  Anura3D - Numerical modelling and simulation of large deformations 
	  !   and soil–water–structure interaction using the material point method (MPM)
      !
	  !	  Copyright (C) 2022  Members of the Anura3D MPM Research Community 
	  !   (See Contributors file "Contributors.txt")
	  !
      !	  This program is free software: you can redistribute it and/or modify
      !	  it under the terms of the GNU Lesser General Public License as published by
      !	  the Free Software Foundation, either version 3 of the License, or
      !	  (at your option) any later version.
	  !
      !	  This program is distributed in the hope that it will be useful,
      !	  but WITHOUT ANY WARRANTY; without even the implied warranty of
      !	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      !	  GNU Lesser General Public License for more details.
	  !
      !	  You should have received a copy of the GNU Lesser General Public License
      !	  along with this program.  If not, see <https://www.gnu.org/licenses/>.
	  !
	  !!*****************************************************************************
	  !
   ! 
   ! module UBCSAND
   ! !**********************************************************************
   ! !   Constitutive model 
   ! !
   ! !
   ! !**********************************************************************
   ! 
   ! 
   ! ! check which of these are needed in this module
   ! ! below are copied from DynmaicImplicitGeneralizedAlphaScheme 
   ! use ModCounters 
   ! use ModReadCalculationData
   ! use ModMPMData
   ! use ModParticle
   ! use ModMeshInfo
   ! use ModWriteTestData
   ! use ModEmptyElements
   ! use ModConvectivePhase
   ! use ModWriteMPMData
   ! use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
   ! use ModLagrangianPhase 
   ! use ModDYNConvectivePhase
   ! use ModConvectivePhase
   ! 
   ! 
   ! 
   ! implicit none 
   ! 
   ! contains 
   ! 
   !
   ! subroutine run 
   ! !**********************************************************************
   ! !
   ! !
   ! !**********************************************************************
   !
   ! implicit none 
   ! 
   ! real(REAL_TYPE), intent(in) :: m_ratmax 
   ! 
   ! ! s is state 
   ! 
   ! call initialize() ! initialize varibales 
   ! 
   ! !plasticity indicator 
   ! plas = 0
   ! 
   ! ! states 
   ! ! shear 
   ! shear_now = 0
   ! shear_past = 0
   ! 
   ! ! tension 
   ! tension_now = 0 
   ! tension_past = 0 
   ! 
   ! ! // this is the first element in the 
   ! ! elasticity matrix (Hooke's law)
   ! m_e1       = m_k + (4.0 * m_g / 3.0)
   ! 
   ! ! // this is the second element in the 
   ! ! elasticity matrix (Hooke's law)
   ! m_e2       = m_k - (2.0 * m_g / 3.0)
   ! 
   ! m_sh2      = 2.0 * m_g
   ! 
   ! if (m_ind .ne. 0.0) then 
   !     m_ind = 2.0
   !     ! // --- get new trial stresses from old, assuming elastic increments ---
   !     cs2 = 1.0 ! I believe this is a directional cosine variable 
   !     alams = 0.0
   !     if (m_ratmax > 0.0) then 
   !        
   !         ! // Trial stress assuming elastic increment
   !         ! strain increment stored in zdeXX
   !         call stnE_(zde11, zde22, zde33, zde12, &
   !                    s, s11, s22, s33, s12)
   !         
   !         ! current stress stored in zsXX
   !         call stnS_(zs11, zs22, zs33, zs12, &
   !                    s, s11, s22, s33, s12)
   !         
   !         ! use elasticity tensor to find 
   !         !  how much would the stress be after
   !         !   elastic stress incdement. 
   !         !// this is sigma_xx put in elasticity matrix
   !         s11i  = zs11 + (zde22 + zde33) * m_e2 + zde11 * m_e1 
   !         !// this is sigma_xx put in elasticity matrix
   !         s22i  = zs22 + (zde11 + zde33) * m_e2 + zde22 * m_e1 
   !         !// fourth row in elasticity tensor
   !         s12i  = zs12 + zde12 * m_sh2
   !         !// assuming Poisson ratio = 0.5
   !         s33i  = .5*(s11i+s22i)
   !         !// calculating stress change
   !         sdif  = s11i - s22i
   !         !// center stress point assuming s11i and s22i are principal
   !         s0    = 0.5 * (s11i + s22i)
   !         !// radius of mohr circle
   !         !// radius = Rm = [((sigma_x - sigma_y)^2)/4 + sigma_xy^2]^0.5
   !         !!// taking half as a common factor
   !         rad   = 0.5 * sqrt (sdif*sdif + 4.0 * s12i*s12i) 
   !         
   !         !//; ----principal stresses ---
   !         si    =  s0 - rad !// s1
   !         sii   =  s0 + rad !// s3
   !         psdif =  si - sii !// principal stress difference
   !         
   !         !//; --- s33 is intermediate (other two cases not allowed) ---
   !         s1    = si
   !         s2    = s33i
   !         s3    = sii
   !         
   !         !//; --- shear yield criterion ---
   !         fs    = (s1 - s3) * m_nphi !// yield surface does not include cohesion
   !         alams = 0.0 !// I don't know what alams is...
   !         !// I guess alam"s" stands for shear
   !         
   !         !//; --- tensile yield criterion ---
   !         !// note positive is tension
   !         !// we are taking it to negative to make it consistent 
   !         !with the yield surface (i.e., if it is less than zero)
   !         ft    = - s3 !// Assumes tensile strength = 0.
   !         alamt = 0.0
   !         !// I guess alam"t" stands for tension
   !         
   !         !//; --- tests for failure ---
   !         !// tension first
   !         if (ft < 0.0) then 
   !         
   !             bisc = sqrt(1.0 + m_nphi * m_nphi) + m_nphi
   !             pdiv = ft + s1*bisc 
   !             
   !             if (pdiv < 0.0) then 
   !                 
   !                 alams = fs / m_x1;
   !                 s1    = s1 - alams * (m_e1 - m_e2 * m_npsi) !// psi = plastic strain increment
   !                 s2    = s2 - alams * m_e2 * (1.0 - m_npsi)  !// psi = plastic strain increment
   !                 s3    = s3 - alams * (m_e2 - m_e1 * m_npsi) !// psi = plastic strain increment
   !                 
   !                 m_ind  = 1.0 !// it looks like there is an indicator for each --> tension and shear
   !                 plas = 1 !// turning on plasticity -> shear failure
   !                 ! I need to think about this
   !                 if ( state_==1 .or. shear_now==1 ) then 
   !                     !s->state_ |= shear_now !// start shearing now
   !                 end if 
   !                 
   !                 
   !             else 
   !                 
   !                 ! //; ---     tension failure ---
   !                 alamt = ft / m_e1 !// dividing by elasticity matrix (first component)
   !                 tco   = alamt * m_e2 !// multiplying by elasticity matrix (second component)
   !                 s1    = s1 + tco !////correcting s1 and s2
   !                 s2    = s2 + tco !//correcting s1 and s2
   !                 s3    = 0.0 !// tension cutoff
   !                 m_ind  = 3.0 !// it looks like there is an indicator for each --> tension only
   !                 plas   = 2 !// -> tension failure
   !                 
   !                 if ( state_==1 .or. tension_now==1 ) then 
   !                     !s->state_ |= tension_now;
   !                 end if 
   !                 
   !                 
   !                 
   !             
   !         
   !             end if   
   !             
   !             
   !         else 
   !             
   !             if (fs < 0.0) then 
   !                 !// shear only 
   !                 !//; ---    shear failure ---
   !                 
   !                 alams = fs / m_x1
   !                 s1    = s1 - alams * (m_e1 - m_e2 * m_npsi) !// psi = plastic strain increment
   !                 s2    = s2 - alams * m_e2 * (1.0 - m_npsi) !// psi = plastic strain increment
   !                 s3    = s3 - alams * (m_e2 - m_e1 * m_npsi) !// psi = plastic strain increment
   !                 m_ind  = 1.0
   !                 plas = 1
   !                 
   !                 if ( state_==1 .or. shear_now==1 ) then 
   !                     !s->state_ |= shear_now
   !                 end if
   !                 
   !             else 
   !                 
   !                 !//; ---    no failure ---
   !                 !// keeping the elastic trail
   !                 
   !                 
   !                 call stnS_(s11i, s22i, s33i, s12i, &
   !                    s, rs11, rs22, zs33, rs33)
   !                 
   !                 
   !                 !s->stnS_.rs11() = s11i
   !                 !s->stnS_.rs22() = s22i
   !                 !s->stnS_.rs33() = s33i
   !                 !s->stnS_.rs12() = s12i
   !                 
   !
   !             end if 
   !         end if 
   !     end if 
   ! end if 
   ! 
   ! 
   ! 
   ! 
   ! !//; --- direction cosines ---
   ! cs2 = 1.0
   ! si2 = 0.0
   ! 
   !
   ! !// if principal stress difference is not equal to zero
   ! !//$psdif =  $si - $sii; // principal stress difference
   ! if ($psdif .ne. 0.0) then !{ // if psdif is not equal to zero
   !      !// $sdif  = $s11i - $s22i;
   !      cs2   = sdif / psdif
   !      si2   = 2.0 * s12i / psdif   
   ! end if
   ! 
   ! 
   ! !//; --- resolve back to global axes ---
   ! if (plas) then 
   !     !// I assume that as long as plas is non-zero, we go into this loop
   !     dc2  = (s1 - s3) * cs2
   !     dss  =  s1 + s3
   !     
   !     
   !     !call stnS_(zs11, zs22, zs33, zs12, & !output 
   !     !               s, s11, s22, s33, s12) !input 
   !   
   !     rs11 = 0.5 * (dss + dc2)
   !     rs22 = 0.5 * (dss - dc2)
   !     rs33 = s2
   !     rs12 = 0.5 * (s1  - s3) * si2
   !     
   !     viscous_ = false
   !     
   ! else 
   !     
   !     viscous_ = true
   ! end if
   ! 
   !     
   !     
   ! 
   ! end if 
   ! 
   !
   ! !//;  -----------------------------------------------------------------------
   ! !//;  -----UBC add on to account for change of elastic and plastic parameters
   ! !//;  -----------------------------------------------------------------------
   ! !//; --------  PLASTIC STRAINS --- 
   ! dZart=s->getSubZoneVolume()
   ! call getSubZoneVolume(dZart, s) !//getting the area of the subzone
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! end subroutine run  
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 
   ! 