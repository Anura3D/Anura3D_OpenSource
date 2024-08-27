	!*****************************************************************************
    !                                       ____  _____  
    !           /\                         |___ \|  __ \ 
    !          /  \   _ __  _   _ _ __ __ _  __) | |  | |
    !         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
    !        / ____ \| | | | |_| | | | (_| |___) | |__| |
    !       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/ 
    !
    !
	!	Anura3D - Numerical modelling and simulation of large deformations 
    !   and soil–water–structure interaction using the material point method (MPM)
    !
    !	Copyright (C) 2023  Members of the Anura3D MPM Research Community 
    !   (See Contributors file "Contributors.txt")
    !
    !	This program is free software: you can redistribute it and/or modify
    !	it under the terms of the GNU Lesser General Public License as published by
    !	the Free Software Foundation, either version 3 of the License, or
    !	(at your option) any later version.
    !
    !	This program is distributed in the hope that it will be useful,
    !	but WITHOUT ANY WARRANTY; without even the implied warranty of
    !	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !	GNU Lesser General Public License for more details.
    !
    !	You should have received a copy of the GNU Lesser General Public License
    !	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	!
    !*****************************************************************************


    module ModIsort
    !**********************************************************************
    !
    !    Function:  
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************
    
    use ModGlobalConstants, only: INTEGER_TYPE, REAL_TYPE
    implicit none
    
    contains
    
    subroutine IHpSort(IA, N)
      !----------------------------------------
      !
      !  function: sort an integer array IA(1:N) into ascending order using the heapsort algorithm
      !            adapted from "1986-92 Numerical Recipes Software"
      !
      !  edit: 2016-03-30 Miriam Mieremet 
      !
      !  IA   I/O   I()   integer array
      !  N    I     I     length of integer array
      !
      !----------------------------------------
      use ModGlobalConstants
      
      implicit none
      
        ! arguments
        integer(INTEGER_TYPE), intent(in) :: N
        integer(INTEGER_TYPE), dimension(N), intent(inout) :: IA
        
        ! local variables
        integer(INTEGER_TYPE) :: IIA
        integer(INTEGER_TYPE) :: L, IR, I, J
        
        if (N < 2) return
        
        L = N / 2 + 1
        IR = N
        
        do
          if (L > 1) then             ! forming heap
            L = L - 1                    
            IIA = IA(L)
          else                           ! sorting heap
            IIA = IA(IR)                 
            IA(IR) = IA(1)
            IR = IR - 1
            if (IR == 1) then
              IA(1) = IIA
              return
            endif
          endif
          
          I = L
          J = L + L
          
          do                            ! sifting down
            if (J > IR) exit
            
            if (J < IR) then
              if (IA(J) < IA(J + 1) ) J = J + 1
            end if
            
            if (IIA < IA(J) ) then
              IA(I) = IA(J)
              I = J
              J = J + J
            else
              J = IR + 1
            endif
          end do
          
          IA(I) = IIA
        end do
        
      end subroutine IHpSort

      subroutine IHpSortKind8(IA, N)
      !----------------------------------------
      !
      !  function: sort an integer array IA(1:N) into ascending order using the heapsort algorithm
      !            adapted from "1986-92 Numerical Recipes Software"
      !
      !  edit: 2016-03-30 Miriam Mieremet 
      !
      !  IA   I/O   I()   integer array (kind = 8)
      !  N    I     I     length of integer array
      !
      !----------------------------------------
      use ModGlobalConstants
      
      implicit none
      
        ! arguments
        integer, intent(in) :: N
        integer(kind = 8), dimension(N), intent(inout) :: IA
        
        ! local variables
        integer(kind = 8) :: IIA
        integer :: L, IR, I, J
        
        if (N < 2) return
        
        L = N / 2 + 1
        IR = N
        
        do
          if (L > 1) then             ! forming heap
            L = L - 1                    
            IIA = IA(L)
          else                           ! sorting heap
            IIA = IA(IR)                 
            IA(IR) = IA(1)
            IR = IR - 1
            if (IR == 1) then
              IA(1) = IIA
              return
            endif
          endif
          
          I = L
          J = L + L
          
          do                            ! sifting down
            if (J > IR) exit
            
            if (J < IR) then
              if (IA(J) < IA(J + 1) ) J = J + 1
            end if
            
            if (IIA < IA(J) ) then
              IA(I) = IA(J)
              I = J
              J = J + J
            else
              J = IR + 1
            endif
          end do
          
          IA(I) = IIA
        end do
        
      end subroutine IHpSortKind8
end module
