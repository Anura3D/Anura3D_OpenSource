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


    module ModInitialiseKernel
    !**********************************************************************
    !
    ! Function : Contains initialisation routines of kernel and project
    !
    !     $Revision: 9707 $
    !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
    !
    !**********************************************************************

    use ModReadCalculationData
    use ModGlobalConstants
      
    implicit none

    contains
    
      subroutine ShowDisclaimer()
      !**********************************************************************
      !
      ! Function : Displays Anura3D disclaimer on screen
      !
      !**********************************************************************

      call GiveMessage('')
      call GiveMessage('Anura3D: Numerical modelling and simulation of large deformations'//&
                           'and soil–water–structure interaction using the material point method (MPM)')
          

      end subroutine ShowDisclaimer
        
      subroutine ShowKernelInformation()
      !**********************************************************************
      !
      ! Function : Displays kernel information on screen
      !
      !**********************************************************************

      implicit none
        
      character(len=1023) :: ExeName, getVersion, getLastChangedDate, getLastCompiledDate
          
      call get_command_argument(0, ExeName)
          
      call GiveMessage('Kernel name       : ' // trim(ExeName))
      call GiveMessage('Kernel version    : ' // trim(getVersion()))
      call GiveMessage('Last update date  : ' // trim(getLastChangedDate()))
      call GiveMessage('Last compiled date: ' // trim(getLastCompiledDate()))
        
      end subroutine ShowKernelInformation
        
      subroutine ReadCommandLineParameters()
      !**********************************************************************
      !
      ! Function : Reads the project name from the command line
      !
      !**********************************************************************

      implicit none
        
      integer(INTEGER_TYPE) :: args ! number of command line arguments
          
#ifdef __INTEL_COMPILER
          args = nargs() 
#else
          args = iargc()
#endif

      call Assert(args > 1, 'Kernel called without arguments. Specify project name.')
      call getarg(1, CalParams%FileNames%ProjectName)
      call GiveMessage(NEW_LINE('A') // 'Project name: ' // trim(CalParams%FileNames%ProjectName) )

      end subroutine ReadCommandLineParameters

    end module ModInitialiseKernel
