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
	
	
	program Anura3D
    
    use ModReadGeometryData

      call InitialiseDimension()
      call Kernel()
  
    end program Anura3D
    
    
      ! DOXYGEN templates for functions and subroutines. Place before the function/subroutine.  
    
      !*************************************************************************************
      !    FUNCTION: NameOfFunction
      ! 
      !    DESCRIPTION:
      !>   Description of function
      !
      !>   @note : Notes
      !
      !>   @param[in] ParameterName : ParameterDescription
      !
      !>   @return NameOfReturnValue : ParameterDescription
      !
      !*************************************************************************************
    
      !*************************************************************************************
      !    SUBROUTINE: NameOfSubroutine
      ! 
      !    DESCRIPTION:
      !>   Description of subroutine
      !
      !>   @note : Notes
      !
      !>   @param[in] ParameterName : ParameterDescription
      !>   @param[out] ParameterName : ParameterDescription
      !>   @param[inout] ParameterName : ParameterDescription
      !
      !*************************************************************************************