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


      module ModFeedback
      !**********************************************************************
      !
      !    Function:  Gives errors or warnings to the user.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants

      implicit none
      private

      public :: WriteMessageInFile, &
                WriteInLogFile

      contains

      subroutine WriteMessageInFile(message, UnitNumber)
      !**********************************************************************
      !
      ! Function:  It writes a message in the file with UnitNumber.
      !
      !**********************************************************************
      implicit none

      ! arguments:
      character*(*), intent(in) :: message
      integer(INTEGER_TYPE), intent(in):: UnitNumber
      character(len = 10) :: IAction
      logical :: isOpened

      ! procedure starts------------------------------------------------
      ! write into file
      inquire(UNIT = UnitNumber, opened = isOpened, action=IAction)
      if (isOpened) then
        if (trim(IAction) == 'WRITE' .or. trim(IAction) == 'READWRITE') then
          write (UnitNumber, fmt='(a)') trim(message)
        endif
      endif

      end subroutine WriteMessageInFile

      subroutine WriteInLogFile(message, FeedbackLevel)
      !**********************************************************************
      !
      ! Function:  It writes a message in the file in LOGunit
      !
      !**********************************************************************
      implicit none
      ! arguments:
      character*(*), intent(in) :: message
      integer(INTEGER_TYPE), optional, intent(in):: FeedbackLevel
      ! local variables
      integer(INTEGER_TYPE) level
      integer(INTEGER_TYPE), external :: GetFeedbackLevel

      if (present(FeedbackLevel)) then
        level = FeedbackLevel
      else
        level = FEEDBACK_LEVEL_DEBUG
      endif

      if (level >= GetFeedbackLevel()) then
        call WriteMessageInFile(message, LOGunit)
      endif

      end subroutine WriteInLogFile

      end module ModFeedback


