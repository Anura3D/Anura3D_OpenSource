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


      module ModErrorHandler
      !**********************************************************************
      !
      !>    Function:  Gives errors or warnings to the user.
      !!    It tries to write the last results into the result files
      !!    and stops the calculation when an error happens.
      !!    If the error happens during resizing some arrays,
      !!    the code is not able to write the results, as the arrays are not compatible anymore.
      !!    In case of warnings or messages, it only displays the message
      !!    on screen and writes it in the log file without stopping the calculation.
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants
      use ModGlobalData
      use ModString
      use ModFeedback

      implicit none
      private

      public :: GiveErrorMessage, &
                GiveNoErrorMessage, &
                GiveWarningMessage, &
                GiveAllocationError, &
                GiveDeallocationError, &
                CheckConditionWithError, &
                CheckConditionWithWarning

      contains

      subroutine GiveNoErrorMessage(message)
      !**********************************************************************
      !
      ! Function:  gives a message without error. It calls GiveErrorOrWarning with proper
      !            error type.
      !
      !**********************************************************************
      implicit none
      ! arguments:
      character*(*), intent(in) :: message

      call GiveErrorOrWarning(message, NO_ERROR)

      end subroutine GiveNoErrorMessage

      subroutine GiveErrorMessage(message)
      !**********************************************************************
      !
      ! Function:  gives error. It calls GiveErrorOrWarning with proper
      !            error type.
      !
      !**********************************************************************
      implicit none
      ! arguments:
      character*(*), intent(in) :: message

      call GiveErrorOrWarning(message, ERROR)

      end subroutine GiveErrorMessage

      subroutine GiveWarningMessage(message)
      !**********************************************************************
      !
      ! Function:  gives warning. It calls GiveErrorOrWarning with proper
      !            error type.
      !
      !**********************************************************************
      implicit none
      ! arguments:
      character*(*), intent(in) :: message

      call GiveErrorOrWarning(message, WARNING)

      end subroutine GiveWarningMessage

      subroutine GiveAllocationError(IErr, ArrayName, Location)
      !**********************************************************************
      !
      ! Function:  checks Allocation errors and gives error message
      !
      !**********************************************************************
      implicit none
      integer(INTEGER_TYPE), intent(in):: IErr
      character*(*), intent(in) :: ArrayName
      character*(*), intent(in) :: Location
      character(len = MessageSize) :: Message

      if (iErr /= 0) then
        write (Message, *) 'Allocation Error: ArrayName: ', trim(ArrayName), &
                                ' in: ', trim(Location)
        call GiveErrorMessage(Message)
      endif

      end subroutine GiveAllocationError

      subroutine GiveDeallocationError(IErr, ArrayName, Location)
      !**********************************************************************
      !
      ! Function:  checks deallocation errors and gives error message
      !
      !**********************************************************************
      implicit none
      integer(INTEGER_TYPE), intent(in):: IErr
      character*(*), intent(in) :: ArrayName
      character*(*), intent(in) :: Location
      character(len = MessageSize) :: Message

      if (iErr /= 0) then
        write (Message, *) 'Deallocation Error: ArrayName: ', trim(ArrayName), &
                                ' in: ', trim(Location)
        call GiveErrorMessage(Message)
      endif

      end subroutine GiveDeallocationError

      subroutine CheckConditionWithWarning(condition, message)
      !**********************************************************************
      !
      ! Function:  checks condition and gives error or warning
      !
      !**********************************************************************
      implicit none

      ! arguments:
      logical, intent(in):: condition
      character*(*), intent(in) :: message

      ! procedure starts------------------------------------------------
      call CheckConditionWithError(condition, message, WARNING)

      end subroutine CheckConditionWithWarning

      subroutine CheckConditionWithError(condition, message, errorTypeInput)
      !**********************************************************************
      !
      ! Function:  checks condition and gives error or warning
      !
      !**********************************************************************
      implicit none

      ! arguments:
      logical, intent(in):: condition
      character*(*), intent(in) :: message
      integer(INTEGER_TYPE), optional, intent(in) :: errorTypeInput
      ! local variables
      integer(INTEGER_TYPE) :: errorType

      ! procedure starts------------------------------------------------
      if (.not.condition) then
        if (present(errorTypeInput)) then
          errorType = errorTypeInput
        else
          errorType = ERROR
        endif
        call giveErrorOrWarning(message, errorType)
      endif

      end subroutine CheckConditionWithError

      subroutine GiveErrorOrWarning(messageInput, errorTypeInput)
      !**********************************************************************
      !
      ! Function:  gives error or warning. It writes the message on screen and
      !            in the out file.
      ! errorType: specifies if it is an error or warning
      !            (default is error)
      !
      !**********************************************************************
      use ModWriteResultData, &
        only: WriteTimeStepResults
      use ModReadCalculationData, &
        only: CalParams
      implicit none

      ! constants:
      character(len = 10), parameter :: ERROR_MSG  = '[ERROR] : '
      character(len = 12), parameter :: WARNING_MSG= '[WARNING] : '
      character(len = 256) :: TimeStepMSG
      ! arguments:
      character*(*), intent(in) :: messageInput
      integer (INTEGER_TYPE), optional, intent(in) :: errorTypeInput
      ! local variables
      integer (INTEGER_TYPE):: errorType

      character(len = MessageSize) :: message

      ! procedure starts------------------------------------------------
      if (present(errorTypeInput)) then
        errorType = errorTypeInput
      else
        errorType = ERROR
      endif

      TimeStepMSG = '(time step: ' // trim(String(CalParams%TimeStep)) // ')'
      if (errorType == WARNING) then
        message = WARNING_MSG // trim(TimeStepMSG) //' ' // trim(messageInput)
      elseif(errorType == ERROR) then
        message = ERROR_MSG // trim(TimeStepMSG) //' ' // trim(messageInput)
      else
        ! just a message
        message = trim(messageInput)
      endif

      ! write on screen
      write (*,fmt='(a)') trim(message)

      ! write into outfile
      call WriteMessageInFile(message, OUTunit)

      ! write into LOGfile
      call WriteMessageInFile(message, LOGunit)

      if (errorType == ERROR) then
        call flush(OUTunit)
        call flush(LOGunit)
        if (CalParams%TimeStep > 1) then
          write (*,fmt='(a)') 'trying to write results...'
          call WriteTimeStepResults(.true.)
        endif
        STOP
      endif

      end subroutine GiveErrorOrWarning

      end module ModErrorHandler


      subroutine GiveMessage(message)
      !**********************************************************************
      !
      !> Gives a message to the user
      !
      !**********************************************************************
      use ModErrorHandler
      implicit none
      ! arguments:
      character*(*), intent(in) :: message

      call GiveNoErrorMessage(message)

      end subroutine GiveMessage


      subroutine GiveError(message)
      !**********************************************************************
      !
      !> Gives an error, tries to write the last results and stops the calculation
      !
      !**********************************************************************      
      use ModErrorHandler
      implicit none
      ! arguments:
      character*(*), intent(in) :: message

      call GiveErrorMessage(message)

      end subroutine GiveError


      subroutine GiveWarning(message)
      !**********************************************************************
      !
      !> Gives a warning to the user with a defined message
      !
      !**********************************************************************       
      use ModErrorHandler
      implicit none
      ! arguments:
      character*(*), intent(in) :: message

      call GiveWarningMessage(message)

      end subroutine GiveWarning


      subroutine AllocationError(IErr, ArrayName, Location)
      !**********************************************************************
      !
      !> Gives an error, in case of allocation error
      !
      !********************************************************************** 
      use ModErrorHandler
      use ModGlobalConstants 
      implicit none
      integer (INTEGER_TYPE), intent(in):: IErr
      character*(*), intent(in) :: ArrayName
      character*(*), intent(in) :: Location

      call GiveAllocationError(IErr, ArrayName, Location)

      end subroutine AllocationError


      subroutine DeallocationError(IErr, ArrayName, Location)
      !**********************************************************************
      !
      !> Gives an error, in case of deallocation error
      !
      !********************************************************************** 
      use ModErrorHandler
      use ModGlobalConstants
      implicit none
      integer (INTEGER_TYPE), intent(in):: IErr
      character*(*), intent(in) :: ArrayName
      character*(*), intent(in) :: Location

      call GiveDeallocationError(IErr, ArrayName, Location)

      end subroutine DeallocationError


      subroutine Assert(condition, message)
      !**********************************************************************
      !
      !> If condition is not true, gives an error
      !
      !********************************************************************** 
      use ModErrorHandler
      implicit none
      ! arguments:
      logical, intent(in):: condition
      character*(*), intent(in) :: message

      if (.not.condition) then
        call CheckConditionWithError(condition, message)
      endif

      end subroutine Assert


      subroutine AssertWarning(condition, message)
      !**********************************************************************
      !
      !> If condition is not true, gives a warning
      !
      !********************************************************************** 
      use ModErrorHandler
      implicit none
      ! arguments:
      logical, intent(in):: condition
      character*(*), intent(in) :: message

      if (.not.condition) then
        call CheckConditionWithWarning(condition, message)
      endif

      end subroutine AssertWarning

