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
	  
	  
	  Module ModTiming
	  !**********************************************************************
      !
      !   Function : Contains routines for timing the code
      !
      !     $Author$
      !     $Revision$
      !     $Date$
      !     $URL$
      !
      !**********************************************************************
	  
	  
      use ModString
      use ModGlobalConstants
      implicit none
      type FunctionTimerType
        logical :: isInitialised
        character(len=50) :: name
        integer :: numberCalls
        integer :: maxWallTimeStep
        integer :: maxCPUTimeStep
        integer :: startWallTime
        real(REAL_TYPE) :: totalWallTime
        real(REAL_TYPE) :: maxWallTime
        real(REAL_TYPE) :: startCPUTime
        real(REAL_TYPE) :: totalCPUTime
        real(REAL_TYPE) :: maxCPUTime
      end type FunctionTimerType

      private
      type(FunctionTimerType), dimension(100) :: functionTimer
      public :: startTimer
      public :: finishTimer
      public :: giveTimerTable
      public :: initialiseTimer
      public :: initialiseTimerItem
      public :: getID
      public :: IDTimerMain

      integer :: IDTimerMain = 0

      contains

      subroutine startTimer(name, ID, messure)
      !**********************************************************************
      !
      !    Function:  Starts timer
      !
      !**********************************************************************
      implicit none
      character(*), intent(in):: name
      integer, intent(out):: ID
      integer :: i, clockRate
      logical, intent(in), optional:: messure
      logical, external :: ConsiderTimingProcedures

      if (present(messure)) then
          if (.not.messure) then
              ID = 0
              RETURN
          endif
      endif

      do i=1,size(functionTimer)      
          if (functionTimer(i)%isInitialised) then
              if (trim(functionTimer(i)%name) == trim(name)) then
                  ID = i
                   functionTimer(i)%numberCalls = functionTimer(i)%numberCalls + 1             
                  call SYSTEM_CLOCK(functionTimer(i)%startWallTime, clockRate)
                  call CPU_TIME(functionTimer(i)%startCPUTime)
                  RETURN
              endif
          else
              functionTimer(i)%name = trim(name)
              ID = i
              functionTimer(i)%numberCalls = functionTimer(i)%numberCalls + 1
              call SYSTEM_CLOCK(functionTimer(i)%startWallTime, clockRate)
              call CPU_TIME(functionTimer(i)%startCPUTime)
              functionTimer(i)%isInitialised = .true.
              RETURN
          endif
      enddo  
           
      ID = 0
            
    end subroutine startTimer

      integer function getID(name) result(res)
      implicit none
      character(*), intent(in):: name
      integer :: i

      res = 0

      do i=1,size(functionTimer)
        if (functionTimer(i)%isInitialised) then
          if (trim(functionTimer(i)%name) == trim(name)) then
            res = i
            RETURN
          endif
        else
          RETURN
        endif
      enddo

      end function getID

      subroutine finishTimer(ID)
      !**********************************************************************
      !
      !    Function:  Ends timer
      !
      !**********************************************************************
      implicit none
      integer, intent(in):: ID

      if (ID <= 0 .or. ID > size(functionTimer)) RETURN

      call finishTimerWall(ID)
      call finishTimerCPU(ID)

      end subroutine finishTimer

      subroutine giveTimerTable(timerName)
      !**********************************************************************
      !
      !    Function:  Gives message regarding timing of each subroutine
      !
      !**********************************************************************
      implicit none
      integer, parameter :: WALL_TYPE = 1
      integer, parameter :: CPT_TYPE  = 2

      character(*), optional, intent(in):: timerName
      integer :: i, IDReference, timerType
      logical, external :: ConsiderTimingProcedures

      if (present(timerName)) then
          if (trim(lower(timerName)) == 'cpu') then
              timerType = CPT_TYPE
          else
            timerType = WALL_TYPE
          endif
      else
          timerType = WALL_TYPE
      endif

      IDReference = max(getID('main'), 1)
      call GiveMessage(trim(String('Item', '(A6)'))         // trim(String(' Name of subroutine/function/procedure')) // &
          trim(String('Times', '(A13)'))       // trim(String(' Time per call (sec)', '(A15)')) // &
          trim(String(' Total time (sec)', '(A15)')) // trim(String(' Percentage', '(A11)')) )
      do i=1,size(functionTimer)
          if (functionTimer(i)%isInitialised) then
               if (timerType == WALL_TYPE) then
                  call giveTimerItemWall(i, IDReference)
              else
                  call giveTimerItemCPU(i, IDReference)
              endif
          else
              RETURN
          endif
      enddo      
  
    end subroutine giveTimerTable

      subroutine giveTimerItemWall(i, IDReference)
      !**********************************************************************
      !
      !    Function:  Gives message regarding timer item 
      !
      !**********************************************************************
      implicit none
      real(REAL_TYPE), parameter :: PERCENT = 100.0
      integer, intent(in) :: i
      integer, intent(in), optional :: IDReference
      integer :: IDRoutine
      real(REAL_TYPE) :: totalTimeRef

      if (present(IDReference)) then
        IDRoutine = IDReference
        if (IDReference > 0) then
          totalTimeRef = functionTimer(IDReference)%totalWallTime
        endif
      else
        IDRoutine = 0
      endif
      
      if (functionTimer(i)%isInitialised) then
        if (IDRoutine == 0 .or. abs(totalTimeRef) < TINY) then
          call GiveMessage(trim(String(i,'(I6)'))       //'| '// &
                           functionTimer(i)%name(1:40)  //'|'// &
                           trim(String(functionTimer(i)%numberCalls,'(I9)')) //'|'// &
                           trim(String(functionTimer(i)%totalWallTime/functionTimer(i)%numberCalls)) //'|'// &
                           trim(String(functionTimer(i)%totalWallTime)) )
        else
          call GiveMessage(trim(String(i,'(I6)'))       //'| '// &
                           functionTimer(i)%name(1:40)  //'|'// &
                           trim(String(functionTimer(i)%numberCalls,'(I9)')) //'|'// &
                           trim(String(functionTimer(i)%totalWallTime/functionTimer(i)%numberCalls)) //'|'// &
                           trim(String(functionTimer(i)%totalWallTime)) //'|'// &
                           trim(String((PERCENT*functionTimer(i)%totalWallTime)/totalTimeRef, '(F6.2)')) //'%')
          
        endif
      endif

      end subroutine giveTimerItemWall

      subroutine giveTimerItemCPU(i, IDReference)
      !**********************************************************************
      !
      !    Function:  Gives message regarding timer item 
      !
      !**********************************************************************
      implicit none
      real(REAL_TYPE), parameter :: PERCENT = 100.0
      integer, intent(in) :: i
      integer, intent(in), optional :: IDReference
      integer :: IDRoutine
      real(REAL_TYPE) :: totalTimeRef

      if (present(IDReference)) then
        IDRoutine = IDReference
        if (IDReference > 0) then
          totalTimeRef = functionTimer(IDReference)%totalCPUTime
        endif
      else
        IDRoutine = 0
      endif

      if (functionTimer(i)%isInitialised) then
        if (IDRoutine == 0 .or. abs(totalTimeRef) < TINY) then
          call GiveMessage(trim(String(i,'(I6)'))       //'| '// &
                           functionTimer(i)%name(1:40)  //'|'// &
                           trim(String(functionTimer(i)%numberCalls,'(I9)')) //'|'// &
                           trim(String(functionTimer(i)%totalCPUTime/functionTimer(i)%numberCalls)) //'|'// &
                           trim(String(functionTimer(i)%totalCPUTime)) )
        else
          call GiveMessage(trim(String(i,'(I6)'))       //'| '// &
                           functionTimer(i)%name(1:40)  //'|'// &
                           trim(String(functionTimer(i)%numberCalls,'(I9)')) //'|'// &
                           trim(String(functionTimer(i)%totalCPUTime/functionTimer(i)%numberCalls)) //'|'// &
                           trim(String(functionTimer(i)%totalCPUTime)) //'|'// &
                           trim(String((PERCENT*functionTimer(i)%totalCPUTime)/totalTimeRef, '(F6.2)')) //'%')
          
        endif
      endif

      end subroutine giveTimerItemCPU

      subroutine initialiseTimerItem(i)
      !**********************************************************************
      !
      !    Function:  Initializes all timer items
      !
      !**********************************************************************
      implicit none
      integer, intent(in) :: i

      functionTimer(i)%isInitialised = .false.
      functionTimer(i)%name = ' '
      functionTimer(i)%totalWallTime = 0
      functionTimer(i)%startWallTime = 0
      functionTimer(i)%totalCPUTime = 0
      functionTimer(i)%startCPUTime = 0
      functionTimer(i)%numberCalls = 0

      functionTimer(i)%maxWallTimeStep = 0
      functionTimer(i)%maxCPUTimeStep = 0
      functionTimer(i)%maxWallTime = 0
      functionTimer(i)%maxCPUTime = 0

      end subroutine initialiseTimerItem

      subroutine initialiseTimer()
      !**********************************************************************
      !
      !    Function:  Initializes timer
      !
      !**********************************************************************
      implicit none
      integer :: i

      do i=1,size(functionTimer)
        call initialiseTimerItem(i)
      enddo

      end subroutine initialiseTimer

      subroutine finishTimerCPU(ID)
      !**********************************************************************
      !
      !    Function:  Ends timer if maxCPU time exceeded
      !
      !**********************************************************************
      implicit none
      integer, intent(in):: ID
      real(REAL_TYPE) :: finishCPUTime, dt

      call CPU_TIME(finishCPUTime)

      dt = finishCPUTime - functionTimer(ID)%startCPUTime
      if (dt > 0.0d0) then
        functionTimer(ID)%totalCPUTime =  functionTimer(ID)%totalCPUTime + dt
        functionTimer(ID)%startCPUTime = huge(functionTimer(ID)%startCPUTime)

        if (dt > functionTimer(ID)%maxCPUTime) then
          functionTimer(ID)%maxCPUTime = dt
          functionTimer(ID)%maxCPUTimeStep = functionTimer(ID)%numberCalls
        endif
      endif

      end subroutine finishTimerCPU

      subroutine finishTimerWall(ID)
      !**********************************************************************
      !
      !    Function:  Ends timer if maximum wall time exceeded
      !
      !**********************************************************************
      implicit none
      integer, intent(in):: ID
      integer :: clockRate, finishTime
      real(REAL_TYPE) :: dt

      call SYSTEM_CLOCK(finishTime, clockRate)

      dt = float(finishTime - functionTimer(ID)%startWallTime) / float(clockRate)
      if (dt > 0.0d0) then
        functionTimer(ID)%totalWallTime =  functionTimer(ID)%totalWallTime + dt
        functionTimer(ID)%startCPUTime = huge(functionTimer(ID)%startCPUTime)

        if (dt > functionTimer(ID)%maxCPUTime) then
          functionTimer(ID)%maxWallTime = dt
          functionTimer(ID)%maxWallTimeStep = functionTimer(ID)%numberCalls
        endif
      endif

      end subroutine finishTimerWall

      end Module ModTiming
