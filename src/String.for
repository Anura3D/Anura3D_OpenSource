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
	  
	  
	  module ModString
      !**********************************************************************
      !
      !    Function:  returns string out of any types of variable.
      !
      !
      !     $Revision: 9707 $
      !     $Date: 2022-04-14 14:56:02 +0200 (do, 14 apr 2022) $
      !
      !**********************************************************************
      use ModGlobalConstants
      
      implicit none

      ! String functions
      !> makes string out of different types of variables
      public:: String 

      !>  in string that are blank or comma separated
      public:: CountItems

      !>  case the text arg
      public:: Lower

      !> returns true if at least one charcter is letter
      public:: ContainLetters

      !> returns a given item in a string
      public:: GetItem

      !> returns location of the first character of the next item
      public:: GetFirstLocationNextItem

      !> returns location of the last character of the next item
      public:: GetLastLocationNextItem

      !> returns true if a key is found in a string
      public:: HasKey


      interface String
        module procedure StringReal
        module procedure StringInt
        module procedure StringInt8
        module procedure StringLogical
        module procedure StringString
        module procedure StringRealArray
        module procedure StringIntArray
      end interface String

      private

    contains
    
        !**********************************************************************
        !    StringReal, StringInt, StringInt8, StringLogical, StringString, 
        !         StringRealArray, StringIntArray
        !
        !    Function:  Make string out of different types of variables
        !
        !**********************************************************************
        
        character(len = 64) function StringReal(Variable, fmt) result(str)
        implicit none
        real(REAL_TYPE), intent(in):: Variable
        character(*), intent(in), optional :: fmt

        if (present(fmt)) then
          write(str, fmt) Variable
        else
          write(str, '(E14.4E3)') Variable
        endif

        end function StringReal

        ! -----------------------------------------------------------------------------------------
        character(len = 16) function StringInt(Variable, fmt) result(str)
        implicit none
        integer(INTEGER_TYPE), intent(in):: Variable
        character(*), intent(in), optional :: fmt

        if (present(fmt)) then
          write(str, fmt) Variable
        else
          if (Variable < 1000000) then
            write(str, '(I6)') Variable
          else
            write(str, '(I12)') Variable
          endif
        endif

        end function StringInt

        ! -----------------------------------------------------------------------------------------
        character(len = 16) function StringLogical(Variable, fmt) result(str)
        implicit none
        logical, intent(in):: Variable
        character(*), intent(in), optional :: fmt

        if (present(fmt)) then
          write(str, fmt) Variable
        else
          write(str, *) Variable
        endif

        end function StringLogical

        ! -----------------------------------------------------------------------------------------
        character(len = 24) function StringInt8(Variable, fmt) result(str)
        implicit none
        integer(8), intent(in):: Variable
        character(*), intent(in), optional :: fmt

        if (present(fmt)) then
          write(str, fmt) Variable
        else
          write(str, '(I14)') Variable
        endif

        end function StringInt8

        ! -----------------------------------------------------------------------------------------
        character(len = 4096) function StringString(Variable, fmt) result(str)
        implicit none
        character(*), intent(in):: Variable
        character(*), intent(in), optional :: fmt

        if (present(fmt)) then
          write(str, fmt) Variable
        else
          write(str, *) Variable
        endif

        end function StringString

        ! -----------------------------------------------------------------------------------------
        character(len = 4096) function StringRealArray(Variable, first, last) result(str)
        implicit none
        real(REAL_TYPE), intent(in):: Variable(*)
        integer, intent(in):: first, last

        write(str, *) Variable(first:last)

        end function StringRealArray

        ! -----------------------------------------------------------------------------------------
        character(len = 4096) function StringIntArray(Variable, first, last) result(str)
        implicit none
        integer(INTEGER_TYPE), intent(in):: Variable(*)
        integer(INTEGER_TYPE), intent(in):: first, last

        write(str, *) Variable(first:last)

        end function StringIntArray


        logical function ContainLetters(s) result(res)
        !**********************************************************************
        !
        !    Function:  Returns true if at least one character is a letter
        !
        !**********************************************************************
        
        implicit none
        character(*),intent(IN) :: s
        integer(INTEGER_TYPE) :: i

        res = .false.
        do i=65,90
          ! Uppercase letters
          res = scan(s, char(i)) > 0
          if (res) RETURN
        enddo

        do i=97,122
          ! Lowercase letters
          res = scan(s, char(i)) > 0
          if (res) RETURN
        enddo

        end function ContainLetters
        
        

        integer(INTEGER_TYPE) function CountItems(s)  
        !**********************************************************************
        !
        !    Function:  Counts items in a string that are blank or comma separated
        !
        !**********************************************************************
        
        implicit none
        character(*), intent(in) :: s
        integer(INTEGER_TYPE) :: i, k

        k = 1
        do i = 1,len_trim(s)-1
          if (.not.IsDelimiter(s(i:i)) .and. IsDelimiter(s(i+1:i+1))) k = k+1
        end do
        i = len_trim(s)
        if (IsDelimiter(s(i:i))) k = k-1 
        
        CountItems = k

        end function CountItems

        subroutine GetItem(s, item, outs)
        !**********************************************************************
        !
        !    Function:  Returns a given item in a string
        !
        !**********************************************************************
        
        implicit none
        character(*), intent(in):: s
        character(*), intent(out):: outs
        integer(INTEGER_TYPE), intent(in) :: item
        integer(INTEGER_TYPE) itemCounter, length, iFirstOfWord, iLastOfWord

        outs = ' '

        ! find location of given item:
        iFirstOfWord = GetFirstLocationNextItem(s)
        iLastOfWord  = GetLastLocationNextItem(s)

        do itemCounter=1,item
          if (itemCounter == item) then
            length = iLastOfWord - iFirstOfWord + 1
            outs(1:length) = s(iFirstOfWord:iLastOfWord)
            RETURN
          else
            iFirstOfWord = GetFirstLocationNextItem(s,iLastOfWord)
            iLastOfWord  = GetLastLocationNextItem(s,iLastOfWord)
          endif
        enddo

        end subroutine GetItem


        function Lower(s1) result (s2)
        !**********************************************************************
        !
        !    Function:  Case the text argument
        !
        !**********************************************************************
        
        implicit none
        character(*)       :: s1
        character(len(s1)) :: s2
        character          :: ch
        integer(INTEGER_TYPE), parameter :: DUC = ichar('A') - ichar('a')
        integer(INTEGER_TYPE)            :: i

        do i = 1,len(s1)
          ch = s1(i:i)
          if (ch >= 'A'.and.ch <= 'Z') ch = char(ichar(ch) - DUC)
          s2(i:i) = ch
        end do
        end function Lower
       
        integer(INTEGER_TYPE) function GetFirstLocationNextItem(s, LocLastCharacter) result(location)
        !**********************************************************************
        !
        !    Function:  Returns location of first character of next item
        !
        !**********************************************************************
        implicit none
        character(*), intent(in):: s
        integer, intent(in), optional :: LocLastCharacter
        integer i, locLast
        
         if (present(LocLastCharacter)) then
          locLast = LocLastCharacter
        else
          locLast = 0
        endif

        ! find location of first character:
        location = len_trim(s)
        do i = locLast+1,len_trim(s)
          if (.not.IsDelimiter(s(i:i))) then
            location = i
            Return
          endif
        enddo

        end function GetFirstLocationNextItem


        integer(INTEGER_TYPE) function GetLastLocationNextItem(s, LocLastCharacter) result(location)
        !**********************************************************************
        !
        !    Function:  Returns location of last character of the next item
        !
        !**********************************************************************
        implicit none
        character(*), intent(in):: s
        integer(INTEGER_TYPE), intent(in), optional :: LocLastCharacter
        integer(INTEGER_TYPE) i, locLast,iFirst

        if (present(LocLastCharacter)) then
          locLast = LocLastCharacter
        else
          locLast = 0
        endif

        ! find location of first character:
        location = len_trim(s)
        do i = locLast+1,len_trim(s)
          if (.not.IsDelimiter(s(i:i))) then
            iFirst = i
            EXIT
          endif
        enddo

        do i = iFirst,len_trim(s)
          if (.not.IsDelimiter(s(i:i)) .and. IsDelimiter(s(i+1:i+1)) ) then
            location = i
            RETURN
          endif
        enddo

        end function GetLastLocationNextItem


        logical function IsDelimiter(c) result(res)
        !**********************************************************************
        !
        !    Function:  Returns true if delimiter is found in string
        !
        !**********************************************************************
        implicit none
        character, intent(in):: c

        res =     c == ' ' &
             .or. c == ',' &
             .or. c == '(' &
             .or. c == ')' &
             .or. c == '[' &
             .or. c == ']' &
             .or. c == '{' &
             .or. c == '}' &
             .or. c == '''' &
             .or. c == '"' &
             .or. c == ';' &
             .or. c == ':'

        end function IsDelimiter


        logical function HasKey(s,key) result(res)
        !**********************************************************************
        !
        !    Function:  Returns true if a key is found in a string
        !
        !**********************************************************************
        implicit none
        character(*), intent(in):: s
        character(*), intent(in):: key
        character(len_trim(s)) :: itemString
        integer(INTEGER_TYPE) :: item

        res = .false.
        do item=1,CountItems(s)
          call GetItem(s, item, itemString)
          if (trim(lower(key)) == trim(lower(itemString))) then
            res = .true.
            RETURN
          endif
        enddo
        end function HasKey

      end module ModString
