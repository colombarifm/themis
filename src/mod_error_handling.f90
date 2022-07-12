!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2022 Themis developers
!
!   This file was written by Asdrubal Lozada-Blanco.
!
!---------------------------------------------------------------------------------------------------
!
!   Themis is free software: you can redistribute it and/or modify it under the terms of the GNU 
!   General Public License as published by the Free Software Foundation, either version 3 of the 
!   License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
!   the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License along with this program. If 
!   not, see <https://www.gnu.org/licenses/>.
!
!---------------------------------------------------------------------------------------------------
!> @file   mod_error_handling.f90
!> @author Asdrubal Lozada-Blanco
!> @brief  This module contains procedures to error handling
!> This is a Standard Fortran 2008 compliance code
!> @date - Nov 2019
!> - module created and incorporated into code
!---------------------------------------------------------------------------------------------------

module mod_error_handling
  use iso_fortran_env , only : output_unit, error_unit

  implicit none
  private

  type, public :: error
    integer                       :: code
    character                     :: type
    character(len=:), allocatable :: message, check, tip
  contains
    procedure                     :: error       => raise_error
    procedure                     :: termination => normal_termination
  end type error

  public :: Normal_termination, Raise_error

contains

  subroutine Raise_error(self, type, code, message, check, tip)
    class(error), intent(inout)        :: self
    character, intent(in)              :: type        
    integer, optional, intent(in)      :: code        
    character(*), optional, intent(in) :: message, check, tip

    self%type = type

    if (present( message )) self % message = message
    if (present( check ))   self % check   = check
    if (present( tip ))     self % tip     = tip
    if (present( code ))    self % code    = code

    select case(type)
      case('w')
        !if(.not.present(message)) write(output_unit,'("Warning: unexpected event.")')     
        if (present( message )) write(output_unit,'(/,T5,"Warning: ",a)') self % message
        if (present( check ))   write(output_unit,'(/,T5,"Warning: ",a)') self % check
        if (present( tip ))     write(output_unit,'(/,T5,"Warning: ",a)') self % tip
        if (present( code ))    write(output_unit,'("Code: ",i0)')        self % code
      case('e')
        !if(.not.present(message)) write(output_unit,'("Error: abnormal condition.")')     
        if (present( message )) write(output_unit,'(/,T5,"Error: abnormal condition ",a)') self % message
        if (present( check ))   write(output_unit,'(/,T5,"Please check ",a)')             self % check
        if (present( tip ))     write(output_unit,'(/,T5,"TIP: ",a)')                     self % tip
        if (present( code ))    write(output_unit,'("Code: ",i0)')                        self % code
    end select

  end subroutine Raise_error   

  subroutine Normal_termination(self,iostat,status)
    class(error), intent(inout) :: self
    integer, intent(in)         :: iostat
    character(*), intent(in)    :: status

    if(iostat == 0) then
      select case(status)
        case('f')
          write(output_unit,'(T5,"Sucessfull termination.")')
          stop
        case('i')
          write(output_unit,'(T5,"Normal termination.")')        
      end select
    end if

  end subroutine Normal_termination

end module mod_error_handling
