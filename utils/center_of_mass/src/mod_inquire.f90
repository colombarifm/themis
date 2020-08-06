!---------------------------------------------------------------------------------------------------
! SAS_GRID: A code to obtain the solvent accessible surface (SAS) around a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Themis developers
!                 Laboratory of Theoretical Chemistry (LQT) - Federal University of São Carlos 
!                 <http://www.lqt.dq.ufscar.br>
!
!   Please cite: 
!
!   This file was written by Felippe M. Colombari and Asdrubal Lozada-Blanco.
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
!> @file   mod_inquire.f90
!> @author Felippe M. Colombari
!> @brief This module performs inquire checks for all files that will be read
!> @date - Nov, 2019                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_inquire
  use mod_error_handling

  implicit none

contains

  subroutine Inquire_file(file_unit,file_name,file_format,file_access)
      
    implicit none

    logical                          :: lexists
    integer                          :: ios
    integer, intent(IN)              :: file_unit             
    character( len = * ), intent(IN) :: file_format, file_access, file_name
    type(error)                      :: err

    ios = 0

    inquire ( file = trim(file_name), exist = lexists, iostat = ios )
      
    if ( ios == 0 ) then 
        
      if ( lexists ) then 

        open( unit = file_unit, file = trim(file_name), status = 'old', &
              form = trim(file_format), access = trim(file_access) ) 

        rewind(file_unit)

      else

        call err%error('e',message="while reading file: "//trim(file_name))
        call err%error('e',check="if file exists in this directory.")

        stop

      endif

    else

      call err%error('e',message="while reading file: "//trim(file_name))
      
      stop

    endif

    return
  end subroutine Inquire_file

end module mod_inquire
