!---------------------------------------------------------------------------------------------------
! COM: A code to calculate the center of mass of a given molecular structure                                                  
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
!   This file was written by Felippe M. Colombari.
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
!> @file   mod_cmd_line.f90
!> @author Felippe M. Colombari
!> @brief  Get command line arguments         
!> @date - Jan, 2020                                                           
!> - independent module created                                                
!> @date - Jan, 2020
!> - update error condition by error_handling module added
!---------------------------------------------------------------------------------------------------

module mod_cmd_line
  use iso_fortran_env    , only : stdout => output_unit
  use mod_constants      , only : DP, char_alphabet, dashline
  use mod_error_handling

  implicit none

  private
  public Parse_Arguments, filename_molecule, filename_com, center_option

  character( len = 64 )                            :: filename_molecule = char(0)
  character( len = 64 )                            :: filename_com      = char(0)
  character( len = 16 )                            :: center_option     = char(0)
  character( len = 64 ), allocatable, dimension(:) :: arg         

  integer                                          :: ierr
  type(error)                                      :: err

contains
  
  !---------------------------------------------------------------------------
  !> @brief Parses the command line arguments
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Parse_arguments

    implicit none

    integer                :: i
    integer                :: ios            = 0
    integer                :: narg           = 0
    integer                :: nochar         = 0
    character( len = 256 ) :: cmd_line       = char(0)
      
    narg = command_argument_count()
    call get_command(cmd_line)

    write(stdout,'(/,T5, A, A)') "COMMAND LINE READ: ", trim(cmd_line)
    write(stdout,'(/,T3, A)') dashline

    if ( narg > 0 ) then
    
      ! to avoid allocation errors if one forget the argument "rad"

      allocate( arg(narg+1), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      arg = char(0)

      do i = 1, narg

        call get_command_argument( i, arg(i) )

        if ( arg(i)(1:2) == '--' ) then

          SELECT CASE( arg(i) )

            CASE( '--help' )

              call Display_help

            CASE( '--license' )

              call Display_license

            CASE( '--version')  

              call display_version

            CASE( '--input' )

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), char_alphabet )

              if ( nochar > 0 ) then

                call err%error('e',message="while reading command line.")
                call err%error('e',check="molecule coordinate file.") 
                call err%error('e',tip="Should be a valid .xyz file.")
               
                stop
                
              else

                read(arg(i+1),*,iostat=ios) filename_molecule

                if ( ios > 0 ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="molecule coordinate file.") 
                  call err%error('e',tip="Should be a valid .xyz file.")
               
                  stop

                endif

              endif

            CASE( '--output' )

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), char_alphabet )

              if ( nochar > 0 ) then

                call err%error('e',message="while reading command line.")
                call err%error('e',check="molecule coordinate file.") 
                call err%error('e',tip="Should be a valid .xyz file.")
               
                stop
                
              else

                read(arg(i+1),*,iostat=ios) filename_com

                if ( ios > 0 ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="molecule coordinate file.") 
                  call err%error('e',tip="Should be a valid .xyz file.")
               
                  stop

                endif

              endif

            CASE( '--center' )

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), char_alphabet )

              if ( nochar > 0 ) then

                call err%error('e',message="while reading command line.")
                call err%error('e',check="centering option.") 
                call err%error('e',tip="Should be TRUE or FALSE.")
               
                stop
                
              else

                read(arg(i+1),*,iostat=ios) center_option

                if ( ios > 0 ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="centering option.") 
                  call err%error('e',tip="Should be TRUE or FALSE.")
               
                  stop

                endif

                if ( ( center_option == 'TRUE' ) .or. ( center_option == 'FALSE' ) ) then

                  continue

                else
                  
                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="centering option.") 
                  call err%error('e',tip="Should be TRUE or FALSE.")
               
                  stop

                endif

              endif

            CASE DEFAULT

              call err%error('e',message="while reading command line.")
              call err%error('e',check="invalid command line flag '"//trim(adjustl(arg(i)))//"'.")
               
              stop

          end SELECT

        else 

          if ( arg(1)(1:2) /= '--' ) then

            call err%error('e',message="while reading command line.")
            call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

            stop

          endif

          if ( ( i > 1 ) .and. ( arg(i-1)(1:2) ) /= '--' ) then

            call err%error('e',message="while reading command line.")
            call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

            stop

          endif

        endif

      enddo

      if ( allocated(arg) ) deallocate(arg)

    else if ( narg == 0 ) then 

      call err%error('e',message="while reading command line.")
      call err%error('e',tip="Command line arguments are missing.")

      stop

    endif

    if ( filename_molecule == char(0) ) then

      call err%error('e',message="while reading command line.")
      call err%error('e',check="molecule coordinate file.") 
      call err%error('e',tip="Should be a valid .xyz file.")

      stop

    else if ( filename_com == char(0) ) then

      call err%error('e',message="while reading command line.")
      call err%error('e',check="output coordinate file.") 
      call err%error('e',tip="Should be a valid .xyz file.")

      stop

    else if ( center_option == char(0) ) then

      call err%error('e',message="while reading command line.")
      call err%error('e',check="centering option.") 
      call err%error('e',tip="Should be TRUE or FALSE.")

      stop

    endif

  end subroutine Parse_arguments

  !---------------------------------------------------------------------------
  !> @brief Displays command line options
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_help
            
    implicit none

    write(stdout,'(/,T20, A)')'Usage:  com --input [FILEIN] --output [FILEOUT] --center [TRUE/FALSE]     '
    write(stdout,'(/,T3, A)') dashline
    write(stdout,'(/,T25, A)')'[FILEIN]     : initial .xyz coordinate file.'
    write(stdout,'(/,T25, A)')'[FILEOUT]    : final .xyz coordinate file.'
    write(stdout,'(/,T25, A)')'[TRUE/FALSE] : place the COM at (0,0,0) ?' 
    write(stdout,'(/,T3, A)') dashline
    
    if ( allocated(arg) ) deallocate(arg)

    call err%termination(0,'f')

  end subroutine Display_help

  !---------------------------------------------------------------------------
  !> @brief Displays the license
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_license

    implicit none

    write(stdout,'(T36, A)')'Copyright 2020 Felippe M. Colombari'
    write(stdout,'(/,T33, A)')'License GPLv3+: GNU GPL version 3 or later' 
    write(stdout,'(/,T6, A)')' This program is free software: you can redistribute it and/or modify it &
                      &under the terms of the'
    write(stdout,'(T5, A)')'GNU General Public License as published by the Free Software Foundation, &
                      &either version 3 of the'
    write(stdout,'(T30, A)')'License, or (at your option) any later version.'
    write(stdout,'(/,T5, A)')'This program is distributed in the hope that it will be useful, but &
                      &WITHOUT ANY WARRANTY; without'
    write(stdout,'(T12, A)')'even the implied warranty of MERCHANTABILITY or FITNESS FOR A &
                      &PARTICULAR PURPOSE.'
    write(stdout,'(T26, A)')'See the GNU General Public License for more details.'
    write(stdout,'(/,T4, A)')'You should have received a copy of the GNU General Public License along & 
                      &with this program. If not,'
    write(stdout,'(T34, A)')'see <https://www.gnu.org/licenses/>.'
    write(stdout,'(/,T36, A)')'E-mail: colombarifm@hotmail.com'
    write(stdout,'(/,T3, A,/)') dashline

    if ( allocated(arg) ) deallocate(arg)

    call err%termination(0,'f')

  end subroutine Display_license

  !---------------------------------------------------------------------------
  !> @brief Displays the version
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_version()

    implicit none
    
    character(len=:), allocatable :: version
    
    ! TODO link to version control 

    version = '1.0.0'

    write(stdout,'("com ",a)') version
     
    call err%termination(0,'f')

  end subroutine Display_version

end module mod_cmd_line
