!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2023 Themis developers
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
!> @file   mod_cmd_line.f90
!> @author Felippe M. Colombari
!> @brief  Get command line arguments         
!> @date - Nov, 2017                                                           
!> - independent module created                                                
!> @date - Jan, 2018                                                           
!> - parsing tests added  
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Jul, 2022 
!> - help, license and citation routines moved to mod_info
!---------------------------------------------------------------------------------------------------

module mod_cmd_line
  use iso_fortran_env    , only : output_unit
  use mod_constants      , only : DP, float_alphabet, char_alphabet, dashline, version
  use mod_error_handling
  use mod_info

  implicit none

  private
  public Parse_Arguments, rad, irun, grid_type, grid_transl, use_plot

  real( kind = DP )                                :: rad         = 0.0_DP
  character( len = 10 )                            :: irun        = char(0)
  character( len = 10 )                            :: grid_type   = char(0)
  character( len = 40 )                            :: grid_transl = char(0) 
  character( len = 20 ), allocatable, dimension(:) :: arg         
  logical                                          :: use_plot    = .false.

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
    integer                :: ios      = 0
    integer                :: narg     = 0
    integer                :: nochar   = 0
    character( len = 256 ) :: cmd_line = char(0)
      
    narg = command_argument_count()
    call get_command(cmd_line)

    write( output_unit, '(/,T5, A, A)' ) "COMMAND LINE READ: ", trim(cmd_line)

    if ( narg > 0 ) then
    
      ! to avoid allocation errors if one forget the argument "rad"

      allocate( arg(narg+1), stat=ierr )
      if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      arg = char(0)

      do i = 1, narg

        call get_command_argument( i, arg(i) )

        if ( arg(i)(1:2) == '--' ) then

          select case( arg(i) )

            case( '--help' )

              call Display_help

              if ( allocated(arg) ) deallocate(arg)

              call err%termination(0,'f')

            case( '--license' )

              call Display_license
              
              if ( allocated(arg) ) deallocate(arg)

              call err%termination(0,'f')

            case( '--citation' )

              call Display_citation
              
              if ( allocated(arg) ) deallocate(arg)

              call err%termination(0,'f')

            case( '--rerun' )
            
              irun = "rerun"

              call Get_command_argument( i+1, arg(i+1) )

              if ( arg(i+1)(1:2) /= '--' ) then

                call err % error( 'e', message = "while reading command line." )
                call err % error( 'e', check   = "'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag." )
               
                write( output_unit, '(/,T3, A)' ) dashline

                stop

              endif

            case( '--run' )

              irun = "run"

              call Get_command_argument( i+1, arg(i+1) )

              if ( len(trim(arg(i+1))) > 0 ) then

                if ( arg(i+1)(1:2) /= '--' ) then

                  call err % error( 'e', message = "while reading command line." )
                  call err % error( 'e', check   = "'"//trim(adjustl(arg(i+1)))//"' argument of '"//&
                                                        trim(adjustl(arg(i)))//"' flag." )

                  write( output_unit, '(/,T3, A)' ) dashline

                  stop

                endif

              endif

            case( '--shell' )

              grid_type = "shell"

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), float_alphabet )

              if ( nochar > 0 ) then

                call err % error( 'e', message = "while reading command line." )
                call err % error( 'e', check   = "spherical grid for translation." ) 
                call err % error( 'e', tip     = "Its radius value (in Angstrom) should be > 1.0." )
               
                write( output_unit, '(/,T3, A)' ) dashline

                stop
                
              else

                read(arg(i+1),*,iostat=ios) rad

                if ( ios > 0 ) then

                  call err % error( 'e', message = "while reading command line." )
                  call err % error( 'e', check   = "spherical grid for translation." ) 
                  call err % error( 'e', tip     = "Its radius value (in Angstrom) should be > 1.0." )
               
                  write( output_unit, '(/,T3, A)' ) dashline

                  stop

                endif

                if ( rad <= 1.0_DP ) then

                  call err % error( 'e', message = "while reading command line." )
                  call err % error( 'e', check   = "spherical grid for translation." ) 
                  call err % error( 'e', tip     = "Its radius value (in Angstrom) should be > 1.0." )
               
                  write( output_unit, '(/,T3, A)' ) dashline

                  stop

                endif

              endif

            case('--user')

              grid_type = "user"

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), char_alphabet )

              if ( nochar > 0 ) then

                call err % error( 'e', message = "while reading command line." )
                call err % error( 'e', check   = "invalid filename '"//trim(adjustl(arg(i+1)))//"'." )
               
                write( output_unit, '(/,T3, A)' ) dashline

                stop

              else

                read(arg(i+1),*,iostat=ios) grid_transl

                if ( ios > 0 ) then

                  call err % error( 'e', message = "while reading command line." )
                  call err % error( 'e', check   = "invalid filename '"//trim(adjustl(arg(i+1)))//"'." )
               
                  write( output_unit, '(/,T3, A)' ) dashline

                  stop

                endif

              endif

              
            case('--plot')

              use_plot = .true.

            case default

              call err % error( 'e', message = "while reading command line." )
              call err % error( 'e', check   = "invalid command line flag '"//trim(adjustl(arg(i)))//"'." )
               
              write( output_unit, '(/,T3, A)' ) dashline

              stop

          end select

        else 

          if ( arg(1)(1:2) /= '--' ) then

            call err % error( 'e', message = "while reading command line." )
            call err % error( 'e', check   = "'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag." )
 
            write( output_unit, '(/,T3, A)' ) dashline
            
            stop

          endif

          if ( ( i > 1 ) .and. ( arg(i-1)(1:2) ) /= '--' ) then

            call err % error( 'e', message = "while reading command line." )
            call err % error( 'e', check   = "'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag." )

            write( output_unit, '(/,T3, A)' ) dashline
            
            stop

          endif

        endif

      enddo

      if ( allocated(arg) ) deallocate(arg)

    else if ( narg == 0 ) then 

      call err % error( 'e', message = "while reading command line." )
      call err % error( 'e', tip     = "Command line arguments are missing. Use 'themis --help' for options." )

      write( output_unit, '(/,T3, A)' ) dashline

      stop

    endif

    if ( irun == char(0) ) then

      call err % error( 'e', message = "while reading command line." )
      call err % error( 'e', check   = "IRUN option." )
      call err % error( 'e', tip     = "options are --run or --rerun." )

      write( output_unit, '(/,T3, A)' ) dashline
      
      stop

    else if ( grid_type == char(0) ) then

      call err % error( 'e', message = "while reading command line." )
      call err % error( 'e', check   = "GRID_TYPE option." )
      call err % error( 'e', tip     = "options are --shell or --user." )

      write( output_unit, '(/,T3, A)' ) dashline
      
      stop

    endif

  end subroutine Parse_arguments

end module mod_cmd_line
