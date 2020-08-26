!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Themis developers
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
!---------------------------------------------------------------------------------------------------

module mod_cmd_line
  use iso_fortran_env    , only: output_unit
  use mod_constants      , only : DP, float_alphabet, char_alphabet, dashline, version
  use mod_error_handling
  use mod_info

  implicit none

  private
  public Parse_Arguments, rad, irun, grid_type, grid_transl

  real( kind = DP )                                :: rad         = 0.0_DP
  character( len = 10 )                            :: irun        = char(0)
  character( len = 10 )                            :: grid_type   = char(0)
  character( len = 40 )                            :: grid_transl = char(0) 
  character( len = 20 ), allocatable, dimension(:) :: arg         

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

    write(output_unit,'(/,T5, A, A)') "COMMAND LINE READ: ", trim(cmd_line)

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

            CASE( '--citation' )

              call Display_citation

            CASE( '--rerun' )
            
              irun = "rerun"

              call Get_command_argument( i+1, arg(i+1) )

              if ( arg(i+1)(1:2) /= '--' ) then

                call err%error('e',message="while reading command line.")
                call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")
               
                write(output_unit,'(/,T3, A)') dashline

                stop

              endif

            CASE( '--run' )

              irun = "run"

              call Get_command_argument( i+1, arg(i+1) )

              if ( len(trim(arg(i+1))) > 0 ) then

                if ( arg(i+1)(1:2) /= '--' ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

                  write(output_unit,'(/,T3, A)') dashline

                  stop

                endif

              endif

            CASE( '--shell' )

              grid_type = "shell"

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), float_alphabet )

              if ( nochar > 0 ) then

                call err%error('e',message="while reading command line.")
                call err%error('e',check="spherical grid for translation.") 
                call err%error('e',tip="Its radius value (in Angstrom) should be > 1.0.")
               
                write(output_unit,'(/,T3, A)') dashline

                stop
                
              else

                read(arg(i+1),*,iostat=ios) rad

                if ( ios > 0 ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="spherical grid for translation.") 
                  call err%error('e',tip="Its radius value (in Angstrom) should be > 1.0.")
               
                  write(output_unit,'(/,T3, A)') dashline

                  stop

                endif

                if ( rad <= 1.0_DP ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="spherical grid for translation.") 
                  call err%error('e',tip="Its radius value (in Angstrom) should be > 1.0.")
               
                  write(output_unit,'(/,T3, A)') dashline

                  stop

                endif

              endif

            CASE('--user')

              grid_type = "user"

              call Get_command_argument( i+1, arg(i+1) )

              nochar = verify( trim( arg(i+1) ), char_alphabet )

              if ( nochar > 0 ) then

                call err%error('e',message="while reading command line.")
                call err%error('e',check="invalid filename '"//trim(adjustl(arg(i+1)))//"'.")
               
                write(output_unit,'(/,T3, A)') dashline

                stop

              else

                read(arg(i+1),*,iostat=ios) grid_transl

                if ( ios > 0 ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="invalid filename '"//trim(adjustl(arg(i+1)))//"'.")
               
                  write(output_unit,'(/,T3, A)') dashline

                  stop

                endif

              endif

            CASE DEFAULT

              call err%error('e',message="while reading command line.")
              call err%error('e',check="invalid command line flag '"//trim(adjustl(arg(i)))//"'.")
               
              write(output_unit,'(/,T3, A)') dashline

              stop

          end SELECT

        else 

          if ( arg(1)(1:2) /= '--' ) then

            call err%error('e',message="while reading command line.")
            call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")
 
            write(output_unit,'(/,T3, A)') dashline
            
            stop

          endif

          if ( ( i > 1 ) .and. ( arg(i-1)(1:2) ) /= '--' ) then

            call err%error('e',message="while reading command line.")
            call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

            write(output_unit,'(/,T3, A)') dashline
            
            stop

          endif

        endif

      enddo

      if ( allocated(arg) ) deallocate(arg)

    else if ( narg == 0 ) then 

      call err%error('e',message="while reading command line.")
      call err%error('e',tip="Command line arguments are missing. Use 'themis --help' for options.")

      write(output_unit,'(/,T3, A)') dashline

      stop

    endif

    if ( irun == char(0) ) then

      call err%error('e',message="while reading command line.")
      call err%error('e',check="IRUN option.")
      call err%error('e',tip="options are --run or --rerun.")

      write(output_unit,'(/,T3, A)') dashline
      
      stop

    else if ( grid_type == char(0) ) then

      call err%error('e',message="while reading command line.")
      call err%error('e',check="GRID_TYPE option.")
      call err%error('e',tip="options are --shell or --user.")

      write(output_unit,'(/,T3, A)') dashline
      
      stop

    endif

  end subroutine Parse_arguments

  !---------------------------------------------------------------------------
  !> @brief Displays command line options
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_help
            
    implicit none

    write(output_unit,'(/,T3,A)')  dashline
    write(output_unit,'(/,T36,A)') 'Usage:  themis [RUNTYPE] [GRID]     '
      
    write(output_unit,'(/,T3,A)')   dashline
    write(output_unit,'(/,T30,A)') '[RUNTYPE] options'       
    write(output_unit,'(/,T8,A18,4x,A)') adjustr('--run'), &
                                        &adjustl('Start new calculation.')
    write(output_unit,'(/,T8,A18,4x,A)') adjustr('--rerun'), &
                                        &adjustl('Calculate properties from:') 
    write(output_unit,'(T10,A10,12x,A)') adjustr(''), &
                                        &adjustl('energy.bin - energies obtained from Themis run;')
    write(output_unit,'(T10,A10,12x,A)') adjustr(''), &
                                        &adjustl('energy.log - energies obtained from external calculation.')
    
    write(output_unit,'(/,T3,A)')  dashline
    write(output_unit,'(/,T30,A)') '[GRID] options'       
    write(output_unit,'(/,T8,A18,4x,A)') adjustr('--shell <radius>'), &
                                        &adjustl('Translation moves will be performed on a spherical shell generated')
    write(output_unit,'(T8,A10,12x,A)')  adjustr(''), &
                                        &adjustl('on the run. The real argument <radius> (in Angstrom) scales its radius.') 

    write(output_unit,'(/,T8,A18,4x,A)') adjustr('--user <file.xyz>'), &
                                        &adjustl('Translation moves will be performed on an user-defined grid read from')
    write(output_unit,'(T8,A10,12x,A)')  adjustr(''), &
                                        &adjustl('<file.xyz> that must be aligned with molecule 1.')
      
    write(output_unit,'(/,T3,A,/)') dashline
    write(output_unit,'(T8,A18,4x,A)') adjustr('--help'), &
                                      &adjustl('Display this help')
    write(output_unit,'(T8,A18,4x,A)') adjustr('--license'), &
                                      &adjustl('Display the license')
    write(output_unit,'(T8,A18,4x,A)') adjustr('--citation'), &
                                      &adjustl('Display Themis papers')
    write(output_unit,'(/,T3,A)') dashline

    write(output_unit,'(/,T3, A)') '   Report bugs to:'
    write(output_unit,'(T10, A)')  '   Felippe M. Colombari   - colombarifm@hotmail.com'
    write(output_unit,'(T10, A)')  '   Asdrubal Lozada-Blanco - aslozada@gmail.com'

    write(output_unit,'(/,T3, A,/)') dashline
    
    if ( allocated(arg) ) deallocate(arg)

    call err%termination(0,'f')

  end subroutine Display_help

  !---------------------------------------------------------------------------
  !> @brief Displays the license
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Display_license

    implicit none
   
    write(output_unit,'(/,T3,A,/)') dashline
    write(output_unit,'(T36,A)')    'Copyright 2020 Felippe M. Colombari'
    write(output_unit,'(/,T33,A)')  'License GPLv3+: GNU GPL version 3 or later' 
    write(output_unit,'(/,T6,A)')   ' This program is free software: you can redistribute it and/or modify it &
                                   &under the terms of the'
    write(output_unit,'(T5,A)')     'GNU General Public License as published by the Free Software Foundation, &
                                   &either version 3 of the'
    write(output_unit,'(T30,A)')   'License, or (at your option) any later version.'
    write(output_unit,'(/,T5,A)')  'This program is distributed in the hope that it will be useful, but &
                                   &WITHOUT ANY WARRANTY; without'
    write(output_unit,'(T12,A)')   'even the implied warranty of MERCHANTABILITY or FITNESS FOR A &
                                   &PARTICULAR PURPOSE.'
    write(output_unit,'(T26,A)')   'See the GNU General Public License for more details.'
    write(output_unit,'(/,T4,A)')  'You should have received a copy of the GNU General Public License along & 
                                   &with this program. If not,'
    write(output_unit,'(T34,A)')   'see <https://www.gnu.org/licenses/>.'
    write(output_unit,'(/,T3,A,/)') dashline

    if ( allocated(arg) ) deallocate(arg)

    call err%termination(0,'f')

  end subroutine Display_license

  subroutine Display_citation

    implicit none

    write(output_unit,'(/,T3, A)') dashline

    write(output_unit,'(/,T5, A)') '***** This version of THEMIS uses:'
    write(output_unit,'(/,T5, A)') 'SPHERE_GRID library (John Burkardt)'
    write(output_unit,'(T5, A)')   'https://people.sc.fsu.edu/~jburkardt/f_src/sphere_grid/sphere_grid.html'
    write(output_unit,'(/,T5, A)') 'XDR library (James Barnett)'
    write(output_unit,'(T5, A)')   'https://github.com/kmtu/xdrfort/blob/master/xdr.F90'

    write(output_unit,'(/,T5, A)') '***** Detailed information about the methodology:'
    write(output_unit,'(/,T5, A)') 'Themis:'
    write(output_unit,'(T5, A)')   'a software to assess association free energies via direct estimative of partition functions.'
    write(output_unit,'(T5, A)')   'ref..., doi....'

    write(output_unit,'(/,T5, A)') '***** Papers that used earlier or adapted Themis versions:'
    write(output_unit,'(/,T5, A)') 'Emergence of complexity in hierarchically organized chiral particles.' 
    write(output_unit,'(T5, A)')   'Science, v. 368, p. 642, 2020, doi: 10.1126/science.aaz7949' 

    write(output_unit,'(/,T5, A)') 'Solvent effect on the regulation of urea hydrolysis reactions by copper complexes.'
    write(output_unit,'(T5, A)')   'Chemistry, v. 2, p. 525, 2020, doi: 10.3390/chemistry2020032'

    write(output_unit,'(/,T5, A)') 'Graphitic Carbon Nitrides as Platforms for Single-Atom Photocatalysis.'
    write(output_unit,'(T5, A)')   'Faraday Discussions, v. xxx, p. xxx, 2020, doi: 10.1039/C9FD00112C'

    write(output_unit,'(/,T5, A)') 'Ion pair free energy surface as a probe of ionic liquid structure.'
    write(output_unit,'(T5, A)')   'The Journal of Chemical Physics, v. 152, p. 014103, 2020, doi: 10.1063/1.5128693.'

    write(output_unit,'(/,T5, A)') 'Low-Temperature Phase Transitions of the Ionic Liquid 1-Ethyl-3-methylimidazolium Dicyanamide'
    write(output_unit,'(T5, A)')   'The Journal of Physical Chemistry B, v. 123, p. 9418, 2019, doi: 10.1021/acs.jpcb.9b07654.'

    write(output_unit,'(/,T5, A)') 'Site-selective photoinduced cleavage and profiling of DNA by chiral semiconductor&
                                    & nanoparticles.'
    write(output_unit,'(T5, A)')   'Nature Chemistry, v. 10, p. 821, 2018, doi: 10.1038/s41557-018-0083-y'
    
    write(output_unit,'(/,T3, A,/)') dashline

    if ( allocated(arg) ) deallocate(arg)

    call err%termination(0,'f')

  end subroutine Display_citation

end module mod_cmd_line
