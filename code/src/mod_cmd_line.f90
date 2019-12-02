!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief This module contains a routine to get command line arguments         
!> @author Felippe M. Colombari                                                
!> - Laboratório de Química Teórica, LQT -- UFSCar                             
!> @date - Jun, 2017                                                           
!> - independent module created                                                
!------------------------------------------------------------------------------

module mod_cmd_line
  use mod_constants, only: dp, float_alphabet, char_alphabet, dashline
  use mod_error_handling

  implicit none

  private
  public Parse_Arguments, rad, irun, grid_type, grid_transl

  real( kind = DP )                                :: rad         = 0.0_dp
  character( len = 10 )                            :: irun        = char(0)
  character( len = 10 )                            :: grid_type   = char(0)
  character( len = 40 )                            :: grid_transl = char(0) 
  character( len = 20 ), allocatable, dimension(:) :: arg         

  integer                                          :: ierr
  type(error)                                      :: err

  contains
  
    !---------------------------------------------------------------------------
    !> @brief This routine parses the command line arguments
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Jun, 2017
    !> - subroutine  created
    !> @date - Dec, 2017
    !> - parsing tests added
    !> @note added error_handling
    !> Asdrubal Lozada-Blanco
    !> @date - Nov 2019
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

      write(*,*)
      write(*,'(T5, A, A)') "COMMAND LINE READ: ", trim(cmd_line)
      write(*,*)
      write(*,'(T3, A)') dashline

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

              CASE( '--rerun' )
            
                irun = "rerun"

                call Get_command_argument( i+1, arg(i+1) )

                if ( arg(i+1)(1:2) /= '--' ) then

                  call err%error('e',message="while reading command line.")
                  call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

                  stop

                endif

              CASE( '--run' )

                irun = "run"

                call Get_command_argument( i+1, arg(i+1) )

                if ( len(trim(arg(i+1))) > 0 ) then

                  if ( arg(i+1)(1:2) /= '--' ) then

                    call err%error('e',message="while reading command line.")
                    call err%error('e',check="'"//trim(adjustl(arg(i+1)))//"' argument of '"//trim(adjustl(arg(i)))//"' flag.")

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
               
                  stop
                
                else

                  read(arg(i+1),*,iostat=ios) rad

                  if ( ios > 0 ) then

                    call err%error('e',message="while reading command line.")
                    call err%error('e',check="spherical grid for translation.") 
                    call err%error('e',tip="Its radius value (in Angstrom) should be > 1.0.")
               
                    stop

                  endif

                  if ( rad <= 1.0_DP ) then

                    call err%error('e',message="while reading command line.")
                    call err%error('e',check="spherical grid for translation.") 
                    call err%error('e',tip="Its radius value (in Angstrom) should be > 1.0.")
               
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
               
                  stop

                else

                  read(arg(i+1),*,iostat=ios) grid_transl

                  if ( ios > 0 ) then

                    call err%error('e',message="while reading command line.")
                    call err%error('e',check="invalid filename '"//trim(adjustl(arg(i+1)))//"'.")
               
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

      if ( irun == char(0) ) then

        call err%error('e',message="while reading command line.")
        call err%error('e',check="IRUN option.")
        call err%error('e',tip="options are --run or --rerun.")

        stop

      else if ( grid_type == char(0) ) then

        call err%error('e',message="while reading command line.")
        call err%error('e',check="GRID_TYPE option.")
        call err%error('e',tip="options are --shell or --user.")

        stop

      endif

    end subroutine Parse_arguments

    subroutine Display_help
            
      implicit none

!      character(10) :: numbers = "1234567890"

      write(*,*)
      write(*,'(T10, A)')'   Usage:  themis [RUNTYPE] [GRIDTYPE] [RADIUS|FILENAME]     '
      
      write(*,*)
      write(*,'(T3, A)') dashline
      write(*,*) 
      write(*,'(T25, A)')'      [RUNTYPE] options'       
      write(*,*) 
      write(*,'(T10, A10, 05x, A)') adjustr('--run'), &
                                   &adjustl('Start new calculation.')
      write(*,*)
      write(*,'(T10, A10, 05x, A)') adjustr('--rerun'), &
                                   &adjustl('Read energies from previous run. User must') 
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('give an energy.bin (from Themis) or an')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('energy.log file (from external programs).')
      
      write(*,*)
      write(*,'(T3, A)') dashline
      write(*,*) 
      write(*,'(T25, A)')'     [GRIDTYPE] options'       
      write(*,*) 
      write(*,'(T10, A10, 05x, A)') adjustr('--shell'), &
                                   &adjustl('Translation moves will be performed on a')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('spherical shell grid with a given RADIUS')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('generated on the run.')

      write(*,*)
      write(*,'(T10, A10, 05x, A)') adjustr('--user'), &
                                   & adjustl('Translation moves will be performed on an')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('user-defined grid read from FILENAME.')
      
      write(*,*)
      write(*,'(T3, A)') dashline
      write(*,*) 
      write(*,'(T25, A)')'       [RADIUS] value'       
      write(*,*) 
      write(*,'(T10, A10, 05x, A)') adjustr('(real)'), &
                                   &adjustl('Scaling factor for the spherical grid')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('radius (in Angstrom)')
      
      write(*,*) 
      write(*,'(T25, A)')'     [FILENAME] value'       
      write(*,*) 
      write(*,'(T10, A10, 05x, A)') adjustr('(char)'), &
                                   &adjustl('XYZ file containing the user-defined')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('translation grid. It must be aligned')
      write(*,'(T10, A10, 05x, A)') adjustr(''), &
                                   &adjustl('with molecule 1')
      
      write(*,*)
      write(*,'(T3, A)') dashline
      write(*,*) 
!      write(*,'(T3, A)') numbers//numbers//numbers//numbers//numbers//numbers//numbers//"1234"
      write(*,'(T25, A)')'        Other options'       
      write(*,*)
      write(*,'(T10, A10, 05x, A)') adjustr('--help'), &
                                   &adjustl('Display this help')
      write(*,*)
      write(*,'(T10, A10, 05x, A)') adjustr('--version'), &
                                   &adjustl('Display the version')
      write(*,*)
      write(*,'(T3, A)') dashline
      write(*,*) 

      if ( allocated(arg) ) deallocate(arg)

      call err%termination(0,'f')

    end subroutine Display_help

    subroutine Display_license

      implicit none
    
      write(*,'(T15, A)')'  Copyright (C) 2019 Felippe Mariano Colombari       '        
      write(*,*) 
      write(*,'(T15, A)')'                 License GPLv3+:                     '
      write(*,'(T15, A)')'           GNU GPL version 3 or later                ' 
      write(*,'(T15, A)')'      see <http://gnu.org/license/gpl.html>          '
      write(*,*)
      write(*,'(T15, A)')'             This is a free software                 '
      write(*,'(T15, A)')' you are free to change it and redistributibe it     '
      write(*,'(T13, A)')' There is NO WARRANTY, to the extent permited by law '
      write(*,*)
      write(*,'(T15, A)')'         Written by Felippe M. Colombari             '
      write(*,'(T15, A)')'        E-mail: colombarifm@hotmail.com              '
      write(*,*)
      write(*,'(T3, A)') dashline
      write(*,*)

      if ( allocated(arg) ) deallocate(arg)

      call err%termination(0,'f')

    end subroutine Display_license

    subroutine display_version()

       implicit none
       character(len=:), allocatable :: version


       ! TODO link to version control 

       version = '1.0.0-beta'

       write(*,'("themis ",a)') version
      
       call err%termination(0,'f')

    end subroutine display_version

end module mod_cmd_line
