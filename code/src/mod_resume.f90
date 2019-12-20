!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief This module contains routines that resume the run.         
!> @author Felippe M. Colombari                                                
!> - Laboratório de Química Teórica, LQT -- UFSCar                             
!> @date - Jun, 2017                                                           
!> - independent module created                                                
!------------------------------------------------------------------------------

module MOD_RESUME
  use MOD_CONSTANTS

  implicit none 

  real( kind = DP ), allocatable, dimension(:) :: A_min(:)
  integer, dimension(:)                        :: pos_A_min(1)

  contains

    !---------------------------------------------------------------------------
    !> @brief This routine writes details of the calculation on screen.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Jun, 2017
    !> - subroutine  created
    !---------------------------------------------------------------------------	

    subroutine ENDING( timet )
      use MOD_CMD_LINE,   only: grid_type, rad, irun
      use MOD_INPUT_READ, only: potential, atom_overlap
      use MOD_READ_GRIDS, only: A
      use MOD_LOOPS,      only: min_ener
      use MOD_SEARCH,     only: n

      implicit none

      real( kind = DP ), intent(IN)                :: timet

      select case( grid_type )

        case( "shell" )

          write(*,'(/, T5, A, T93, f6.2, A)') "Shell mapping done for r: ", rad, " A" 

        case ( "user" )

          write(*,'(/, T5, A, T93, f6.2, A)') "Grid mapping done."

        case default

          continue

      end select

      if ( irun == "run" ) then

        write(*,'(/,T5,A,T92,i9,/)') "Skipped configurations:", count(atom_overlap)

      else

        write(*,*)

      endif

      if ( potential /= "none" ) then

        write(*,'(T5,a34,T88,es13.5e3,/)') "Lowest energy structure (kJ/mol): ", min_ener
        write(*,'(T5,a39,T94, i7,/)') "Number of structures within 0.5 * kBT: ", n

        write(*,'(T5,a30,T96,i5,/)') "Lowest free energy grid point: ", minloc(A)
        write(*,'(T9,a28,T88,es13.5E3,/)') "Free energy value (kJ/mol): ", minval(A)
        write(*,'(T5,a31,T96,i5,/)') "Highest free energy grid point: ", maxloc(A)
        write(*,'(T9,a28,T88,es13.5E3,/)') "Free energy value (kJ/mol): ", maxval(A)

      endif

      write(*,'(T3, A,/)') dashline
      write(*,'(T5,a20,T83,f10.2,a8,/)') "Total running time: ", timet, " seconds."
      write(*,'(T3, A,/)') dashline

      return
    end subroutine ENDING

    !---------------------------------------------------------------------------
    !> @brief This routine sorts central grid points according to their free 
    !>  energy values.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Aug, 2017
    !> - subroutine  created
    !> @note \n 
    !> - The method used to sort the free energy array is non-standard, but is 
    !>   much more easier to implement and faster to run: \n
    !>   i)   The lowest free energy value is found on array A and saved at first 
    !>        position of array B;
    !>   ii)  In array A, such position is replaced by INFINITY; \n
    !>   iii) Now the lowest free energy value from array A will be the second 
    !>        value from array B and so on...
    !---------------------------------------------------------------------------	
  
    subroutine SORT_OUTPUT
      use MOD_READ_GRIDS
      use MOD_LOOPS, only: ATOTAL, mTSTOTAL, ETOTALavg, Ztrans, kBT, min_ener
      use mod_error_handling
    
      implicit none

      integer                           :: t
      real( kind = DP )                 :: log_Z
      character( len = : ), allocatable :: write_fmt

      integer                           :: ierr
      type(error)                       :: err

      open(unit=21,file='output-sort.log',position='append',status='replace')

      write_fmt = '(a1,2x,3(a11,1x),3x,a5,3x,4(a13,3x))'

      write(21,*) grid_trans % numpoint

      write(21, write_fmt ) " ", "X (A)", "Y (A)", "Z (A)", "point", "PROB", &
        &"A (kJ/mol)", "-TS (kJ/mol)", "E (kJ/mol)"

      allocate( A_min( grid_trans % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      A_min = 0.0_dp

      write_fmt = '(A,2x,3(f11.5,1x),3x,i5,3x,4(es13.5E3,3x))'
      
      do t = 1, grid_trans % numpoint

        A_min(t) = minval(A) 

        pos_A_min = minloc(A)

        A(pos_A_min) = FPINF

        write( 21, write_fmt ) "X", grid_trans % points(pos_A_min) % grid_xyz(1), &
                                    grid_trans % points(pos_A_min) % grid_xyz(2), & 
                                    grid_trans % points(pos_A_min) % grid_xyz(3), &
                                    pos_A_min, probT(pos_A_min), A_min(t), mTS(pos_A_min), Eavg(pos_A_min)

      enddo

      write(21,'(a)') dashline

      log_Z = dlog(Ztrans) - min_ener / kBT

      write_fmt = '(9x,"TOTAL OVER TRANSLATIONAL GRID",12x,4(es13.5E3,3x))'
      
      write(21, write_fmt) sum(probT), ATOTAL, mTSTOTAL, ETOTALavg

      close(21)

      deallocate( A_min )
      deallocate( A )
      deallocate( mTS )
      deallocate( Eavg )

      return
    end subroutine SORT_OUTPUT

end module MOD_RESUME
