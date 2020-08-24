!---------------------------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition function estimation                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Themis developers
!                 Laboratory of Theoretical Chemistry (LQT) - Federal University of São Carlos 
!                 <http://www.lqt.dq.ufscar.br>
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
!> @file   mod_resume.f90
!> @author Felippe M. Colombari
!> @brief  This module contains routines that resume the run.         
!> @date - Jun, 2017                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_resume
  use iso_fortran_env, only: output_unit
  use mod_constants, only: DP, FPINF, DASHLINE

  implicit none 

  real( kind = DP ), allocatable, dimension(:) :: A_min(:)
  integer, dimension(:)                        :: pos_A_min(1)

contains

  !---------------------------------------------------------------------------
  !> @brief This routine writes details of the calculation on screen.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Ending( timet )
    use mod_cmd_line,          only: grid_type, rad, irun
    use mod_input_read,        only: potential, atom_overlap
    use mod_grids,             only: A
    use mod_loops,             only: min_ener
    use mod_search_structures, only: n

    implicit none

    real( kind = DP ), intent(IN)                :: timet

    select case( grid_type )

      case( "shell" )

        write(output_unit,'(/, T5, A, T93, f6.2, A)') "Shell mapping done for r: ", rad, " A" 

      case ( "user" )

        write(output_unit,'(/, T5, A, T93, f6.2, A)') "Grid mapping done."

      case default

        continue

    end select

    if ( irun == "run" ) then

      write(output_unit,'(/,T5,A,T92,i9,/)') "Skipped configurations:", count(atom_overlap)

    else

      write(output_unit,*)

    endif

    if ( potential /= "none" ) then

      write(output_unit,'(T5,a34,T88,es13.5e3,/)') "Lowest energy structure (kJ/mol): ", min_ener
      write(output_unit,'(T5,a39,T94, i7,/)') "Number of structures within 0.5 * kBT: ", n

      write(output_unit,'(T5,a30,T96,i5,/)') "Lowest free energy grid point: ", minloc(A)
      write(output_unit,'(T9,a28,T88,es13.5E3,/)') "Free energy value (kJ/mol): ", minval(A)
      write(output_unit,'(T5,a31,T96,i5,/)') "Highest free energy grid point: ", maxloc(A)
      write(output_unit,'(T9,a28,T88,es13.5E3,/)') "Free energy value (kJ/mol): ", maxval(A)

    endif

    write(output_unit,'(T3, A,/)') dashline
    write(output_unit,'(T5,a20,T83,f10.2,a8,/)') "Total running time: ", timet, " seconds."

    return
  end subroutine Ending

  !---------------------------------------------------------------------------
  !> @brief This routine "marotamente" sorts central grid points according to their free energy values.
  !> @author Felippe M. Colombari
  !> @note \n 
  !> - The method used to sort the free energy array is non-standard, but is much more easier to 
  !> implement and faster to run: \n
  !>   i)   The lowest free energy value is found on array A and saved at first position of array B;
  !>   ii)  In array A, such position is replaced by INFINITY; \n
  !>   iii) Now the lowest free energy value from array A will be the second value from array B and so on...
  !---------------------------------------------------------------------------	 
  subroutine Sort_output
    use mod_grids
    use mod_loops, only: ATOTAL, mTSTOTAL, ETOTALavg, Ztrans, kBT, min_ener
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

    write(21, write_fmt ) " ", "X (A)", "Y (A)", "Z (A)", "point", "PROB", "A (kJ/mol)", "-TS (kJ/mol)", "E (kJ/mol)"

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
  end subroutine Sort_output

end module Mod_resume
