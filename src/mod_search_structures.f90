!---------------------------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition function estimation                                                  
!---------------------------------------------------------------------------------------------------
!   Copyright 2020 Felippe M. Colombari
!
!   This program is free software: you can redistribute it and/or modify it under the terms of the 
!   GNU General Public License as published by the Free Software Foundation, either version 3 of the 
!   License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
!   the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License along with this program. If 
!   not, see <https://www.gnu.org/licenses/>.
!---------------------------------------------------------------------------------------------------
!> @file   mod_search_structures.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  This module contains routines to search and write the lowest energy structures for the run
!> @date - Jun, 2017                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_search_structures
  use iso_fortran_env, only: output_unit
  use mod_constants, only: DP, FPINF

  implicit none

  integer                                      :: pos_min_energy(3), n
  real( kind = DP )                            :: lowest
  real( kind = DP ), allocatable, dimension(:) :: energy_ordered

contains

  !---------------------------------------------------------------------------
  !> @brief This routine counts the number of structures with energy closer 
  !>  (within 0.5 * kBT) to the most stable one.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Count_structures
    use mod_input_read, only: rot2_factor, inter_energy 
    use mod_grids,      only: grid_trans, grid_rot1
    use mod_loops,      only: min_ener, kBT

    implicit none

    integer :: t, r1, r2

    n = 0

    do t = 1, grid_trans % numpoint

      do r1 = 1, grid_rot1 % numpoint

        do r2 = 1, rot2_factor

          if ( dabs( min_ener - inter_energy(r2,r1,t) ) <= kBT/2 ) then

            n = n + 1

          endif

        enddo

      enddo

    enddo

    ! to write out all structures within this energy range, please uncomment:
    !nstruc = n

    return
  end subroutine Count_structures

  !---------------------------------------------------------------------------
  !> @brief This routine "marotamente" sorts the energy array and write its first n values 
  !>  to energy-sort.log
  !> @author Felippe M. Colombari
  !> @note \n 
  !> - The method used to sort the free energy array is non-standard, but is much more easier to 
  !> implement and faster to run: \n
  !>   i)   The lowest free energy value is found on array A and saved at first position of array B;
  !>   ii)  In array A, such position is replaced by INFINITY; \n
  !>   iii) Now the lowest free energy value from array A will be the second value from array B and so on...
  !---------------------------------------------------------------------------	 
  subroutine Sort_energy
    use mod_input_read, only: nstruc, inter_energy
    use mod_loops,      only: kBT, min_ener, ztrans
    use mod_error_handling

    implicit none
  
    integer                                      :: i
    real( kind = DP )                            :: inv_Ztrans
    real( kind = DP ), allocatable, dimension(:) :: prob
    integer                                      :: file_unit   = 17        
    character( len = 9 )                         :: file_format = "formatted"
    character( len = 7 )                         :: file_status = "unknown"
    character( len = 15 )                        :: file_name   = "energy-sort.log"
    integer                                      :: ierr
    type(error)                                  :: err

    open( unit = file_unit, file = file_name, form = file_format, status = file_status )

    write(file_unit,'(A)') '# inter_energy(g,r,t) # r2  #   r1  #   t   #    prob.   # sum prob.'

    allocate( energy_ordered ( nstruc ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate( prob ( nstruc ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    energy_ordered = 0.0_dp

    prob = 0.0_dp

    inv_Ztrans = 1.0_DP / Ztrans 

    do i = 1, nstruc

      energy_ordered(i) = minval(inter_energy)

      prob(i) = inv_Ztrans * dexp( ( -energy_ordered(i) + min_ener ) / kBT )

      pos_min_energy = minloc(inter_energy)

      inter_energy(pos_min_energy(1),pos_min_energy(2),pos_min_energy(3)) = FPINF

      write(file_unit,'(4x,es13.5E3,3x,3(i5,3x),2(2x,es10.3E3))') energy_ordered(i), pos_min_energy, prob(i), sum(prob)

    enddo

    close(file_unit)

    if ( allocated(prob) ) deallocate( prob )
      
    return
  end subroutine Sort_energy

  !---------------------------------------------------------------------------
  !> @brief This routine generates the configurations for structures written at energy-sort.log.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Search_structures
    use mod_inquire        , only: Inquire_file
    use mod_input_read     , only: nstruc, vector1, ref1, vector2, ref2
    use mod_read_molecules
    use mod_grids
    use mod_loops
    use mod_error_handling

    implicit none

    integer                                      :: n, r2, r1, t
    integer, allocatable, dimension(:)           :: r2_position, r1_position, t_position
    character( len = 4 )                         :: nfrm
    integer                                      :: file_unit   = 17        
    character( len = 9 )                         :: file_format = "formatted"
    character( len = 10 )                        :: file_access = "sequential"
    character( len = 15 )                        :: file_name   = "energy-sort.log"
    integer                                      :: ierr
    type(error)                                  :: err

    call Inquire_file( file_unit, file_name, file_format, file_access )

    write(output_unit,'(/,T5,a23)',advance="no") "SEARCHING STRUCTURES..."   

    allocate( r2_position( nstruc ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate( r1_position( nstruc ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate( t_position( nstruc ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    read(file_unit,*)

    costhetar = -1.0_DP
    sinthetar =  0.0_DP
    cosphir   =  0.0_DP
    sinphir   =  1.0_DP

    call mol1 % Read_molecule( "conf1.xyz" )
    call mol1 % Translate_molecule( ref1 )
    call mol1 % Align_molecule( vector1, ref1 )
    call mol1 % Rotate_molecule( 1 ) 

    call mol2 % Read_molecule( "conf2.xyz" )
    call mol2 % Translate_molecule( ref2 )
    call mol2 % Align_molecule( vector2, ref2 )
 
    do n = 1, nstruc

      read(file_unit,*) energy_ordered(n), r2_position(n), r1_position(n), t_position(n)

      t = t_position(n)

      call mol2 % Translate_molecule ( ref2 )

      grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) * mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) * mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) * mol2 % rot_vector

      r1 = r1_position(n)
            
      costhetar = dcos(thetar(r1))
      sinthetar = dsin(thetar(r1))
      cosphir   = dcos(phir(r1))
      sinphir   = dsin(phir(r1))

      r2 = r2_position(n)

      call mol2 % Rotate_molecule( r2 )

      mol2 % atoms(:) % xyz(1) = mol2 % atoms(:) % xyz(1) + grid_trans % points(t) % grid_xyz(1)
      mol2 % atoms(:) % xyz(2) = mol2 % atoms(:) % xyz(2) + grid_trans % points(t) % grid_xyz(2)
      mol2 % atoms(:) % xyz(3) = mol2 % atoms(:) % xyz(3) + grid_trans % points(t) % grid_xyz(3)

      write(nfrm,'(I4.4)') n

      lowest = energy_ordered(n)

      open( unit = 66, file = 'lowest_'//nfrm//'.xyz', status = 'unknown' )
 
      call dimers % Build_dimer

      call dimers % Write_xyz( lowest )
 
      close(66)
                
      mol2 % atoms(:) % xyz(1) = mol2 % atoms(:) % xyz_old(1)
      mol2 % atoms(:) % xyz(2) = mol2 % atoms(:) % xyz_old(2)
      mol2 % atoms(:) % xyz(3) = mol2 % atoms(:) % xyz_old(3)

      grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) / mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) / mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) / mol2 % rot_vector

    enddo

    if ( allocated(r2_position) ) deallocate( r2_position )
    if ( allocated(r1_position) ) deallocate( r1_position )
    if ( allocated(t_position) ) deallocate( t_position )
    if ( allocated(energy_ordered) ) deallocate( energy_ordered )

    close(file_unit)

    write(output_unit,'(T70,A)') "DONE"
    write(output_unit,'(/,T3,A)') dashline

    return
  end subroutine Search_structures

end module mod_search_structures
