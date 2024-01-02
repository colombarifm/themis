!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2024 Themis developers
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
!> @file   mod_search_structures.f90
!> @author Felippe M. Colombari
!> @brief  This module contains routines to search and write the lowest energy structures for the run
!> @date - Jun, 2017                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Apr, 2020
!> - major revision
!> @date - May, 2021
!> - added support for ensemble of structures of molecule 2
!> @date - Jul, 2022
!> - added support to PDB files (read and write)
!---------------------------------------------------------------------------------------------------

module mod_search_structures
  use iso_fortran_env , only : output_unit
  use mod_constants   , only : DP, FPINF

  implicit none

  integer                                      :: pos_min_energy(4), n
  real( kind = DP )                            :: lowest
  real( kind = DP ), allocatable, dimension(:) :: energy_ordered

contains

  !---------------------------------------------------------------------------
  !> @brief This routine counts the number of structures with energy closer 
  !>  (within 0.5 * kBT) to the most stable one.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Count_structures
    use mod_input_read , only : axis_rot_moves, nconf2, inter_energy 
    use mod_grids      , only : grid_trans, grid_rot1
    use mod_loops      , only : min_ener, kBT

    implicit none

    integer :: n_trans, n_conf, n_rot1, n_rot2

    n = 0

    do n_trans = 1, grid_trans % numpoint

      do n_conf = 1, nconf2

        do n_rot1 = 1, grid_rot1 % numpoint

          do n_rot2 = 1, axis_rot_moves

            if ( dabs( min_ener - inter_energy( n_rot2, n_rot1, n_conf, n_trans ) ) <= kBT/2 ) then

              n = n + 1

            endif

          enddo

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
    use mod_inquire        , only: Inquire_file, Get_new_unit
    use mod_input_read     , only: nstruc, inter_energy
    use mod_loops          , only: kBT, min_ener, ztrans
    use mod_error_handling

    implicit none
  
    integer                                      :: n
    real( kind = DP )                            :: inv_Ztrans, delta_E
    real( kind = DP ), allocatable, dimension(:) :: prob
    integer                                      :: file_unit          
    character( len = * ), parameter              :: file_access = "sequential"
    character( len = * ), parameter              :: file_format = "formatted"
    character( len = * ), parameter              :: file_status = "unknown"
    character( len = * ), parameter              :: file_name   = "energy-sort.log"
    integer                                      :: ierr
    type(error)                                  :: err

    file_unit = Get_new_unit(10)

    open( unit = file_unit, file = trim(file_name), status = file_status, &
          form = trim(file_format), access = trim(file_access) ) 

    write( file_unit, '(A)' ) '# energy(r2,r1,c,t) #   r2 #   r1 #  c #     t #      delta_E #        prob. #    sum_prob.'

    allocate( energy_ordered ( nstruc ), stat=ierr )
    if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( prob ( nstruc ), stat=ierr )
    if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    energy_ordered = 0.0_dp

    prob = 0.0_dp

    inv_Ztrans = 1.0_DP / Ztrans 

    do n = 1, nstruc

      energy_ordered(n) = minval(inter_energy)

      delta_E = energy_ordered(n) - energy_ordered(1)

      prob(n) = inv_Ztrans * dexp( ( -energy_ordered(n) + min_ener ) / kBT )

      pos_min_energy = minloc(inter_energy)

      inter_energy( pos_min_energy(1), pos_min_energy(2), pos_min_energy(3), pos_min_energy(4) ) = FPINF

      write( file_unit, '(es16.5E3,3x,2i7,i5,i8,3es15.5E3)' ) energy_ordered(n), pos_min_energy, delta_E, prob(n), sum(prob)

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
    use mod_inquire        , only: Inquire_file, Get_new_unit
    use mod_input_read     , only: nstruc, vector1, ref1, vector2, ref2, nconf2, file_type
    use mod_read_molecules
    use mod_grids
    use mod_loops
    use mod_error_handling

    implicit none

    integer                                      :: n, n_rot2, n_rot1, n_conf, n_trans
    integer, allocatable, dimension(:)           :: rot2_position, rot1_position, conf_position, trans_position
    character( len = 4 )                         :: nfrm
    integer                                      :: file_unit          
    character( len = * ), parameter              :: file_status = "old"
    character( len = * ), parameter              :: file_format = "formatted"
    character( len = * ), parameter              :: file_access = "sequential"
    character( len = * ), parameter              :: file_name   = "energy-sort.log"
    integer                                      :: ierr
    type( error )                                :: err

    file_unit = Get_new_unit(10)

    call Inquire_file( file_unit, file_name, file_status, file_format, file_access )

    write( output_unit, '(/,T5,a23)', advance = "no" ) "SEARCHING STRUCTURES..."   

    allocate( rot2_position( nstruc ), stat=ierr )
    if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( rot1_position( nstruc ), stat=ierr )
    if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( conf_position( nstruc ), stat=ierr )
    if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( trans_position( nstruc ), stat=ierr )
    if( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    read( file_unit, * )

    costhetar = -1.0_DP
    sinthetar =  0.0_DP
    cosphir   =  0.0_DP
    sinphir   =  1.0_DP

    call mol1 % Read_molecule( "conf1."//file_type, 1, file_type )
    call mol2 % Read_molecule( "conf2."//file_type, nconf2, file_type )

    call mol1 % Translate_molecule( ref1, 1 )
    call mol1 % Align_molecule( vector1, ref1, 1 )
    call mol1 % Rotate_molecule( 1, 1 ) 
    
    do n_conf = 1, nconf2
    
      call mol2 % Translate_molecule( ref2, n_conf )
      call mol2 % Align_molecule( vector2, ref2, n_conf )

    enddo
 
    do n = 1, nstruc

      read( file_unit, * ) energy_ordered(n), rot2_position(n), rot1_position(n), conf_position(n), trans_position(n)

      n_trans = trans_position(n)

      n_conf = conf_position(n)

      call mol2 % Translate_molecule ( ref2, n_conf )

      grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) * mol2 % rot_vector( n_conf )
      grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) * mol2 % rot_vector( n_conf )
      grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) * mol2 % rot_vector( n_conf )

      n_rot1 = rot1_position(n)
            
      costhetar = dcos(thetar( n_rot1 ))
      sinthetar = dsin(thetar( n_rot1 ))
      cosphir   = dcos(phir( n_rot1 ))
      sinphir   = dsin(phir( n_rot1 ))

      n_rot2 = rot2_position(n)

      call mol2 % Rotate_molecule( n_rot2, n_conf )

      mol2 % atoms( n_conf, : ) % xyz(1) = mol2 % atoms( n_conf, : ) % xyz(1) + grid_trans % points( n_trans ) % grid_xyz(1)
      mol2 % atoms( n_conf, : ) % xyz(2) = mol2 % atoms( n_conf, : ) % xyz(2) + grid_trans % points( n_trans ) % grid_xyz(2)
      mol2 % atoms( n_conf, : ) % xyz(3) = mol2 % atoms( n_conf, : ) % xyz(3) + grid_trans % points( n_trans ) % grid_xyz(3)

      write( nfrm, '(I4.4)' ) n

      lowest = energy_ordered(n)

      if ( file_type == "xyz" ) then

        open( unit = 66, file = 'lowest_'//nfrm//'.xyz', status = 'unknown' )

        call dimers % Build_dimer( n_conf, file_type )

        call dimers % Write_xyz( lowest, n_conf )
 
        close( 66 )

      else if ( file_type == "pdb" ) then

        open( unit = 66, file = 'lowest_'//nfrm//'.pdb', status = 'unknown' )

        call dimers % Build_dimer( n_conf, file_type )

        call dimers % Write_pdb( n_conf )
 
        close( 66 )

      endif
                
      mol2 % atoms( n_conf, : ) % xyz(1) = mol2 % atoms( n_conf, : ) % xyz_old(1)
      mol2 % atoms( n_conf, : ) % xyz(2) = mol2 % atoms( n_conf, : ) % xyz_old(2)
      mol2 % atoms( n_conf, : ) % xyz(3) = mol2 % atoms( n_conf, : ) % xyz_old(3)

      grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) / mol2 % rot_vector( n_conf )
      grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) / mol2 % rot_vector( n_conf )
      grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) / mol2 % rot_vector( n_conf )

    enddo

    if ( allocated(rot2_position) )  deallocate( rot2_position )
    if ( allocated(rot1_position) )  deallocate( rot1_position )
    if ( allocated(conf_position) )  deallocate( conf_position )
    if ( allocated(trans_position) ) deallocate( trans_position )
    if ( allocated(energy_ordered) ) deallocate( energy_ordered )

    close(file_unit)

    write( output_unit, '(T70,A)' ) "DONE"
    write( output_unit, '(/,T3,A)' ) dashline

    return
  end subroutine Search_structures

end module mod_search_structures
