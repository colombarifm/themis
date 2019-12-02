!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief This module contains routines to search and write the lowest energy 
!>  structures for the run.         
!> @author Felippe M. Colombari                                                
!> - Laboratório de Química Teórica, LQT -- UFSCar                             
!> @date - Sep, 2017                                                           
!> - independent module created                                                
!> @date - Jan, 2018                                                           
!> - documentation and support added                                           
!------------------------------------------------------------------------------

module MOD_SEARCH
  use MOD_CONSTANTS

  implicit none

  integer                                      :: pos_min_energy(3), n
  real( kind = DP )                            :: lowest
  real( kind = DP ), allocatable, dimension(:) :: energy_ordered

  contains

    !---------------------------------------------------------------------------
    !> @author Felippe M. Colombari
    !> Laboratório de Química Teórica, LQT -- UFSCar
    !> @brief This routine counts the number of structures with energy closer 
    !>  (within 0.5 * kBT) to the most stable one.
    !---------------------------------------------------------------------------	

    subroutine COUNT_STRUCTURES
      use MOD_INPUT_READ, only: gyr_factor, inter_energy 
      use MOD_READ_GRIDS, only: grid_trans, grid_reo
      use MOD_LOOPS,      only: min_ener, kBT

      implicit none

      integer :: t, r, g

      n = 0

      do t = 1, grid_trans % numpoint

        do r = 1, grid_reo % numpoint

          do g = 1, gyr_factor

            if ( dabs( min_ener - inter_energy(g,r,t) ) <= kBT/2 ) then

              n = n + 1

            endif

          enddo

        enddo

      enddo

!      nstruc = n

      return
    end subroutine COUNT_STRUCTURES

    !---------------------------------------------------------------------------
    !> @brief This routine sorts the energy array and write its first n values 
    !>  to energy-sort.log
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Sep, 2017
    !> - subroutine  created
    !> @note \n 
    !> - The method used to sort the energy array is non-standard, but is much 
    !>   more easier to implement and faster to run: \n
    !>   i)   The lowest energy value is found on array A and saved at first 
    !>        position of array B;
    !>   ii)  In array A, such position is replaced by INFINITY; \n
    !>   iii) Now the lowest energy value from array A will be the second value 
    !>        from array B and so on...
    !---------------------------------------------------------------------------	

    subroutine SORT_ENERGY
      use MOD_INPUT_READ, only: nstruc, inter_energy
      use MOD_LOOPS,      only: kBT, min_ener, ztrans
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

      write(file_unit,'(A)') '# inter_energy(g,r,t) # g   #   r   #   t   #    prob.   # sum prob.'

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
    end subroutine SORT_ENERGY

    !---------------------------------------------------------------------------
    !> @brief This routine generates the configurations for structures written 
    !>  at energy-sort.log.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Sep, 2017
    !> - subroutine  created
    !---------------------------------------------------------------------------	
  
    subroutine SEARCH_STRUCTURES
      use mod_inquire        , only: Inquire_file
      use MOD_INPUT_READ     , only: nstruc, vector1, ref1, vector2, ref2
      use MOD_READ_MOLECULES
      use MOD_READ_GRIDS
      use MOD_LOOPS
      use mod_error_handling

      implicit none

      integer                                      :: n, g, r, t
      integer, allocatable, dimension(:)           :: g_position, r_position, t_position
      character( len = 4 )                         :: nfrm
      integer                                      :: file_unit   = 17        
      character( len = 9 )                         :: file_format = "formatted"
      character( len = 10 )                        :: file_access = "sequential"
      character( len = 15 )                        :: file_name   = "energy-sort.log"
      integer                                      :: ierr
      type(error)                                  :: err

      call Inquire_file( file_unit, file_name, file_format, file_access )

      write(*,*)
      write(*,'(T5,a23)',advance="no") "SEARCHING STRUCTURES..."   

      allocate( g_position( nstruc ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      allocate( r_position( nstruc ), stat=ierr )
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

        read(file_unit,*) energy_ordered(n), g_position(n), r_position(n), t_position(n)

        t = t_position(n)

        call mol2 % Translate_molecule ( ref2 )

        grid_reo % points(:) % grid_xyz(1) = grid_reo % points(:) % grid_xyz(1) * mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(2) = grid_reo % points(:) % grid_xyz(2) * mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(3) = grid_reo % points(:) % grid_xyz(3) * mol2 % rot_vector

        r = r_position(n)
            
        costhetar = dcos(thetar(r))
        sinthetar = dsin(thetar(r))
        cosphir   = dcos(phir(r))
        sinphir   = dsin(phir(r))

        g = g_position(n)

        call mol2 % Rotate_molecule( g )

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

        grid_reo % points(:) % grid_xyz(1) = grid_reo % points(:) % grid_xyz(1) / mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(2) = grid_reo % points(:) % grid_xyz(2) / mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(3) = grid_reo % points(:) % grid_xyz(3) / mol2 % rot_vector

      enddo

      if ( allocated(g_position) ) deallocate( g_position )
      if ( allocated(r_position) ) deallocate( r_position )
      if ( allocated(t_position) ) deallocate( t_position )
      if ( allocated(energy_ordered) ) deallocate( energy_ordered )

      close(file_unit)

      write(*,'(T44,A)') "DONE"
      write(*,*)
      write(*,'(T3,A)') repeat('-',74)

      return
    end subroutine SEARCH_STRUCTURES

end module MOD_SEARCH
