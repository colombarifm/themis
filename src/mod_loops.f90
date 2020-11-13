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
!> @file   mod_loops.f90
!> @author Felippe M. Colombari
!> @brief  This module contains the main program loops for Z(N,V,T) calculation
!> @date - Oct, 2017                                                           
!> - independent module created     
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_loops
  use iso_fortran_env, only: output_unit
  use mod_constants

  implicit none

  character( len = 4 )                         :: ntra, nrot1, nrot2
  real( kind = DP )                            :: t0, time, kBT, min_ener
  real( kind = DP )                            :: ATOTAL, mTSTOTAL, ETOTALavg
  real( kind = DP )                            :: Ztrans, sumVexpVtrans
  real( kind = DP ), allocatable, dimension(:) :: min_ener_t
  integer                                      :: tini, tend, rate

contains

  !---------------------------------------------------------------------------
  !> @brief This routine performs the configurational sampling and energy calculation
  !> @author Felippe M. Colombari
  !> @date - Jun, 2017
  !> - subroutine  created
  !> - energy values written in binary file
  !> @date - Dec, 2017
  !> - lowest energy values are subtracted from each value to avoid problems
  !>   due to exponential explosion
  !> @date - Sep, 2019                                                           
  !> - pointers for potential selection addded
  !---------------------------------------------------------------------------	
  subroutine Run_loops
    use mod_cmd_line
    use mod_input_read
    use mod_read_molecules
    use mod_grids
    use mod_pot_ljc
    use mod_pot_bhc
    use xdr,           only: xtcfile
    use mod_error_handling

    implicit none

    type( xtcfile )                   :: xtc_out
    real( kind = SP ), dimension(:,:) :: pos(3,mol1 % num_atoms + mol2 % num_atoms)
    real( kind = SP ), dimension(3,3) :: box
    real( kind = SP )                 :: frametime, prec
    integer                           :: frame, i, j, r2, r1, t

    real( kind = DP )                 :: lowest
    character( len = 20 )             :: prefix

    type( ljc_dimer ), target         :: ljc_target
    type( bhc_dimer ), target         :: bhc_target
    type( dimer ), target             :: none_target
    class( dimer ), pointer           :: potential_pointer

    integer                                      :: ierr
    type(error)                                  :: err

    allocate( inter_energy( rot2_factor, grid_rot1 % numpoint, grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate( atom_overlap( rot2_factor, grid_rot1 % numpoint, grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(        Zrot( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(  min_ener_t( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate( sumVexpVrot( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    frame         = 1
    min_ener      = 0.0_DP
    min_ener_t    = 0.0_DP
    Zrot          = 0.0_DP
    sumVexpVrot   = 0.0_DP
    t0            = 0.0_DP
    inter_energy  = 0.0_DP
    atom_overlap  = .false.
    kBT           = KB * temp
      
    SELECT CASE (wrtxtc)

      CASE ("yes", "true", "T")

        box      = 0.0_SP
        prec     = 1000.0_SP
          
        call xtc_out % init( "full_ensemble.xtc", "w" )

    end SELECT
  
    open( unit = 79, file = 'grid_log.log', status = 'replace' )

    write(79,'("t point   rejected structures   time (s)")') 
        
    flush(79)
      
    ! WRITE ENERGY VALUES ON BINARY FILE (FASTER WRITING AND +80% SPACE SAVING!):
      
    open( unit = 271, file = 'energy.bin', form = 'unformatted', access = 'stream', status = 'replace' )

    SELECT CASE (potential)

      CASE ("lj-coul")

        potential_pointer => ljc_target

        call mol1_ljc % Read_LJC_params( "parameters1", mol1 % num_atoms )
        call mol2_ljc % Read_LJC_params( "parameters2", mol2 % num_atoms )

        call mol1_ljc % Check_LJC_params( 1, mol1 % num_atoms, "parameters1" )
        call mol2_ljc % Check_LJC_params( 2, mol2 % num_atoms, "parameters2" )

      CASE ("bh-coul")

        potential_pointer => bhc_target

        call mol1_bhc % Read_BHC_params( "parameters1", mol1 % num_atoms )
        call mol2_bhc % Read_BHC_params( "parameters2", mol2 % num_atoms )

        call mol1_bhc % Check_BHC_params( 1, mol1 % num_atoms, "parameters1" )
        call mol2_bhc % Check_BHC_params( 2, mol2 % num_atoms, "parameters2" )

      CASE ("none")

        potential_pointer => none_target

    end SELECT

    call Check_moves
    call potential_pointer % Calc_cross

    write(output_unit,'(/,T5,a17)',advance='no') "ENTERING LOOPS..."
    
    call system_clock(tini, rate)

    tlp: do t = 1, grid_trans % numpoint

      call mol2 % Translate_molecule( ref2 )
      
      grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) * mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) * mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) * mol2 % rot_vector
          
      r1lp: do r1 = 1, grid_rot1 % numpoint

        costhetar = dcos(thetar(r1))
        sinthetar = dsin(thetar(r1))
        cosphir   = dcos(phir(r1))
        sinphir   = dsin(phir(r1))

        r2lp: do r2 = 1, rot2_factor

          call mol2 % Rotate_molecule( r2 ) 

          do j = 1, mol2 % num_atoms

            mol2 % atoms(j) % xyz(:) = mol2 % atoms(j) % xyz(:) + grid_trans % points(t) % grid_xyz(:)

          enddo

          call potential_pointer % Calc_energy( r2, r1, t )

          write(271) inter_energy( r2, r1, t )

          if ( atom_overlap( r2, r1, t ) .eqv. .false. ) then
            
            write(prefix,'("point_",I4.4,"_",I4.4,"_",I4.4)') t, r1, r2 
              
            SELECT CASE (writeframe)

              CASE ("XYZ", "xyz", "Xyz")

                open( unit = 66, file = prefix//'.xyz', status = 'unknown' )

                lowest = 0.0_DP

                call dimers % Build_dimer
             
                call dimers % Write_xyz( lowest )
                            
                close(66)

              CASE ("MOP", "mop", "Mop")
  
                open( unit = 66, file = prefix//'.mop', status = 'unknown' )
        
                call dimers % Build_dimer
          
                call dimers % Write_mop( mopac_head )
                            
                close(66)
              
            END SELECT

          endif 

          SELECT CASE (wrtxtc)

            CASE ( "yes", "true", "T" )

              do i = 1, mol1 % num_atoms

                pos(:,i) = mol1 % atoms(i) % xyz(:) / 10.0_SP

              enddo

              do j = 1, mol2 % num_atoms

                pos(:,j+mol1 % num_atoms) = mol2 % atoms(j) % xyz(:) / 10.0_SP

              enddo

              frametime = frame / 1000.0_SP

              call xtc_out % write( mol1 % num_atoms + mol2 % num_atoms, frame, frametime, box, pos, prec )

              frame = frame + 1

          end SELECT

          mol2 % atoms(:) % xyz(1) = mol2 % atoms(:) % xyz_old(1)
          mol2 % atoms(:) % xyz(2) = mol2 % atoms(:) % xyz_old(2)
          mol2 % atoms(:) % xyz(3) = mol2 % atoms(:) % xyz_old(3)

        enddo r2lp

      enddo r1lp

      grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) / mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) / mol2 % rot_vector
      grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) / mol2 % rot_vector
          
      call system_clock(tend)

      time = real( tend - tini, DP ) / real( rate, DP ) 

      write(79,'(i7,3x,i7," of ", i7,1x,f10.3)') t, count( atom_overlap(:,:,t) ), rot_total, time-t0

      flush(79)

      t0 = time

      min_ener_t(t) = minval( inter_energy(:,:,t) )

    enddo tlp

    SELECT CASE (wrtxtc)

      CASE ("yes", "true", "T")

        call xtc_out % close

    end SELECT

    close(79)
    close(271)

    min_ener = minval( inter_energy )

    do t = 1, grid_trans % numpoint

      SELECT CASE (potential)

        CASE ("lj-coul", "bh-coul")

          Zrot(t) = sum( dexp( -( inter_energy(:,:,t) - min_ener_t(t) ) / kBT ) )
            
          sumVexpVrot(t) = sum( ( inter_energy(:,:,t) - min_ener_t(t) ) * &
            dexp( -( inter_energy(:,:,t) - min_ener_t(t) ) / kBT ) ) 

      end SELECT

    enddo 

    write(output_unit,'(T76,A)') "DONE"

    return
  end subroutine Run_loops

  !---------------------------------------------------------------------------
  !> @brief This routine performs the main loop and reads already calculated
  !>  energy values from energy file.
  !> @author Felippe M. Colombari
  !> @date - Jun, 2017
  !> - subroutine  created
  !> - energy values read from binary file or formatted file from external programs
  !> @date - Dec, 2017
  !> - lowest energy values are subtracted from each value to avoid problems
  !>   in exponential calculation
  !---------------------------------------------------------------------------	
  subroutine Rerun_loops
    use mod_input_read, only: rot2_factor, temp, potential, inter_energy
    use mod_grids,      only: Zrot, sumVexpVrot, &
                              Check_Moves, grid_trans, grid_rot1
    use mod_inquire,    only: Inquire_file, Get_new_unit       
    use mod_error_handling

    implicit none

    integer                           :: r2, r1, t
    integer                           :: ios         = 0
    integer                           :: file_unit   
    character( len = : ), allocatable :: file_format 
    character( len = : ), allocatable :: file_access 
    character( len = : ), allocatable :: file_name   
    character( len = : ), allocatable :: file_status 
    character( len = 64 )             :: line1
    integer                           :: ierr
    type(error)                       :: err

    allocate(        Zrot( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(  min_ener_t( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate( sumVexpVrot( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate( inter_energy( rot2_factor, grid_rot1 % numpoint , grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    kBT = KB * temp
      
    ios         = 0
    Zrot        = 0.0_DP
    min_ener    = 0.0_DP
    sumVexpVrot = 0.0_DP
    inter_energy  = 0.0_DP

    if ( potential == "none" ) then

      file_unit   = Get_new_unit(10)
      file_format = "formatted"
      file_access = "sequential"
      file_name   = "energy.log"
      file_status = "old"

      call Inquire_file( file_unit , file_name , file_status, file_format , file_access )

      read(file_unit,*,iostat=ios) line1
  
      if ( ios /= 0 ) then

        call err%error('e',message="Error while reading energy.log file reader.")

        stop

      endif

      do t = 1, grid_trans % numpoint

        do r1 = 1, grid_rot1 % numpoint

          do r2 = 1, rot2_factor

            read(file_unit,*,iostat=ios) inter_energy( r2, r1, t )

            if ( ios /= 0 ) then

              call err%error('e',message="Error while reading energy.log file.")

              stop

            endif

          enddo 

        enddo 

        min_ener_t(t) = minval( inter_energy(:,:,t) )

      enddo 

      close(file_unit)

    else

      file_unit   = Get_new_unit(10)
      file_format = "unformatted"
      file_access = "stream"
      file_name   = "energy.bin"
      file_status = "old"

      call Inquire_file( file_unit , file_name , file_status, file_format , file_access )
        
      do t = 1, grid_trans % numpoint

        do r1 = 1, grid_rot1 % numpoint

          do r2 = 1, rot2_factor

            read(file_unit, iostat=ios) inter_energy( r2, r1, t )

            if ( ios /= 0 ) then

              call err%error('e',message="Error while reading energy.bin file.")

              stop

            endif

          enddo 

        enddo 

        min_ener_t(t) = minval( inter_energy(:,:,t) )

      enddo 

      close(file_unit)

    endif

    min_ener = minval( inter_energy )

    do t = 1, grid_trans % numpoint

      Zrot(t) = sum( dexp( -( inter_energy(:,:,t) - min_ener_t(t) ) / kBT ) )
            
      sumVexpVrot(t) = sum( ( inter_energy(:,:,t) - min_ener_t(t) ) * dexp( -( inter_energy(:,:,t) - min_ener_t(t) ) / kBT ) )

    enddo 
      
    return
  end subroutine Rerun_loops

  !---------------------------------------------------------------------------
  !> @brief This routine performs partition function calculation for each grid
  !>  point and also for the whole grid.
  !> @author Felippe M. Colombari
  !> @date - Jun, 2017
  !> - subroutine  created
  !> @date - Dec, 2017
  !> - absolute energy values are recovered and thermodynamic properties are
  !>   obtained from them
  !> @date - Mar, 2019
  !> - -TdS from ideal gas is subtracted from free energy
  !---------------------------------------------------------------------------	
  subroutine Calculate_ztotal
    use mod_input_read, only: temp, rot2_factor, inter_energy
    use mod_grids,      only: probT, Zrot, A, Eavg, mTS, sumVexpVrot, grid_trans, grid_rot1
    use mod_error_handling

    implicit none

    integer                       :: t
    character(len=:), allocatable :: write_fmt
    real(kind=DP)                 :: mTS_ideal, mTS_ideal_t
    integer                       :: ierr
    type(error)                   :: err

    kBT = KB * temp

    allocate(     A( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate(   mTS( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate(  Eavg( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate( probT( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    A             = 0.0_DP
    Eavg          = 0.0_DP
    mTS           = 0.0_DP
    probT         = 0.0_DP
    mTS_ideal     = -kBT * log ( real(grid_trans % numpoint * grid_rot1 % numpoint * rot2_factor) )
    mTS_ideal_t   = -kBT * log ( real(grid_rot1 % numpoint * rot2_factor) )

    open( unit = 20, file = 'output.log', position = 'append', status = 'replace' )

    write_fmt = '(a1,2x,3(a11,1x),3x,a5,3x,4(a13,3x))'

    write(20,*) grid_trans % numpoint
    write(20, write_fmt ) " ", "X (A)", "Y (A)", "Z (A)", "point", "PROB", "A (kJ/mol)", "-TS (kJ/mol)", "E (kJ/mol)"

    write(output_unit,*)
    write(output_unit,'(T5,a23)',advance="no") "CALCULATING Z(N,V,T)..."   

    Ztrans = sum( dexp( -( inter_energy - min_ener ) / kBT ) ) 

    sumVexpVtrans = sum( ( inter_energy - min_ener ) * dexp( -( inter_energy - min_ener ) / kBT ) ) 

    write_fmt = '(A,2x,3(f11.5,1x),3x,i5,3x,4(es13.5E3,3x))'

    do t = 1, grid_trans % numpoint

      probT(t) = Zrot(t) / Ztrans * dexp( ( min_ener - min_ener_t(t) ) / kBT ) 

      A(t) = -kBT * dlog( Zrot(t) ) + min_ener_t(t) - mTS_ideal_t

      Eavg(t) = sumVexpVrot(t) / Zrot(t) + min_ener_t(t)

      mTS(t) = A(t) - Eavg(t)

      write( 20, write_fmt ) "X", grid_trans % points(t) % grid_xyz(:), t, probT(t), A(t), mTS(t), Eavg(t)

    enddo 

    ATOTAL = -kBT * dlog( Ztrans ) + min_ener - mTS_ideal

    ETOTALavg = sumVexpVtrans / Ztrans + min_ener

    mTSTOTAL = ATOTAL - ETOTALavg

    write_fmt = '(9x,"TOTAL OVER TRANSLATIONAL GRID",12x,5(es13.5E3,3x))'

    write(20,'(a)') repeat('-',152)
    write(20, write_fmt) sum(probT), ATOTAL, mTSTOTAL, ETOTALavg

    close(20)

    write(output_unit,'(T70,A)') "DONE"

    return
  end subroutine Calculate_ztotal

end module mod_loops
