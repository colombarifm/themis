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
!> @file   mod_loops.f90
!> @author Felippe M. Colombari
!> @brief  This module contains the main program loops for Z(N,V,T) calculation
!> @date - Oct, 2017                                                           
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

module mod_loops
  use iso_fortran_env , only: output_unit
  use mod_constants

  implicit none

  integer                                      :: tini, tend, rate
  real( kind = DP )                            :: t0, time, kBT, min_ener
  real( kind = DP )                            :: ATOTAL, mTSTOTAL, ETOTALavg
  real( kind = DP )                            :: Ztrans, sumVexpVtrans
  real( kind = DP ), allocatable, dimension(:) :: min_ener_t

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
    !use mod_pot_bhc
    use xdr                , only : xtcfile
    use mod_error_handling

    implicit none

    type( xtcfile )                   :: xtc_out
    real( kind = SP ), dimension(:,:) :: pos(3,mol1 % num_atoms + mol2 % num_atoms)
    real( kind = SP ), dimension(3,3) :: box
    real( kind = SP )                 :: frametime, prec
    integer                           :: frame, n_atom_1, n_atom_2, n_rot2, n_rot1, n_conf, n_trans

    real( kind = DP )                 :: lowest
    character( len = 24 )             :: prefix

    type( ljc_dimer ), target         :: ljc_target
    !type( bhc_dimer ), target         :: bhc_target
    type( dimer ), target             :: none_target
    class( dimer ), pointer           :: potential_pointer

    integer                           :: ierr
    type(error)                       :: err

    allocate( inter_energy( axis_rot_moves, grid_rot1 % numpoint, nconf2, grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err % error( 'e', message = "abnormal memory allocation" )
      
    allocate( atom_overlap( axis_rot_moves, grid_rot1 % numpoint, nconf2, grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( Zrot( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( min_ener_t( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err % error( 'e', message = "abnormal memory allocation" )
      
    allocate( sumVexpVrot( grid_trans % numpoint ), stat=ierr )
    if(ierr/=0) call err % error( 'e', message = "abnormal memory allocation" )

    frame        = 1
    min_ener     = 0.0_DP
    min_ener_t   = 0.0_DP
    Zrot         = 0.0_DP
    sumVexpVrot  = 0.0_DP
    t0           = 0.0_DP
    inter_energy = 0.0_DP
    atom_overlap = .false.
    kBT          = KB * temp
      
    select case ( wrtxtc )

      case ( "yes", "true", "T" )

        box  = 0.0_SP
        prec = 1000.0_SP
          
        call xtc_out % init( "full_ensemble.xtc", "w" )
  
    end select
  
    open( unit = 79, file = 'grid_log.log', status = 'replace' )

    write( 79,'("#n_trans  rejected structures   time (s)")' ) 
        
    flush( 79 )
      
    ! WRITE ENERGY VALUES ON BINARY FILE (FASTER WRITING AND +80% SPACE SAVING!):
      
    open( unit = 271, file = 'energy.bin', form = 'unformatted', access = 'stream', status = 'replace' )

    select case ( potential )

      case ( "lj-coul" )

        potential_pointer => ljc_target

        call mol1_ljc % Read_LJC_params( "parameters1", mol1 % num_atoms )
        call mol2_ljc % Read_LJC_params( "parameters2", mol2 % num_atoms )

        call mol1_ljc % Check_LJC_params( 1, mol1 % num_atoms, "parameters1" )
        call mol2_ljc % Check_LJC_params( 2, mol2 % num_atoms, "parameters2" )

      !TO DO
      !CASE ("lj")
      !CASE ("lj-coul_AB")
      !CASE ("lj_AB")

      !case ( "bh-coul" )

        !potential_pointer => bhc_target

        !call mol1_bhc % Read_BHC_params( "parameters1", mol1 % num_atoms )
        !call mol2_bhc % Read_BHC_params( "parameters2", mol2 % num_atoms )

        !call mol1_bhc % Check_BHC_params( 1, mol1 % num_atoms, "parameters1" )
        !call mol2_bhc % Check_BHC_params( 2, mol2 % num_atoms, "parameters2" )

      case ( "none" )

        potential_pointer => none_target

    end select

    call Check_moves
    call potential_pointer % Calc_cross

    write( output_unit, '(/,T5,a17)', advance = 'no' ) "ENTERING LOOPS..."
    
    call system_clock( tini, rate )

    tlp: do n_trans = 1, grid_trans % numpoint

      clp: do n_conf = 1, nconf2  ! mol2 % nconf

        call mol2 % Translate_molecule( ref2, n_conf )
      
        grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) * mol2 % rot_vector( n_conf )
        grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) * mol2 % rot_vector( n_conf )
        grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) * mol2 % rot_vector( n_conf )
          
        r1lp: do n_rot1 = 1, grid_rot1 % numpoint

          costhetar = dcos(thetar( n_rot1 ))
          sinthetar = dsin(thetar( n_rot1 ))
          cosphir   = dcos(phir( n_rot1 ))
          sinphir   = dsin(phir( n_rot1 ))

          r2lp: do n_rot2 = 1, axis_rot_moves

            call mol2 % Rotate_molecule( n_rot2, n_conf ) 

            do n_atom_2 = 1, mol2 % num_atoms

              mol2 % atoms( n_conf, n_atom_2 ) % xyz(:) = mol2 % atoms( n_conf, n_atom_2 ) % xyz(:) &
                                                        + grid_trans % points( n_trans ) % grid_xyz(:)

            enddo

            call potential_pointer % Calc_energy( n_rot2, n_rot1, n_conf, n_trans )

            inter_energy( n_rot2, n_rot1, n_conf, n_trans ) = inter_energy( n_rot2, n_rot1, n_conf, n_trans ) &
                                                            + mol2 % conf_energy( n_conf )

            write( 271 ) inter_energy( n_rot2, n_rot1, n_conf, n_trans )

            if ( atom_overlap( n_rot2, n_rot1, n_conf, n_trans ) .eqv. .false. ) then
            
              write( prefix, '("point_",I5.5,"_",I2.2,"_",I4.4,"_",I4.4)' ) n_trans, n_conf, n_rot1, n_rot2 
              
              select case ( writeframe )

                case ( "xyz" )

                  open( unit = 66, file = prefix//'.xyz', status = 'unknown' )

                  lowest = 0.0_DP

                  call dimers % Build_dimer( n_conf, file_type )
             
                  call dimers % Write_xyz( lowest, n_conf )
                            
                  close( 66 )

                case ( "pdb" )

                  open( unit = 66, file = prefix//'.pdb', status = 'unknown' )

                  call dimers % Build_dimer( n_conf, file_type )
             
                  call dimers % Write_pdb( n_conf )
                           
                  close( 66 )

                case ( "mop" )
  
                  open( unit = 66, file = prefix//'.mop', status = 'unknown' )
        
                  call dimers % Build_dimer( n_conf, "xyz" )
          
                  call dimers % Write_mop( mopac_head, n_conf )
                            
                  close( 66 )
              
              end select

            endif 

            select case ( wrtxtc )

              case ( "yes", "true", "T" )

                do n_atom_1 = 1, mol1 % num_atoms

                  pos( :, n_atom_1 ) = mol1 % atoms( 1, n_atom_1 ) % xyz(:) / 10.0_SP

                enddo

                do n_atom_2 = 1, mol2 % num_atoms

                  pos( :, n_atom_2+mol1 % num_atoms ) = mol2 % atoms( n_conf, n_atom_2 ) % xyz(:) / 10.0_SP

                enddo

                frametime = frame / 1000.0_SP

                call xtc_out % write( mol1 % num_atoms + mol2 % num_atoms, frame, frametime, box, pos, prec )

                frame = frame + 1

            end select

            mol2 % atoms(n_conf,:) % xyz(1) = mol2 % atoms(n_conf,:) % xyz_old(1)
            mol2 % atoms(n_conf,:) % xyz(2) = mol2 % atoms(n_conf,:) % xyz_old(2)
            mol2 % atoms(n_conf,:) % xyz(3) = mol2 % atoms(n_conf,:) % xyz_old(3)

          enddo r2lp

        enddo r1lp

        grid_rot1 % points(:) % grid_xyz(1) = grid_rot1 % points(:) % grid_xyz(1) / mol2 % rot_vector( n_conf )
        grid_rot1 % points(:) % grid_xyz(2) = grid_rot1 % points(:) % grid_xyz(2) / mol2 % rot_vector( n_conf )
        grid_rot1 % points(:) % grid_xyz(3) = grid_rot1 % points(:) % grid_xyz(3) / mol2 % rot_vector( n_conf )
        
      enddo clp
      
      call system_clock(tend)

      time = real( tend - tini, DP ) / real( rate, DP ) 

      write( 79,'(i7,3x,i7," of ", i7,1x,f10.3)' ) n_trans, count( atom_overlap(:,:,:,n_trans) ), conf_total, time-t0

      flush( 79 )

      t0 = time

      min_ener_t( n_trans ) = minval( inter_energy(:,:,:,n_trans) )

    enddo tlp

    select case ( wrtxtc )

      case ( "yes", "true", "T" )

        call xtc_out % close

    end select

    close( 79 )
    close( 271 )

    min_ener = minval( inter_energy )

    do n_trans = 1, grid_trans % numpoint

      select case ( potential )

        case ( "lj-coul", "bh-coul" )

          Zrot( n_trans ) = sum( dexp( -( inter_energy(:,:,:,n_trans) - min_ener_t(n_trans) ) / kBT ) )
            
          sumVexpVrot( n_trans ) = sum( ( inter_energy(:,:,:,n_trans) - min_ener_t(n_trans) ) * &
            dexp( -( inter_energy(:,:,:,n_trans) - min_ener_t(n_trans) ) / kBT ) ) 

      end select

    enddo 

    write( output_unit,'(T76,A)' ) "DONE"

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
    use mod_input_read     , only : axis_rot_moves, nconf2, temp, potential, inter_energy
    use mod_grids          , only : Zrot, sumVexpVrot, Check_Moves, grid_trans, grid_rot1
    use mod_inquire        , only : Inquire_file, Get_new_unit       
    use mod_error_handling

    implicit none

    integer                           :: n_rot2, n_rot1, n_conf, n_trans
    integer                           :: ios         
    integer                           :: file_unit   
    character( len = : ), allocatable :: file_format 
    character( len = : ), allocatable :: file_access 
    character( len = : ), allocatable :: file_name   
    character( len = : ), allocatable :: file_status 
    character( len = 64 )             :: line1
    integer                           :: ierr
    type(error)                       :: err

    allocate( Zrot( grid_trans % numpoint ), stat=ierr )
    if ( ierr/= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( min_ener_t( grid_trans % numpoint ), stat=ierr )
    if ( ierr/= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( sumVexpVrot( grid_trans % numpoint ), stat=ierr )
    if ( ierr/= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    allocate( inter_energy( axis_rot_moves, grid_rot1 % numpoint , nconf2, grid_trans % numpoint ), stat=ierr )
    if ( ierr/= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    kBT = KB * temp
      
    ios           = 0
    Zrot          = 0.0_DP
    min_ener      = 0.0_DP
    sumVexpVrot   = 0.0_DP
    inter_energy  = 0.0_DP

    if ( potential == "none" ) then

      file_unit   = Get_new_unit(10)
      file_format = "formatted"
      file_access = "sequential"
      file_name   = "energy.log"
      file_status = "old"

      call Inquire_file( file_unit , file_name , file_status, file_format , file_access )

      read( file_unit, *, iostat = ios ) line1
  
      if ( ios /= 0 ) then

        call err % error( 'e', message = "Error while reading energy.log file reader." )

        stop

      endif

      do n_trans = 1, grid_trans % numpoint

        do n_conf = 1, nconf2

          do n_rot1 = 1, grid_rot1 % numpoint

            do n_rot2 = 1, axis_rot_moves

              read( file_unit, *, iostat = ios ) inter_energy( n_rot2, n_rot1, n_conf, n_trans )

              if ( ios /= 0 ) then

                call err % error( 'e', message = "Error while reading energy.log file." )

                stop

              endif

            enddo 

          enddo 

        enddo

        min_ener_t(n_trans) = minval( inter_energy(:,:,:,n_trans) )

      enddo 

      close(file_unit)

    else

      file_unit   = Get_new_unit(10)
      file_format = "unformatted"
      file_access = "stream"
      file_name   = "energy.bin"
      file_status = "old"

      call Inquire_file( file_unit , file_name , file_status, file_format , file_access )
        
      do n_trans = 1, grid_trans % numpoint

        do n_conf = 1, nconf2

          do n_rot1 = 1, grid_rot1 % numpoint

            do n_rot2 = 1, axis_rot_moves

              read( file_unit, iostat = ios ) inter_energy( n_rot2, n_rot1, n_conf, n_trans )

              if ( ios /= 0 ) then

                call err % error( 'e', message = "Error while reading energy.bin file." )

                stop

              endif

            enddo 

          enddo

        enddo 

        min_ener_t(n_trans) = minval( inter_energy(:,:,:,n_trans) )

      enddo 

      close( file_unit )

    endif

    min_ener = minval( inter_energy )

    do n_trans = 1, grid_trans % numpoint

      Zrot(n_trans) = sum( dexp( -( inter_energy(:,:,:,n_trans) - min_ener_t(n_trans) ) / kBT ) )
            
      sumVexpVrot(n_trans) = sum( ( inter_energy(:,:,:,n_trans) - min_ener_t(n_trans) ) &
        * dexp( -( inter_energy(:,:,:,n_trans) - min_ener_t(n_trans) ) / kBT ) )

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
    use mod_input_read     , only : temp, axis_rot_moves, nconf2, inter_energy
    use mod_grids          , only : probT, Zrot, A, Eavg, mTS, sumVexpVrot, grid_trans, grid_rot1
    use mod_error_handling

    implicit none

    integer                           :: n_trans
    character( len = : ), allocatable :: write_fmt
    real( kind = DP )                 :: mTS_ideal, mTS_ideal_t
    integer                           :: ierr
    type( error )                     :: err

    kBT = KB * temp

    allocate( A( grid_trans % numpoint ), stat=ierr )
    if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
      
    allocate( mTS( grid_trans % numpoint ), stat=ierr )
    if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
      
    allocate( Eavg( grid_trans % numpoint ), stat=ierr )
    if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
      
    allocate( probT( grid_trans % numpoint ), stat=ierr )
    if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    A           = 0.0_DP
    Eavg        = 0.0_DP
    mTS         = 0.0_DP
    probT       = 0.0_DP
    mTS_ideal   = -kBT * log ( real(grid_trans % numpoint * nconf2 * grid_rot1 % numpoint * axis_rot_moves) )
    mTS_ideal_t = -kBT * log ( real(nconf2 * grid_rot1 % numpoint * axis_rot_moves) )

    open( unit = 20, file = 'output.log', position = 'append', status = 'replace' )

    write_fmt = '(a1,2x,3(a11,1x),3x,a5,3x,4(a13,3x))'

    write( 20, * ) grid_trans % numpoint
    write( 20, write_fmt ) " ", "X (A)", "Y (A)", "Z (A)", "point", "PROB", "A (kJ/mol)", "-TS (kJ/mol)", "E (kJ/mol)"

    write( output_unit, * )
    write( output_unit, '(T5,a23)', advance = "no" ) "CALCULATING Z(N,V,T)..."   

    Ztrans = sum( dexp( -( inter_energy - min_ener ) / kBT ) ) 

    sumVexpVtrans = sum( ( inter_energy - min_ener ) * dexp( -( inter_energy - min_ener ) / kBT ) ) 

    write_fmt = '(A,2x,3(f11.5,1x),3x,i5,3x,4(es13.5E3,3x))'

    do n_trans = 1, grid_trans % numpoint

      probT( n_trans ) = Zrot( n_trans ) / Ztrans * dexp( ( min_ener - min_ener_t( n_trans ) ) / kBT ) 

      A( n_trans ) = -kBT * dlog( Zrot( n_trans ) ) + min_ener_t( n_trans ) - mTS_ideal_t

      Eavg( n_trans ) = sumVexpVrot( n_trans ) / Zrot( n_trans ) + min_ener_t( n_trans )

      mTS( n_trans ) = A( n_trans ) - Eavg( n_trans )

      write( 20, write_fmt ) "X", grid_trans % points( n_trans ) % grid_xyz(:), n_trans, probT( n_trans ), &
        A( n_trans ), mTS( n_trans ), Eavg( n_trans )

    enddo 

    ATOTAL = -kBT * dlog( Ztrans ) + min_ener - mTS_ideal

    ETOTALavg = sumVexpVtrans / Ztrans + min_ener

    mTSTOTAL = ATOTAL - ETOTALavg

    write_fmt = '(9x,"TOTAL OVER TRANSLATIONAL GRID",12x,5(es13.5E3,3x))'

    write( 20, '(a)' ) repeat('-',152)
    write( 20, write_fmt ) sum(probT), ATOTAL, mTSTOTAL, ETOTALavg

    close( 20 )

    write( output_unit, '(T70,A)' ) "DONE"

    return
  end subroutine Calculate_ztotal

end module mod_loops
