!------------------------------------------------------------------------------
! THEMIS: 
! Copyright (C) 2017 Felippe M. Colombari
!------------------------------------------------------------------------------
!> @brief This module contains the main program loop for Z(N,V,T) calculation
!> @author Felippe M. Colombari 
!> - Laboratório de Química Teórica, LQT -- UFSCar
!> @date - Oct, 2017 
!> - module created
!------------------------------------------------------------------------------

module MOD_LOOPS
  use MOD_CONSTANTS

  implicit none

  character( len = 4 )                         :: ntra, nreo, ngyr
  real( kind = DP )                            :: t0, time, kBT, min_ener
  real( kind = DP )                            :: ATOTAL, mTSTOTAL, ETOTALavg
  real( kind = DP )                            :: Ztrans, sumVexpVtrans
  real( kind = DP ), allocatable, dimension(:) :: min_ener_t
  integer                                      :: tini, tend, rate

  contains

    !---------------------------------------------------------------------------
    !> @brief This routine performs the main loop and calls subroutine s to 
    !>  perform the configurational sampling and energy calculations.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Jun, 2017
    !> - subroutine  created
    !> - energy values written in binary file
    !> @date - Dec, 2017
    !> - lowest energy values are subtracted from each value to avoid problems
    !>   due to exponential explosion
    !---------------------------------------------------------------------------	

    subroutine CALC
      use MOD_CMD_LINE
      use MOD_INPUT_READ
      use MOD_READ_MOLECULES
      use MOD_READ_GRIDS
      use MOD_POT_LJC
      use MOD_POT_BHC
      use MOD_POT_ljc_pair
      use XDR,           only: xtcfile
      use mod_error_handling

      implicit none

      type( xtcfile )                   :: xtc_out
      real( kind = SP ), dimension(:,:) :: pos(3,mol1 % num_atoms + mol2 % num_atoms)
      real( kind = SP ), dimension(3,3) :: box
      real( kind = SP )                 :: frametime, prec
      integer                           :: frame, i, j, g, r, t

      real( kind = DP )                 :: lowest
      character( len = 20 )             :: prefix

      type( ljc_dimer ), target         :: ljc_target
      type( bhc_dimer ), target         :: bhc_target
      type( ljc_pair_dimer ), target    :: ljc_pair_target
      type( dimer ), target             :: none_target
      class( dimer ), pointer           :: potential_pointer

      integer                                      :: ierr
      type(error)                                  :: err

      allocate( inter_energy( gyr_factor, grid_reo % numpoint, grid_trans % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
      allocate( atom_overlap( gyr_factor, grid_reo % numpoint, grid_trans % numpoint ), stat=ierr )
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

      kBT = kB * TEMP
      
      SELECT CASE (wrtxtc)

        CASE ("yes", "true", "T")

          box(1,1) = 0.0_SP 
          box(2,2) = 0.0_SP
          box(3,3) = 0.0_SP

          box(2,1) = 0.0_SP 
          box(3,1) = 0.0_SP 

          box(1,2) = 0.0_SP 
          box(3,2) = 0.0_SP

          box(1,3) = 0.0_SP
          box(2,3) = 0.0_SP

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

        CASE ("ljc_pair")

          potential_pointer => ljc_pair_target

          call ljc_pair_dimers % Read_ljc_pair_params ( "parameters" )

        CASE ("none")

          potential_pointer => none_target

      end SELECT

      call Check_moves
      call potential_pointer % Calc_cross

      write(*,'(/,T5,a17)',advance='no') "ENTERING LOOPS..."
    
      call system_clock(tini, rate)

      tlp: do t = 1, grid_trans % numpoint

        call mol2 % Translate_molecule( ref2 )
      
        grid_reo % points(:) % grid_xyz(1) = grid_reo % points(:) % grid_xyz(1) * mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(2) = grid_reo % points(:) % grid_xyz(2) * mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(3) = grid_reo % points(:) % grid_xyz(3) * mol2 % rot_vector
          
        rlp: do r = 1, grid_reo % numpoint

          costhetar = dcos(thetar(r))
          sinthetar = dsin(thetar(r))
          cosphir   = dcos(phir(r))
          sinphir   = dsin(phir(r))

          glp: do g = 1, gyr_factor

            call mol2 % Rotate_molecule( g ) 

            do j = 1, mol2 % num_atoms

              mol2 % atoms(j) % xyz(:) = mol2 % atoms(j) % xyz(:) + grid_trans % points(t) % grid_xyz(:)

            enddo

            call potential_pointer % Calc_energy( g, r, t )

            write(271) inter_energy( g, r, t )

            if ( atom_overlap( g, r, t ) .eqv. .false. ) then
            
              write(prefix,'("point_",I4.4,"_",I4.4,"_",I4.4)') t, r, g 
              
              SELECT CASE (writeframe)

                CASE ("XYZ")

                  open( unit = 66, file = prefix//'.xyz', status = 'unknown' )

                  lowest = 0.0_DP

                  call dimers % Build_dimer
             
                  call dimers % Write_xyz( lowest )
                            
                  close(66)

                CASE ("MOP")
  
                  open( unit = 66, file = prefix//'.mop', status = 'unknown' )
        
                  call dimers % Build_dimer
          
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

          enddo glp

        enddo rlp

        grid_reo % points(:) % grid_xyz(1) = grid_reo % points(:) % grid_xyz(1) / mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(2) = grid_reo % points(:) % grid_xyz(2) / mol2 % rot_vector
        grid_reo % points(:) % grid_xyz(3) = grid_reo % points(:) % grid_xyz(3) / mol2 % rot_vector
          
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

          CASE ("lj-coul", "bh-coul", "ljc_pair")

            Zrot(t) = sum( dexp( -( inter_energy(:,:,t) - min_ener_t(t) ) / kBT ) )
            
            sumVexpVrot(t) = sum( ( inter_energy(:,:,t) - min_ener_t(t) ) * &
              dexp( -( inter_energy(:,:,t) - min_ener_t(t) ) / kBT ) ) 

        end SELECT

      enddo 

      write(*,'(T50,A)') "DONE"

      return
    end subroutine CALC

    !---------------------------------------------------------------------------
    !> @brief This routine performs the main loop and reads already calculated
    !>  energy values from energy file.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Jun, 2017
    !> - subroutine  created
    !> - energy values read from binary file
    !> @date - Dec, 2017
    !> - lowest energy values are subtracted from each value to avoid problems
    !>   in exponential calculation
    !---------------------------------------------------------------------------	

    subroutine RECALC
      use MOD_INPUT_READ, only: gyr_factor, TEMP, potential, inter_energy, scale_factor
      use MOD_READ_GRIDS, only: Zrot, sumVexpVrot, &
                                Check_Moves, grid_trans, grid_reo
      use MOD_INQUIRE,    only: Inquire_file             
      use mod_error_handling

      implicit none

      integer                           :: g, r, t
      integer                           :: ios         = 0
      integer                           :: file_unit   
      character( len = : ), allocatable :: file_format 
      character( len = : ), allocatable :: file_access 
      character( len = : ), allocatable :: file_name   
      character( len = 64 )             :: line1
      integer                           :: ierr
      type(error)                       :: err

      allocate(        Zrot( grid_trans % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      allocate(  min_ener_t( grid_trans % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      allocate( sumVexpVrot( grid_trans % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      allocate( inter_energy( gyr_factor, grid_reo % numpoint , grid_trans % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      kBT = kB * TEMP
      
      ios         = 0
      Zrot        = 0.0_DP
      min_ener    = 0.0_DP
      sumVexpVrot = 0.0_DP

      inter_energy  = 0.0_DP

      if ( potential == "none" ) then

        file_unit   = 16
        file_format = "formatted"
        file_access = "sequential"
        file_name   = "energy.log"

        call Inquire_file( file_unit, file_name, file_format, file_access )

        read(file_unit,*,iostat=ios) line1
  
        if ( ios /= 0 ) then

          call err%error('e',message="Error while reading energy.log file reader.")

          stop

        endif

        do t = 1, grid_trans % numpoint

          do r = 1, grid_reo % numpoint

            do g = 1, gyr_factor

              read(file_unit,*,iostat=ios) inter_energy( g, r, t )

              inter_energy( g, r, t ) = inter_energy( g, r, t ) / scale_factor

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

        file_unit   = 271
        file_format = "unformatted"
        file_access = "stream"
        file_name   = "energy.bin"

        call Inquire_file( file_unit, file_name, file_format, file_access )
        
        do t = 1, grid_trans % numpoint

          do r = 1, grid_reo % numpoint

            do g = 1, gyr_factor

              read(file_unit, iostat=ios) inter_energy( g, r, t )

              inter_energy( g, r, t ) = inter_energy( g, r, t ) / scale_factor
              
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
    end subroutine RECALC

    !---------------------------------------------------------------------------
    !> @brief This routine performs partition function calculation for each grid
    !>  point and also for the whole grid.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Jun, 2017
    !> - subroutine  created
    !> @date - Dec, 2017
    !> - absolute energy values are recovered and thermodynamic properties are
    !>   obtained from them
    !> @date - Mar, 2019
    !> - -TdS from ideal gas is subtracted from free energy
    !---------------------------------------------------------------------------	

    subroutine CALC_ZTOTAL
      use MOD_INPUT_READ, only: TEMP, gyr_factor, inter_energy
      use MOD_READ_GRIDS, only: probT, Zrot, A, Eavg, mTS, sumVexpVrot, grid_trans, grid_reo
      use mod_error_handling

      implicit none

      integer                       :: t
      character(len=:), allocatable :: write_fmt
      real(kind=DP)                 :: mTS_ideal, mTS_ideal_t
      integer                       :: ierr
      type(error)                   :: err

      kBT = kB * TEMP

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
      mTS_ideal     = -kBT * log ( real(grid_trans % numpoint * grid_reo % numpoint * gyr_factor) )
      mTS_ideal_t   = -kBT * log ( real(grid_reo % numpoint * gyr_factor) )

      open( unit = 20, file = 'output.log', position = 'append', status = 'replace' )

      write_fmt = '(a1,2x,3(a11,1x),3x,a5,3x,4(a13,3x))'

      write(20,*) grid_trans % numpoint
      write(20, write_fmt ) " ", "X (A)", "Y (A)", "Z (A)", "point", "PROB", &
        &"A (kJ/mol)", "-TS (kJ/mol)", "E (kJ/mol)"

      write(*,*)
      write(*,'(T5,a23)',advance="no") "CALCULATING Z(N,V,T)..."   

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

      write(*,'(T44,A)') "DONE"

      return
    end subroutine CALC_ZTOTAL

end module MOD_LOOPS
