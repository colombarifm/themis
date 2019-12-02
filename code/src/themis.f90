!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief Main module of THEMIS 
!> @author Felippe M. Colombari 
!> - Laboratório de Química Teórica, LQT -- UFSCar
!> @date - Jun, 2017
!> - initial test version
!> @date - Oct, 2017 
!> - modular version
!------------------------------------------------------------------------------

program themis
  use mod_constants      , only: DP, dashline
  use mod_cmd_line       , only: Parse_arguments, irun, grid_type, rad, grid_transl
  use mod_input_read     , only: Read_input_file, ref1, vector1, ref2, vector2, reo_factor, trans_factor, potential
  use MOD_READ_MOLECULES , only: mol1, mol2
  use MOD_READ_GRIDS     , only: grid_trans, grid_reo
  use MOD_LOOPS          , only: Calc, Recalc, Calc_Ztotal
  use MOD_RESUME         , only: Ending, Sort_output
  use MOD_VMD            , only: Write_VMD_files
  use MOD_SEARCH         , only: Count_structures, Sort_energy, Search_structures
  use MOD_DEALLOCATE_ALL , only: Deallocate_arrays

  implicit none
  
  character( len = 16 ), parameter :: version = 'beta'
  integer,dimension(8)             :: values
  real( kind = DP )                :: timet
  integer                          :: finish, start, rate2
  
  call Date_and_time( VALUES = values )

  write(*,'(/, T3, A)') dashline
  write(*,'(/, T25, A, A)') " Program THEMIS version ", trim(version)
  write(*,'(/, T5, "STARTED AT: ", i2.2, "/", i2.2, "/", i4, " - ", &
                               &i2.2, ":", i2.2, ":", i2.2)')    &
                               &values(3), values(2), values(1), &
                               &values(5), values(6), values(7)

  call System_clock( start, rate2 )

  call Parse_arguments
  call Read_input_file

  SELECT CASE ( irun )

    CASE ( "RUN", "run" )

      call mol1 % Read_molecule ( "conf1.xyz" )
      call mol1 % Check_molecule ( ref1, vector1, 1 )
      call mol1 % Translate_molecule( ref1 )
      call mol1 % Align_molecule( vector1, ref1 )
      call mol1 % Rotate_molecule( 1 )

      SELECT CASE (grid_type)

        CASE ( "shell" )

          call grid_trans % Build_translation_sphere( trans_factor, rad )
        
        CASE ( "user" ) 

          call grid_trans % Read_grid ( grid_transl )
          call grid_trans % Translate_grid 
          call grid_trans % Align_grid 
          call grid_trans % Rotate_grid 

      end SELECT

      call mol2 % Read_molecule ( "conf2.xyz" )
      call mol2 % Check_molecule ( ref2, vector2, 2 )
      call mol2 % Translate_molecule( ref2 )
      call mol2 % Align_molecule( vector2, ref2 )

      call grid_reo % Build_reorientation_sphere( reo_factor )
      call Calc

      SELECT CASE ( potential )

        CASE ( "lj-coul", "bh-coul", "ljc_pair" )

          call Calc_Ztotal

          call System_clock( finish )
  
          timet = real(finish - start) / real(rate2)

          call WRITE_VMD_FILES
          call COUNT_STRUCTURES
          call SORT_ENERGY
          call SEARCH_STRUCTURES
          call ENDING( timet )
          call SORT_OUTPUT
          call DEALLOCATE_ARRAYS

        CASE ( "none" )

          call ENDING( timet )

      end SELECT

    CASE ( "RERUN", "rerun" )

      call mol1 % Read_molecule ( "conf1.xyz" )
      call mol1 % Translate_molecule( ref1 )
      call mol1 % Align_molecule( vector1, ref1 )
      call mol1 % Rotate_molecule( 1 )
    
      SELECT CASE ( grid_type )

        CASE ( "shell" )

          call grid_trans % Build_translation_sphere( trans_factor, rad )

        CASE ( "user" ) 

          call grid_trans % Read_grid ( grid_transl )
          call grid_trans % Translate_grid 
          call grid_trans % Align_grid 
          call grid_trans % Rotate_grid 

      end SELECT
      
      call mol2 % Read_molecule ( "conf2.xyz" )
      call mol2 % Translate_molecule( ref2 )
      call mol2 % Align_molecule( vector2, ref2 )

      call grid_reo % Build_reorientation_sphere( reo_factor )
      call RECALC
      call CALC_ZTOTAL

      call System_clock( finish )
  
      timet = real(finish - start) / real(rate2)

      call WRITE_VMD_FILES
      call SORT_ENERGY
      call SEARCH_STRUCTURES
      call ENDING( timet )
      call SORT_OUTPUT
      call DEALLOCATE_ARRAYS

  end SELECT

  call Date_and_time( VALUES = values )

  write(*,'(T5, "FINISHED AT: ", i2.2, "/", i2.2, "/", i4, " - ", &
                                &i2.2, ":", i2.2, ":", i2.2, /)')    &
                                &values(3), values(2), values(1), &
                                &values(5), values(6), values(7)
  
  write(*,'(T3, A)') dashline

end program
