!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2021 Themis developers
!
!   This file was written by Felippe M. Colombari.
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
!> @file   themis.f90
!> @author Felippe M. Colombari
!> @brief  Main module of THEMIS 
!> @date - Jun, 2017                                                           
!> - initial test version (dimer.x)
!> @date - Oct, 2017
!> - modular version
!> @date - Jan, 2020
!> - module reorganization
!> @date - Apr, 2020
!> - major revision
!> @date - May, 2021
!> - added support for ensemble of structures of molecule 2
!> @date - Jul, 2022
!> - added support to PDB files (read and write)
!---------------------------------------------------------------------------------------------------

program themis
  use mod_info              , only : Display_header, Display_date_time
  use mod_constants         , only : DP
  use mod_error_handling    , only : Normal_termination, Raise_error, error
  use mod_cmd_line          , only : Parse_arguments, irun, grid_type, rad, grid_transl, use_plot
  use mod_input_read        , only : Read_input_file, ref1, vector1, ref2, vector2, nconf2, &
                                     point_rot_factor, trans_factor, potential, file_type
  use mod_read_molecules    , only : mol1, mol2
  use mod_grids             , only : grid_trans, grid_rot1
  use mod_loops             , only : Run_loops, Rerun_loops, Calculate_ztotal
  use mod_resume            , only : Ending, Sort_output
  use mod_write_vmd         , only : Write_VMD_files
  use mod_search_structures , only : Count_structures, Sort_energy, Search_structures
  use mod_deallocate_all    , only : Deallocate_arrays
  use mod_plot             

  implicit none
  
  real( kind = DP )                :: timet
  integer                          :: finish, start, rate2, n_conf
  type( error )                    :: err
  
  call display_header()
  call Parse_arguments
  call Display_date_time( "BEGAN AT: " )
  call System_clock( start, rate2 )

  call Read_input_file

  SELECT CASE ( irun )

    CASE ( "run" )

      call mol1 % Read_molecule ( "conf1."//file_type, 1, file_type )
      call mol2 % Read_molecule ( "conf2."//file_type, nconf2, file_type )
      
      call mol1 % Check_molecule ( ref1, vector1, 1 )
      call mol1 % Translate_molecule( ref1, 1 )
      call mol1 % Align_molecule( vector1, ref1, 1 )
      call mol1 % Rotate_molecule( 1, 1 )

      SELECT CASE (grid_type)

        CASE ( "shell" )

          call grid_trans % Build_translation_sphere( trans_factor, rad )
        
        CASE ( "user" ) 

          call grid_trans % Read_grid ( grid_transl )
          call grid_trans % Translate_grid 
          call grid_trans % Align_grid 
          call grid_trans % Rotate_grid 

      end SELECT

      call mol2 % Check_molecule ( ref2, vector2, 2 )

      do n_conf = 1, nconf2

        call mol2 % Translate_molecule( ref2, n_conf )
        call mol2 % Align_molecule( vector2, ref2, n_conf )

      enddo

      call grid_rot1 % Build_rotation1_sphere( point_rot_factor )
      call Run_loops

      SELECT CASE ( potential )

        CASE ( "lj-coul", "bh-coul" )

          call Calculate_ztotal

          call System_clock( finish )
  
          timet = real(finish - start) / real(rate2)

          call Write_vmd_files
          call Count_structures
          call Sort_energy
          call Search_structures
          call Ending( timet )
          call Sort_output
          call Deallocate_arrays

        CASE ( "none" )

          call Ending( timet )

      END SELECT

    CASE ( "rerun" )

      call mol1 % Read_molecule ( "conf1."//file_type, 1, file_type )
      call mol2 % Read_molecule ( "conf2."//file_type, nconf2, file_type )
      
      call mol1 % Translate_molecule( ref1, 1 )
      call mol1 % Align_molecule( vector1, ref1, 1 )
      call mol1 % Rotate_molecule( 1, 1 )
    
      SELECT CASE ( grid_type )

        CASE ( "shell" )

          call grid_trans % Build_translation_sphere( trans_factor, rad )

        CASE ( "user" ) 

          call grid_trans % Read_grid ( grid_transl )
          call grid_trans % Translate_grid 
          call grid_trans % Align_grid 
          call grid_trans % Rotate_grid 

      END SELECT
      
      call mol2 % Translate_molecule( ref2, nconf2 )
      
      do n_conf = 1, nconf2

        call mol2 % Align_molecule( vector2, ref2, nconf2 )

      enddo

      call grid_rot1 % Build_rotation1_sphere( point_rot_factor )
      call Rerun_loops
      call Calculate_ztotal

      call System_clock( finish )
  
      timet = real(finish - start) / real(rate2)

      call Write_vmd_files
      call Sort_energy
      call Search_structures
      call Ending( timet )
      call Sort_output
      call Deallocate_arrays

  END SELECT

  if(use_plot) Call Call_plot()

  call Display_date_time( "FINISHED AT: " )           
  call err % termination(0,'f')

end program
