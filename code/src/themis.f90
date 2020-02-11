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
!> @file   themis.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  Main module of THEMIS 
!> @date - Jun, 2017                                                           
!> - initial test version (dimer.x)
!> @date - Oct, 2017
!> - modular version
!---------------------------------------------------------------------------------------------------

program themis
  use mod_info              , only : Display_header, Display_date_time
  use mod_constants         , only : DP, dashline
  use mod_cmd_line          , only : Parse_arguments, irun, grid_type, rad, grid_transl
  use mod_input_read        , only : Read_input_file, ref1, vector1, ref2, vector2, reo_factor, trans_factor, potential
  use mod_read_molecules    , only : mol1, mol2
  use mod_grids             , only : grid_trans, grid_reo
  use mod_loops             , only : Run_loops, Rerun_loops, Calculate_ztotal
  use mod_resume            , only : Ending, Sort_output
  use mod_write_vmd         , only : Write_VMD_files
  use mod_search_structures , only : Count_structures, Sort_energy, Search_structures
  use mod_deallocate_all    , only : Deallocate_arrays

  implicit none
  
  real( kind = DP )                :: timet
  integer                          :: finish, start, rate2
  
  call display_header()
  call Display_date_time( "BEGAN AT: " )
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
      call Run_loops

      SELECT CASE ( potential )

        CASE ( "lj-coul", "bh-coul", "ljc_pair" )

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

  end SELECT

  call Display_date_time( "FINISHED AT: " )

  write(*,'(/, T3, A)') dashline

end program
