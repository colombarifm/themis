!---------------------------------------------------------------------------------------------------
! COM: A code to calculate the center of mass of a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Themis developers
!                 Laboratory of Theoretical Chemistry (LQT) - Federal University of São Carlos 
!                 <http://www.lqt.dq.ufscar.br>
!
!   Please cite: 
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
!> @file   com.f90
!> @author Felippe M. Colombari
!> @brief  This program calculates the center of mass of a given molecular structure
!> @date - Jan, 2020                                                           
!> - independent module created                                                
!> @note 
!> - reads "filein.xyz" and writes "fileout.xyz" with an extra XX site corresponding to the center 
!>   of mass. XX coordinates can ba placed at the origin.
!---------------------------------------------------------------------------------------------------

program com

  use mod_info              , only : Display_header, Display_date_time
  use mod_cmd_line          , only : Parse_arguments, filename_molecule, filename_com, center_option
  use mod_read_molecule     , only : mol
  use mod_calc_com          , only : mol_mass, Assign_mass, Calc_com
  use mod_deallocate_all    , only : Deallocate_arrays

  implicit none

  call display_header()

  call Parse_arguments

  call mol % Read_molecule( filename_molecule )

  call mol_mass % Assign_mass

  call mol_mass % Calc_com( filename_com, center_option )

  !  call Deallocate_arrays

end program com
