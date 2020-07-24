!---------------------------------------------------------------------------------------------------
! SAS_GRID: A code to obtain the solvent accessible surface (SAS) around a given molecular structure                                                  
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
!> @file   sas_grid.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  Mais module of sas_grid 
!> @date - ~2018
!> - first version
!> @date - Nov, 2019                                                           
!> - modular version                                                
!---------------------------------------------------------------------------------------------------

program sas_grid

  use mod_info              , only : Display_header, Display_date_time
  use mod_constants         , only : dashline
  use mod_cmd_line          , only : Parse_arguments, filename, factor, radius
  use mod_read_molecule     , only : mol 
  use mod_grids             , only : grid_sphere, grid_sas 
  use mod_deallocate_all    , only : Deallocate_arrays

  implicit none

  call display_header()

  call Parse_arguments

  call mol % Read_molecule( filename )

  call mol % Read_vdw_radii

  call grid_sphere % Build_sphere_grid( factor, radius )

  call grid_sas % Build_sas_grid

  call Deallocate_arrays

  write(*,'(/, T3, A)') dashline

end program sas_grid
