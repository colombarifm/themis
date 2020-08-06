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
!> @file   mod_deallocate_all.f90
!> @author Felippe M. Colombari
!> @brief  Deallocates all arrays prior to program termination
!> @date - Jan, 2020                                                           
!> - independent module created                                                
!> @note
!> - to check for memory leaks and/or final status of allocatable arrays please use: \n 
!>   valgrind --leak-check=full --show-leak-kinds=all -v ./com [ options ]
!> @note
!> - to check memory usage along execution, please use: \n
!>   valgrind --tool=massif ./com [options]
!---------------------------------------------------------------------------------------------------

module mod_deallocate_all

  implicit none

  private
  public Deallocate_arrays

contains

  subroutine Deallocate_arrays
    use mod_read_molecule , only : mol, mass_found

    implicit none

    !!!!! ARRAYS FROM MOD_READ_XYZ !!!!!

    deallocate( mol % atoms )
    deallocate( mass_found )

    return
  end subroutine Deallocate_arrays

end module mod_deallocate_all
