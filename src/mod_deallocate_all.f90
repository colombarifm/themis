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
!> @file   mod_deallocate_all.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  Deallocates all arrays prior to program termination
!> @date - Dec, 2017                                                           
!> - independent module created                                                
!> @note
!> - to check for memory leaks and/or final status of allocatable arrays please use: \n 
!>   valgrind --leak-check=full --show-leak-kinds=all -v ./themis [ options ]
!> @note
!> - to check memory usage along execution, please use: \n
!>   valgrind --tool=massif ./themis [options]
!---------------------------------------------------------------------------------------------------

module mod_deallocate_all

  implicit none

  private
  public Deallocate_arrays

contains

  subroutine Deallocate_arrays
    use MOD_CMD_LINE       , only : irun
    use MOD_INPUT_READ     , only : potential, atom_overlap, inter_energy
    use MOD_READ_MOLECULES , only : mol1, mol2
    use MOD_GRIDS          , only : grid_trans, grid_rot1, rr, thetar, phir, Zrot, sumVexpVrot, probT
    use MOD_POT_LJC        , only : mol1_ljc, mol2_ljc, ljc_dimer
    use MOD_POT_BHC        , only : mol1_bhc, mol2_bhc, bhc_dimer
    use MOD_LOOPS          , only : min_ener_t

    implicit none

    !!!!! ARRAYS FROM MOD_READ_XYZ !!!!!

    deallocate( mol1 % atoms )
    deallocate( mol2 % atoms )

    !!!!! ARRAYS FROM MOD_INPUT_READ !!!!!

    if ( allocated(atom_overlap) ) deallocate( atom_overlap )
    deallocate( inter_energy )

    !!!!! ARRAYS FROM MOD_READ_GRIDS

    deallocate( grid_trans % points )
    deallocate( grid_rot1 % points )

    deallocate( rr )
    deallocate( thetar )
    deallocate( phir )
    deallocate( Zrot )
    deallocate( sumVexpVrot )
    deallocate( probT )
    
    SELECT CASE (potential)

      CASE ("lj-coul")
  
        !!!!! ARRAYS FROM MOD_POT_LJC.F90 !!!!!

        if ( irun .eq. "run" ) then

          deallocate( mol1_ljc % ljc_atoms )
          deallocate( mol2_ljc % ljc_atoms )

        endif

      CASE ("bh-coul")
  
        !!!!! ARRAYS FROM MOD_POT_BHC.F90 !!!!!

        if ( irun .eq. "run" ) then

          deallocate( mol1_bhc % bhc_atoms )
          deallocate( mol2_bhc % bhc_atoms )

        endif

      CASE ("none")

        CONTINUE

    end SELECT

    !!!!! ARRAY FROM MOD_LOOPS.F90 !!!!!

    deallocate( min_ener_t )

    return
  end subroutine Deallocate_arrays

end module mod_deallocate_all
