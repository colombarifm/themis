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
!> @file   mod_write_vmd.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  This module contains instructions to write VMD scripts for visualization
!> @date - Oct, 2017                                                           
!> - independent module created                                                
!---------------------------------------------------------------------------------------------------

module mod_write_vmd
  
  implicit none

contains

  !---------------------------------------------------------------------------
  !> @brief This routine generates .vmd script files to visualize A, E and -TS
  !>  along central grid.
  !> @author Felippe M. Colombari, Asdrubal Lozada-Blanco
  !> @date - Jun, 2017
  !> - subroutine  created
  !> @date - Dec, 2018
  !> - no PDB files needed.... \o/
  !---------------------------------------------------------------------------	
  subroutine Write_vmd_files
    use mod_constants, only: DP
    use mod_grids,     only: A, mTS, Eavg

    implicit none

    integer                             :: i
    character( len = 18 ), dimension(3) :: property
    real( kind = DP ), dimension(3)     :: min_val, max_val

    property(1) = "free-energy"
    property(2) = "entropic-penalty"
    property(3) = "energy"

    min_val(1) = minval(A)
    max_val(1) = maxval(A)
    min_val(2) = minval(mTS)
    max_val(2) = maxval(mTS)
    min_val(3) = minval(Eavg)
    max_val(3) = maxval(Eavg)

    do i = 1, 3

      open( unit = 666, file = "surf_"//trim(property(i))//".vmd", status = "replace" )

      write(666,'("# Generates the surface representation of ", A, " along translation grid", /)') trim(property(i))
      write(666,'("set min_val ", es13.5E3)') min_val(i)
      write(666,'("set max_val ", es13.5E3)') max_val(i)
      write(666,'(/, "display projection orthographic", /)')
      write(666,'("mol new output.log type xyz")')
      write(666,'("mol modcolor 0 top User")')
      write(666,'("mol colupdate 0 top 1")')
      write(666,'("mol scaleminmax top 0 $min_val $max_val")')
      write(666,'("mol modstyle 0 0 CPK 0.7 0.0 50 0", /)')
      write(666,'("set numatoms [ molinfo top get numatoms ]")')
      write(666,'("set data [ open ""output.log"" r ]")')
      write(666,'("set dummy1 [ gets $data ]")')
      write(666,'("set dummy2 [ gets $data ]", /)')
      write(666,'("for { set j 0 } { $j < ( $numatoms ) } { incr j } {")')
      write(666,'("  set ke [ gets $data ]")')
      write(666,'("  set value [ lindex $ke ", i1, " ]")') i+5
      write(666,'("  set atomsel [ atomselect top ""index $j"" ]")')
      write(666,'("  $atomsel set user $value")')
      write(666,'("  $atomsel delete")')
      write(666,'("}", /)')
      write(666,'("mol new lowest_0001.xyz type xyz waitfor all")')
      write(666,'("mol delrep 0 top")')
      write(666,'("mol modcolor 0 top Name")')
      write(666,'("mol representation VDW 1.0 50")')
      write(666,'("mol addrep top")')

      close(666)

    enddo

    return
  end subroutine Write_vmd_files

end module mod_write_vmd
