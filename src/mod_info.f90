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
!> @file   mod_info.f90
!> @author Asdrubal Lozada-Blanco, Felippe Mariano Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  aslozada@gmail.com, colombarifm@hotmail.com
!> @brief  This module contains the program header
!> @date - Dec, 2019 
!> - module created and incorporated into code
!---------------------------------------------------------------------------------------------------

module mod_info
  use iso_fortran_env, only: output_unit
  use mod_constants, only: dashline, DP, version, revision
  ! TODO: include in makefile 
  !   include 'revision.inc'     

  implicit none
  private

  public :: Display_header, display_date_time

contains

  subroutine Display_header()

    implicit none

    write(output_unit,'(T3, A)') dashline
    write(output_unit,'(T50, A, A)') "THEMIS"
    write(output_unit,'(/,T8, A)')" A Software to Assess Association Free Energies via Direct Estimative of Partion Function"
    write(output_unit,'(/,T40, A, A)') "Author: Felippe M. Colombari"
    write(output_unit,'(T50, A, A)') "Brazil"
    write(output_unit,'(/,T42, A, A)') "Program version: ",trim(version)
    write(output_unit,'(T44, A, A)') "Revision: ", revision
    write(output_unit,'(T3, A)') dashline
    write(output_unit,'(/,T5, A)') "Contributions" 
    write(output_unit,'(T8, A, A)') "Asdrubal Lozada: ", "Error handling revision"
    write(output_unit,'(/,T5, A)') "This version of THEMIS uses:"
    write(output_unit,'(T8, A)') "SPHERE_GRID library (John Burkardt)"
    write(output_unit,'(T8, A)') "XDR library (James Barnett)"
    write(output_unit,'(/,T8, A)') "For citations, please refer to:" 
    write(output_unit,'(T10,A)') "https://people.sc.fsu.edu/~jburkardt/f_src/sphere_grid/sphere_grid.html"
    write(output_unit,'(T10,A)') "https://github.com/kmtu/xdrfort/blob/master/xdr.F90"
    write(output_unit,'(/,T3, A)') dashline
      
  end subroutine Display_header

  subroutine Display_date_time( string )

    implicit none

    character( len = * ), intent(IN) :: string
    integer,dimension(8)             :: values
    
    write(output_unit,'(T3, A)') dashline

    call Date_and_time( VALUES = values )
    write(output_unit,'(/, T5, A, i2.2, "/", i2.2, "/", i4, " - ", i2.2, ":", i2.2, ":", i2.2)') string,   &
                              &values(3), values(2), values(1), values(5), values(6), values(7)

    write(output_unit,'(/,T3, A)') dashline

  end subroutine Display_date_time
         
end module mod_info
