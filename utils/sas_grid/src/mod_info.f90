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
  use iso_fortran_env , only : stdout => output_unit
  use mod_constants   , only : dashline, DP
  ! TODO: include in makefile 
  !   include 'revision.inc'     

  implicit none
  private

  public :: Display_header, display_date_time

contains

  subroutine Display_header()

    implicit none

    character( len = 16 ), parameter :: version = 'beta'
    character(len=7), parameter      :: revision = 'ae2e8ac'
      
    write(stdout,'(T3, A)') dashline
    write(stdout,'(T50, A, A)') "SAS_GRID"
    write(stdout,'(/,T8, A)')" A code to obtain the solvent accessible surface (SAS) around a given molecular structure" 
    write(stdout,'(/,T40, A, A)') "Author: Felippe M. Colombari"
    write(stdout,'(T50, A, A)') "Brazil"
    write(stdout,'(/,T42, A, A)') "Program version: ",trim(version)
    write(stdout,'(T44, A, A)') "Revision: ", revision
    write(stdout,'(T3, A)') dashline
    write(stdout,'(/,T5, A)') "Contributions" 
    write(stdout,'(T8, A, A)') "Asdrubal Lozada: ", "Error handling revision"
    write(stdout,'(/,T5, A)') "This version of SAS_GRID uses:"
    write(stdout,'(T8, A)') "SPHERE_GRID library (John Burkardt)"
    write(stdout,'(/,T3, A)') dashline
      
  end subroutine Display_header

  subroutine Display_date_time( string )

    implicit none

    character( len = * ), intent(IN) :: string
    integer,dimension(8)             :: values
    
    call Date_and_time( VALUES = values )
    write(stdout,'(/, T5, A, i2.2, "/", i2.2, "/", i4, " - ", i2.2, ":", i2.2, ":", i2.2)') string,   &
                              &values(3), values(2), values(1), values(5), values(6), values(7)


    end subroutine Display_date_time
         
end module mod_info
