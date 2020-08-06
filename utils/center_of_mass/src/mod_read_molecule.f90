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
!> @file   mod_read_molecule.f90
!> @author Felippe M. Colombari
!> @brief  This module reads coordinates for a given .xyz file and sets the mass for each atom
!> @date - Jan, 2020                                                           
!> - independent module created                                                
!> @date - Jan 2020
!> - update error condition by error_handling module added 
!---------------------------------------------------------------------------------------------------

module mod_read_molecule
  use mod_constants

  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type atom 
    character(len=3)          :: label
    real(kind=dp)             :: xyz(3)
  end type atom

  type molecule
    type(atom), allocatable, dimension(:)  :: atoms
    integer                                :: num_atoms
  contains
    procedure, pass                        :: Read_molecule
  end type molecule

  type( molecule )                         :: mol

  logical, allocatable, dimension(:) :: mass_found

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  !---------------------------------------------------------------------------
  !> @brief This routine reads coordinates for xyz_files.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Read_molecule ( this, molecule_filename )
    use mod_inquire, only: Inquire_file
    use mod_error_handling

    implicit none

    class( molecule ), intent(inout) :: this
    character( len = * ), intent(in) :: molecule_filename
    integer                          :: i
    integer                          :: ios         = 0
    integer                          :: file_unit   = 10        
    character( len = 15 )            :: file_format = "formatted"
    character( len = 15 )            :: file_access = "sequential"

    integer                          :: ierr
    type(error)                      :: err

    call Inquire_file( file_unit, molecule_filename, file_format, file_access )

    read(file_unit,*,iostat=ios) this % num_atoms
    read(file_unit,*)

    if ( allocated ( this % atoms ) ) deallocate ( this % atoms )
    allocate( this % atoms( this % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    do i = 1, this % num_atoms

      read(file_unit,*,iostat=ios) this % atoms(i) % label, this % atoms(i) % xyz(:)

    enddo

    write(666,*) this % num_atoms
    write(666,*)

    do i = 1, this % num_atoms

      write(666,'(a3,3(f12.7))') this % atoms(i) % label, this % atoms(i) % xyz(:)

    enddo

    close(file_unit)

    return
  end subroutine Read_molecule

end module mod_read_molecule
