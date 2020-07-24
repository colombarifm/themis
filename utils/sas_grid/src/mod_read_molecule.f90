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
!> @file   mod_read_molecule.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  This module reads coordinates for a given .xyz file and sets the vdW radii for each atom
!> @date - Nov, 2019                                                           
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_read_molecule
  use mod_constants, only: DP, PI

  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type atom 
    real( kind = dp )                       :: xyz(3), vdw
    character( len = 2 )                    :: label
  end type atom

  type molecule
    type( atom ), allocatable,dimension(:) :: atoms
    integer                                :: num_atoms
  contains
    procedure, pass                        :: Read_molecule
    procedure, pass                        :: Read_vdw_radii
  end type molecule

  type( molecule )                         :: mol

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type atom_entry
    character( len = 2 )                   :: symbol 
    real( kind = DP )                      :: radii 
  end type atom_entry

  type atom_list
    integer                                     :: n_list
    type(atom_entry), allocatable, dimension(:) :: listed_atoms
  end type atom_list

  type(atom_list)   :: list

  logical, allocatable, dimension(:) :: radii_found

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

    close(file_unit)

    return
  end subroutine Read_molecule

  subroutine Read_vdw_radii( this )
    use mod_cmd_line , only : radius

    implicit none

    class( molecule ), intent(inout)     :: this
    integer                              :: i, j

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                  !!
  !!   READ VDW RADII TABLE GIVEN IN vdw-radii.txt    !!
  !!                                                  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  inquire( file = 'vdw-radii.txt' )

  open(10,file='vdw-radii.txt',status='old')
  
  read(10,*) list % n_list
 
  allocate( list % listed_atoms( list % n_list ) )

  do i = 1, list % n_list
  
    read(10,*) list % listed_atoms(i) % symbol, list % listed_atoms(i) % radii

  enddo
      
  close (10)
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                  !!
  !!     READ COORDINATES FROM REFERENCE MOLECULE     !!
  !! AND COMPARE SYMBOLS WITH THOSE IN vdw-radii.txt. !!
  !!   THEN ASSOCIATE EACH ATOM WITH ITS RESPECTIVE   !!
  !!                    VDW RADII.                    !!
  !!                                                  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate( radii_found( mol % num_atoms ) )

  radii_found = .false. 

  do i = 1, mol % num_atoms

    do j = 1, list % n_list

      if ( mol % atoms(i) % label .eq. list % listed_atoms(j) % symbol ) then
        
        this % atoms(i) % vdw = list % listed_atoms(j) % radii + radius

        radii_found(i) = .true.

        write(*,'(T5, "Atom ", i4.4,": ", a2, " found!", f7.3)') i, this % atoms(i) % label, &
        list % listed_atoms(j) % radii

      endif

    enddo
    
    if ( radii_found(i) .eqv. .false. ) then

      write(*,'(T5, "Atom ", i4.4,": ", a2, " not found!")') i, this % atoms(i) % label
      write(*,*)
      write(*,'(T5, "Check your vdw-radii.txt file for the required entry.")')

      stop

    endif

  enddo
      
  end subroutine Read_vdw_radii

end module mod_read_molecule
