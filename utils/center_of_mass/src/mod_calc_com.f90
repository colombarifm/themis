!---------------------------------------------------------------------------------------------------
! COM: A code to calculate the center of mass of a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!   Copyright 2019 Felippe M. Colombari
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
!> @file   mod_calc_com.f90
!> @author Felippe M. Colombari
!> @brief  This program calculates the center of mass of a given molecular structure
!> @date - Jan, 2020                                                           
!> - independent module created                                                
!> @note 
!> - reads "filein.xyz" and writes "fileout.xyz" with an extra XX site corresponding to the center 
!>   of mass. XX coordinates can ba placed at the origin.
!> - usage: ./com <filein.xyz> <fileout.xyz> center/nocenter
!---------------------------------------------------------------------------------------------------

module mod_calc_com
  use mod_constants
  use mod_read_molecule

  implicit none

  type, extends( atom ), private   :: mass_atom
    real(kind=dp)                  :: mass  
  end type

  type, extends( molecule ), private               :: mass_molecules
    real( kind = DP )                              :: total_mass, mass_components(3), com(3)
    type( mass_atom ), allocatable, dimension(:)   :: mass_atoms
  contains 
    procedure, pass                                :: Assign_mass
    procedure, pass                                :: Calc_com
  end type mass_molecules

  type( mass_molecules )                           :: mol_mass

contains

  subroutine Assign_mass( this )
    use mod_error_handling

    implicit none

    class( mass_molecules ), intent(inout) :: this
    integer                                :: i
    integer                                :: ierr
    type(error)                            :: err

    this % total_mass      = 0.0_dp
    this % mass_components = 0.0_dp
    
    if ( allocated ( this % mass_atoms ) ) deallocate ( this % mass_atoms )
    allocate( this % mass_atoms( mol % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
    
    do i = 1, mol % num_atoms
 
      if ( mol % atoms(i) % label == 'H' ) then

        this % mass_atoms(i) % mass = mass_H

    else if ( mol % atoms(i) % label == 'He' ) then
   
      this % mass_atoms(i) % mass = mass_He
     
    else if ( mol % atoms(i) % label == 'Be' ) then
   
      this % mass_atoms(i) % mass = mass_Be
     
    else if ( mol % atoms(i) % label == 'B' ) then
   
      this % mass_atoms(i) % mass = mass_B
     
    else if ( mol % atoms(i) % label == 'C' ) then
   
      this % mass_atoms(i) % mass = mass_C
     
    else if ( mol % atoms(i) % label == 'N' ) then
   
      this % mass_atoms(i) % mass = mass_N
     
    else if ( mol % atoms(i) % label == 'O' ) then
   
      this % mass_atoms(i) % mass = mass_O
     
    else if ( mol % atoms(i) % label == 'F' ) then
   
      this % mass_atoms(i) % mass = mass_F
     
    else if ( mol % atoms(i) % label == 'Ne' ) then
   
      this % mass_atoms(i) % mass = mass_Ne
     
    else if ( mol % atoms(i) % label == 'Na' ) then
   
      this % mass_atoms(i) % mass = mass_Na
    
    else if ( mol % atoms(i) % label == 'Mg' ) then
   
      this % mass_atoms(i) % mass = mass_Mg
    
    else if ( mol % atoms(i) % label == 'Al' ) then
   
      this % mass_atoms(i) % mass = mass_Al
    
    else if ( mol % atoms(i) % label == 'Si' ) then
   
      this % mass_atoms(i) % mass = mass_Si
    
    else if ( mol % atoms(i) % label == 'P' ) then
   
      this % mass_atoms(i) % mass = mass_P
    
    else if ( mol % atoms(i) % label == 'S' ) then
   
      this % mass_atoms(i) % mass = mass_S
    
    else if ( mol % atoms(i) % label == 'Cl' ) then
  
      this % mass_atoms(i) % mass = mass_Cl
    
    else if ( mol % atoms(i) % label == 'Ar' ) then
   
      this % mass_atoms(i) % mass = mass_Ar
    
    else if ( mol % atoms(i) % label == 'K' ) then
   
      this % mass_atoms(i) % mass = mass_K
    
    else if ( mol % atoms(i) % label == 'Ca' ) then
   
      this % mass_atoms(i) % mass = mass_Ca
    
    else if ( mol % atoms(i) % label == 'Sc' ) then
   
      this % mass_atoms(i) % mass = mass_Sc
    
    else if ( mol % atoms(i) % label == 'Ti' ) then
   
      this % mass_atoms(i) % mass = mass_Ti
    
    else if ( mol % atoms(i) % label == 'V' ) then
   
      this % mass_atoms(i) % mass = mass_V
    
    else if ( mol % atoms(i) % label == 'Cr' ) then
   
      this % mass_atoms(i) % mass = mass_Cr
    
    else if ( mol % atoms(i) % label == 'Mn' ) then
   
      this % mass_atoms(i) % mass = mass_Mn
    
    else if ( mol % atoms(i) % label == 'Fe' ) then
   
      this % mass_atoms(i) % mass = mass_Fe
    
    else if ( mol % atoms(i) % label == 'Co' ) then
   
      this % mass_atoms(i) % mass = mass_Co
    
    else if ( mol % atoms(i) % label == 'Ni' ) then
   
      this % mass_atoms(i) % mass = mass_Ni
    
    else if ( mol % atoms(i) % label == 'Cu' ) then
   
      this % mass_atoms(i) % mass = mass_Cu
    
    else if ( mol % atoms(i) % label == 'Zn' ) then
   
      this % mass_atoms(i) % mass = mass_Zn
    
    else if ( mol % atoms(i) % label == 'Ga' ) then
   
      this % mass_atoms(i) % mass = mass_Ga
    
    else if ( mol % atoms(i) % label == 'Ge' ) then
   
      this % mass_atoms(i) % mass = mass_Ge
    
    else if ( mol % atoms(i) % label == 'As' ) then
   
      this % mass_atoms(i) % mass = mass_As
    
    else if ( mol % atoms(i) % label == 'Se' ) then
   
      this % mass_atoms(i) % mass = mass_Se
    
    else if ( mol % atoms(i) % label == 'Br' ) then
   
      this % mass_atoms(i) % mass = mass_Br
    
    else if ( mol % atoms(i) % label == 'Kr' ) then
   
      this % mass_atoms(i) % mass = mass_Kr
    
    else if ( mol % atoms(i) % label == 'Rb' ) then
   
      this % mass_atoms(i) % mass = mass_Rb
    
    else if ( mol % atoms(i) % label == 'Sr' ) then
   
      this % mass_atoms(i) % mass = mass_Sr
    
    else if ( mol % atoms(i) % label == 'Y' ) then
   
      this % mass_atoms(i) % mass = mass_Y
    
    else if ( mol % atoms(i) % label == 'Zr' ) then
   
      this % mass_atoms(i) % mass = mass_Zr
    
    else if ( mol % atoms(i) % label == 'Nb' ) then
   
      this % mass_atoms(i) % mass = mass_Nb
    
    else if ( mol % atoms(i) % label == 'Mo' ) then
   
      this % mass_atoms(i) % mass = mass_Mo

    else if ( mol % atoms(i) % label == 'Tc' ) then
   
      this % mass_atoms(i) % mass = mass_Tc

    else if ( mol % atoms(i) % label == 'Ru' ) then
   
      this % mass_atoms(i) % mass = mass_Ru
     
    else if ( mol % atoms(i) % label == 'Rh' ) then
   
      this % mass_atoms(i) % mass = mass_Rh
     
    else if ( mol % atoms(i) % label == 'Pd' ) then
   
      this % mass_atoms(i) % mass = mass_Pd
     
    else if ( mol % atoms(i) % label == 'Ag' ) then
   
      this % mass_atoms(i) % mass = mass_Ag
     
    else if ( mol % atoms(i) % label == 'Cd' ) then
   
      this % mass_atoms(i) % mass = mass_Cd
     
    else if ( mol % atoms(i) % label == 'In' ) then
   
      this % mass_atoms(i) % mass = mass_In
     
    else if ( mol % atoms(i) % label == 'Sn' ) then
   
      this % mass_atoms(i) % mass = mass_Sn
     
    else if ( mol % atoms(i) % label == 'Sb' ) then
   
      this % mass_atoms(i) % mass = mass_Sb
     
    else if ( mol % atoms(i) % label == 'Te' ) then
   
      this % mass_atoms(i) % mass = mass_Te
     
    else if ( mol % atoms(i) % label == 'I' ) then
   
      this % mass_atoms(i) % mass = mass_I
     
    else if ( mol % atoms(i) % label == 'Xe' ) then
   
      this % mass_atoms(i) % mass = mass_Xe
     
    else if ( mol % atoms(i) % label == 'Cs' ) then
   
      this % mass_atoms(i) % mass = mass_Cs

    else if ( mol % atoms(i) % label == 'Ba' ) then
   
      this % mass_atoms(i) % mass = mass_Ba
     
    else if ( mol % atoms(i) % label == 'Hf' ) then
   
      this % mass_atoms(i) % mass = mass_Hf
     
    else if ( mol % atoms(i) % label == 'Ta' ) then
   
      this % mass_atoms(i) % mass = mass_Ta
     
    else if ( mol % atoms(i) % label == 'W' ) then
   
      this % mass_atoms(i) % mass = mass_W
     
    else if ( mol % atoms(i) % label == 'Re' ) then
   
      this % mass_atoms(i) % mass = mass_Re
     
    else if ( mol % atoms(i) % label == 'Os' ) then
   
      this % mass_atoms(i) % mass = mass_Os
     
    else if ( mol % atoms(i) % label == 'Ir' ) then
   
      this % mass_atoms(i) % mass = mass_Ir
     
    else if ( mol % atoms(i) % label == 'Pt' ) then
   
      this % mass_atoms(i) % mass = mass_Pt
     
    else if ( mol % atoms(i) % label == 'Au' ) then
   
      this % mass_atoms(i) % mass = mass_Au
     
    else if ( mol % atoms(i) % label == 'Hg' ) then
   
      this % mass_atoms(i) % mass = mass_Hg
     
    else if ( mol % atoms(i) % label == 'Tl' ) then
   
      this % mass_atoms(i) % mass = mass_Tl
     
    else if ( mol % atoms(i) % label == 'Pb' ) then
   
      this % mass_atoms(i) % mass = mass_Pb

    else if ( mol % atoms(i) % label == 'Bi' ) then
   
      this % mass_atoms(i) % mass = mass_Bi
     
    else

      write(*,*) "missing atomic mass for element ", this % mass_atoms(i) % label
      write(*,*) "enter mass value"
      read(*,*) this % mass_atoms(i) % mass
  
    endif

    this % total_mass = this % total_mass + this % mass_atoms(i) % mass

    this % mass_components(1) = this % mass_components(1) + &
      & this % mass_atoms(i) % mass * mol % atoms(i) % xyz(1)

    this % mass_components(2) = this % mass_components(2) + &
      & this % mass_atoms(i) % mass * mol % atoms(i) % xyz(2)
   
    this % mass_components(3) = this % mass_components(3) + &
      & this % mass_atoms(i) % mass * mol % atoms(i) % xyz(3)

  enddo

  end subroutine Assign_mass

  subroutine Calc_com( this, com_filename, option_center )

    implicit none

    class( mass_molecules ), intent(inout) :: this
    character( len = * ), intent(in)       :: com_filename, option_center
    integer                                :: i
    integer                                :: file_unit   = 11        
    character( len = 15 )                  :: file_format = "formatted"
    character( len = 15 )                  :: file_access = "sequential"
    
    this % com(:) = this % mass_components(:) / this % total_mass

    open( unit = file_unit, file = trim(com_filename), status = 'unknown', &
              form = trim(file_format), access = trim(file_access) ) 

    write(file_unit,'(i9)') mol % num_atoms + 1
    write(file_unit,*)

    if ( option_center == "TRUE" ) then

      do i = 1, mol % num_atoms
 
        write(file_unit,'(2x,a2,3(5x,f12.6))') mol % atoms(i) % label, &
          mol % atoms(i) % xyz(:) - this % com(:)
   
      enddo
      
      write(file_unit,'(2x,a2,3(5x,f12.6))') 'X', this % com(:) - this % com(:)

    else if ( option_center == "FALSE" ) then

      do i = 1, mol % num_atoms
 
        write(file_unit,'(2x,a2,3(5x,f12.6))') mol % atoms(i) % label, &
          mol % atoms(i) % xyz(:) 
   
      enddo
      
      write(file_unit,'(2x,a2,3(5x,f12.6))') 'X', this % com(:)

    endif

    close(file_unit)

  end subroutine Calc_com

end module mod_calc_com
