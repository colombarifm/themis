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
!> @file   com.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  This program calculates the center of mass of a given molecular structure
!> @date - Jan, 2019                                                           
!> - independent module created                                                
!> @note 
!> - reads "filein.xyz" and writes "fileout.xyz" with an extra XX site corresponding to the center 
!>   of mass. XX coordinates can ba placed at the origin.
!> - usage: ./com <filein.xyz> <fileout.xyz> center/nocenter
!---------------------------------------------------------------------------------------------------

program com

  implicit none

  integer, parameter                :: dp = selected_real_kind(14)
  
  !! masses from http://www.chemicalelements.com/show/mass.html

  real(kind=dp), parameter          :: mass_H  =   1.007940_dp
  real(kind=dp), parameter          :: mass_He =   4.002602_dp
  real(kind=dp), parameter          :: mass_Li =   6.941000_dp
  real(kind=dp), parameter          :: mass_Be =   9.012182_dp
  real(kind=dp), parameter          :: mass_B  =  10.811000_dp
  real(kind=dp), parameter          :: mass_C  =  12.010700_dp
  real(kind=dp), parameter          :: mass_N  =  14.006740_dp
  real(kind=dp), parameter          :: mass_O  =  15.999400_dp
  real(kind=dp), parameter          :: mass_F  =  18.998400_dp
  real(kind=dp), parameter          :: mass_Ne =  20.179700_dp
  real(kind=dp), parameter          :: mass_Na =  22.989770_dp
  real(kind=dp), parameter          :: mass_Mg =  24.305000_dp
  real(kind=dp), parameter          :: mass_Al =  26.981539_dp
  real(kind=dp), parameter          :: mass_Si =  28.085500_dp
  real(kind=dp), parameter          :: mass_P  =  30.973760_dp
  real(kind=dp), parameter          :: mass_S  =  32.066600_dp
  real(kind=dp), parameter          :: mass_Cl =  35.452700_dp
  real(kind=dp), parameter          :: mass_Ar =  39.948000_dp
  real(kind=dp), parameter          :: mass_K  =  39.098300_dp
  real(kind=dp), parameter          :: mass_Ca =  40.078000_dp
  real(kind=dp), parameter          :: mass_Sc =  44.955910_dp
  real(kind=dp), parameter          :: mass_Ti =  47.867000_dp
  real(kind=dp), parameter          :: mass_V  =  50.941500_dp
  real(kind=dp), parameter          :: mass_Cr =  51.996100_dp
  real(kind=dp), parameter          :: mass_Mn =  54.938049_dp
  real(kind=dp), parameter          :: mass_Fe =  55.845000_dp
  real(kind=dp), parameter          :: mass_Co =  58.933200_dp
  real(kind=dp), parameter          :: mass_Ni =  58.693400_dp
  real(kind=dp), parameter          :: mass_Cu =  63.546000_dp
  real(kind=dp), parameter          :: mass_Zn =  65.390000_dp
  real(kind=dp), parameter          :: mass_Ga =  69.723000_dp
  real(kind=dp), parameter          :: mass_Ge =  72.610000_dp
  real(kind=dp), parameter          :: mass_As =  74.921600_dp
  real(kind=dp), parameter          :: mass_Se =  78.960000_dp
  real(kind=dp), parameter          :: mass_Br =  79.904000_dp
  real(kind=dp), parameter          :: mass_Kr =  83.800000_dp
  real(kind=dp), parameter          :: mass_Rb =  85.467800_dp
  real(kind=dp), parameter          :: mass_Sr =  87.620000_dp
  real(kind=dp), parameter          :: mass_Y  =  88.905850_dp
  real(kind=dp), parameter          :: mass_Zr =  91.224000_dp
  real(kind=dp), parameter          :: mass_Nb =  92.906380_dp
  real(kind=dp), parameter          :: mass_Mo =  95.950000_dp
  real(kind=dp), parameter          :: mass_Tc =  98.000000_dp
  real(kind=dp), parameter          :: mass_Ru = 101.070000_dp
  real(kind=dp), parameter          :: mass_Rh = 102.905500_dp
  real(kind=dp), parameter          :: mass_Pd = 106.420000_dp
  real(kind=dp), parameter          :: mass_Ag = 107.868200_dp
  real(kind=dp), parameter          :: mass_Cd = 112.411000_dp
  real(kind=dp), parameter          :: mass_In = 114.818000_dp
  real(kind=dp), parameter          :: mass_Sn = 118.710000_dp
  real(kind=dp), parameter          :: mass_Sb = 121.760000_dp
  real(kind=dp), parameter          :: mass_Te = 127.600000_dp
  real(kind=dp), parameter          :: mass_I  = 126.904470_dp
  real(kind=dp), parameter          :: mass_Xe = 131.290000_dp
  real(kind=dp), parameter          :: mass_Cs = 132.905450_dp
  real(kind=dp), parameter          :: mass_Ba = 137.327000_dp
  real(kind=dp), parameter          :: mass_Hf = 178.490000_dp
  real(kind=dp), parameter          :: mass_Ta = 180.947900_dp
  real(kind=dp), parameter          :: mass_W  = 183.840000_dp
  real(kind=dp), parameter          :: mass_Re = 186.207000_dp
  real(kind=dp), parameter          :: mass_Os = 190.230000_dp
  real(kind=dp), parameter          :: mass_Ir = 192.217000_dp
  real(kind=dp), parameter          :: mass_Pt = 195.084900_dp
  real(kind=dp), parameter          :: mass_Au = 196.966550_dp
  real(kind=dp), parameter          :: mass_Hg = 200.590000_dp
  real(kind=dp), parameter          :: mass_Tl = 204.383300_dp
  real(kind=dp), parameter          :: mass_Pb = 207.200000_dp
  real(kind=dp), parameter          :: mass_Bi = 208.980380_dp

  integer                     :: i

  type atom
    character(len=3)          :: label
    real(kind=dp)             :: mass_atom  
    real(kind=dp)             :: xyz(3)
  end type atom

  type molecule
    type(atom), allocatable, dimension(:)  :: atoms
    integer                                :: numat
    real(kind=dp)                          :: com(3)
    real(kind=dp)                          :: mass_molecule
    real(kind=dp)                          :: mass_components(3)
  end type molecule

  type(molecule)  :: mol_1
  
  character(30)                            :: namein, nameout, option

  call get_command_argument( 1, namein )
  call get_command_argument( 2, nameout )
  call get_command_argument( 3, option )

  option = trim(option)

  open(unit=7,file=namein,status='old')
  open(unit=8,file=nameout,status='replace')

  read(7,*) mol_1 % numat
  read(7,*)

  allocate( mol_1 % atoms ( mol_1 % numat ) )
  
  mol_1 % mass_molecule   = 0.0_dp
  mol_1 % mass_components = 0.0_dp

  do i = 1, mol_1 % numat
 
    read(7,*) mol_1 % atoms(i) % label, mol_1 % atoms(i) % xyz(:)

    if ( mol_1 % atoms(i) % label == 'H' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_H
     
    else if ( mol_1 % atoms(i) % label == 'He' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_He
     
    else if ( mol_1 % atoms(i) % label == 'Be' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Be
     
    else if ( mol_1 % atoms(i) % label == 'B' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_B
     
    else if ( mol_1 % atoms(i) % label == 'C' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_C
     
    else if ( mol_1 % atoms(i) % label == 'N' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_N
     
    else if ( mol_1 % atoms(i) % label == 'O' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_O
     
    else if ( mol_1 % atoms(i) % label == 'F' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_F
     
    else if ( mol_1 % atoms(i) % label == 'Ne' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ne
     
    else if ( mol_1 % atoms(i) % label == 'Na' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Na
    
    else if ( mol_1 % atoms(i) % label == 'Mg' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Mg
    
    else if ( mol_1 % atoms(i) % label == 'Al' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Al
    
    else if ( mol_1 % atoms(i) % label == 'Si' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Si
    
    else if ( mol_1 % atoms(i) % label == 'P' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_P
    
    else if ( mol_1 % atoms(i) % label == 'S' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_S
    
    else if ( mol_1 % atoms(i) % label == 'Cl' ) then
  
      mol_1 % atoms(i) % mass_atom = mass_Cl
    
    else if ( mol_1 % atoms(i) % label == 'Ar' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ar
    
    else if ( mol_1 % atoms(i) % label == 'K' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_K
    
    else if ( mol_1 % atoms(i) % label == 'Ca' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ca
    
    else if ( mol_1 % atoms(i) % label == 'Sc' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Sc
    
    else if ( mol_1 % atoms(i) % label == 'Ti' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ti
    
    else if ( mol_1 % atoms(i) % label == 'V' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_V
    
    else if ( mol_1 % atoms(i) % label == 'Cr' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Cr
    
    else if ( mol_1 % atoms(i) % label == 'Mn' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Mn
    
    else if ( mol_1 % atoms(i) % label == 'Fe' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Fe
    
    else if ( mol_1 % atoms(i) % label == 'Co' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Co
    
    else if ( mol_1 % atoms(i) % label == 'Ni' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ni
    
    else if ( mol_1 % atoms(i) % label == 'Cu' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Cu
    
    else if ( mol_1 % atoms(i) % label == 'Zn' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Zn
    
    else if ( mol_1 % atoms(i) % label == 'Ga' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ga
    
    else if ( mol_1 % atoms(i) % label == 'Ge' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ge
    
    else if ( mol_1 % atoms(i) % label == 'As' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_As
    
    else if ( mol_1 % atoms(i) % label == 'Se' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Se
    
    else if ( mol_1 % atoms(i) % label == 'Br' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Br
    
    else if ( mol_1 % atoms(i) % label == 'Kr' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Kr
    
    else if ( mol_1 % atoms(i) % label == 'Rb' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Rb
    
    else if ( mol_1 % atoms(i) % label == 'Sr' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Sr
    
    else if ( mol_1 % atoms(i) % label == 'Y' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Y
    
    else if ( mol_1 % atoms(i) % label == 'Zr' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Zr
    
    else if ( mol_1 % atoms(i) % label == 'Nb' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Nb
    
    else if ( mol_1 % atoms(i) % label == 'Mo' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Mo

    else if ( mol_1 % atoms(i) % label == 'Tc' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Tc

    else if ( mol_1 % atoms(i) % label == 'Ru' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ru
     
    else if ( mol_1 % atoms(i) % label == 'Rh' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Rh
     
    else if ( mol_1 % atoms(i) % label == 'Pd' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Pd
     
    else if ( mol_1 % atoms(i) % label == 'Ag' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ag
     
    else if ( mol_1 % atoms(i) % label == 'Cd' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Cd
     
    else if ( mol_1 % atoms(i) % label == 'In' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_In
     
    else if ( mol_1 % atoms(i) % label == 'Sn' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Sn
     
    else if ( mol_1 % atoms(i) % label == 'Sb' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Sb
     
    else if ( mol_1 % atoms(i) % label == 'Te' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Te
     
    else if ( mol_1 % atoms(i) % label == 'I' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_I
     
    else if ( mol_1 % atoms(i) % label == 'Xe' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Xe
     
    else if ( mol_1 % atoms(i) % label == 'Cs' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Cs

    else if ( mol_1 % atoms(i) % label == 'Ba' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ba
     
    else if ( mol_1 % atoms(i) % label == 'Hf' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Hf
     
    else if ( mol_1 % atoms(i) % label == 'Ta' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ta
     
    else if ( mol_1 % atoms(i) % label == 'W' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_W
     
    else if ( mol_1 % atoms(i) % label == 'Re' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Re
     
    else if ( mol_1 % atoms(i) % label == 'Os' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Os
     
    else if ( mol_1 % atoms(i) % label == 'Ir' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Ir
     
    else if ( mol_1 % atoms(i) % label == 'Pt' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Pt
     
    else if ( mol_1 % atoms(i) % label == 'Au' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Au
     
    else if ( mol_1 % atoms(i) % label == 'Hg' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Hg
     
    else if ( mol_1 % atoms(i) % label == 'Tl' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Tl
     
    else if ( mol_1 % atoms(i) % label == 'Pb' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Pb

    else if ( mol_1 % atoms(i) % label == 'Bi' ) then
   
      mol_1 % atoms(i) % mass_atom = mass_Bi
     
    else

      write(*,*) "missing atomic mass for element ", mol_1 % atoms(i) % label
      write(*,*) "enter mass value"
      read(*,*) mol_1 % atoms(i) % mass_atom
  
    endif

    mol_1 % mass_molecule = mol_1 % mass_molecule  + mol_1 % atoms(i) % mass_atom

    mol_1 % mass_components(1) = mol_1 % mass_components(1) + &
      & mol_1 % atoms(i) % mass_atom * mol_1 % atoms(i) % xyz(1)

    mol_1 % mass_components(2) = mol_1 % mass_components(2) + &
      & mol_1 % atoms(i) % mass_atom * mol_1 % atoms(i) % xyz(2)
   
    mol_1 % mass_components(3) = mol_1 % mass_components(3) + &
      & mol_1 % atoms(i) % mass_atom * mol_1 % atoms(i) % xyz(3)

  enddo

  mol_1 % com(:) = mol_1 % mass_components(:) / mol_1 % mass_molecule

  write(8,*) mol_1 % numat + 1
  write(8,*)

  write(*,*) mol_1 % com(:), mol_1 % mass_components(:)

  if ( option == "center" ) then

    do i = 1, mol_1 % numat
 
      write(8,'(2x,a2,3(5x,f12.6))') mol_1 % atoms(i) % label, mol_1 % atoms(i) % xyz(:) - mol_1 % com(:)
   
    enddo
      
    write(8,'(2x,a2,3(5x,f12.6))') 'X', mol_1 % com(:) - mol_1 % com(:)

  else if ( option == "nocenter" ) then

    do i = 1, mol_1 % numat
 
      write(8,'(2x,a2,3(5x,f12.6))') mol_1 % atoms(i) % label, mol_1 % atoms(i) % xyz(:) 
   
    enddo
      
    write(8,'(2x,a2,3(5x,f12.6))') 'X', mol_1 % com(:)

  else 

    write(*,*) "invalid option for centering."
    write(*,*) "usage is: ./com <filein> <fileout> center/nocenter."

  endif

  deallocate(mol_1 % atoms)

  close(7)
  close(8)

end program

