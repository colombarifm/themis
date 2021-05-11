!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2021 Themis developers
!
!   This file was written by Felippe M. Colombari and Asdrubal Lozada-Blanco.
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
!> @file   mod_pot_bhc.f90
!> @author Felippe M. Colombari
!> @brief  This module performs BH + COULOMBIC potential calculations.
!> @date - Oct, 2017                                                           
!> - independent module created                                                
!> @date - Jan, 2018                                                           
!> - support added  
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @note Buckingham potential used is the one by MATSUI (1991)
!------------------------------------------------------------------------------

module mod_pot_bhc
  use iso_fortran_env, only: output_unit
  use mod_constants      , only: CCON, DP
  use mod_read_molecules , only: mol1, mol2, atom, molecule, dimer

  implicit none

  type, extends( atom ), public                    :: bhc_atom
    real( kind = DP )                              :: q, A_bh, B_bh, C_bh
  end type

  type, extends( molecule ), public                :: bhc_molecule
    type( bhc_atom ), allocatable, dimension(:)    :: bhc_atoms
  contains  
    procedure, pass                              :: Read_BHC_params
    procedure, pass                              :: Check_BHC_params
  end type bhc_molecule

  type( bhc_molecule )                             :: mol1_bhc, mol2_bhc
 
  type, extends( dimer ), public                   :: bhc_dimer
    real( kind = DP )                              :: pot_bh, pot_bh_rep, pot_bh_disp, pot_coul
    real( kind = DP ), allocatable, dimension(:,:) :: A_bh_ij, B_bh_ij, C_bh_ij, q_ij 
  contains
    procedure, pass, public                      :: Calc_cross  => Calc_BHC_cross
    procedure, pass, public                      :: Calc_energy => Calc_BHC_energy 
  end type bhc_dimer

  type( bhc_dimer )                                :: bhc_dimers

contains

  !---------------------------------------------------------------------------
  !> @brief This routine reads all necessary parameters from files.
  !> @author Felippe M. Colombari
  !> @date - Jun, 2017
  !> - subroutine  created
  !> @note Files read are:\n
  !> - parameters1: contains labels, atomic charges, A, B and C values for molecule 1
  !> - parameters2: contains labels, atomic charges, A, B and C values for molecule 2
  !---------------------------------------------------------------------------	
  subroutine Read_BHC_Params( this, bhc_filename, numat )
    use mod_constants      , only: dashline
    use mod_inquire        , only: Inquire_file, Get_new_unit
    use mod_error_handling

    implicit none

    class( bhc_molecule ), intent(INOUT) :: this
    character( len = * ), intent(IN)     :: bhc_filename
    integer, intent(IN)                  :: numat
    integer                              :: i
    integer                              :: file_unit           
    integer                              :: ios         = 0
    character( len = * ), parameter      :: file_status = "old"
    character( len = * ), parameter      :: file_format = "formatted"
    character( len = * ), parameter      :: file_access = "sequential"
    character( len = 2 )                 :: dummy
    character( len = 10 )                :: line_number
    integer                              :: ierr
    type(error)                          :: err

    file_unit = Get_new_unit(10)

    call Inquire_file( file_unit, bhc_filename, file_status, file_format, file_access )

    allocate( this % bhc_atoms( numat ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
  
    read( file_unit, * )

    do i = 1, numat

      read( file_unit, *, iostat = ios ) dummy, this % bhc_atoms(i) % q, &
                                                this % bhc_atoms(i) % A_bh, &
                                                this % bhc_atoms(i) % B_bh, &
                                                this % bhc_atoms(i) % C_bh

      if ( ios /= 0 ) then

        write(line_number,'(i10)') i+1
        call err%error('e',message="while reading file "//trim(adjustl(bhc_filename))//".")
        call err%error('e',check="line "//trim(adjustl(line_number))//".")
        write(output_unit,'( /, T3, A )') dashline

        stop

      endif

    enddo

    close( file_unit )
      
    return
  end subroutine Read_BHC_Params
      
  !------------------------------------------------------------------------------
  !> @brief This routine checks for inconsistent entries on parameter files.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine Check_bhc_params( this, mol, numat, bhc_filename )
    use mod_constants, only: dashline
    use mod_error_handling

    implicit none

    class( bhc_molecule ), intent(inout) :: this
    integer, intent(IN)                  :: numat, mol
    character( len = * ), intent(IN)     :: bhc_filename
    type(error)                          :: err

    write(output_unit,'( /, T3, A )')                             dashline
    write(output_unit,'( T5, "Checking molecule ", i1, " ..." )') mol 
    write(output_unit,'( T3, A )')                                dashline
    write(output_unit,'( /, T5, "Number of atoms", T65, i10   )') numat
    write(output_unit,'( /, T5, "Total charge", T65, f10.5, / )') sum(this % bhc_atoms % q)

    if ( ANY( this % bhc_atoms % A_bh < 0.0_DP ) ) then

      call err%error('e',message="while reading file "//trim(adjustl(bhc_filename))//".")
      call err%error('e',check="for A values < 0.")
      write(output_unit,'( /, T3, A )') dashline

      stop

    else if ( ANY( this % bhc_atoms % B_bh < 0.0_DP ) ) then

      call err%error('e',message="while reading file "//trim(adjustl(bhc_filename))//".")
      call err%error('e',check="for B values < 0.")
      write(output_unit,'( /, T3, A )') dashline

      stop

    else if ( ANY( this % bhc_atoms % C_bh < 0.0_DP ) ) then

      call err%error('e',message="while reading file "//trim(adjustl(bhc_filename))//".")
      call err%error('e',check="for C values < 0.")
      write(output_unit,'( /, T3, A )') dashline

      stop

    else

      write(output_unit,'( T5, "BH parameters:", T73, "OK" )') 

    endif

    return 
  end subroutine Check_BHC_params

  !------------------------------------------------------------------------------
  !> @brief This routine calculates crossing values for BH + coulomb potential before the loops.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine Calc_bhc_cross( this )
    use mod_error_handling

    implicit none

    class( bhc_dimer ), intent(INOUT) :: this
    integer                           :: i, j
    integer                           :: ierr
    type(error)                       :: err

    allocate(     this % q_ij( mol1 % num_atoms, mol2 % num_atoms ),stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate(  this % A_bh_ij( mol1 % num_atoms, mol2 % num_atoms ),stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate(  this % B_bh_ij( mol1 % num_atoms, mol2 % num_atoms ),stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate(  this % C_bh_ij( mol1 % num_atoms, mol2 % num_atoms ),stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    this % A_bh_ij = 0.0_DP
    this % B_bh_ij = 0.0_DP
    this % C_bh_ij = 0.0_DP
    this % q_ij    = 0.0_DP

    do j = 1, mol2 % num_atoms
        
      do i = 1, mol1 % num_atoms
     
        this % A_bh_ij(i,j) =        mol1_bhc % bhc_atoms(i) % A_bh + mol2_bhc % bhc_atoms(j) % A_bh
        this % B_bh_ij(i,j) =        mol1_bhc % bhc_atoms(i) % B_bh + mol2_bhc % bhc_atoms(j) % B_bh
        this % C_bh_ij(i,j) =        mol1_bhc % bhc_atoms(i) % C_bh * mol2_bhc % bhc_atoms(j) % C_bh
        this % q_ij(i,j)    = CCON * mol1_bhc % bhc_atoms(i) % q    * mol2_bhc % bhc_atoms(j) % q

      enddo 

    enddo 

    return
  end subroutine Calc_bhc_cross

  !---------------------------------------------------------------------------
  !> @brief This routine calculates the interaction energy 
  !> @author Felippe M. Colombari
  !> @date - Oct, 2017
  !> - modular subroutine  created
  !> @date - Aug, 2018
  !> - dummy sites (X*) are skipped
  !---------------------------------------------------------------------------	
  subroutine Calc_bhC_energy( this, r2, r1, t )
    use mod_input_read, only: rcut_sqr, atom_overlap, inter_energy

    implicit none

    class( bhc_dimer ), intent(INOUT) :: this
    integer, intent(IN)               :: r2, r1, t
    integer                           :: i, j
    real( kind = DP )                 :: rij, rijrij

    this % pot_coul    = 0.0_DP
    this % pot_bh      = 0.0_DP
    this % pot_bh_rep  = 0.0_DP
    this % pot_bh_disp = 0.0_DP
    rijrij             = 0.0_DP
    rij                = 0.0_DP

    jlp: do j = 1, mol2 % num_atoms

      ilp: do i = 1, mol1 % num_atoms

        rijrij = sum( ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) * &
                      ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) )

        if ( ( rijrij < rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) /= 'X' ) ) then

          atom_overlap( r2, r1, t ) = .true.

          this % pot_bh       = 1.0E10_DP

          this % pot_coul     = 0.0_DP

          exit jlp

        else if ( ( rijrij < rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) == 'X' ) ) then

          continue

        else if ( ( rijrij > rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) == 'X' ) ) then

          continue

        else if ( ( rijrij > rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) /= 'X' ) ) then

          rij = dsqrt(rijrij)
            
          this % pot_bh_disp = -this % C_bh_ij(i,j) / ( rij ** 6 )
            
          this % pot_bh_rep  = 4.184_DP * this % B_bh_ij(i,j) * dexp( ( this % A_bh_ij(i,j) - rij ) / this % B_bh_ij(i,j) )
            
          this % pot_coul    = this % pot_coul + this % q_ij(i,j) / rij 
            
          this % pot_bh      = this % pot_bh + this % pot_bh_disp + this % pot_bh_rep

        endif 

      enddo ilp

    enddo jlp

    inter_energy( r2, r1, t ) = ( this % pot_bh + this % pot_coul ) 
  
    return
  end subroutine Calc_bhc_energy

end module mod_pot_bhc
