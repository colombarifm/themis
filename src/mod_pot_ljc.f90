!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Themis developers
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
!> @file   mod_pot_ljc.f90
!> @author Felippe M. Colombari
!> @brief  This module performs LJ + coulombic potential calculations
!> @date - Oct, 2017                                                           
!> - independent module created                                                
!> @date - Jan, 2018                                                           
!> - support added  
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_pot_ljc
  use iso_fortran_env, only: output_unit
  use mod_constants      , only: CCON, DP
  use mod_read_molecules , only: mol1, mol2, atom, molecule, dimer

  implicit none

  type, extends( atom ), private                     :: ljc_atom
    real( kind = DP )                                :: q, sig, eps
  end type

  type, extends( molecule ), private                 :: ljc_molecules
    type( ljc_atom ), allocatable, dimension(:)      :: ljc_atoms
  contains 
    procedure, pass                                :: Read_LJC_params
    procedure, pass                                :: Check_LJC_params
  end type ljc_molecules

  type( ljc_molecules )                              :: mol1_ljc, mol2_ljc
 
  type, extends( dimer ), public                     :: ljc_dimer
    real( kind = DP )                                :: pot_lj, pot_coul
    real( kind = DP ), allocatable, dimension(:,:)   :: eps_ij, sig_ij, q_ij
    logical, allocatable, dimension(:,:)             :: dummy
  contains
    procedure, pass, public                        :: Calc_cross  => Calc_LJC_cross
    procedure, pass, public                        :: Calc_energy => Calc_LJC_energy 
  end type ljc_dimer

  type( ljc_dimer )                                  :: ljc_dimers

contains

  !---------------------------------------------------------------------------
  !> @brief This routine reads all necessary parameters from files.
  !> @author Felippe M. Colombari
  !> @date - Jun, 2017
  !> - subroutine  created
  !> @note Files read are:\n
  !> - parameters1: contains labels, atomic charges, sigma (A) and epsilon (kJ/mol) for molecule 1
  !> - parameters2: contains labels, atomic charges, sigma (A) and epsilon (kJ/mol) for molecule 2
  !---------------------------------------------------------------------------	

  subroutine  Read_ljc_params( this, ljc_filename, numat )
    use mod_constants      , only: dashline
    use mod_inquire        , only: Inquire_file, Get_new_unit
    use mod_error_handling

    implicit none

    class( ljc_molecules ), intent(INOUT) :: this
    character( len = * ), intent(IN)     :: ljc_filename
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

    call Inquire_file( file_unit, ljc_filename, file_status, file_format, file_access )

    allocate( this % ljc_atoms( numat ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    read( file_unit, * )

    do i = 1, numat
        
      read( file_unit, *, iostat = ios ) dummy, this % ljc_atoms(i) % q, &
                                                this % ljc_atoms(i) % sig, &
                                                this % ljc_atoms(i) % eps

      if ( ios /= 0 ) then

        write(line_number,'(i10)') i+1
        call err%error('e',message="while reading file "//trim(adjustl(ljc_filename))//".")
        call err%error('e',check="line "//trim(adjustl(line_number))//".")
        write(output_unit,'( /, T3, A )') dashline

        stop

      endif

    enddo

    close( file_unit )
      
    return
  end subroutine  Read_ljc_params

  !------------------------------------------------------------------------------
  !> @brief This routine checks for inconsistent entries on parameter files.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine  Check_ljc_params( this, mol, numat, ljc_filename )
    use mod_constants, only: dashline
    use mod_error_handling

    implicit none

    class( ljc_molecules ), intent(inout) :: this
    integer, intent(IN)                   :: numat, mol
    character( len = * ), intent(IN)      :: ljc_filename
    type(error)                           :: err

    write(output_unit,'( /, T3, A )')                             dashline
    write(output_unit,'( T5, "Checking molecule ", i1, " ..." )') mol
    write(output_unit,'( T3, A )')                                dashline
    write(output_unit,'( /, T5, "Number of atoms", T91, i10   )') numat
    write(output_unit,'( /, T5, "Total charge", T91, f10.5, / )') sum(this % ljc_atoms % q)

    if ( ANY( this % ljc_atoms % sig < 0.0_DP ) ) then

      call err%error('e',message="while reading file "//trim(adjustl(ljc_filename))//".")
      call err%error('e',check="for sigma values < 0.")
      write(output_unit,'( /, T3, A )') dashline

      stop

    else if ( ANY( this % ljc_atoms % eps < 0.0_DP ) ) then

      call err%error('e',message="while reading file "//trim(adjustl(ljc_filename))//".")
      call err%error('e',check="for epsilon values < 0.")
      write(output_unit,'( /, T3, A )') dashline
        
      stop

    else

      write(output_unit,'( T5, "LJ parameters:", T99, "OK" )') 

    endif

    return 
  end subroutine Check_ljc_params

  !------------------------------------------------------------------------------
  !> @brief This routine calculates crossing values for LJ + coulomb potential before the loops.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine  Calc_ljc_cross( this )
    use mod_error_handling

    implicit none

    class( ljc_dimer ), intent(INOUT) :: this
    integer                           :: i, j
    integer                           :: ierr
    type(error)                       :: err

    allocate(       this % q_ij( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(     this % sig_ij( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(     this % eps_ij( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(     this % dummy( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    this % eps_ij     = 0.0_DP
    this % sig_ij     = 0.0_DP
    this % q_ij       = 0.0_DP
    
    this % dummy      = .false.

    do j = 1, mol2 % num_atoms
        
      do i = 1, mol1 % num_atoms
      
        this % eps_ij(i,j)     = 4.0_DP * dsqrt( ( mol1_ljc % ljc_atoms(i) % eps * mol2_ljc % ljc_atoms(j) % eps ) )
        this % sig_ij(i,j)     =                   mol1_ljc % ljc_atoms(i) % sig * mol2_ljc % ljc_atoms(j) % sig
        this % q_ij(i,j)       =            CCON * mol1_ljc % ljc_atoms(i) % q   * mol2_ljc % ljc_atoms(j) % q

        if ( ( mol1 % atoms(i) % symbol(1:1) == 'X') .or. ( mol2 % atoms(j) % symbol(1:1) == 'X' ) ) then

          this % dummy(i,j) = .true.

        endif

      enddo 

    enddo 

    return
  end subroutine Calc_ljc_cross

  !------------------------------------------------------------------------------
  !> @brief This routine calculates LJ + coulomb energy for each valid configuration.
  !> @author Felippe M. Colombari
  !> @date - Oct, 2017
  !> - modular subroutine  created
  !> @date - Aug, 2018
  !> - dummy sites (X*) are skipped
  !------------------------------------------------------------------------------
  subroutine  Calc_ljc_energy( this, r2, r1, t )
    use mod_input_read, only: rcut_sqr, atom_overlap, inter_energy

    implicit none

    class( ljc_dimer ), intent(INOUT) :: this
    integer, intent(IN)               :: r2, r1, t
    integer                           :: i, j
    real( kind = DP )                 :: lj_factor_cube, rij, rijrij

    this % pot_coul = 0.0_DP
    this % pot_lj   = 0.0_DP
    rijrij          = 0.0_DP
    rij             = 0.0_DP

     jlp: do j = 1, mol2 % num_atoms

       ilp: do i = 1, mol1 % num_atoms

        rijrij = sum( ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) * &
                      ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) )

        if ( ( rijrij <= rcut_sqr ) .and. ( this % dummy(i,j) .eqv. .false. ) ) then 

          atom_overlap( r2, r1, t ) = .true.

          this % pot_lj       = 1.0E10_DP

          this % pot_coul     = 0.0_DP

          exit jlp

        else if ( ( rijrij > rcut_sqr ) .and. ( this % dummy(i,j) .eqv. .false. ) ) then

          lj_factor_cube = ( this % sig_ij(i,j) / rijrij ) ** 3

          this % pot_lj = this % pot_lj + this % eps_ij(i,j) * ( lj_factor_cube * ( lj_factor_cube - 1 ) )

          rij = dsqrt(rijrij)
            
          this % pot_coul = this % pot_coul + this % q_ij(i,j) / rij

        else if ( ( rijrij <= rcut_sqr ) .and. ( this % dummy(i,j) .eqv. .true. ) ) then

          cycle

        else if ( ( rijrij > rcut_sqr ) .and. ( this % dummy(i,j) .eqv. .true. ) ) then

          cycle

        endif

      enddo ilp

    enddo jlp

    inter_energy( r2, r1, t ) = ( this % pot_coul + this % pot_lj ) 
  
    return
  end subroutine Calc_ljc_energy

end module mod_pot_ljc
