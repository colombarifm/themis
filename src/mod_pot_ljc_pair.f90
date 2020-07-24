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
!> @file   mod_pot_ljc_pair.f90
!> @author Felippe M. Colombari
!>         Laboratory of Theoretical Chemistry - LQT
!>         Federal University of São Carlos
!>         <http://www.lqt.dq.ufscar.br>
!> @email  colombarifm@hotmail.com
!> @brief  This module performs martini-like LJ + coulombic potential calculations
!> @date - Sep, 2019
!> - independent module created                                                
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_pot_ljc_pair
  use iso_fortran_env, only: output_unit
  use mod_constantS      , only: CCON, DP
  use mod_read_molecules , only: mol1, mol2, atom, molecule, dimer

  implicit none

  type, extends( atom ), private                     :: ljc_pair_atom
    real( kind = DP )                                :: q, sig, eps
  end type

  type, extends( molecule ), private                 :: ljc_pair_molecule
    type( ljc_pair_atom ), allocatable, dimension(:)  :: ljc_pair_atoms
  end type ljc_pair_molecule

  type( ljc_pair_molecule )                           :: mol1_ljc_pair, mol2_ljc_pair
 
  type, extends( dimer ), public                     :: ljc_pair_dimer
    integer                                          :: ntypes1, ntypes2
    real( kind = DP )                                :: pot_lj, pot_coul
    real( kind = DP ), allocatable, dimension(:)     :: c6, c12, qiqj
    real( kind = DP ), allocatable, dimension(:,:)   :: c6_ij, c12_ij, q_ij
    character( len = 2), allocatable, dimension(:,:) :: pair
  contains
    procedure, pass, public                        :: Read_ljc_pair_params
    procedure, pass, public                        :: Calc_cross  => Calc_ljc_pair_cross
    procedure, pass, public                        :: Calc_energy => Calc_ljc_pair_energy 
  end type ljc_pair_dimer

  type( ljc_pair_dimer )                              :: ljc_pair_dimers

contains

  !------------------------------------------------------------------------------
  !> @brief This routine reads atom i, atom j, C6, C12 and qi*qj values from parameter file.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine Read_ljc_pair_Params( this, ljc_pair_filename )
    use mod_constants, only: dashline
    use mod_inquire, only: Inquire_file
    use mod_error_handling

    implicit none

    class( ljc_pair_dimer ), intent(INOUT) :: this
    character( len = * ), intent(IN)     :: ljc_pair_filename
    integer                              :: i
    integer                              :: ios         = 0
    integer, parameter                   :: file_unit   = 12            
    character( len = * ), parameter      :: file_format = "formatted"
    character( len = * ), parameter      :: file_access = "sequential"
    character( len = 10 )                :: line_number
    integer                              :: ierr
    type(error)                          :: err

    call Inquire_file( file_unit, ljc_pair_filename, file_format, file_access )

    read( file_unit, * )
    read( file_unit, * ) this % ntypes1, this % ntypes2

    allocate(   this % c6(    this % ntypes1 * this % ntypes2 ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(  this % c12(    this % ntypes1 * this % ntypes2 ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate( this % qiqj(    this % ntypes1 * this % ntypes2 ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate( this % pair( 2, this % ntypes1 * this % ntypes2 ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    do i = 1, this % ntypes1 * this % ntypes2

      read( file_unit, *, iostat = ios ) this % pair(1,i), this % pair(2,i), this % c6(i), this % c12(i), this % qiqj(i)
        
      if ( ios /= 0 ) then

        write(line_number,'(i10)') i+1
        call err%error('e',message="while reading file "//trim(adjustl(ljc_pair_filename))//".")
        call err%error('e',check="line "//trim(adjustl(line_number))//".")
        write(output_unit,'( /, T3, A )') dashline

        stop

      endif

    enddo

    close( file_unit )
      
    return
  end subroutine Read_ljc_pair_params

  !------------------------------------------------------------------------------
  !> @brief This routine calculates crossing values for LJ + coulomb potential before the loops.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine Calc_ljc_pair_cross( this )
    use mod_error_handling

    implicit none

    class( ljc_pair_dimer ), intent( INOUT ) :: this
    integer                                 :: i, j, k
    integer                                 :: ierr
    type(error)                             :: err

    allocate( this % c12_ij( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate(  this % c6_ij( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    allocate(   this % q_ij( mol1 % num_atoms, mol2 % num_atoms ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
    do k = 1, ljc_pair_dimers % ntypes1 * ljc_pair_dimers % ntypes2

      do j = 1, mol2 % num_atoms

        do i = 1, mol1 % num_atoms

          if ( ( mol2 % atoms(j) % symbol == ljc_pair_dimers % pair(2,k) ) .and. &
               ( mol1 % atoms(i) % symbol == ljc_pair_dimers % pair(1,k) ) ) then

            this % c12_ij(i,j) = ljc_pair_dimers % c12(k)
            this % c6_ij(i,j)  = ljc_pair_dimers % c6(k)
            this % q_ij(i,j)   = CCON * ljc_pair_dimers % qiqj(k)

          endif

        enddo

      enddo

    enddo

    return
  end subroutine Calc_ljc_pair_cross

  !------------------------------------------------------------------------------
  !> @brief This routine calculates LJ + coulomb energy for each valid configuration.
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine Calc_ljc_pair_energy( this, r2, r1, t )
    use mod_input_read, only: rcut_sqr, cutoff_sqr, atom_overlap, inter_energy, scale_factor

    implicit none

    class( ljc_pair_dimer ), intent(INOUT) :: this
    integer, intent(IN)               :: r2, r1, t
    integer                           :: i, j
    real( kind = DP )                 :: rij, rijrij

    this % pot_coul = 0.0_DP
    this % pot_lj   = 0.0_DP
    rijrij          = 0.0_DP
    rij             = 0.0_DP

    jlp: do j = 1, mol2 % num_atoms

      ilp: do i = 1, mol1 % num_atoms

        rijrij = sum( ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) * &
                      ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) )

        if ( ( rijrij < rcut_sqr ) .and. ( this % c6_ij(i,j) > 1.0d-6 ) ) then

          atom_overlap( r2, r1 ,t ) = .true.

          this % pot_lj       = 1.0E10_DP

          this % pot_coul     = 0.0_DP

          exit jlp

        else if ( ( rijrij < rcut_sqr ) .and. ( this % c6_ij(i,j) < 1.0d-6 ) ) then

          continue

        else if ( ( rijrij > rcut_sqr ) .and. ( this % c6_ij(i,j) < 1.0d-6 ) ) then

          continue

        else if ( ( rijrij > rcut_sqr ) .and. ( rijrij < cutoff_sqr ) .and. ( this % c6_ij(i,j) > 1.0d-6 ) ) then

          this % pot_lj = this % pot_lj + this % c12_ij(i,j) / rijrij**6 - this % c6_ij(i,j) / rijrij**3

        endif

      enddo ilp

    enddo jlp

    inter_energy( r2, r1, t ) = this % pot_lj / scale_factor 
  
    return
  end subroutine calc_ljc_pair_energy

end module mod_pot_ljc_pair
