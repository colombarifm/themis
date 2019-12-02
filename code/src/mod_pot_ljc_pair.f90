!-----------------------------------------------------------------------------
! THEMIS: 
! Copyright (C) 2019 Felippe M. Colombari
!------------------------------------------------------------------------------
!> @brief This module performs LJ + coulombic potential calculations
!> @author Felippe M. Colombari
!> - Laboratório de Química Teórica, LQT -- UFSCar
!> @date - Oct, 2017 
!> - module created
!> @note added error_handling
!  Asdrubal Lozada-Blanco
!> @date - Nov 2019
!------------------------------------------------------------------------------

module mod_pot_ljc_pair
  use mod_constantS      , only: CCON, DP
  use mod_read_molecules , only: mol1, mol2, atom, molecule, dimer

  implicit none

  type, extends( atom ), private                     :: ljc_pair_atom
    real( kind = DP )                                :: q, sig, eps
  end type

  type, extends( molecule ), private                 :: ljc_pair_molecule
    type( ljc_pair_atom ), allocatable, dimension(:)  :: ljc_pair_atoms
    contains 
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
    !> @brief This routine reads atom names, charges, sigma and epsilon values from parameter files.
    !> @author Felippe M. Colombari

    subroutine  Read_ljc_pair_Params( this, ljc_pair_filename )
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
          write(*,'( /, T3, A )') dashline

          stop

        endif

      enddo

      close( file_unit )
      
      return
    end subroutine  Read_ljc_pair_params

    !------------------------------------------------------------------------------
    !> @brief This routine calculates crossing values for LJ + coulomb potential before the loops.
    !> @author Felippe M. Colombari

    subroutine  Calc_ljc_pair_cross( this )
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
    end subroutine  Calc_ljc_pair_cross

    !------------------------------------------------------------------------------
    !> @brief This routine calculates LJ + coulomb energy for each valid configuration.
    !> @author Felippe M. Colombari
    
    subroutine  Calc_ljc_pair_energy( this, g, r, t )
      use MOD_INPUT_READ, only: rcut_sqr, cutoff_sqr, atom_overlap, inter_energy, scale_factor

      implicit none

      class( ljc_pair_dimer ), intent(INOUT) :: this
      integer, intent(IN)               :: g, r, t
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

            atom_overlap( g, r ,t ) = .true.

            this % pot_lj       = 1.0E10_DP

            this % pot_coul     = 0.0_DP

            exit jlp

          else if ( ( rijrij < rcut_sqr ) .and. ( this % c6_ij(i,j) < 1.0d-6 ) ) then

            continue

          else if ( ( rijrij > rcut_sqr ) .and. ( this % c6_ij(i,j) < 1.0d-6 ) ) then

            continue

          else if ( ( rijrij > rcut_sqr ) .and. ( rijrij < cutoff_sqr ) .and. ( this % c6_ij(i,j) > 1.0d-6 ) ) then

            this % pot_lj = this % pot_lj + this % c12_ij(i,j) / rijrij**6 - this % c6_ij(i,j) / rijrij**3
!            this % pot_lj = 0.0_DP 

          endif

!          if ( dabs(this % q_ij(i,j)) < 1.0d-6 ) then

!            continue

!          else if ( ( dabs(this % q_ij(i,j)) > 1.0d-6 ) .and. ( rijrij < cutoff_sqr ) ) then

!            rij = dsqrt(rijrij)
            
!            this % pot_coul = this % pot_coul + this % q_ij(i,j) / rij

!          endi

        enddo ilp

      enddo jlp

      inter_energy( g, r, t ) = this % pot_lj / scale_factor 
  
      return
    end subroutine  CALC_ljc_pair_ENERGY

end module MOD_POT_ljc_pair

