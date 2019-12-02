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

module MOD_POT_LJC
  use MOD_CONSTANTS      , only: CCON, DP
  use MOD_READ_MOLECULES , only: mol1, mol2, atom, molecule, dimer

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
    contains
      procedure, pass, public                        :: Calc_cross  => Calc_LJC_cross
      procedure, pass, public                        :: Calc_energy => Calc_LJC_energy 
  end type ljc_dimer

  type( ljc_dimer )                                  :: ljc_dimers

  contains

    !------------------------------------------------------------------------------
    !> @brief This routine reads atom names, charges, sigma and epsilon values from parameter files.
    !> @author Felippe M. Colombari

    subroutine  Read_LJC_Params( this, ljc_filename, numat )
      use MOD_CONSTANTS, only: dashline
      use mod_inquire, only: Inquire_file
      use mod_error_handling

      implicit none

      class( ljc_molecules ), intent(INOUT) :: this
      character( len = * ), intent(IN)     :: ljc_filename
      integer, intent(IN)                  :: numat
      integer                              :: i
      integer                              :: ios         = 0
      integer, parameter                   :: file_unit   = 12            
      character( len = * ), parameter      :: file_format = "formatted"
      character( len = * ), parameter      :: file_access = "sequential"
      character( len = 2 )                 :: dummy
      character( len = 10 )                :: line_number
      integer                              :: ierr
      type(error)                          :: err


      call Inquire_file( file_unit, ljc_filename, file_format, file_access )

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
          write(*,'( /, T3, A )') dashline

          stop

        endif

      enddo

      close( file_unit )
      
      return
    end subroutine  Read_LJC_params

    !------------------------------------------------------------------------------
    !> @brief This routine checks for inconsistent entries on parameter files.
    !> @author Felippe M. Colombari

    subroutine  Check_LJC_Params( this, mol, numat, ljc_filename )
      use MOD_CONSTANTS, only: dashline
      use mod_error_handling

      implicit none

      class( ljc_molecules ), intent(inout) :: this
      integer, intent(IN)                   :: numat, mol
      character( len = * ), intent(IN)      :: ljc_filename
      type(error)                           :: err

      write(*,'( /, T3, A )')                             dashline
      write(*,'( T5, "Checking molecule ", i1, " ..." )') mol
      write(*,'( T3, A )')                                dashline
      write(*,'( /, T5, "Number of atoms", T65, i10   )') numat
      write(*,'( /, T5, "Total charge", T65, f10.5, / )') sum(this % ljc_atoms % q)

      if ( ANY( this % ljc_atoms % sig < 0.0_DP ) ) then

        call err%error('e',message="while reading file "//trim(adjustl(ljc_filename))//".")
        call err%error('e',check="for sigma values < 0.")
        write(*,'( /, T3, A )') dashline

        stop

      else if ( ANY( this % ljc_atoms % eps < 0.0_DP ) ) then

        call err%error('e',message="while reading file "//trim(adjustl(ljc_filename))//".")
        call err%error('e',check="for epsilon values < 0.")
        write(*,'( /, T3, A )') dashline
        
        stop

      else

        write(*,'( T5, "LJ parameters:", T73, "OK" )') 

      endif

      return 
    end subroutine  Check_LJC_params

    !------------------------------------------------------------------------------
    !> @brief This routine calculates crossing values for LJ + coulomb potential before the loops.
    !> @author Felippe M. Colombari

    subroutine  Calc_LJC_cross( this )
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

      this % eps_ij     = 0.0_DP
      this % sig_ij     = 0.0_DP
      this % q_ij       = 0.0_DP
      
      do j = 1, mol2 % num_atoms
        
        do i = 1, mol1 % num_atoms
      
          this % eps_ij(i,j)     = 4.0_DP * dsqrt( ( mol1_ljc % ljc_atoms(i) % eps * mol2_ljc % ljc_atoms(j) % eps ) )
          this % sig_ij(i,j)     =                   mol1_ljc % ljc_atoms(i) % sig * mol2_ljc % ljc_atoms(j) % sig
          this % q_ij(i,j)       =            CCON * mol1_ljc % ljc_atoms(i) % q   * mol2_ljc % ljc_atoms(j) % q

        enddo 

      enddo 

      return
    end subroutine  Calc_LJC_cross

    !------------------------------------------------------------------------------
    !> @brief This routine calculates LJ + coulomb energy for each valid configuration.
    !> @author Felippe M. Colombari
    
    subroutine  Calc_LJC_energy( this, g, r, t )
      use MOD_INPUT_READ, only: rcut_sqr, atom_overlap, inter_energy, scale_factor

      implicit none

      class( ljc_dimer ), intent(INOUT) :: this
      integer, intent(IN)               :: g, r, t
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

          if ( ( rijrij < rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) /= 'X' ) ) then

            atom_overlap( g, r, t ) = .true.

            this % pot_lj       = 1.0E10_DP

            this % pot_coul     = 0.0_DP

            exit jlp

          else if ( ( rijrij < rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) == 'X' ) ) then

            continue

          else if ( ( rijrij < rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) == 'X' ) ) then

            continue

          else if ( ( rijrij > rcut_sqr ) .and. ( mol2 % atoms(j) % symbol(1:1) /= 'X' ) ) then

            lj_factor_cube = ( this % sig_ij(i,j) / rijrij ) ** 3

            this % pot_lj = this % pot_lj + this % eps_ij(i,j) * ( lj_factor_cube * ( lj_factor_cube - 1 ) )

          endif

          if ( dabs(this % q_ij(i,j)) < 1.0d-6 ) then

            continue

          else

            rij = dsqrt(rijrij)
            
            this % pot_coul = this % pot_coul + this % q_ij(i,j) / rij

          endif

        enddo ilp

      enddo jlp

      inter_energy( g, r, t ) = ( this % pot_coul + this % pot_lj ) / scale_factor
  
      return
    end subroutine  CALC_LJC_ENERGY

end module MOD_POT_LJC
