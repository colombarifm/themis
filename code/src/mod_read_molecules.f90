!------------------------------------------------------------------------------
! THEMIS: A program for 
! Copyright (C) 2017 Felippe M. Colombari
!------------------------------------------------------------------------------
!> @brief This module reads coordinate files for xyz_files 1 and 2
!> @author Felippe M. Colombari
!> Laboratório de Química Teórica, LQT -- UFSCar
!> @note Revision history
!> @date Jan 2018: module created
!> @date Jul 2019: OO routine
!------------------------------------------------------------------------------

module MOD_READ_MOLECULES
  use MOD_CONSTANTS, only: DP, PI
  use MOD_INPUT_READ, only: rcut_sqr

  implicit none
  
  type atom 
    real( kind = dp )                       :: xyz(3), xyz_old(3), xyz_rot(3)
    character( len = 4 )                    :: symbol
  end type atom

  type molecule
    type( atom ), allocatable,dimension(:) :: atoms
    real( kind = dp )                      :: atom_ref1(3), atom_ref2(3), rot_vector
    integer                                :: num_atoms
    contains
      procedure, pass                      :: Read_molecule
      procedure, pass                      :: Check_molecule
      procedure, pass                      :: Translate_molecule
      procedure, pass                      :: Align_molecule
      procedure, pass                      :: Rotate_molecule
  end type molecule

  type( molecule )                         :: mol1, mol2

  type dimer
    type( molecule ), dimension(2)         :: molecules
    real( kind = DP )                      :: pot_none
    contains
      procedure, pass                      :: Build_dimer
      procedure, pass                      :: Write_xyz
      procedure, pass                      :: Write_mop
      procedure, pass, public              :: Calc_cross  => Calc_none_cross
      procedure, pass, public              :: Calc_energy => Calc_none_energy
  end type dimer

  type( dimer )                            :: dimers

  real( kind = DP )                        :: dtheta, dphi, dalpha
  real( kind = DP )                        :: cosphi, sinphi, costheta, sintheta, cosalpha, sinalpha
  real( kind = DP )                        :: cosphir   =  0.0_DP
  real( kind = DP )                        :: sinphir   =  1.0_DP
  real( kind = DP )                        :: costhetar = -1.0_DP
  real( kind = DP )                        :: sinthetar =  0.0_DP

  contains

    subroutine Calc_none_cross( this )
      class( dimer ), intent(INOUT) :: this

      this % pot_none = 0.0_DP

      return
    end subroutine  calc_none_cross

    !------------------------------------------------------------------------------
    !> @brief This routine calculates interatomic distances and 
    !> @author Felippe M. Colombari
    
    Subroutine Calc_none_energy( this, g, r, t )
      use MOD_INPUT_READ, only: atom_overlap, inter_energy

      class( dimer ), intent(INOUT) :: this
      integer, intent(IN)           :: g, r, t
      integer                       :: i, j
      real( kind=DP )               :: rijrij
 
      jlp: do j = 1, mol2 % num_atoms

        ilp: do i = 1, mol1 % num_atoms

          rijrij = sum( ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) * &
                        ( mol1 % atoms(i) % xyz(:) - mol2 % atoms(j) % xyz(:) ) )

          if ( rijrij < rcut_sqr ) then

            if ( ( mol2 % atoms(j) % symbol(1:1) == 'X' ) .or. ( mol1 % atoms(i) % symbol(1:1) == 'X' ) ) then

              continue
              
            else
              
              atom_overlap( g, r, t ) = .true.

              this % pot_none = 1.0E10_DP

              exit jlp

            endif

          else

            atom_overlap( g, r, t ) = .false.

            this % pot_none = 1.0E10_DP

          endif

        enddo ilp

      enddo jlp

      inter_energy( g, r, t ) = this % pot_none

      return
    End subroutine  Calc_none_energy

    !---------------------------------------------------------------------------
    !> @author Felippe M. Colombari
    !> Laboratório de Química Teórica, LQT -- UFSCar
    !> @brief This routine reads coordinates for xyz_files 1 and 2.
    !> @note Revision history
    !> @date Jun 2017: first version
    !> @date Oct 2017: modular routine
    !> @date Jul 2019: OO routine
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

        read(file_unit,*,iostat=ios) this % atoms(i) % symbol, this % atoms(i) % xyz(:)

      enddo

      close(file_unit)

      return
    end subroutine Read_molecule

    !---------------------------------------------------------------------------
    !> @author Felippe M. Colombari
    !> Laboratório de Química Teórica, LQT -- UFSCar
    !> @brief This routine checks the reference sites molecules 1 and 2.
    !> @note Revision history
    !> @date Jul 2019: routine created
    !---------------------------------------------------------------------------	
    
    subroutine Check_molecule( this, reference, vector, mol_number )
      use mod_constants, only: dashline
      use mod_error_handling

      implicit none

      class( molecule ), intent(inout) :: this
      integer, intent(in)              :: reference, vector, mol_number
      character( len = 1 )             :: numol
      character( len = 6 )             :: numat
      type(error)                      :: err

      write(numat,'(i6)') this % num_atoms
      write(numol,'(i1)') mol_number

      if ( ( reference < 1 ) .or. ( reference > this % num_atoms ) ) then

        write(*,'(/, T3, A)') dashline
        call err%error('e',message="while reading INPUT file.")
        call err%error('e',check="translation point for molecule "//numol)
        call err%error('e',tip="use an integer ( n > 0 and n <= "//trim(adjustl(numat))//" ) &
          &to specify the index of translation point.")
        write(*,'(/, T3, A)') dashline
        
        stop 

      endif

      if ( ( vector < 1 ) .or. ( vector > this % num_atoms ) ) then

        write(*,'(/, T3, A)') dashline
        call err%error('e',message="while reading INPUT file.")
        call err%error('e',check="translation point for molecule "//numol)
        call err%error('e',tip="use an integer ( n > 0 and n <= "//trim(adjustl(numat))//" ) &
          &to specify the index of rotation point.")
        write(*,'(/, T3, A)') dashline
        
        
        stop 

      endif

      if ( reference == vector ) then

        write(*,'(/, T3, A)') dashline
        call err%error('e',message="while reading INPUT file.")
        call err%error('e',check="translation and rotation points for molecule "//numol)
        call err%error('e',tip="use different indexes for translation point and &
          &rotation point.")
        write(*,'(/, T3, A)') dashline
        
        stop 

      endif

      return
    end subroutine Check_molecule

    !---------------------------------------------------------------------------
    !> @author Felippe M. Colombari
    !> Laboratório de Química Teórica, LQT -- UFSCar
    !> @brief This routine translates the molecule 
    !---------------------------------------------------------------------------	

    subroutine Translate_molecule( this, reference )

      implicit none

      class( molecule ), intent(inout) :: this
      integer, intent(in)              :: reference
      integer                          :: i
 
      this % atom_ref1(:) = this % atoms(reference) % xyz(:)

      do i = 1, this % num_atoms

        this % atoms(i) % xyz(:) = this % atoms(i) % xyz(:) - this % atom_ref1(:)

      enddo

      return
    end subroutine Translate_molecule

    !---------------------------------------------------------------------------
    !> @author Felippe M. Colombari
    !> Laboratório de Química Teórica, LQT -- UFSCar
    !> @brief This routine alings the molecule rotation vector along Z axis
    !---------------------------------------------------------------------------	

    subroutine Align_molecule( this, vector, reference )

      implicit none

      class( molecule ), intent(inout) :: this
      integer, intent(IN)              :: vector, reference
      integer                          :: j

      this % atom_ref2(:) = this % atoms(vector) % xyz(:)

      this % rot_vector = dsqrt( sum( ( this % atom_ref2 - this % atoms(reference) % xyz(:) )**2 ) )
      
      !! CHANGES REFERENCE AXIS TO ALIGN MOL
    
      dphi   = ( PI / 2.0_DP ) - datan2( this % atoms(vector) % xyz(2) , this % atoms(vector) % xyz(1) ) 
      dtheta = dacos( this % atoms(vector) % xyz(3) / this % rot_vector ) 

      !! ROTATE AROUND z AXIS: PLACE ROTATION VECTOR AT yz PLANE

      cosphi = dcos(dphi)
      sinphi = dsin(dphi)

      costheta = dcos(dtheta)
      sintheta = dsin(dtheta)

      do j = 1, this % num_atoms

        this % atoms(j) % xyz_rot(1) = this % atoms(j) % xyz(1) * cosphi - this % atoms(j) % xyz(2) * sinphi
        this % atoms(j) % xyz_rot(2) = this % atoms(j) % xyz(1) * sinphi + this % atoms(j) % xyz(2) * cosphi
        this % atoms(j) % xyz_rot(3) = this % atoms(j) % xyz(3)

        this % atoms(j) % xyz(:)     = this % atoms(j) % xyz_rot(:)

        this % atoms(j) % xyz_old(:) = this % atoms(j) % xyz(:)

      enddo

      return
    end subroutine Align_molecule

    !---------------------------------------------------------------------------
    !> @author Felippe M. Colombari
    !> Laboratório de Química Teórica, LQT -- UFSCar
    !> @brief This routine performs rotations of mol2 around mol1
    !---------------------------------------------------------------------------	

    subroutine Rotate_molecule( this, gyration )
      use MOD_CONSTANTS, only: deg2rad
      use MOD_INPUT_READ, only: max_gyr, gyr_factor

      implicit none

      class( molecule ), intent(inout) :: this
      integer                          :: j
      integer, intent(IN)              :: gyration
 
      dalpha = ( max_gyr / deg2rad ) * ( gyration - 1 ) / gyr_factor

      cosalpha = dcos(dalpha)
      sinalpha = dsin(dalpha)

      do j = 1, this % num_atoms
   
        !! ROTATE AROUND x AXIS: PLACE ROTATION VECTOR AT z AXIS

        this % atoms(j) % xyz_rot(1) = this % atoms(j) % xyz(1)
        this % atoms(j) % xyz_rot(2) = this % atoms(j) % xyz(2) * costheta - this % atoms(j) % xyz(3) * sintheta
        this % atoms(j) % xyz_rot(3) = this % atoms(j) % xyz(2) * sintheta + this % atoms(j) % xyz(3) * costheta

        this % atoms(j) % xyz(:) = this % atoms(j) % xyz_rot(:)

        !! DONE: ROTATION VECTOR IS AT ORIGIN !! PRECESSION MOVES AROUND z AXIS

        this % atoms(j) % xyz_rot(1) = this % atoms(j) % xyz(1) * cosalpha - this % atoms(j) % xyz(2) * sinalpha
        this % atoms(j) % xyz_rot(2) = this % atoms(j) % xyz(1) * sinalpha + this % atoms(j) % xyz(2) * cosalpha
        this % atoms(j) % xyz_rot(3) = this % atoms(j) % xyz(3)

        this % atoms(j) % xyz(:) = this % atoms(j) % xyz_rot(:)

        !! ROTATION MOVE TO rth POINT OF THE SPHERICAL GRID !! ROTATION AROUND y AXIS

        this % atoms(j) % xyz_rot(1) =  this % atoms(j) % xyz(1) * costhetar + this % atoms(j) % xyz(3) * sinthetar
        this % atoms(j) % xyz_rot(2) =  this % atoms(j) % xyz(2)
        this % atoms(j) % xyz_rot(3) = -this % atoms(j) % xyz(1) * sinthetar + this % atoms(j) % xyz(3) * costhetar

        this % atoms(j) % xyz(:) = this % atoms(j) % xyz_rot(:)

        !! ROTATION AROUND z AXIS TO PLACE ROTATION VECTOR AT THE rth SPHERE POINT !!

        this % atoms(j) % xyz_rot(1) = this % atoms(j) % xyz(1) * cosphir - this % atoms(j) % xyz(2) * sinphir
        this % atoms(j) % xyz_rot(2) = this % atoms(j) % xyz(1) * sinphir + this % atoms(j) % xyz(2) * cosphir
        this % atoms(j) % xyz_rot(3) = this % atoms(j) % xyz(3) 

        this % atoms(j) % xyz(:) = this % atoms(j) % xyz_rot(:)

      enddo

      return
    end subroutine Rotate_molecule

    subroutine Build_dimer( this )
       use mod_error_handling

      implicit none

      class( dimer ), intent(inout)     :: this
      integer                           :: i, j

      integer                           :: ierr
      type(error)                       :: err

      allocate( this % molecules(1) % atoms( mol1 % num_atoms ),stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      allocate( this % molecules(2) % atoms( mol2 % num_atoms ),stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      do i = 1, mol1 % num_atoms

        this % molecules(1) % atoms(i) % xyz(:) = mol1 % atoms(i) % xyz(:)
        this % molecules(1) % atoms(i) % symbol = mol1 % atoms(i) % symbol

      enddo

      do j = 1, mol2 % num_atoms

        this % molecules(2) % atoms(j) % xyz(:) = mol2 % atoms(j) % xyz(:)
        this % molecules(2) % atoms(j) % symbol = mol2 % atoms(j) % symbol

      enddo

      return
    end subroutine Build_dimer

    subroutine Write_xyz( this, lowest )

      implicit none

      class( dimer ), intent(inout)     :: this
      real( kind = DP ), intent(IN)     :: lowest
      integer                           :: i, j

      write(66,*) mol1 % num_atoms + mol2 % num_atoms

      write(66,'("Energy = ",es15.7E2)') lowest

      do i = 1, mol1 % num_atoms

        write(66,'(a4,1x,3(f10.4,1x))') this % molecules(1) % atoms(i) % symbol, &
          this % molecules(1) % atoms(i) % xyz(:)

      enddo

      do j = 1, mol2 % num_atoms

        write(66,'(a4,1x,3(f10.4,1x))') this % molecules(2) % atoms(j) % symbol, &
          this % molecules(2) % atoms(j) % xyz(:)

      enddo

      deallocate( this % molecules(1) % atoms )
      deallocate( this % molecules(2) % atoms )

      return
    end subroutine Write_xyz

    subroutine Write_mop( this, header )

      implicit none

      class( dimer ), intent(inout)      :: this
      character( len = 128 ), intent(in) :: header
      integer                            :: i, j

      write(66,'(A)') header
      write(66,*)
      write(66,*)

      do i = 1, mol1 % num_atoms

        write(66,'(a4,1x,3(f10.4,2x,a1,1x))') this % molecules(1) % atoms(i) % symbol, &
          this % molecules(1) % atoms(i) % xyz(1), "0", &
          this % molecules(1) % atoms(i) % xyz(2), "0", &
          this % molecules(1) % atoms(i) % xyz(3), "0"

      enddo

      do j = 1, mol2 % num_atoms

        write(66,'(a4,1x,3(f10.4,2x,a1,1x))') this % molecules(2) % atoms(j) % symbol, &
          this % molecules(2) % atoms(j) % xyz(1), "0", &
          this % molecules(2) % atoms(j) % xyz(2), "0", &
          this % molecules(2) % atoms(j) % xyz(3), "0"

      enddo

      deallocate( this % molecules(1) % atoms )
      deallocate( this % molecules(2) % atoms )

      return
    end subroutine Write_mop

end module MOD_READ_MOLECULES
