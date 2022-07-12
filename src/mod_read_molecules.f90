!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2022 Themis developers
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
!> @file   mod_read_molecules.f90
!> @author Felippe M. Colombari
!> @brief  This module reads coordinate files for xyz_files 1 and 2
!> @date - Jan, 2018                                                           
!> - independent module created                                                
!> @date - Nov, 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Apr, 2020
!> - major revision
!> @date - May, 2021
!> - added support for ensemble of structures of molecule 2
!> @date - Jul, 2022
!> - internal energies of molecule 2 are read from independent file
!> - added support to PDB files (read and write)
!> - PDB format http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
!---------------------------------------------------------------------------------------------------

module mod_read_molecules
  use iso_fortran_env , only : output_unit
  use mod_constants   , only : DP, PI
  use mod_input_read  , only : rcut_sqr, nconf2!, dphi, sinphi, cosphi, dtheta, sintheta, costheta

  implicit none
  
  type atom 
    real( kind = dp )                              :: xyz(3), xyz_old(3), xyz_rot(3)
    character( len = 4 )                           :: symbol
    ! for PDB
    real( kind = DP )                              :: occup
    real( kind = DP )                              :: t_fact
    integer                                        :: atom_num
    integer                                        :: res_num
    character( len = 1 )                           :: chain_id, ins_code, alt_loc
    character( len = 2 )                           :: elem_symbol
    character( len = 3 )                           :: res_name
    character( len = 6 )                           :: rec_type
  end type atom

  type molecule
    type( atom ), allocatable,dimension(:,:)       :: atoms
    real( kind = dp ), allocatable, dimension(:,:) :: atom_ref1, atom_ref2
    real( kind = DP ), allocatable, dimension(:)   :: conf_energy
    real( kind = dp ), allocatable, dimension(:)   :: rot_vector
    real( kind = DP ), allocatable, dimension(:)   :: dtheta, dphi
    real( kind = DP ), allocatable, dimension(:)   :: cosphi, sinphi, costheta, sintheta
    integer                                        :: num_atoms
  contains
    procedure, pass                                :: Read_molecule
    procedure, pass                                :: Check_molecule
    procedure, pass                                :: Translate_molecule
    procedure, pass                                :: Align_molecule
    procedure, pass                                :: Rotate_molecule
  end type molecule

  type( molecule )                                 :: mol1, mol2

  type dimer
    type( molecule ), dimension(2)                 :: molecules
    real( kind = DP )                              :: pot_none
  contains
    procedure, pass                                :: Build_dimer
    procedure, pass                                :: Write_xyz
    procedure, pass                                :: Write_mop
    procedure, pass                                :: Write_pdb
    procedure, pass, public                        :: Calc_cross  => Calc_none_cross
    procedure, pass, public                        :: Calc_energy => Calc_none_energy
  end type dimer

  type( dimer )                                    :: dimers

  real( kind = DP )                                :: dalpha
  real( kind = DP )                                :: cosalpha, sinalpha
  real( kind = DP )                                :: cosphir   =  0.0_DP
  real( kind = DP )                                :: sinphir   =  1.0_DP
  real( kind = DP )                                :: costhetar = -1.0_DP
  real( kind = DP )                                :: sinthetar =  0.0_DP

contains

  !------------------------------------------------------------------------------
  !> @brief This routine does nothing :( 
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  subroutine Calc_none_cross( this )
    class( dimer ), intent(INOUT) :: this

    this % pot_none = 0.0_DP

    return
  end subroutine Calc_none_cross

  !------------------------------------------------------------------------------
  !> @brief This routine calculates interatomic distances and marks forbbiden microstates 
  !> @author Felippe M. Colombari
  !------------------------------------------------------------------------------
  Subroutine Calc_none_energy( this, n_rot2, n_rot1, n_conf, n_trans )
    use mod_input_read , only : atom_overlap, inter_energy

    implicit none

    class( dimer ), intent(INOUT) :: this
    integer, intent(IN)           :: n_rot2, n_rot1, n_conf, n_trans
    integer                       :: n_atom_1, n_atom_2
    real( kind=DP )               :: rijrij
 
    jlp: do n_atom_2 = 1, mol2 % num_atoms

      ilp: do n_atom_1 = 1, mol1 % num_atoms

        rijrij = sum( ( mol1 % atoms( 1, n_atom_1 ) % xyz(:) - mol2 % atoms( n_conf, n_atom_2 ) % xyz(:) ) * &
                      ( mol1 % atoms( 1, n_atom_1 ) % xyz(:) - mol2 % atoms( n_conf, n_atom_2 ) % xyz(:) ) )

        if ( rijrij < rcut_sqr ) then

          if ( ( mol2 % atoms( n_conf, n_atom_2 ) % symbol(1:1) == 'X' ) &
            .or. ( mol1 % atoms( 1, n_atom_1 ) % symbol(1:1) == 'X' ) ) then

            continue
              
          else
              
            atom_overlap( n_rot2, n_rot1, n_conf, n_trans ) = .true.

            this % pot_none = 1.0E10_DP

            exit jlp

          endif

        else

          atom_overlap( n_rot2, n_rot1, n_conf, n_trans ) = .false.

          this % pot_none = 1.0E10_DP

        endif

      enddo ilp

    enddo jlp

    inter_energy( n_rot2, n_rot1, n_conf, n_trans ) = this % pot_none

    return
  End subroutine  Calc_none_energy

  !---------------------------------------------------------------------------
  !> @brief This routine reads coordinates for xyz_files 1 and 2.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Read_molecule ( this, molecule_filename, num_conformations, structure_format ) ! pass file format
    use mod_inquire        , only: Inquire_file, Get_new_unit
    use mod_error_handling

    implicit none

    class( molecule ), intent(inout) :: this
    character( len = * ), intent(in) :: molecule_filename, structure_format
    integer, intent(in)              :: num_conformations
    integer                          :: n_conf, c_natom
    integer                          :: file_unit   
    integer                          :: ios         = 0
    character( len = * ), parameter  :: file_format = "formatted"
    character( len = * ), parameter  :: file_access = "sequential"
    character( len = * ), parameter  :: file_status = "old"

    real(kind = DP )                 :: min_conf
    character( len = 16 )            :: line_energy

    integer                          :: ierr
    type(error)                      :: err

    character(len=80)                :: line_buffer
    
    file_unit = Get_new_unit(10)

    if ( structure_format == "XYZ" ) then

      call Inquire_file( file_unit, molecule_filename, file_status, file_format, file_access )

      read(file_unit,*,iostat=ios) this % num_atoms

      rewind( file_unit ) !

      if ( allocated ( this % atoms ) ) deallocate ( this % atoms )
      allocate( this % atoms( num_conformations, this % num_atoms ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      if ( allocated ( this % conf_energy ) ) deallocate ( this % conf_energy )
      allocate( this % conf_energy( num_conformations ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      this % conf_energy = 0.0_DP

      do n_conf = 1, num_conformations !

        read( file_unit, *, iostat=ios ) this % num_atoms 
        read( file_unit, '(A)', iostat=ios ) line_energy
    
        if (line_energy .eq. '') then
    
          line_energy="0"
      
        else
      
          read(line_energy,*) this % conf_energy( n_conf )
      
        endif

        if ( ios/=0 ) then
        
          call err % error( 'e', message = "error while reading "//molecule_filename )

          call err % error( 'e', tip = "check number of conformations set in INPUT file.")

          stop

        endif

        do c_natom = 1, this % num_atoms

          read(file_unit,*,iostat=ios) this % atoms( n_conf, c_natom ) % symbol, this % atoms( n_conf, c_natom ) % xyz(:)

        enddo

      enddo

      min_conf = minval( this % conf_energy )

      do n_conf = 1, num_conformations

        this % conf_energy( n_conf ) = this % conf_energy( n_conf ) - min_conf

      enddo

    else if ( structure_format == "PDB" ) then

      call Inquire_file( file_unit, molecule_filename, file_status, file_format, file_access )

      ! get number of atoms
      
      ios = 0
      this % num_atoms = 0

      do
    
        read( file_unit, '(A)', iostat = ios ) line_buffer 

        if ( ios /= 0 ) exit

        if ( ( line_buffer(1:6) == "HETATM" ) .or. ( line_buffer(1:4) == "ATOM" ) ) then

          this % num_atoms = this % num_atoms + 1

        endif

        if ( line_buffer(1:6) == "ENDMDL" ) then

          exit

        endif

      enddo

      rewind( file_unit )

      if ( allocated ( this % atoms ) ) deallocate ( this % atoms )
      allocate( this % atoms( num_conformations, this % num_atoms ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      if ( allocated ( this % conf_energy ) ) deallocate ( this % conf_energy )
      allocate( this % conf_energy( num_conformations ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      this % conf_energy = 0.0_DP

      ios = 0
      c_natom = 0
      n_conf = 1

      do
    
        read( file_unit, '(A)', iostat = ios ) line_buffer 

        if ( ios /= 0 ) exit

        if ( ( line_buffer(1:6) == "HETATM" ) .or. ( line_buffer(1:4) == "ATOM" ) ) then

          c_natom = c_natom + 1

          read(line_buffer(1:6)  , * ) this % atoms( n_conf, c_natom ) % rec_type  
          read(line_buffer(7:11) , * ) this % atoms( n_conf, c_natom ) % atom_num
          read(line_buffer(14:16), * ) this % atoms( n_conf, c_natom ) % symbol         
          
          !read(line_buffer(17:17), * ) this % atoms( n_conf, c_natom ) % alt_loc        not needed?
          !this % atoms( n_conf, c_natom ) % alt_loc = "A"
          
          read(line_buffer(18:20), * ) this % atoms( n_conf, c_natom ) % res_name
          
          !read(line_buffer(22:22), * ) this % atoms( n_conf, c_natom ) % chain_id       often missing?
          
          read(line_buffer(23:26), * ) this % atoms( n_conf, c_natom ) % res_num
          
          !read(line_buffer(27:27), * ) this % atoms( n_conf, c_natom ) % ins_code not needed?
          !this % atoms( n_conf, c_natom ) % ins_code = "A"

          read(line_buffer(31:38), * ) this % atoms( n_conf, c_natom ) % xyz(1)
          read(line_buffer(39:46), * ) this % atoms( n_conf, c_natom ) % xyz(2)             
          read(line_buffer(47:54), * ) this % atoms( n_conf, c_natom ) % xyz(3)             
          
          !read(line_buffer(55:60), '(A)') this % atoms( n_conf, c_natom ) % occup   often missing?
          this % atoms( n_conf, c_natom ) % occup = 1.00 
          
          !read(line_buffer(61:66), '(A)') this % atoms( n_conf, c_natom ) % t_fact    often missing?
          this % atoms( n_conf, c_natom ) % t_fact = 0.00 

          !read(line_buffer(77:78), '(A)') this % atoms( n_conf, c_natom ) % elem_symbol often missing?
          this % atoms( n_conf, c_natom ) % elem_symbol = this % atoms( n_conf, c_natom ) % symbol(1:1) 

        else if ( line_buffer(1:6) == "ENDMDL" ) then

          c_natom = 0

          n_conf = n_conf + 1

        endif

      enddo

    endif

    close(file_unit)

    return
  end subroutine Read_molecule

  !---------------------------------------------------------------------------
  !> @brief This routine checks the reference sites molecules 1 and 2.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Check_molecule( this, reference, vector, mol_number )
    use mod_constants      , only: dashline
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

      write(output_unit,'(/, T3, A)') dashline

      call err%error( 'e', message = "while reading INPUT file." )
      call err%error( 'e', check   = "translation point for molecule "//numol )
      call err%error( 'e', tip     = "use an integer ( n > 0 and n <= "//trim(adjustl(numat))//" ) &
                                     &to specify the index of translation point." )

      write(output_unit,'(/, T3, A)') dashline
        
      stop 

    endif

    if ( ( vector < 1 ) .or. ( vector > this % num_atoms ) ) then

      write(output_unit,'(/, T3, A)') dashline
      
      call err%error( 'e', message = "while reading INPUT file." )
      call err%error( 'e', check   = "translation point for molecule "//numol )
      call err%error( 'e', tip     = "use an integer ( n > 0 and n <= "//trim(adjustl(numat))//" ) &
                                     &to specify the index of rotation point." )

      write(output_unit,'(/, T3, A)') dashline
        
      stop 

    endif

    if ( reference == vector ) then

      write(output_unit,'(/, T3, A)') dashline

      call err%error( 'e', message = "while reading INPUT file." )
      call err%error( 'e', check   = "translation and rotation points for molecule "//numol )
      call err%error( 'e', tip     = "use different indexes for translation point and  rotation point." )

      write(output_unit,'(/, T3, A)') dashline
        
      stop 

    endif

    return
  end subroutine Check_molecule

  !---------------------------------------------------------------------------
  !> @brief This routine translates the molecule 
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Translate_molecule( this, reference, n_conf )

    implicit none

    class( molecule ), intent(inout) :: this
    integer, intent(in)              :: reference, n_conf
    integer                          :: c_natom

    if ( .not. allocated(this % atom_ref1) ) allocate( this % atom_ref1( nconf2, 3 ) )
    if ( .not. allocated(this % atom_ref2) ) allocate( this % atom_ref2( nconf2, 3 ) )

    this % atom_ref1( n_conf, : ) = this % atoms( n_conf, reference ) % xyz(:)

    do c_natom = 1, this % num_atoms

      this % atoms( n_conf, c_natom ) % xyz(:) = this % atoms( n_conf, c_natom ) % xyz(:) - this % atom_ref1( n_conf, : )

    enddo

    return
  end subroutine Translate_molecule

  !---------------------------------------------------------------------------
  !> @brief This routine alings the molecule rotation vector along Z axis
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Align_molecule( this, vector, reference, n_conf )

    implicit none

    class( molecule ), intent(inout) :: this
    integer, intent(IN)              :: vector, reference!, nconf2
    integer, intent(IN)              :: n_conf
    integer                          :: c_natom

    if ( .not. allocated(this % rot_vector) ) allocate( this % rot_vector( nconf2 ) )
    if ( .not. allocated(this % dphi) )       allocate( this % dphi( nconf2 ) )
    if ( .not. allocated(this % dtheta) )     allocate( this % dtheta( nconf2 ) )
    if ( .not. allocated(this % cosphi) )     allocate( this % cosphi( nconf2 ) )
    if ( .not. allocated(this % sinphi) )     allocate( this % sinphi( nconf2 ) )
    if ( .not. allocated(this % costheta) )   allocate( this % costheta( nconf2 ) )
    if ( .not. allocated(this % sintheta) )   allocate( this % sintheta( nconf2 ) )

    this % atom_ref2( n_conf, : ) = this % atoms( n_conf, vector ) % xyz(:)

    this % rot_vector( n_conf ) = dsqrt( sum( ( this % atom_ref2( n_conf, : ) &
                                              - this % atoms( n_conf, reference ) % xyz(:) )**2 ) )
    
    !! CHANGES REFERENCE AXIS TO ALIGN MOL
    
    this % dphi( n_conf )   = ( PI / 2.0_DP ) - datan2( this % atoms( n_conf, vector ) % xyz(2) , &
                                                        this % atoms( n_conf, vector ) % xyz(1) ) 

    this % dtheta( n_conf ) = dacos( this % atoms( n_conf , vector ) % xyz(3) / this % rot_vector( n_conf ) ) 

    !! ROTATE AROUND z AXIS: PLACE ROTATION VECTOR AT yz PLANE

    this % cosphi( n_conf ) = dcos( this % dphi( n_conf ) )
    this % sinphi( n_conf ) = dsin( this % dphi( n_conf ) )

    this % costheta( n_conf ) = dcos( this % dtheta( n_conf ) )
    this % sintheta( n_conf ) = dsin( this % dtheta( n_conf ) )

    do c_natom = 1, this % num_atoms

      this % atoms( n_conf, c_natom ) % xyz_rot(1) = this % atoms( n_conf, c_natom ) % xyz(1) * this % cosphi( n_conf ) &
                                                   - this % atoms( n_conf, c_natom ) % xyz(2) * this % sinphi( n_conf )

      this % atoms( n_conf, c_natom ) % xyz_rot(2) = this % atoms( n_conf, c_natom ) % xyz(1) * this % sinphi( n_conf ) &
                                                   + this % atoms( n_conf, c_natom ) % xyz(2) * this % cosphi( n_conf )

      this % atoms( n_conf, c_natom ) % xyz_rot(3) = this % atoms( n_conf, c_natom ) % xyz(3)

      this % atoms( n_conf, c_natom ) % xyz(:)     = this % atoms( n_conf, c_natom ) % xyz_rot(:)

      this % atoms( n_conf, c_natom ) % xyz_old(:) = this % atoms( n_conf, c_natom ) % xyz(:)

    enddo

    !deallocate( dphi )
    !deallocate( dtheta )
    !deallocate( cosphi )
    !deallocate( sinphi )
    !deallocate( costheta )
    !deallocate( sintheta )

    return
  end subroutine Align_molecule

  !---------------------------------------------------------------------------
  !> @brief This routine performs rotations of mol2 around mol1
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Rotate_molecule( this, rotation2, n_conf )
    use mod_constants  , only : deg2rad
    use mod_input_read , only : axis_rot_range, axis_rot_moves

    implicit none

    class( molecule ), intent(inout) :: this
    integer                          :: c_natom
    integer, intent(IN)              :: n_conf
    integer, intent(IN)              :: rotation2!, num_conformations
 
    dalpha = ( axis_rot_range / deg2rad ) * ( rotation2 - 1 ) / axis_rot_moves

    cosalpha = dcos(dalpha)
    sinalpha = dsin(dalpha)

    do c_natom = 1, this % num_atoms
   
      !! ROTATE AROUND x AXIS: PLACE ROTATION VECTOR AT z AXIS

      this % atoms( n_conf, c_natom ) % xyz_rot(1) = this % atoms( n_conf, c_natom ) % xyz(1)

      this % atoms( n_conf, c_natom ) % xyz_rot(2) = this % atoms( n_conf, c_natom ) % xyz(2) * this % costheta( n_conf ) &
                                                   - this % atoms( n_conf, c_natom ) % xyz(3) * this % sintheta( n_conf )

      this % atoms( n_conf, c_natom ) % xyz_rot(3) = this % atoms( n_conf, c_natom ) % xyz(2) * this % sintheta( n_conf ) &
                                                   + this % atoms( n_conf, c_natom ) % xyz(3) * this % costheta( n_conf )

      this % atoms( n_conf, c_natom ) % xyz(:) = this % atoms( n_conf, c_natom ) % xyz_rot(:)

      !! DONE: ROTATION VECTOR IS AT ORIGIN !! PRECESSION MOVES AROUND z AXIS

      this % atoms( n_conf, c_natom ) % xyz_rot(1) = this % atoms( n_conf, c_natom ) % xyz(1) * cosalpha &
                                                   - this % atoms( n_conf, c_natom ) % xyz(2) * sinalpha

      this % atoms( n_conf, c_natom ) % xyz_rot(2) = this % atoms( n_conf, c_natom ) % xyz(1) * sinalpha &
                                                   + this % atoms( n_conf, c_natom ) % xyz(2) * cosalpha

      this % atoms( n_conf, c_natom ) % xyz_rot(3) = this % atoms( n_conf, c_natom ) % xyz(3)

      this % atoms( n_conf, c_natom ) % xyz(:) = this % atoms( n_conf, c_natom ) % xyz_rot(:)

      !! ROTATION MOVE TO rth POINT OF THE SPHERICAL GRID !! ROTATION AROUND y AXIS

      this % atoms( n_conf, c_natom ) % xyz_rot(1) = this % atoms( n_conf, c_natom ) % xyz(1) * costhetar &
                                                   + this % atoms( n_conf, c_natom ) % xyz(3) * sinthetar

      this % atoms( n_conf, c_natom ) % xyz_rot(2) = this % atoms( n_conf, c_natom ) % xyz(2)

      this % atoms( n_conf, c_natom ) % xyz_rot(3) = -this % atoms( n_conf, c_natom ) % xyz(1) * sinthetar &
                                                   +  this % atoms( n_conf, c_natom ) % xyz(3) * costhetar

      this % atoms( n_conf, c_natom ) % xyz(:) = this % atoms( n_conf, c_natom ) % xyz_rot(:)

      !! ROTATION AROUND z AXIS TO PLACE ROTATION VECTOR AT THE rth SPHERE POINT !!

      this % atoms( n_conf, c_natom ) % xyz_rot(1) = this % atoms( n_conf, c_natom ) % xyz(1) * cosphir &
                                                   - this % atoms( n_conf, c_natom ) % xyz(2) * sinphir

      this % atoms( n_conf, c_natom ) % xyz_rot(2) = this % atoms( n_conf, c_natom ) % xyz(1) * sinphir &
                                                   + this % atoms( n_conf, c_natom ) % xyz(2) * cosphir

      this % atoms( n_conf, c_natom ) % xyz_rot(3) = this % atoms( n_conf, c_natom ) % xyz(3) 

      this % atoms( n_conf, c_natom ) % xyz(:) = this % atoms( n_conf, c_natom ) % xyz_rot(:)

    enddo

    return
  end subroutine Rotate_molecule

  subroutine Build_dimer( this, n_conf, structure_format )
    use mod_error_handling

    implicit none

    class( dimer ), intent(inout)     :: this
    character( len = * ), intent(in)  :: structure_format
    integer, intent(in)               :: n_conf
    integer                           :: n_atom_1, n_atom_2
    integer                           :: ierr
    type(error)                       :: err
    character( len = 53 )             :: pdb_fmt

    allocate( this % molecules(1) % atoms( 1 , mol1 % num_atoms ),stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")
    allocate( this % molecules(2) % atoms( nconf2, mol2 % num_atoms ),stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    if ( structure_format == "XYZ" ) then

      do n_atom_1 = 1, mol1 % num_atoms

        this % molecules(1) % atoms( 1, n_atom_1 ) % xyz(:) = mol1 % atoms( 1, n_atom_1 ) % xyz(:)
        this % molecules(1) % atoms( 1, n_atom_1 ) % symbol = mol1 % atoms( 1, n_atom_1 ) % symbol

      enddo

      do n_atom_2 = 1, mol2 % num_atoms

        this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(:) = mol2 % atoms( n_conf, n_atom_2 ) % xyz(:)
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % symbol = mol2 % atoms( n_conf, n_atom_2 ) % symbol

      enddo

    else if ( structure_format == "PDB" ) then
      
      pdb_fmt='(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,10x,a2)'

      do n_atom_1 = 1, mol1 % num_atoms

        this % molecules(1) % atoms( 1, n_atom_1 ) % rec_type    = mol1 % atoms( 1, n_atom_1 ) % rec_type
        this % molecules(1) % atoms( 1, n_atom_1 ) % atom_num    = mol1 % atoms( 1, n_atom_1 ) % atom_num
        this % molecules(1) % atoms( 1, n_atom_1 ) % symbol      = mol1 % atoms( 1, n_atom_1 ) % symbol
        this % molecules(1) % atoms( 1, n_atom_1 ) % alt_loc     = "A" 
        this % molecules(1) % atoms( 1, n_atom_1 ) % res_name    = mol1 % atoms( 1, n_atom_1 ) % res_name
        this % molecules(1) % atoms( 1, n_atom_1 ) % chain_id    = "A"
        this % molecules(1) % atoms( 1, n_atom_1 ) % res_num     = mol1 % atoms( 1, n_atom_1 ) % res_num
        this % molecules(1) % atoms( 1, n_atom_1 ) % ins_code    = "A" 
        this % molecules(1) % atoms( 1, n_atom_1 ) % xyz(:)      = mol1 % atoms( 1, n_atom_1 ) % xyz(:)
        this % molecules(1) % atoms( 1, n_atom_1 ) % occup       = mol1 % atoms( 1, n_atom_1 ) % occup
        this % molecules(1) % atoms( 1, n_atom_1 ) % t_fact      = mol1 % atoms( 1, n_atom_1 ) % t_fact
        this % molecules(1) % atoms( 1, n_atom_1 ) % elem_symbol = mol1 % atoms( 1, n_atom_1 ) % elem_symbol

      enddo

      do n_atom_2 = 1, mol2 % num_atoms

        this % molecules(2) % atoms( n_conf, n_atom_2 ) % rec_type    = mol2 % atoms( n_conf, n_atom_2 ) % rec_type
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % atom_num    = mol2 % atoms( n_conf, n_atom_2 ) % atom_num
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % symbol      = mol2 % atoms( n_conf, n_atom_2 ) % symbol
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % alt_loc     = "B" 
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % res_name    = mol2 % atoms( n_conf, n_atom_2 ) % res_name
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % chain_id    = "B"
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % res_num     = mol2 % atoms( n_conf, n_atom_2 ) % res_num
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % ins_code    = "B" 
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(:)      = mol2 % atoms( n_conf, n_atom_2 ) % xyz(:)
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % occup       = mol2 % atoms( n_conf, n_atom_2 ) % occup
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % t_fact      = mol2 % atoms( n_conf, n_atom_2 ) % t_fact
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % elem_symbol = mol2 % atoms( n_conf, n_atom_2 ) % elem_symbol

      enddo

    endif

    return
  end subroutine Build_dimer

  subroutine Write_xyz( this, lowest, n_conf )

    implicit none

    class( dimer ), intent(inout)     :: this
    real( kind = DP ), intent(IN)     :: lowest
    integer                           :: n_atom_1, n_atom_2
    integer, intent(IN)               :: n_conf

    write(66,*) mol1 % num_atoms + mol2 % num_atoms

    write(66,'("Energy = ",es15.7E2)') lowest

    do n_atom_1 = 1, mol1 % num_atoms

      write(66,'(a4,1x,3(f10.4,1x))') this % molecules(1) % atoms( 1, n_atom_1 ) % symbol, &
                                      this % molecules(1) % atoms( 1, n_atom_1 ) % xyz(:)

    enddo

    do n_atom_2 = 1, mol2 % num_atoms

      write(66,'(a4,1x,3(f10.4,1x))') this % molecules(2) % atoms( n_conf, n_atom_2 ) % symbol, &
                                      this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(:)

    enddo

    deallocate( this % molecules(1) % atoms )
    deallocate( this % molecules(2) % atoms )

    return
  end subroutine Write_xyz

  subroutine Write_pdb( this, n_conf )

    implicit none

    class( dimer ), intent(inout)       :: this
    integer                             :: n_atom_1, n_atom_2
    integer, intent(IN)                 :: n_conf
    character( len = 53 )               :: pdb_fmt

    pdb_fmt='(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,10x,a2)'

    write(66,'(a5,i9)') "MODEL", 1

    do n_atom_1 = 1, mol1 % num_atoms

      write(66,pdb_fmt) this % molecules(1) % atoms( 1, n_atom_1 ) % rec_type, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % atom_num, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % symbol, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % alt_loc, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % res_name, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % chain_id, & 
                        this % molecules(1) % atoms( 1, n_atom_1 ) % res_num, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % ins_code, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % xyz(:), &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % occup, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % t_fact, &
                        this % molecules(1) % atoms( 1, n_atom_1 ) % elem_symbol

    enddo

    do n_atom_2 = 1, mol2 % num_atoms

      write(66,pdb_fmt) this % molecules(2) % atoms( n_conf, n_atom_2 ) % rec_type, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % atom_num, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % symbol, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % alt_loc, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % res_name, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % chain_id, & 
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % res_num, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % ins_code, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(:), &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % occup, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % t_fact, &
                        this % molecules(2) % atoms( n_conf, n_atom_2 ) % elem_symbol

    enddo

    write(66,'(a3)') "TER"
    write(66,'(a6)') "ENDMDL"

    deallocate( this % molecules(1) % atoms )
    deallocate( this % molecules(2) % atoms )

    return
  end subroutine Write_pdb

  subroutine Write_mop( this, header, n_conf )

    implicit none

    class( dimer ), intent(inout)      :: this
    character( len = 128 ), intent(in) :: header
    integer                            :: n_atom_1, n_atom_2
    integer, intent(IN)                :: n_conf

    write(66,'(A)') header
    write(66,*)
    write(66,*)

    do n_atom_1 = 1, mol1 % num_atoms

      write(66,'(a4,1x,3(f10.4,2x,a1,1x))') this % molecules(1) % atoms( n_conf, n_atom_1 ) % symbol, &
        this % molecules(1) % atoms( n_conf, n_atom_1 ) % xyz(1), "0", &
        this % molecules(1) % atoms( n_conf, n_atom_1 ) % xyz(2), "0", &
        this % molecules(1) % atoms( n_conf, n_atom_1 ) % xyz(3), "0"

    enddo

    do n_atom_2 = 1, mol2 % num_atoms

      write(66,'(a4,1x,3(f10.4,2x,a1,1x))') this % molecules(2) % atoms( n_conf, n_atom_2 ) % symbol, &
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(1), "0", &
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(2), "0", &
        this % molecules(2) % atoms( n_conf, n_atom_2 ) % xyz(3), "0"

    enddo

    deallocate( this % molecules(1) % atoms )
    deallocate( this % molecules(2) % atoms )

    return
  end subroutine Write_mop

end module mod_read_molecules
