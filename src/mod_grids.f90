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
!> @file   mod_grids.f90
!> @author Felippe M. Colombari
!> @brief  This module manipulates the coordinates of translation and rotation grids
!> @date - Jan, 2018                                                           
!> - independent module created                                                
!> @date - Jun, 2019                                                           
!> - subroutine for spherical grids generation was added
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!> @date - Jul 2022
!> - accept trans_factor = 0 and point_rot_factor = 0 to generate single configurations
!---------------------------------------------------------------------------------------------------

module mod_grids
  use iso_fortran_env , only : output_unit
  use mod_constants   , only : DP

  implicit none

  type point
    real( kind = DP )                              :: grid_xyz(3), grid_xyz_old(3), grid_xyz_rot(3)
    character( len = 3 )                           :: grid_symbol
  end type point

  type grid
    type( point ), allocatable, dimension(:)       :: points
    integer                                        :: numpoint
  contains
    procedure, pass                                :: Read_grid  
    procedure, pass                                :: Translate_grid
    procedure, pass                                :: Align_grid
    procedure, pass                                :: Rotate_grid
    procedure, pass                                :: Build_translation_sphere
    procedure, pass                                :: Build_rotation1_sphere
  end type grid

  type( grid )                                     :: grid_trans, grid_rot1

  real( kind = DP ), allocatable, dimension(:)     :: A, mTS, Eavg
  real( kind = DP ), allocatable, dimension(:)     :: Zrot, sumVexpVrot, probT
  real( kind = DP ), allocatable, dimension(:)     :: rr, thetar, phir
  integer                                          :: ntotal, conf_total

contains

  !---------------------------------------------------------------------------
  !> @brief This routine reads coordinates for user-defined translation grid.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Read_grid( this, grid_filename )
    use mod_inquire        , only : Inquire_file, Get_new_unit
    use mod_error_handling

    implicit none

    class( grid ), intent( inout )     :: this
    character( len = * ), intent( in ) :: grid_filename
    integer                            :: n_trans
    integer                            :: file_unit        
    integer                            :: ios         = 0
    character( len = * ), parameter    :: file_status = "old"
    character( len = * ), parameter    :: file_format = "formatted"
    character( len = * ), parameter    :: file_access = "sequential"
    integer                            :: ierr
    type(error)                        :: err

    file_unit = Get_new_unit(10)

    call Inquire_file( file_unit , grid_filename , file_status, file_format , file_access )

    read( file_unit, *, iostat = ios ) this % numpoint
    read( file_unit, * )

    if ( allocated ( this % points ) ) deallocate ( this % points )
    allocate( this % points( this % numpoint ), stat=ierr )
    if (ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

    do n_trans = 1, this % numpoint

      read( file_unit, *, iostat = ios ) this % points( n_trans ) % grid_symbol, this % points( n_trans ) % grid_xyz(:)

    enddo

    close( file_unit )

    return
  end subroutine Read_grid

  !---------------------------------------------------------------------------
  !> @brief This routine centers user-defined translation grid altogether with mol1.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Translate_grid( this )
    use mod_read_molecules !, only: mol1

    implicit none

    class( grid ), intent( inout ) :: this
    integer                        :: n_trans

    do n_trans = 1, this % numpoint

      this % points( n_trans ) % grid_xyz(:) = this % points( n_trans ) % grid_xyz(:) - mol1 % atom_ref1(1,:)

    enddo

    return
  end subroutine Translate_grid

  !---------------------------------------------------------------------------
  !> @brief This routine align user-defined translation grid along Z altogether with mol1.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Align_grid( this ) 
    use mod_read_molecules !, only: sinphi, cosphi

    implicit none

    class( grid ), intent( inout ) :: this
    integer                        :: n_trans

    do n_trans = 1, this % numpoint

      this % points( n_trans ) % grid_xyz_rot(1) = this % points( n_trans ) % grid_xyz(1) * mol1 % cosphi( 1 ) - &
                                                   this % points( n_trans ) % grid_xyz(2) * mol1 % sinphi( 1 )

      this % points( n_trans ) % grid_xyz_rot(2) = this % points( n_trans ) % grid_xyz(1) * mol1 % sinphi( 1 ) + &
                                                   this % points( n_trans ) % grid_xyz(2) * mol1 % cosphi( 1 )

      this % points( n_trans ) % grid_xyz_rot(3) = this % points( n_trans ) % grid_xyz(3)

      this % points( n_trans ) % grid_xyz(:)     = this % points( n_trans ) % grid_xyz_rot(:)

      this % points( n_trans ) % grid_xyz_old(:) = this % points( n_trans ) % grid_xyz(:)

    enddo

    return
  end subroutine Align_grid

  !---------------------------------------------------------------------------
  !> @brief This routine rotate user-defined translation grid along Z altogether with mol1.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Rotate_grid( this ) 
    use mod_read_molecules  !, only: costheta, sintheta, cosalpha, sinalpha, costhetar, sinthetar, cosphir, sinphir

    implicit none

    class( grid ), intent( inout ) :: this
    integer                        :: n_trans

    do n_trans = 1, this % numpoint
   
      !! ROTATE AROUND x AXIS: PLACE ROTATION VECTOR AT z AXIS

      this % points( n_trans ) % grid_xyz_rot(1) = this % points( n_trans ) % grid_xyz(1)

      this % points( n_trans ) % grid_xyz_rot(2) = this % points( n_trans ) % grid_xyz(2) * mol1 % costheta( 1 ) - &
                                                   this % points( n_trans ) % grid_xyz(3) * mol1 % sintheta( 1 )

      this % points( n_trans ) % grid_xyz_rot(3) = this % points( n_trans ) % grid_xyz(2) * mol1 % sintheta( 1 ) + &
                                                   this % points( n_trans ) % grid_xyz(3) * mol1 % costheta( 1 )

      this % points( n_trans ) % grid_xyz(:)     = this % points( n_trans ) % grid_xyz_rot(:)

      !! DONE: ROTATION VECTOR IS AT ORIGIN !! PRECESSION MOVES AROUND z AXIS

      this % points( n_trans ) % grid_xyz_rot(1) = this % points( n_trans ) % grid_xyz(1) * cosalpha - &
                                                   this % points( n_trans ) % grid_xyz(2) * sinalpha

      this % points( n_trans ) % grid_xyz_rot(2) = this % points( n_trans ) % grid_xyz(1) * sinalpha + &
                                                   this % points( n_trans ) % grid_xyz(2) * cosalpha

      this % points( n_trans ) % grid_xyz_rot(3) = this % points( n_trans ) % grid_xyz(3)

      this % points( n_trans ) % grid_xyz(:)     = this % points( n_trans ) % grid_xyz_rot(:)

      !! ROTATION MOVE TO rth POINT OF THE SPHERICAL GRID !! ROTATION AROUND y AXIS

      this % points( n_trans ) % grid_xyz_rot(1) =  this % points( n_trans ) % grid_xyz(1) * costhetar + &
                                                    this % points( n_trans ) % grid_xyz(3) * sinthetar

      this % points( n_trans ) % grid_xyz_rot(2) =  this % points( n_trans ) % grid_xyz(2)

      this % points( n_trans ) % grid_xyz_rot(3) = -this % points( n_trans ) % grid_xyz(1) * sinthetar + &
                                                    this % points( n_trans ) % grid_xyz(3) * costhetar

      this % points( n_trans ) % grid_xyz(:)     = this % points( n_trans ) % grid_xyz_rot(:)

      !! ROTATION AROUND z AXIS TO PLACE ROTATION VECTOR AT THE rth SPHERE POINT !!

      this % points( n_trans ) % grid_xyz_rot(1) = this % points( n_trans ) % grid_xyz(1) * cosphir - &
                                                   this % points( n_trans ) % grid_xyz(2) * sinphir

      this % points( n_trans ) % grid_xyz_rot(2) = this % points( n_trans ) % grid_xyz(1) * sinphir + &
                                                   this % points( n_trans ) % grid_xyz(2) * cosphir

      this % points( n_trans ) % grid_xyz_rot(3) = this % points( n_trans ) % grid_xyz(3) 

      this % points( n_trans ) % grid_xyz(:)     = this % points( n_trans ) % grid_xyz_rot(:)

    enddo

    return
  end subroutine Rotate_grid

  !---------------------------------------------------------------------------
  !> @brief This routine generates the translation sphere by tesselation.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Build_translation_sphere( this, trans_factor, rad )
    use mod_spherical_grids
    use mod_error_handling

    implicit none

    class( grid ), intent( inout )                   :: this
    integer, intent( in )                            :: trans_factor
    real( kind = DP ), intent( in )                  :: rad
    real( kind = DP ), allocatable, dimension(:,:)   :: node_xyz
    integer                                          :: node_num
    integer                                          :: ierr
    type(error)                                      :: err

    ! if parameter is zero, do not perform molecular translations, just place it at (0,0,rad)!
    
    if ( trans_factor == 0 ) then

      node_num = 1

      this % numpoint = node_num

      allocate( this % points( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      this % points(:) % grid_xyz(1) = 0.0_DP 
      this % points(:) % grid_xyz(2) = 0.0_DP 
      this % points(:) % grid_xyz(3) = 1.0_DP * rad

    else

      node_num = 12 + 10 * 3 * ( trans_factor - 1 ) + 10 * ( trans_factor - 2 ) * ( trans_factor - 1 )
      !edge_num = 30 * factor * factor
      !face_num = 20 * factor * factor

      this % numpoint = node_num

      allocate( node_xyz(3,node_num), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      allocate( this % points( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      call sphere_icos1_points ( trans_factor, node_num, node_xyz )

      this % points(:) % grid_xyz(1) = node_xyz(1,:) * rad
      this % points(:) % grid_xyz(2) = node_xyz(2,:) * rad
      this % points(:) % grid_xyz(3) = node_xyz(3,:) * rad

      if ( allocated(node_xyz) ) deallocate(node_xyz)
      
    endif

    return
  end subroutine Build_translation_sphere

  subroutine Build_rotation1_sphere( this, point_rot_factor )
    use mod_read_molecules
    use mod_spherical_grids
    use mod_error_handling
      
    implicit none

    class( grid ), intent( inout )                   :: this
    integer, intent( in )                            :: point_rot_factor
    integer                                          :: n_rot1
    integer                                          :: node_num
    real( kind = DP ), allocatable, dimension(:,:)   :: node_xyz
    integer                                          :: ierr
    type(error)                                      :: err

    ! if parameter is zero, do not perform molecular reorientations, just align it along Z !

    if ( point_rot_factor == 0 ) then

      node_num = 1

      this % numpoint = node_num

      allocate( this % points( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
        
      allocate(     rr( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
        
      allocate(   phir( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
        
      allocate( thetar( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      phir(1)   = PI / 2.0_DP
      thetar(1) = PI

    else

      node_num = 12 + 10 * 3 * ( point_rot_factor - 1 ) + 10 * ( point_rot_factor - 2 ) * ( point_rot_factor - 1 )
      ! edge_num = 30 * factor * factor
      ! face_num = 20 * factor * factor

      this % numpoint = node_num

      allocate( node_xyz(3,node_num), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
        
      allocate( this % points( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      call sphere_icos1_points ( point_rot_factor, node_num, node_xyz )

      this % points(:) % grid_xyz(1) = node_xyz(1,:) 
      this % points(:) % grid_xyz(2) = node_xyz(2,:)
      this % points(:) % grid_xyz(3) = node_xyz(3,:)

      if ( allocated(node_xyz) ) deallocate(node_xyz)

      allocate( rr( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      allocate( phir( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )
        
      allocate( thetar( node_num ), stat=ierr )
      if ( ierr /= 0 ) call err % error( 'e', message = "abnormal memory allocation" )

      !! CONVERTS ROTATION GRID CARTESIAN COORDS DO SPHERICAL COORDS

      do n_rot1 = 1, node_num

        rr( n_rot1 ) = dsqrt( this % points( n_rot1 ) % grid_xyz(1)**2 &
                            + this % points( n_rot1 ) % grid_xyz(2)**2 &
                            + this % points( n_rot1 ) % grid_xyz(3)**2 )

        phir( n_rot1 )   = datan2( this % points( n_rot1 ) % grid_xyz(2) , this % points( n_rot1 ) % grid_xyz(1) )
        thetar( n_rot1 ) = dacos( this % points( n_rot1 ) % grid_xyz(3) / rr( n_rot1 ) )

      enddo

    endif

    return
  end subroutine Build_rotation1_sphere

  !---------------------------------------------------------------------------
  !> @brief This routine prints some checking information on screen.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Check_moves
    use mod_constants  , only : dashline
    use mod_input_read , only : nconf2, axis_rot_moves, axis_rot_range

    implicit none

    character( len = 10 ) :: var_nt, var_nconf, var_nr1, var_nr2, var_max_r2

    write(var_nt,    '(i10)' ) grid_trans % numpoint
    write(var_nconf, '(i10)' ) nconf2
    write(var_nr1,   '(i10)' ) grid_rot1 % numpoint
    write(var_nr2,   '(i10)' ) axis_rot_moves
    write(var_max_r2,'(f6.2)') axis_rot_range

    conf_total = nconf2 * grid_rot1 % numpoint * axis_rot_moves
    ntotal     = conf_total * grid_trans % numpoint 

    write(output_unit,'(/, T3, A)') dashline      
    write(output_unit,'(T5, "Phase space discretization"                     )') 
    write(output_unit,'(T3, A)') dashline      
    write(output_unit,'(/, T5, "Translation grid points:", T91, A)') trim(var_nt)
    write(output_unit,'(/, T5, "Molecule 2 conformations:", T91, A)') trim(var_nconf)
    write(output_unit,'(/, T5, "Reorientation grid points:", T91, A)') trim(var_nr1)
    write(output_unit,'(/, T5, "Rotation moves around mol2 axis (from 0 to ", A, " degree):", T91,  A)') &
             & trim(var_max_r2), trim(var_nr2)
    write(output_unit,'(/, T9, "Total number of rotations:", T91, i10)') conf_total
    write(output_unit,'(/, T5, "TOTAL NUMBER OF CONFIGURATIONS:", T91, i10 )') ntotal
    write(output_unit,'(/, T3, A)') dashline

    return
  end subroutine Check_moves

end module mod_grids
