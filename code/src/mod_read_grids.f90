!------------------------------------------------------------------------------
! THEMIS: A program for 
! Copyright (C) 2017 Felippe M. Colombari
!------------------------------------------------------------------------------
!> @brief This module reads coordinate files for translation and rotation grids
!> @author Felippe M. Colombari
!> - Laboratório de Química Teórica, LQT -- UFSCar
!> @date - Jan, 2018 
!> - independent module created
!> @date - Jan, 2018 
!> - documentation and support added
!> @note
!> - grids are standard XYZ coordinate files.
!> @note added error_handling
!> Asdrubal Lozada-Blanco
!> @date - Nov 2019
!------------------------------------------------------------------------------

module MOD_READ_GRIDS
  use MOD_CONSTANTS, only: DP

  implicit none

  type point
    real( kind = DP )                              :: grid_xyz(3), grid_xyz_old(3), grid_xyz_rot(3)
    character( len = 3 )                           :: grid_symbol
  end type point

  type grid
    type( point ), allocatable, dimension(:)       :: points
    integer                                        :: numpoint
    contains
      procedure, pass                              :: Read_grid  
      procedure, pass                              :: Translate_grid
      procedure, pass                              :: Align_grid
      procedure, pass                              :: Rotate_grid
      procedure, pass                              :: Build_translation_sphere
      procedure, pass                              :: Build_reorientation_sphere
  end type grid

  type( grid )                                     :: grid_trans, grid_reo

  real( kind = DP ), allocatable, dimension(:)     :: A, mTS, Eavg
  real( kind = DP ), allocatable, dimension(:)     :: Zrot, sumVexpVrot, probT
  real( kind = DP ), allocatable, dimension(:)     :: rr, thetar, phir
  integer                                          :: ntotal, rot_total

  contains

    !---------------------------------------------------------------------------
    !> @brief This routine reads coordinates for user-defined translation grid.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Jun, 2017
    !> - first version 
    !> @date - Dec, 2017
    !> - modular subroutine  created
    !> @date - Jan, 2018
    !> - final version, documentation and support 
    !---------------------------------------------------------------------------	

    subroutine Read_grid( this, grid_filename )
      use mod_inquire, only: Inquire_file
      use mod_error_handling

      implicit none

      class( grid ), intent(inout) :: this
      character( len = * ), intent(in) :: grid_filename
      integer                          :: t
      integer                          :: ios         = 0
      integer                          :: file_unit   = 10        
      character( len = 15 )            :: file_format = "formatted"
      character( len = 15 )            :: file_access = "sequential"

      integer                          :: ierr

      type(error) :: err

      call Inquire_file( file_unit, grid_filename, file_format, file_access )

      read(file_unit,*,iostat=ios) this % numpoint
      read(file_unit,*)

      if ( allocated ( this % points ) ) deallocate ( this % points )
      allocate( this % points( this % numpoint ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      do t = 1, this % numpoint

        read(file_unit,*,iostat=ios) this % points(t) % grid_symbol, this % points(t) % grid_xyz(:)

      enddo

      close(file_unit)

      return
    end subroutine Read_grid

    !---------------------------------------------------------------------------
    !> @brief This routine centers user-defined translation grid altogether with mol1.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !---------------------------------------------------------------------------	

    subroutine Translate_grid( this )
      use MOD_READ_MOLECULES !, only: mol1

      implicit none

      class( grid ), intent(inout) :: this
      integer                          :: t

      do t = 1, this % numpoint

        this % points(t) % grid_xyz(:) = this % points(t) % grid_xyz(:) - mol1 % atom_ref1(:)

      enddo

      return
    end subroutine Translate_grid

    !---------------------------------------------------------------------------
    !> @brief This routine align user-defined translation grid along Z altogether with mol1.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !---------------------------------------------------------------------------	

    subroutine Align_grid( this ) 
      use MOD_READ_MOLECULES , only: sinphi, cosphi

      implicit none

      class( grid ), intent(inout) :: this
      integer                          :: t


      do t = 1, this % numpoint

        this % points(t) % grid_xyz_rot(1) = this % points(t) % grid_xyz(1) * cosphi - &
                                             this % points(t) % grid_xyz(2) * sinphi

        this % points(t) % grid_xyz_rot(2) = this % points(t) % grid_xyz(1) * sinphi + &
                                             this % points(t) % grid_xyz(2) * cosphi

        this % points(t) % grid_xyz_rot(3) = this % points(t) % grid_xyz(3)

        this % points(t) % grid_xyz(:)     = this % points(t) % grid_xyz_rot(:)

        this % points(t) % grid_xyz_old(:) = this % points(t) % grid_xyz(:)

      enddo

      return
    end subroutine Align_grid

    subroutine Rotate_grid( this ) 
      use MOD_READ_MOLECULES  , only: costheta, sintheta, cosalpha, sinalpha, costhetar, sinthetar, cosphir, sinphir

      implicit none

      class( grid ), intent(inout) :: this
      integer                      :: t

      do t = 1, this % numpoint
   
        !! ROTATE AROUND x AXIS: PLACE ROTATION VECTOR AT z AXIS

        this % points(t) % grid_xyz_rot(1) = this % points(t) % grid_xyz(1)

        this % points(t) % grid_xyz_rot(2) = this % points(t) % grid_xyz(2) * costheta - &
                                             this % points(t) % grid_xyz(3) * sintheta

        this % points(t) % grid_xyz_rot(3) = this % points(t) % grid_xyz(2) * sintheta + &
                                             this % points(t) % grid_xyz(3) * costheta

        this % points(t) % grid_xyz(:)     = this % points(t) % grid_xyz_rot(:)

        !! DONE: ROTATION VECTOR IS AT ORIGIN !! PRECESSION MOVES AROUND z AXIS

        this % points(t) % grid_xyz_rot(1) = this % points(t) % grid_xyz(1) * cosalpha - &
                                             this % points(t) % grid_xyz(2) * sinalpha

        this % points(t) % grid_xyz_rot(2) = this % points(t) % grid_xyz(1) * sinalpha + &
                                             this % points(t) % grid_xyz(2) * cosalpha

        this % points(t) % grid_xyz_rot(3) = this % points(t) % grid_xyz(3)

        this % points(t) % grid_xyz(:)     = this % points(t) % grid_xyz_rot(:)

        !! ROTATION MOVE TO rth POINT OF THE SPHERICAL GRID !! ROTATION AROUND y AXIS

        this % points(t) % grid_xyz_rot(1) =  this % points(t) % grid_xyz(1) * costhetar + &
                                              this % points(t) % grid_xyz(3) * sinthetar

        this % points(t) % grid_xyz_rot(2) =  this % points(t) % grid_xyz(2)

        this % points(t) % grid_xyz_rot(3) = -this % points(t) % grid_xyz(1) * sinthetar + &
                                              this % points(t) % grid_xyz(3) * costhetar

        this % points(t) % grid_xyz(:)     = this % points(t) % grid_xyz_rot(:)

        !! ROTATION AROUND z AXIS TO PLACE ROTATION VECTOR AT THE rth SPHERE POINT !!

        this % points(t) % grid_xyz_rot(1) = this % points(t) % grid_xyz(1) * cosphir - &
                                             this % points(t) % grid_xyz(2) * sinphir

        this % points(t) % grid_xyz_rot(2) = this % points(t) % grid_xyz(1) * sinphir + &
                                             this % points(t) % grid_xyz(2) * cosphir

        this % points(t) % grid_xyz_rot(3) = this % points(t) % grid_xyz(3) 

        this % points(t) % grid_xyz(:)     = this % points(t) % grid_xyz_rot(:)

      enddo

      return
    end subroutine Rotate_grid

    subroutine Build_translation_sphere( this, trans_factor, rad )
      use MOD_SPHERICAL_GRIDS
      use mod_error_handling

      implicit none

      class( grid ), intent(inout)                   :: this
      integer, intent(IN)                            :: trans_factor
      real( kind = DP ), intent(IN)                  :: rad
      integer                                        :: node_num
      real( kind = DP ), allocatable, dimension(:,:) :: node_xyz
      integer                                        :: ierr
      type(error)                                    :: err

      node_num = 12 + 10 * 3 * ( trans_factor - 1 ) + 10 * ( trans_factor - 2 ) * ( trans_factor - 1 )
!      edge_num = 30 * factor * factor
!      face_num = 20 * factor * factor

      this % numpoint = node_num

      allocate( node_xyz(3,node_num), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      allocate( this % points( node_num ), stat=ierr )
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      call sphere_icos1_points ( trans_factor, node_num, node_xyz )

      this % points(:) % grid_xyz(1) = node_xyz(1,:) * rad
      this % points(:) % grid_xyz(2) = node_xyz(2,:) * rad
      this % points(:) % grid_xyz(3) = node_xyz(3,:) * rad

      if ( allocated(node_xyz) ) deallocate(node_xyz)
      
      return
    end subroutine Build_translation_sphere

    subroutine Build_reorientation_sphere( this, reo_factor )
      use MOD_READ_MOLECULES
      use MOD_SPHERICAL_GRIDS
      use mod_error_handling
      
      implicit none

      class( grid ), intent(inout)                   :: this
      integer, intent(IN)                            :: reo_factor
      integer                                        :: r
      integer                                        :: node_num
      real( kind = DP ), allocatable, dimension(:,:) :: node_xyz
      integer                                        :: ierr
      type(error)                                    :: err

      ! if parameter is zero, do not perform molecular reorientations, just align it along Z !

      if ( reo_factor == 0 ) then

        node_num = 1

        this % numpoint = node_num

        allocate( this % points( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")
        
        allocate(     rr( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")
        
        allocate(   phir( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")
        
        allocate( thetar( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")

        phir(1) = PI / 2
        thetar(1) = PI

      else

        node_num = 12 + 10 * 3 * ( reo_factor - 1 ) + 10 * ( reo_factor - 2 ) * ( reo_factor - 1 )
        ! edge_num = 30 * factor * factor
        ! face_num = 20 * factor * factor

        this % numpoint = node_num

        allocate( node_xyz(3,node_num), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")
        
        allocate( this % points( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")

        call sphere_icos1_points ( reo_factor, node_num, node_xyz )

        this % points(:) % grid_xyz(1) = node_xyz(1,:) 
        this % points(:) % grid_xyz(2) = node_xyz(2,:)
        this % points(:) % grid_xyz(3) = node_xyz(3,:)

        if ( allocated(node_xyz) ) deallocate(node_xyz)

        allocate(     rr( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")

        allocate(   phir( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")
        
        allocate( thetar( node_num ), stat=ierr )
        if(ierr/=0) call err%error('e',message="abnormal memory allocation")

        !! CONVERTS ROTATION GRID CARTESIAN COORDS DO SPHERICAL COORDS

        do r = 1, node_num

          rr(r) = dsqrt( this % points(r) % grid_xyz(1)**2 &
                       + this % points(r) % grid_xyz(2)**2 &
                       + this % points(r) % grid_xyz(3)**2 )

          phir(r) = datan2( this % points(r) % grid_xyz(2) , this % points(r) % grid_xyz(1) )
          thetar(r) = dacos( this % points(r) % grid_xyz(3) / rr(r) )

        enddo

      endif

      return
    end subroutine Build_reorientation_sphere

    !---------------------------------------------------------------------------
    !> @brief This routine prints some checking information on screen.
    !> @author Felippe M. Colombari
    !> - Laboratório de Química Teórica, LQT -- UFSCar
    !> @date - Dec, 2017
    !> - subroutine  created and documented
    !---------------------------------------------------------------------------	
    
    subroutine CHECK_MOVES
      use MOD_CONSTANTS,  only: dashline
      use MOD_INPUT_READ, only: gyr_factor, max_gyr

      implicit none

      character( len = 10 ) :: var_nt, var_nr, var_ng, var_max_gyr

      write(var_nt,'(i10)') grid_trans % numpoint
      write(var_nr,'(i10)') grid_reo % numpoint
      write(var_ng,'(i10)') gyr_factor
      write(var_max_gyr,'(f6.2)') max_gyr

      rot_total = grid_reo % numpoint * gyr_factor
      ntotal    = rot_total * grid_trans % numpoint 

      write(*,'(/, T3, A)') dashline      
      write(*,'(T5, "Phase space discretization"                     )') 
      write(*,'(T3, A)') dashline      
      write(*,'(/, T5, "Translation grid points:", T91, A)') trim(var_nt)
      write(*,'(/, T5, "Reorientation grid points:", T91, A)') trim(var_nr)
      write(*,'(/, T5, "Gyration points around mol2 axis (from 0 to ", A, " degree):", T91,  A)') &
               & trim(var_max_gyr), trim(var_ng)
      write(*,'(/, T9, "Total number of rotations:", T91, i10)') rot_total
      write(*,'(/, T5, "TOTAL NUMBER OF CONFIGURATIONS:", T91, i10 )') ntotal
      write(*,'(/, T3, A)') dashline

      return
    end subroutine CHECK_MOVES

end module MOD_READ_GRIDS
