!---------------------------------------------------------------------------------------------------
! SAS_GRID: A code to obtain the solvent accessible surface (SAS) around a given molecular structure                                                  
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
!> @file   mod_grids.f90
!> @author Felippe M. Colombari
!> @email  colombarifm@hotmail.com
!> @brief  This module manipulates the coordinates of translation and rotation grids
!> @date - Nov, 2019                                                           
!> - independent module created                                                
!> @date - Nov, 2019                                                           
!> - spherical grids generation subroutine added
!> @date - Nov, 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_grids
  use mod_constants , only : DP

  implicit none

  type point
    real( kind = DP )                              :: grid_xyz(3)
    character( len = 3 )                           :: grid_symbol
  end type point

  type sphere_grid
    type( point ), allocatable, dimension(:)       :: points_sphere
    integer                                        :: numpoints_sphere
  contains
    procedure, pass                                :: Build_sphere_grid
  end type sphere_grid

  type( sphere_grid )                              :: grid_sphere

  type sas_grid
    integer                                        :: numpoints_sas
  contains
    procedure, pass                                :: Build_sas_grid
  end type sas_grid

  type( sas_grid )                              :: grid_sas

  real(kind=dp)                 :: x1_sphere, y1_sphere, z1_sphere
  real(kind=dp)                 :: x2_sphere, y2_sphere, z2_sphere
  real(kind=dp)                 :: dx, dy, dz, r_sqr

  logical, allocatable, dimension(:,:)    :: check

contains

  !---------------------------------------------------------------------------
  !> @brief This routine generates the translation sphere by tesselation.
  !> @author Felippe M. Colombari
  !---------------------------------------------------------------------------	
  subroutine Build_sphere_grid( this, tessellation_factor, sphere_radius )
    use mod_spherical_grids
    use mod_error_handling

    implicit none

    class( sphere_grid ), intent(inout)            :: this
    integer, intent(IN)                            :: tessellation_factor
    real( kind = DP ), intent(IN)                  :: sphere_radius
    integer                                        :: node_num
    real( kind = DP ), allocatable, dimension(:,:) :: node_xyz
    integer                                        :: ierr, i
    type(error)                                    :: err

    node_num = 12 + 10 * 3 * ( tessellation_factor - 1 ) + 10 * ( tessellation_factor - 2 ) * ( tessellation_factor - 1 )
    !edge_num = 30 * factor * factor
    !face_num = 20 * factor * factor

    this % numpoints_sphere = node_num

    allocate( node_xyz(3,node_num), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    allocate( this % points_sphere( node_num ), stat=ierr )
    if(ierr/=0) call err%error('e',message="abnormal memory allocation")

    call sphere_icos1_points ( tessellation_factor, node_num, node_xyz )

    this % points_sphere(:) % grid_xyz(1) = node_xyz(1,:)
    this % points_sphere(:) % grid_xyz(2) = node_xyz(2,:)
    this % points_sphere(:) % grid_xyz(3) = node_xyz(3,:)

!    do i = 1, this % numpoints_sphere

!      write(666,'(a2,3(2x,f12.6))') 'XX', this % points_sphere(i) % grid_xyz(:)

!    enddo

    if ( allocated(node_xyz) ) deallocate(node_xyz)
      
    return
  end subroutine Build_sphere_grid

  subroutine Build_sas_grid( this )
    use mod_read_molecule
    use mod_error_handling

    implicit none

    class( sas_grid ), intent(inout)              :: this

    character( len = 4 )                          :: char_numpoints_sphere

    integer                                       :: i, j, k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                                                  !!
    !!    READ SPHERICAL GRID THAT WILL BE PLACED AT    !!
    !!  EVERY ATOMIC CENTER FOR THE REFERENCE MOLECULE. !!
    !!  THEN PERFORM SOME DISTANCE CHECKS TO BUILD THE  !!
    !!        SAS GRID AROUND REFERENCE MOLECULE.       !!
    !!                                                  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    allocate( check ( mol % num_atoms, grid_sphere % numpoints_sphere ) )

    check = .true.

    do i = 1, mol % num_atoms

      do j = 1, grid_sphere % numpoints_sphere

        x1_sphere = mol % atoms(i) % xyz(1) + &
                  & mol % atoms(i) % vdw * &
                  & grid_sphere % points_sphere(j) % grid_xyz(1)

        y1_sphere = mol % atoms(i) % xyz(2) + &
                  & mol % atoms(i) % vdw * &
                  & grid_sphere % points_sphere(j) % grid_xyz(2)

        z1_sphere = mol % atoms(i) % xyz(3) + &
                  & mol % atoms(i) % vdw * &
                  & grid_sphere % points_sphere(j) % grid_xyz(3)

        do k = 1, mol % num_atoms

          if ( i /= k ) then

            x2_sphere = mol % atoms(k) % xyz(1)

            y2_sphere = mol % atoms(k) % xyz(2)
          
            z2_sphere = mol % atoms(k) % xyz(3)
            
            dx = x2_sphere - x1_sphere
            
            dy = y2_sphere - y1_sphere
            
            dz = z2_sphere - z1_sphere
            
            r_sqr = dsqrt( dx*dx + dy*dy + dz*dz )
            
            if ( r_sqr < mol % atoms(k) % vdw ) then
            
              check(i,j) = .false.
              
            endif
            
          endif
        
        enddo
        
      enddo
      
    enddo

    grid_sas % numpoints_sas = 0

    do i = 1, mol % num_atoms
    
      do j = 1, grid_sphere % numpoints_sphere
          
        if ( check(i,j) .eqv. .true. ) then
        
          this % numpoints_sas = this % numpoints_sas + 1
          
        endif
        
      enddo
      
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                                                  !!
    !!      WRITE RESULTING SAS GRID TO A XYZ FILE.     !!
    !!                                                  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(char_numpoints_sphere,'(i4.4)') grid_sphere % numpoints_sphere

    open( unit = 40, file = 'sas_'//char_numpoints_sphere//'.xyz', status = 'unknown' )
      
    write (40,*) this % numpoints_sas
    
    write (40,'(a)')
      
    do i = 1, mol % num_atoms

      do j = 1, grid_sphere % numpoints_sphere
      
        if ( check(i,j) .eqv. .true. ) then
        
          x1_sphere = mol % atoms(i) % xyz(1) + &
                    & mol % atoms(i) % vdw * &
                    & grid_sphere % points_sphere(j) % grid_xyz(1)

          y1_sphere = mol % atoms(i) % xyz(2) + &
                    & mol % atoms(i) % vdw * &
                    & grid_sphere % points_sphere(j) % grid_xyz(2)

          z1_sphere = mol % atoms(i) % xyz(3) + &
                    & mol % atoms(i) % vdw * &
                    & grid_sphere % points_sphere(j) % grid_xyz(3)

          write (40,'(a2,3(3x,f12.8))') 'XX', x1_sphere, y1_sphere, z1_sphere
          
        endif
        
      enddo
      
    enddo
      
    close(40)
     
  end subroutine Build_sas_grid 

end module mod_grids
