!> @note added error handling
! Asdrubal Lozada-Blanco
!> @date - Nov 2019

module MOD_SPHERICAL_GRIDS
  use MOD_CONSTANTS, only: DP, PI

  implicit none

  contains

    subroutine ICOS_SHAPE(point_num,edge_num,face_num,face_order_max,point_coord,&
        edge_point,face_order,face_point )

      !*****************************************************************************80
      !
      !! ICOS_SHAPE describes an icosahedron.
      !
      !  Discussion:
      !
      !    The input data required for this routine can be retrieved from ICOS_SIZE.
      !
      !    The vertices lie on the unit sphere.
      !
      !    The dual of an icosahedron is a dodecahedron.
      !
      !    The data has been rearranged from a previous assignment.  
      !    The STRIPACK program refuses to triangulate data if the first
      !    three nodes are "collinear" on the sphere.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    22 July 2007
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) POINT_NUM, the number of points (12).
      !
      !    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges (30).
      !
      !    Input, integer ( kind = 4 ) FACE_NUM, the number of faces (20).
      !
      !    Input, integer ( kind = 4 ) FACE_ORDER_MAX, the maximum number of 
      !    vertices per face (3).
      !
      !    Output, real(kind=DP) POINT_COORD(3,POINT_NUM), the points.
      !
      !    Output, integer ( kind = 4 ) EDGE_POINT(2,EDGE_NUM), the points that 
      !    make up each edge, listed in ascending order of their indexes.
      !
      !    Output, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
      !    per face.
      !
      !    Output, integer ( kind = 4 ) FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
      !    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
      !    points are listed in the counter clockwise direction defined
      !    by the outward normal at the face.  The nodes of each face are ordered 
      !    so that the lowest index occurs first.  The faces are then sorted by
      !    nodes.
      !

      implicit none

      integer(kind=4)             ::  edge_num
      integer(kind=4), parameter  ::  edge_order = 2
      integer(kind=4)             ::  face_num
      integer(kind=4)             ::  face_order_max
      integer(kind=4)             ::  point_num

      real(kind=DP)               ::  a
      real(kind=DP)               ::  b
      integer(kind=4)             ::  edge_point(edge_order,edge_num)
      integer(kind=4)             ::  face_order(face_num)
      integer(kind=4)             ::  face_point(face_order_max,face_num)
      real(kind=DP)               ::  phi
      real(kind=DP)               ::  point_coord(3,point_num)
      real(kind=DP)               ::  z


      !  Set the point coordinates.

      phi = 0.5_DP * ( sqrt ( 5.0_DP ) + 1.0_DP )

      a = phi / sqrt ( 1.0_DP + phi * phi )
      b = 1.0_DP / sqrt ( 1.0_DP + phi * phi )
      z = 0.0_DP

      !  A*A + B*B + Z*Z = 1.

      point_coord(1:3,1:point_num) = reshape ( (/ &
        a,  b,  z, &
        a, -b,  z, &
        b,  z,  a, &
        b,  z, -a, &
        z,  a,  b, &
        z,  a, -b, &
        z, -a,  b, &
        z, -a, -b, &
       -b,  z,  a, &
       -b,  z, -a, &
       -a,  b,  z, &
       -a, -b,  z /), (/ 3, point_num /) )
  
      !  Set the edges.

      edge_point(1:edge_order,1:edge_num) = reshape ( (/ &
        1,  2, &
        1,  3, &
        1,  4, &
        1,  5, &
        1,  6, &
        2,  3, &
        2,  4, &
        2,  7, &
        2,  8, &
        3,  5, &
        3,  7, &
        3,  9, &
        4,  6, &
        4,  8, &
        4, 10, &
        5,  6, &
        5,  9, &
        5, 11, &
        6, 10, &
        6, 11, &
        7,  8, &
        7,  9, &
        7, 12, &
        8, 10, &
        8, 12, &
        9, 11, &
        9, 12, &
       10, 11, &
       10, 12, &
       11, 12 /), (/ edge_order, edge_num /) )

      !  Set the face orders.

      face_order(1:face_num) = (/ &
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)

      !  Set the faces.

      face_point(1:face_order_max,1:face_num) = reshape ( (/ &
        1,  2,  4, &
        1,  3,  2, &
        1,  4,  6, &
        1,  5,  3, &
        1,  6,  5, &
        2,  3,  7, &
        2,  7,  8, &
        2,  8,  4, &
        3,  5,  9, &
        3,  9,  7, &
        4,  8, 10, &
        4, 10,  6, &
        5,  6, 11, &
        5, 11,  9, &
        6, 10, 11, &
        7,  9, 12, &
        7, 12,  8, &
        8, 12, 10, &
        9, 11, 12, &
       10, 12, 11 /), (/ face_order_max, face_num /) )

      return

    end subroutine ICOS_SHAPE

    subroutine ICOS_NUM(point_num,edge_num,face_num,face_order_max)

      !*****************************************************************************80
      !
      !! ICOS_NUM gives "sizes" for an icosahedron in 3D.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    19 July 2007
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Output, integer ( kind = 4 ) POINT_NUM, the number of points.
      !
      !    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
      !
      !    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
      !
      !    Output, integer ( kind = 4 ) FACE_ORDER_MAX, the maximum order of any face.
      !
    
      implicit none

      integer(kind=4)       ::  edge_num
      integer(kind=4)       ::  face_num
      integer(kind=4)       ::  face_order_max
      integer(kind=4)       ::  point_num

      point_num = 12
      edge_num = 30
      face_num = 20
      face_order_max = 3
  
      return

    end subroutine ICOS_NUM

    FUNCTION R8VEC_NORM(n,a)

      !*****************************************************************************80
      !
      !! R8VEC_NORM returns the L2 norm of an R8VEC.
      !
      !  Discussion:
      !
      !    An R8VEC is a vector of R8's.
      !
      !    The vector L2 norm is defined as:
      !
      !      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    21 August 2010
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the number of entries in A.
      !
      !    Input, real(kind=DP) A(N), the vector whose L2 norm is desired.
      !
      !    Output, real(kind=DP) R8VEC_NORM, the L2 norm of A.
      !

      implicit none

      integer(kind=4)       ::  n

      real(kind=DP)         ::  a(n)
      real(kind=DP)         ::  r8vec_norm

      r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )
  
      return
    
    end FUNCTION R8VEC_NORM

    subroutine SPHERE_ICOS1_POINTS(factor,node_num,node_xyz)
       use mod_error_handling
      !*****************************************************************************80
      !
      !! SPHERE_ICOS1_POINTS returns icosahedral grid points on a sphere.
      !
      !  Discussion:
      !
      !    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
      !
      !    With FACTOR = 2, each triangle of the icosahedron is subdivided into
      !    2x2 subtriangles, resulting in 80 faces and 
      !    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
      !
      !    With FACTOR = 3, each triangle of the icosahedron is subdivided into
      !    3x3 subtriangles, resulting in 180 faces and 
      !    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
      !
      !    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
      !    resulting in 20 * FACTOR * FACTOR faces and
      !      12 
      !    + 20 * 3          * (FACTOR-1) / 2 
      !    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
      !
      !    This routine uses a simple, but only approximate, method of
      !    carrying out the subdivision.  For each spherical triangle of the
      !    face, we actually work in the planar triangle defined by those
      !    three points.  All subdivisions are done on that planar triangle,
      !    and the resulting points are then projected onto the sphere.
      !    While these points are equally spaced on the planar triangle,
      !    that is only approximately true on the sphere.
      !
      !    See SPHERE_ICOS2_POINTS for a more accurate method of subdivision
      !    that works on the spherical triangle itself.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    22 July 2007
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) FACTOR, the subdivision factor, which must
      !    be at least 1.
      !
      !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes, as reported
      !    by SPHERE_GRID_ICOS_SIZE.
      !
      !    Output, real(kind=DP) NODE_XYZ(3,NODE_NUM), the node coordinates.
      !
      !  Local Parameters:
      !
      !    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters 
      !    associated with the icosahedron, and POINT_COORD, EDGE_POINT, 
      !    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
      !    We need to refer to this data to generate the grid.
      !
      !    NODE counts the number of nodes we have generated so far.  At the
      !    end of the routine, it should be equal to NODE_NUM.
      !
      
      implicit none

      integer(kind=4)       ::  node_num
      integer(kind=4)       ::  a
      integer(kind=4)       ::  b
      integer(kind=4)       ::  c
      integer(kind=4)       ::  edge
      integer(kind=4)       ::  edge_num
      integer(kind=4), allocatable, dimension(:,:)  ::  edge_point
      integer(kind=4)       ::  f
      integer(kind=4)       ::  f1
      integer(kind=4)       ::  f2
      integer(kind=4)       ::  face
      integer(kind=4)       ::  face_num
      integer(kind=4), allocatable, dimension(:) ::   face_order
      integer(kind=4), allocatable, dimension ( :, : ) :: face_point
      integer(kind=4)       ::  face_order_max
      integer(kind=4)       ::  factor
      integer(kind=4)       ::  node
      real(kind=DP)         ::  node_norm
      real(kind=DP)         ::  node_xyz(3,node_num)
      real(kind=DP), allocatable, dimension(:,:)  ::  point_coord
      integer(kind=4)       ::  point_num

      integer               :: ierr
      type(error)           :: err

      !  Size the icosahedron.

      call icos_num ( point_num, edge_num, face_num, face_order_max )

      !  Set the icosahedron.

      allocate(point_coord(1:3,1:point_num),stat=ierr)
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")
      
      allocate(edge_point(1:2,1:edge_num),stat=ierr)  
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")
    
      allocate(face_order(1:face_num),stat=ierr)
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      allocate(face_point(1:face_order_max,1:face_num),stat=ierr)
      if(ierr/=0) call err%error('e',message="abnormal memory allocation")

      call icos_shape(point_num,edge_num,face_num,face_order_max,point_coord, &
        edge_point,face_order,face_point)

      !  Generate the point coordinates.
      
      !  A.  Points that are the icosahedral vertices.

      node = 0
      node_xyz(1:3,1:point_num) = point_coord(1:3,1:point_num)
  
      !  B. Points in the icosahedral edges, at 
      !  1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.

      node = 12

      do edge = 1, edge_num

        a = edge_point(1,edge)
        b = edge_point(2,edge)

        do f = 1, factor - 1

          node = node + 1

          node_xyz(1:3,node) = &
            ( real ( factor - f, kind = 8 ) * point_coord(1:3,a)   &
            + real (          f, kind = 8 ) * point_coord(1:3,b) ) &
            / real ( factor,     kind = 8 )

          node_norm = r8vec_norm ( 3, node_xyz(1:3,node) )

          node_xyz(1:3,node) = node_xyz(1:3,node) / node_norm

        end do

      end do

      !  C.  Points in the icosahedral faces.

      do face = 1, face_num

        a = face_point(1,face)
        b = face_point(2,face)
        c = face_point(3,face)

        do f1 = 1, factor - 1
      
          do f2 = 1, factor - f1 - 1
          
            node = node + 1

            node_xyz(1:3,node) = &
              ( real ( factor - f1 - f2, kind=DP ) * point_coord(1:3,a)   &
              + real (          f1,      kind=DP ) * point_coord(1:3,b)   &
              + real (               f2, kind=DP ) * point_coord(1:3,c) ) &
              / real ( factor,           kind=DP )

            node_norm = r8vec_norm ( 3, node_xyz(1:3,node) )

            node_xyz(1:3,node) = node_xyz(1:3,node) / node_norm

          end do
  
        end do

      end do

      if ( allocated( edge_point) ) deallocate(edge_point)
      if ( allocated( face_order) ) deallocate(face_order)
      if ( allocated( face_point) ) deallocate(face_point)
      if ( allocated( point_coord ) ) deallocate(point_coord)
  
      return

    end subroutine SPHERE_ICOS1_POINTS

end module MOD_SPHERICAL_GRIDS
