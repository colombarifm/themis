!---------------------------------------------------------------------------------------------------
!> @file   xdr.f90
!> @author James W. Barnett 
!>         https://github.com/wesbarnett/
!> @email  jbarnet4@tulane.edu 
!> @brief XDR Fortran Interface with Wrappers
!> @date - 2014                                                           
!> - module created                                                
!> @date - Jan, 2018                                                           
!> - .trr functions removed; read options removed (Felippe M. Colombari)
!---------------------------------------------------------------------------------------------------

module xdr
  use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_DOUBLE, C_INT
#ifdef XDR_DOUBLE
#define F_REAL KIND(0.0D0)
! TODO: use iso_fortran_env, only: REAL64
#define C_REAL C_DOUBLE
#else
#define F_REAL KIND(0.0)
! TODO: use iso_fortran_env, only: REAL_KINDS
#define C_REAL C_FLOAT
#endif
  implicit none
  private

  type, abstract :: trjfile
    type(xdrfile), pointer :: xd 
    integer(C_INT) :: natoms, step, stat
    character(len=1) :: mode
  contains
    procedure :: init => init_xdr
    procedure :: close => close_xdr
  end type

  ! *** xtcfile type
  ! box     - triclinic pbc box of the configuration
  ! natoms  - number of atoms in the configuration.
  ! pos     - positions read in (3, natoms)
  ! prec    - precision of the coordinates read in
  ! step    - step number of configuration.
  ! stat    - status of operation. 0 = good
  ! time    - time of the configuration
  ! xd      - pointer from libxdrfile.

  ! Should always call init first. Then call read in a loop and do your
  ! calculations. After the loops call close.
  type, extends(trjfile), public :: xtcfile
    real(C_FLOAT) :: box(3, 3), time
    real(C_FLOAT), allocatable :: pos(:, :)
    real(C_FLOAT) :: prec
  contains
    procedure :: write => write_xtcfile
  end type

  ! the data type located in libxdrfile
  type, bind(C) :: xdrfile
    type(C_PTR) :: fp, xdr
    character(kind=C_CHAR) :: mode
    type(C_PTR) :: buf1, buf2
    integer(C_INT) :: buf1size, buf2size
  end type xdrfile

  ! interface with libxdrfile
  interface 
    type(C_PTR) function xdrfile_open(filename, mode) bind(C, name='xdrfile_open')
      import
      character(kind=C_CHAR), intent(in) :: filename(*), mode(*)
    end function

    integer(C_INT) function xdrfile_close(xd) bind(C, name='xdrfile_close')
      import
      type(xdrfile), intent(in) :: xd
    end function

    integer(C_INT) function write_xtc(xd, natoms, step, time, box, x, prec) bind(C)
      import
      type(xdrfile), intent(in) :: xd
      integer(C_INT), intent(in), value :: natoms, step
      real(C_FLOAT), intent(in), value :: time, prec
      real(C_FLOAT), intent(in) :: box(*), x(*)
    end function

  end interface

contains
  ! our wrappers for the trjfile class
  subroutine  init_xdr(trj, filename_in, mode_opt)
    use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_CHAR, c_f_pointer
    implicit none
    class(trjfile), intent(inout) :: trj
    type(C_PTR) :: xd_c
    character (len=*), intent(in) :: filename_in
    character (len=206) :: filename
    character(len=1), optional, intent(in) :: mode_opt

    trj%mode = mode_opt 

    ! Set the file name to be read in for C.
    filename = trim(filename_in)//C_NULL_CHAR

    ! Open the file for reading or writing. Convert C pointer to Fortran pointer.
    xd_c = xdrfile_open(filename, trj%mode)
    call c_f_pointer(xd_c, trj%xd)
  end subroutine  init_xdr

  subroutine  write_xtcfile(xtc, natoms, step, time, box, pos, prec)
    implicit none
    class(xtcfile), intent(inout) :: xtc
    integer, intent(in) :: natoms, step
    real, intent(in) :: time, box(3, 3), pos(:, :), prec

    xtc%stat = write_xtc(xtc%xd, natoms, step, time, box, pos, prec)
  end subroutine  write_xtcfile

  subroutine  close_xdr(trj)
    implicit none
    class(trjfile), intent(inout) :: trj

    trj%stat = xdrfile_close(trj%xd)
  end subroutine  close_xdr

end module xdr
