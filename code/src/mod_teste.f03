MODULE mod_teste

  USE MOD_CONSTANTS
  USE MOD_POT_LJC
  USE MOD_POT_BHC

  implicit none

  integer :: n

  interface

    subroutine calc_LJC_energy(p, r, t)
    end subroutine calc_LJC_energy(p, r, t)
    
    subroutine calc_BHC_energy(p, r, t)
    end subroutine calc_BHC_energy(p, r, t)

  end interface

  procedure(forced), pointer:: calc_energy => NULL()

  select case( potential )

    case( "lj-coul" )
      calc_energy => ljc_dimer % calc_LJC_energy(p, r, t)

    case( "bh-coul" )
      calc_energy => bhc_dimer % calc_BHC_energy(p, r, t)

  end select

  call calc_energy

  end module
