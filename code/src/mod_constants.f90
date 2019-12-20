!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief This module defines a set constants.
!> @author Felippe M. Colombari                                                
!> - Laboratório de Química Teórica, LQT -- UFSCar                             
!> @date - Dec, 2017                                                           
!> - module created                                                
!------------------------------------------------------------------------------

module mod_constants

  implicit none

  integer, public, parameter           :: DP = selected_real_kind(15, 307) !  double precision constant for portability
  integer, public, parameter           :: SP = selected_real_kind(6, 37)   !  single precision constant for portability
  
  real( kind = DP ), public, parameter :: PI      = 3.14159265358979_DP    !         PI constant with 14 decimal places
  real( kind = DP ), public, parameter :: DEG2RAD = 180.0_DP / PI          !                         degrees to radians
  real( kind = DP ), public, parameter :: kB      = 0.0083144621_DP        !             boltzmann constant in kJ/mol/K
  real( kind = DP ), public, parameter :: CCON    = 1389.354578_DP         ! coulombic energy to kJ/mol in kJ.A/mol/e^2
  real( kind = DP ), public, parameter :: FPZERO  = tiny(1.0_DP)           !              define machine-precision ZERO
  real( kind = DP ), public, parameter :: FPINF   = huge(1.0_DP)           !          define machine-precision INFINITY
  real( kind = DP ), public, parameter :: ms      = 0.001_DP               !                            convert ms to s
  
  character( len = 11 ), parameter     :: int_alphabet   = '1234567890'    !             allowed character for integers
  character( len = 12 ), parameter     :: float_alphabet = '.-1234567890'  !               allowed character for floats
  character( len = 66 ), parameter     :: char_alphabet  = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ._-1234567890 '   !              allowed character for strings
  character( len = 100 ), parameter     :: dashline = repeat('-',100)        !                            just a dashline

end module mod_constants
