!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2024 Themis developers
!
!   This file was written by Felippe M. Colombari.
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
!> @file   mod_constants.f90
!> @author Felippe M. Colombari
!> @brief  Defines a set constants.
!> @date - Dec, 2017                                                           
!> - independent module created                                                
!> @date - Jan, 2023                                                           
!> - add function to convert uppercase strings to lowercase strings                   
!---------------------------------------------------------------------------------------------------

module mod_constants

  implicit none

  character( len = 8 ), parameter      :: version = '3.0.0'
  character( len = 8 ), parameter      :: revision = 'ae2e8ac'

  integer, public, parameter           :: DP = selected_real_kind(15, 307) !  double precision constant for portability
  integer, public, parameter           :: SP = selected_real_kind(6, 37)   !  single precision constant for portability
  
  real( kind = DP ), public, parameter :: PI      = 3.14159265358979_DP    !         PI constant with 14 decimal places
  real( kind = DP ), public, parameter :: DEG2RAD = 180.0_DP / PI          !                         degrees to radians
  real( kind = DP ), public, parameter :: KB      = 0.0083144621_DP        !             boltzmann constant in kJ/mol/K
  real( kind = DP ), public, parameter :: CCON    = 1389.354578_DP         ! coulombic energy to kJ/mol in kJ.A/mol/e^2
  real( kind = DP ), public, parameter :: FPZERO  = tiny(1.0_DP)           !              define machine-precision ZERO
  real( kind = DP ), public, parameter :: FPINF   = huge(1.0_DP)           !          define machine-precision INFINITY
  real( kind = DP ), public, parameter :: MS      = 0.001_DP               !                            convert ms to s
  
  character( len = 100 ), public, parameter  :: DASHLINE = repeat('-',100)      !                       just a dashline
  character( len = 11 ), public, parameter   :: INT_ALPHABET   = '1234567890'   !        allowed character for integers
  character( len = 12 ), public, parameter   :: FLOAT_ALPHABET = '.-1234567890' !          allowed character for floats
  character( len = 66 ), public, parameter   :: CHAR_ALPHABET  = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ._-1234567890 '        !         allowed character for strings

  character( len = 26 ), public, parameter           :: UPPERCASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character( len = 26 ), public, parameter           :: LOWERCASE = 'abcdefghijklmnopqrstuvwxyz'

contains

  Pure Function to_lower (input_string) Result (output_string)

    implicit none
    character( len = * ), intent(IN)          :: input_string
    character( len = len(input_string) )      :: output_string

    integer :: char_index, char_counter

    output_string = input_string
    
    do char_counter = 1, len_trim(input_string)
      char_index = index(UPPERCASE, input_string(char_counter:char_counter))
      if (char_index > 0) output_string(char_counter:char_counter) = LOWERCASE(char_index:char_index)
    end do

  end function to_lower

end module mod_constants
