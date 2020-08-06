!---------------------------------------------------------------------------------------------------
! COM: A code to calculate the center of mass of a given molecular structure                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2020 Themis developers
!                 Laboratory of Theoretical Chemistry (LQT) - Federal University of São Carlos 
!                 <http://www.lqt.dq.ufscar.br>
!
!   Please cite: 
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
!> @file   com.f90
!> @author Felippe M. Colombari
!> @brief  Defines a set constants.
!> @date - Jan, 2020                                                           
!> - independent module created                                                
!> @note 
!> - masses from http://www.chemicalelements.com/show/mass.html
!---------------------------------------------------------------------------------------------------

module mod_constants

  implicit none

  integer, public, parameter           :: DP = selected_real_kind(15, 307) !  double precision constant for portability

  real( kind = DP ), parameter          :: mass_H  =   1.007940_DP
  real( kind = DP ), parameter          :: mass_He =   4.002602_DP
  real( kind = DP ), parameter          :: mass_Li =   6.941000_DP
  real( kind = DP ), parameter          :: mass_Be =   9.012182_DP
  real( kind = DP ), parameter          :: mass_B  =  10.811000_DP
  real( kind = DP ), parameter          :: mass_C  =  12.010700_DP
  real( kind = DP ), parameter          :: mass_N  =  14.006740_DP
  real( kind = DP ), parameter          :: mass_O  =  15.999400_DP
  real( kind = DP ), parameter          :: mass_F  =  18.998400_DP
  real( kind = DP ), parameter          :: mass_Ne =  20.179700_DP
  real( kind = DP ), parameter          :: mass_Na =  22.989770_DP
  real( kind = DP ), parameter          :: mass_Mg =  24.305000_DP
  real( kind = DP ), parameter          :: mass_Al =  26.981539_DP
  real( kind = DP ), parameter          :: mass_Si =  28.085500_DP
  real( kind = DP ), parameter          :: mass_P  =  30.973760_DP
  real( kind = DP ), parameter          :: mass_S  =  32.066600_DP
  real( kind = DP ), parameter          :: mass_Cl =  35.452700_DP
  real( kind = DP ), parameter          :: mass_Ar =  39.948000_DP
  real( kind = DP ), parameter          :: mass_K  =  39.098300_DP
  real( kind = DP ), parameter          :: mass_Ca =  40.078000_DP
  real( kind = DP ), parameter          :: mass_Sc =  44.955910_DP
  real( kind = DP ), parameter          :: mass_Ti =  47.867000_DP
  real( kind = DP ), parameter          :: mass_V  =  50.941500_DP
  real( kind = DP ), parameter          :: mass_Cr =  51.996100_DP
  real( kind = DP ), parameter          :: mass_Mn =  54.938049_DP
  real( kind = DP ), parameter          :: mass_Fe =  55.845000_DP
  real( kind = DP ), parameter          :: mass_Co =  58.933200_DP
  real( kind = DP ), parameter          :: mass_Ni =  58.693400_DP
  real( kind = DP ), parameter          :: mass_Cu =  63.546000_DP
  real( kind = DP ), parameter          :: mass_Zn =  65.390000_DP
  real( kind = DP ), parameter          :: mass_Ga =  69.723000_DP
  real( kind = DP ), parameter          :: mass_Ge =  72.610000_DP
  real( kind = DP ), parameter          :: mass_As =  74.921600_DP
  real( kind = DP ), parameter          :: mass_Se =  78.960000_DP
  real( kind = DP ), parameter          :: mass_Br =  79.904000_DP
  real( kind = DP ), parameter          :: mass_Kr =  83.800000_DP
  real( kind = DP ), parameter          :: mass_Rb =  85.467800_DP
  real( kind = DP ), parameter          :: mass_Sr =  87.620000_DP
  real( kind = DP ), parameter          :: mass_Y  =  88.905850_DP
  real( kind = DP ), parameter          :: mass_Zr =  91.224000_DP
  real( kind = DP ), parameter          :: mass_Nb =  92.906380_DP
  real( kind = DP ), parameter          :: mass_Mo =  95.950000_DP
  real( kind = DP ), parameter          :: mass_Tc =  98.000000_DP
  real( kind = DP ), parameter          :: mass_Ru = 101.070000_DP
  real( kind = DP ), parameter          :: mass_Rh = 102.905500_DP
  real( kind = DP ), parameter          :: mass_Pd = 106.420000_DP
  real( kind = DP ), parameter          :: mass_Ag = 107.868200_DP
  real( kind = DP ), parameter          :: mass_Cd = 112.411000_DP
  real( kind = DP ), parameter          :: mass_In = 114.818000_DP
  real( kind = DP ), parameter          :: mass_Sn = 118.710000_DP
  real( kind = DP ), parameter          :: mass_Sb = 121.760000_DP
  real( kind = DP ), parameter          :: mass_Te = 127.600000_DP
  real( kind = DP ), parameter          :: mass_I  = 126.904470_DP
  real( kind = DP ), parameter          :: mass_Xe = 131.290000_DP
  real( kind = DP ), parameter          :: mass_Cs = 132.905450_DP
  real( kind = DP ), parameter          :: mass_Ba = 137.327000_DP
  real( kind = DP ), parameter          :: mass_Hf = 178.490000_DP
  real( kind = DP ), parameter          :: mass_Ta = 180.947900_DP
  real( kind = DP ), parameter          :: mass_W  = 183.840000_DP
  real( kind = DP ), parameter          :: mass_Re = 186.207000_DP
  real( kind = DP ), parameter          :: mass_Os = 190.230000_DP
  real( kind = DP ), parameter          :: mass_Ir = 192.217000_DP
  real( kind = DP ), parameter          :: mass_Pt = 195.084900_DP
  real( kind = DP ), parameter          :: mass_Au = 196.966550_DP
  real( kind = DP ), parameter          :: mass_Hg = 200.590000_DP
  real( kind = DP ), parameter          :: mass_Tl = 204.383300_DP
  real( kind = DP ), parameter          :: mass_Pb = 207.200000_DP
  real( kind = DP ), parameter          :: mass_Bi = 208.980380_DP

  character( len = 66 ), public, parameter   :: CHAR_ALPHABET  = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ._-1234567890 '   !              allowed character for strings
  character( len = 100 ), public, parameter  :: DASHLINE = repeat('-',100)      !                       just a dashline

end module mod_constants
