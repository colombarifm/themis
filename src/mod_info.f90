!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2023 Themis developers
!
!   This file was written by Asdrubal Lozada-Blanco and Felippe M. Colombari.
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
!> @file   mod_info.f90
!> @author Asdrubal Lozada-Blanco, Felippe Mariano Colombari
!> @brief  This module contains the program header
!> @date - Dec, 2019 
!> - module created and incorporated into code
!> @date - Jul, 2022 
!> - help, license and citation routines moved to this module
!---------------------------------------------------------------------------------------------------

module mod_info
  use iso_fortran_env , only : output_unit
  use mod_constants   , only : dashline, DP, version, revision
  ! TODO: include in makefile 
  !   include 'revision.inc'     
  implicit none
  private

  public :: Display_header, Display_date_time, Display_help, Display_license, Display_citation

contains

  subroutine Display_header()
    implicit none

    write(output_unit,'(T3, A)') dashline
    write(output_unit,'(T50, A, A)'  ) "THEMIS"
    write(output_unit,'(/,T7, A)'    ) " A software to assess association free energies via direct estimative of&
                                        &partition functions"
    write(output_unit,'(/,T40, A, A)') "Author: Felippe M. Colombari"
    write(output_unit,'(/,T36, A, A)') "Contributions: Asdrubal Lozada-Blanco"
    write(output_unit,'(/,T43, A, A)') "Program version: ",trim(version)
    write(output_unit,'(T46, A, A)'  ) "Revision ", revision
    write(output_unit,'(T3, A)') dashline

  end subroutine Display_header

  subroutine Display_date_time( string )
    implicit none

    character( len = * ), intent(IN) :: string
    integer,dimension(8)             :: values
    
    write(output_unit,'(T3, A)') dashline

    call Date_and_time( VALUES = values )
    write(output_unit,'(/, T5, A, i2.2, "/", i2.2, "/", i4, " - ", i2.2, ":", i2.2, ":", i2.2)') string,   &
                              &values(3), values(2), values(1), values(5), values(6), values(7)

    write(output_unit,'(/,T3, A)') dashline

  end subroutine Display_date_time
  
  subroutine Display_help
    implicit none

    write( output_unit, '(/,T3,A)' )  dashline
    write( output_unit, '(/,T36,A)' ) 'Usage:  themis [RUNTYPE] [GRID]     '
      
    write( output_unit, '(/,T3,A)' )   dashline
    write( output_unit, '(/,T30,A)' ) '[RUNTYPE] options'       
    write( output_unit, '(/,T8,A18,4x,A)' ) adjustr('--run'), &
                                        &adjustl('Start new calculation.')
    write( output_unit, '(/,T8,A18,4x,A)' ) adjustr('--rerun'), &
                                        &adjustl('Calculate properties from:') 
    write( output_unit, '(T10,A10,12x,A)' ) adjustr(''), &
                                        &adjustl('energy.bin - energies obtained from Themis run;')
    write( output_unit, '(T10,A10,12x,A)' ) adjustr(''), &
                                        &adjustl('energy.log - energies obtained from external calculation.')
    
    write( output_unit, '(/,T3,A)' )  dashline
    write( output_unit, '(/,T30,A)' ) '[GRID] options'       
    write( output_unit, '(/,T8,A18,4x,A)' ) adjustr('--shell <radius>'), &
                                        &adjustl('Translation moves will be performed on a spherical shell generated')
    write( output_unit, '(T8,A10,12x,A)' )  adjustr(''), &
                                        &adjustl('on the run. The real argument <radius> (in Angstrom) scales its radius.') 

    write( output_unit, '(/,T8,A18,4x,A)' ) adjustr('--user <file.xyz>'), &
                                        &adjustl('Translation moves will be performed on an user-defined grid read from')
    write( output_unit, '(T8,A10,12x,A)' )  adjustr(''), &
                                        &adjustl('<file.xyz> that must be aligned with molecule 1.')
      
    write( output_unit, '(/,T3,A,/)' ) dashline
    write( output_unit, '(T8,A18,4x,A)' ) adjustr('--help'), &
                                      &adjustl('Display this help')
    write( output_unit, '(T8,A18,4x,A)' ) adjustr('--license'), &
                                      &adjustl('Display the license')
    write( output_unit, '(T8,A18,4x,A)' ) adjustr('--citation'), &
                                      &adjustl('Display Themis papers')
    write( output_unit, '(/,T3,A)' ) dashline

    write( output_unit, '(/,T3, A)' ) '   Report bugs to:'
    write( output_unit, '(T10, A)' )  '   Felippe M. Colombari   - colombarifm@hotmail.com'

    write( output_unit, '(/,T3, A,/)' ) dashline
    
  end subroutine Display_help

  subroutine Display_license
    implicit none
   
    write( output_unit, '(/,T3,A,/)' ) dashline
    
    write( output_unit, '(T37,A)' ) 'Copyright 2021 Felippe M. Colombari'
    
    write( output_unit, * ) 
    write( output_unit, '(T33,A)' ) 'License GPLv3+: GNU GPL version 3 or later' 
    write( output_unit, * ) 
    write( output_unit, '(T6,A)' ) 'This program is free software: you can redistribute it and/or modify it &
                                  &under the terms of the'
    write( output_unit, '(T5,A)' ) 'GNU General Public License as published by the Free Software Foundation, &
                                 &either version 3 of the'
    write( output_unit, '(T31,A)' ) 'License, or (at your option) any later version.'
    write( output_unit, * ) 
    write( output_unit, '(T5,A)' ) 'This program is distributed in the hope that it will be useful, but &
                                 &WITHOUT ANY WARRANTY; without'
    write( output_unit, '(T12,A)' ) 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A &
                                 &PARTICULAR PURPOSE.'
    write( output_unit, '(T28,A)' ) 'See the GNU General Public License for more details.'
    write( output_unit, * ) 
    write( output_unit, '(T4,A)' ) 'You should have received a copy of the GNU General Public License along & 
                                 &with this program. If not,'
    write( output_unit, '(T36,A)' ) 'see <https://www.gnu.org/licenses/>.'

    write( output_unit, '(/,T3,A,/)' ) dashline

  end subroutine Display_license

  subroutine Display_citation
    implicit none

    write( output_unit, '(/,T3, A,/)' ) dashline

    write( output_unit, '(T5, A)' ) '***** This version of THEMIS uses:'
    
    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'SPHERE_GRID library (John Burkardt)'
    write( output_unit, '(T5, A)' )   'https://people.sc.fsu.edu/~jburkardt/f_src/sphere_grid/sphere_grid.html'
    
    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'XDR library (James Barnett)'
    write( output_unit, '(T5, A)' ) 'https://github.com/kmtu/xdrfort/blob/master/xdr.F90'

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) '***** Detailed information about the methodology:'
     
    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Themis: A Software to Assess Association Freenergies via Direct Estimative of &
                                   &Partition Functions.' 
    write( output_unit, '(T5, A)' )   'ChemRxiv, 2020, doi: 10.26434/chemrxiv.12925550.v2' 

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) '***** Papers that used earlier or adapted Themis versions:'

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Graphitic Carbon Nitrides as Platforms for Single-Atom Photocatalysis.'
    write( output_unit, '(T5, A)' ) 'Faraday Discussions, v. 227, p. 306, 2021, doi: 10.1039/C9FD00112C'

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Emergence of complexity in hierarchically organized chiral particles.' 
    write( output_unit, '(T5, A)' ) 'Science, v. 368, p. 642, 2020, doi: 10.1126/science.aaz7949' 

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Solvent effect on the regulation of urea hydrolysis reactions by copper complexes.'
    write( output_unit, '(T5, A)' ) 'Chemistry, v. 2, p. 525, 2020, doi: 10.3390/chemistry2020032'

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Ion pair free energy surface as a probe of ionic liquid structure.'
    write( output_unit, '(T5, A)' ) 'The Journal of Chemical Physics, v. 152, p. 014103, 2020, doi: 10.1063/1.5128693.'

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Low-Temperature Phase Transitions of the Ionic Liquid 1-Ethyl-3-methylimidazolium Dicyanamide'
    write( output_unit, '(T5, A)' ) 'The Journal of Physical Chemistry B, v. 123, p. 9418, 2019, doi: 10.1021/acs.jpcb.9b07654.'

    write( output_unit, * ) 
    write( output_unit, '(T5, A)' ) 'Site-selective photoinduced cleavage and profiling of DNA by chiral semiconductor&
                                & nanoparticles.'
    write( output_unit, '(T5, A)' ) 'Nature Chemistry, v. 10, p. 821, 2018, doi: 10.1038/s41557-018-0083-y'
    
    write( output_unit, '(/,T3, A,/)' ) dashline

  end subroutine Display_citation

end module mod_info
