!---------------------------------------------------------------------------------------------------
! THEMIS: A software to assess association free energies via direct estimative of partition functions                                                  
!---------------------------------------------------------------------------------------------------
!
!   Free software, licensed under GNU GPL v3
!
!   Copyright (c) 2017 - 2023 Themis developers
!
!   This file was written by Asdrubal Lozada-Blanco
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
!> @file   mod_plot.f90
!> @author Asdrubal Lozada-Blanco
!> @brief  Call plot subroutine.
!> @date - Oct, 2023                                                         
!> - independent module created                                                
!---------------------------------------------------------------------------------------------------

module mod_plot

contains
  !TODO
  ! Add quatities
  subroutine Call_plot()
    character(len=10), parameter :: data_file = 'output.log'
    character(len=8), parameter :: plot_file =  'plot.plt'
    logical :: exists

    exists = .true.

    inquire(file=plot_file,exist=exists)  
    if(.not.exists) then
      open(99,file=plot_file,status='unknown')
      write(99,*) '# plot.plt'
      write(99,*) 'set term x11'
      write(99,*) 'set title "Quantities"'
      write(99,*) 'set xlabel "A"'
      write(99,*) 'set ylabel "(kJ/mol)"'
      write(99,*) 'plot "'//trim(adjustl(data_file))//'" using 2:7 with lines ls 7 notitle'
    end if

    close(99)

    !TODO
    ! check if plotting program exists

    call execute_command_line('gnuplot -p '//plot_file)

  end subroutine call_plot        
end module mod_plot
