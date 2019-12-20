!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief This module contains a header
!> @note added Asdrubal Lozada-Blanco
!> - Laboratório de Química Teórica, LQT -- UFSCar                             
!> @date - Dec, 2019
!------------------------------------------------------------------------------
module mod_info
   use iso_fortran_env, only: output_unit
   use mod_constants, only: dashline, DP
   ! TODO: include in makefile 
   !   include 'revision.inc'     

   implicit none
   private

   public :: display_header

contains
   subroutine display_header()
      character( len = 16 ), parameter :: version = 'beta'
      integer,dimension(8)             :: values
      real( kind = DP )                :: timet
      integer                          :: finish, start, rate2
      character(len=7), parameter      :: revision = 'ae2e8ac'
      
!      character(len=10)                :: numbers = '1234567890'

      write(output_unit,'(T3, A)') dashline
!      write(output_unit,*) " ", numbers, numbers, numbers, numbers, numbers, numbers, numbers, numbers, numbers, numbers
      write(output_unit,'(T50, A, A)') "THEMIS"
      write(output_unit,'(/,T8, A)')" A Software to Assess Association Free Energies via Direct Estimative of Partion Function"
      write(output_unit,'(/,T40, A, A)') "Author: Felippe Colombari"
      write(output_unit,'(T50, A, A)') "Brazil"
      write(output_unit,'(/,T42, A, A)') "Program version: ",trim(version)
      write(output_unit,'(T44, A, A)') "Revision: ", revision
      write(output_unit,'(T3, A)') dashline
      write(output_unit,*)
      write(output_unit,'(T5, A)') "Contributions" 
      write(output_unit,'(T8, A, A)') "Asdrubal Lozada: ", "Error handling revision"
      write(output_unit,*)
      write(output_unit,'(T5, A)') "This version THEMIS uses:"
      write(output_unit,'(T8, A)') "BURKARD library:"
      write(output_unit,*)
      write(output_unit,'(T3, A)') dashline
      
     call Date_and_time( VALUES = values )
     write(output_unit,'(/, T5, "STARTED AT: ", i2.2, "/", i2.2, "/", i4, " - ", &
                               &i2.2, ":", i2.2, ":", i2.2)')    &
                               &values(3), values(2), values(1), &
                               &values(5), values(6), values(7)
         
   end subroutine display_header

end module mod_info
