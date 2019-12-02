!------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition      
!         function estimation                                                  
!                                                                                   
! Copyright (C) 2017 Felippe M. Colombari                                      
!------------------------------------------------------------------------------
!> @brief This module contains procedures to error handling
!> This is a Standard Fortran 2008 compliance code
!> @note Added by Asdrubal Lozada-Blanco                                                
!> - Laboratório de Química Teórica, LQT -- UFSCar                             
!> @date - Nov, 2019                                                         
!> - independent module created
!------------------------------------------------------------------------------
module mod_error_handling
   use iso_fortran_env, only: output_unit

   implicit none
   private

   type, public :: error
      integer                       :: code
      character                     :: type
      character(len=:), allocatable :: message, check, tip
   contains
      procedure :: error => raise_error
      procedure :: termination => normal_termination
   end type error

   public :: normal_termination

contains
   subroutine raise_error(self, type, code, message, check, tip)
      class(error), intent(inout)            :: self
      character,  intent(in)                 :: type        
      integer, optional,  intent(in)         :: code        
      character(*), optional, intent(in)     :: message, check, tip

      self%type = type

      if(present(code)) self%code = code
      if(present(message)) self%message = message
      if(present(check)) self%check = check
      if(present(tip)) self%tip = tip

      select case(type)
         case('w')
         !   if(.not.present(message)) write(output_unit,'("Warning: unexpected event.")')     
            if(present(message)) write(output_unit,'(/,T5,"Warning: ",a)') self%message
            if(present(check)) write(output_unit,'(/,T5,"Warning: ",a)') self%check
            if(present(tip)) write(output_unit,'(/,T5,"Warning: ",a)') self%tip
            if(present(code)) write(output_unit,'("Code: ",i0)') self%code
         case('e')
         !   if(.not.present(message)) write(output_unit,'("Error: abnormal condition.")')     
            if(present(message)) write(output_unit,'(/,T5,"Error: abnormal condition ",a)') self%message
            if(present(check)) write(output_unit,'(/,T5,"Please check ",a)') self%check
            if(present(tip)) write(output_unit,'(/,T5,"TIP: ",a)') self%tip
            if(present(code)) write(output_unit,'("Code: ",i0)') code
      end select

   end subroutine raise_error   

   subroutine normal_termination(self,iostat,status)
      class(error), intent(inout)            :: self
      integer, intent(in)      :: iostat
      character(*), intent(in) :: status

      if(iostat == 0) then
        select case(status)
           case('f')
               write(output_unit,'("Sucessfull termination.")')
               stop
           case('i')
               write(output_unit,'("Normal termination.")')        
        end select
      end if

   end subroutine normal_termination


end module mod_error_handling
