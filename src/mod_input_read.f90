!---------------------------------------------------------------------------------------------------
! THEMIS: A code to study intermolecular recognition via direct partition function estimation                                                  
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
!   This file was written by Felippe M. Colombari and Asdrubal Lozada-Blanco.
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
!> @file   mod_input_read.f90
!> @author Felippe M. Colombari
!> @brief  This module contains a routine for INPUT reading and checking.
!> @date - Jun, 2017                                                           
!> - independent module created                                                
!> @date - Jan, 2018                                                           
!> - improved read of keywords and attributes. error messages added.
!> @date - Nov 2019
!> - update error condition by error_handling module added by Asdrubal Lozada-Blanco
!---------------------------------------------------------------------------------------------------

module mod_input_read
  use iso_fortran_env, only: output_unit
  use mod_info       , only: Display_date_time
  use mod_constants  , only: dp, int_alphabet, float_alphabet, char_alphabet, dashline, fpinf
  use mod_error_handling

  implicit none

  private 
  public Read_input_file, potential, writeframe, wrtxtc, temp, rcut_sqr, cutoff_sqr, rot1_factor, rot2_factor, &
    trans_factor, max_rot2, ref1, ref2, vector1, vector2, nstruc, atom_overlap, inter_energy, scale_factor

  ! Attributes in keyword
    
  character( len = 250 ) :: msg_line     = char(0) 
  character( len = 240 ) :: mopac_head   = char(0)
  character( len = 10 )  :: potential    = char(0) 
  character( len = 10 )  :: writeframe   = char(0) 
  character( len = 5 )   :: wrtxtc       = char(0)
  real( kind = DP )      :: temp         = 0.0
  real( kind = DP )      :: rcut         = 0.0
  real( kind = DP )      :: rcut_sqr     = 0.0
  real( kind = DP )      :: cutoff       = 0.0
  real( kind = DP )      :: cutoff_sqr   = 0.0
  real( kind = DP )      :: max_rot2     = 0.0
  integer                :: scale_factor = 0
  integer                :: trans_factor = 0
  integer                :: rot1_factor  = 0
  integer                :: rot2_factor  = 0
  integer                :: ref1         = 0
  integer                :: ref2         = 0
  integer                :: vector1      = 0
  integer                :: vector2      = 0
  integer                :: nstruc       = 0

  logical                :: key_translation_factor   = .false.
  logical                :: key_rot1_factor          = .false.
  logical                :: key_rot2_range           = .false.
  logical                :: key_rot2_factor          = .false.
  logical                :: key_potential            = .false.
  logical                :: key_temperature          = .false.
  logical                :: key_scale_factor         = .false.
  logical                :: key_write_frames         = .false.
  logical                :: key_ref_mol1             = .false.
  logical                :: key_ref_mol2             = .false.
  logical                :: key_rot_ref_mol1         = .false.
  logical                :: key_rot_ref_mol2         = .false.
  logical                :: key_shortest_distance    = .false.
  logical                :: key_cutoff_distance      = .false.
  logical                :: key_write_xtc            = .false.
  logical                :: key_lowest_structures    = .false.
  logical                :: key_mopac_job            = .false.

  logical, allocatable, dimension(:,:,:)           :: atom_overlap 
  real( kind = DP ), allocatable, dimension(:,:,:) :: inter_energy

  type(error) :: err

contains

  subroutine Read_input_file
    use mod_inquire, only: Inquire_file, Get_new_unit
    use mod_error_handling

    implicit none

    integer                          :: file_unit      
    integer                          :: nochar      = 0
    integer                          :: line        = 0
    integer                          :: ios         = 0      
    character( len = 15 ), parameter :: file_access = "sequential"
    character( len = 15 ), parameter :: file_format = "formatted"
    character( len = 40 ), parameter :: file_name   = "INPUT"
    integer                          :: c_pos
    character( len = 250)            :: buffer, keyword, attribute
    character( len = 10 )            :: line_number
    type(error)                      :: err
  
    file_unit = Get_new_unit(10)

    call Inquire_file( file_unit , file_name , file_format , file_access )

    do while ( ios == 0 )

      read( file_unit, '(A)', iostat=ios ) buffer

      buffer = adjustl(buffer)
      buffer = trim(buffer)

      if ( ios == 0 ) then

        line = line + 1

        write(line_number,'(i3)') line

        c_pos     = scan( buffer , ':' )

        keyword   = buffer(1:c_pos-2) 
        keyword   = adjustl( keyword )
        keyword   = trim( keyword )

        attribute = buffer(c_pos+2:32)
        attribute = adjustl( attribute )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ( buffer(1:19) == 'translation_factor ' ) then

          key_translation_factor = .true.

          attribute = trim( attribute )

          nochar = verify( trim(attribute), int_alphabet)

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) .or. ( attribute == '0' ) ) then

            msg_line = "Please use an integer ( > 0 ) to define the dodecahedron tesselation level."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )

            stop

          else

            read( attribute, *, iostat=ios ) trans_factor

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:12) == 'rot1_factor ' ) then

          key_rot1_factor = .true.

          nochar = verify( trim(attribute),int_alphabet)

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use an integer ( >= 0 ) to define the dodecahedron tesselation level."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )

            stop

          else

            read( attribute, *, iostat=ios ) rot1_factor

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:11) == 'rot2_range ' ) then

          key_rot2_range = .true.

          nochar = verify( trim( attribute ), float_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use a float to specify a maximum value for &
                       &rotation around mol2 axis (in degree)."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )

            stop

          else

            read(attribute, *, iostat=ios) max_rot2

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:12) == 'rot2_factor ' ) then

          key_rot2_factor = .true.

          nochar = verify( trim( attribute ), int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use an integer to specify the number of &
                       &rotation points around mol2 axis."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )

            stop

          else

            read(attribute, *, iostat=ios) rot2_factor

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:15) == 'scaling_factor ' ) then

          key_scale_factor = .true.

          nochar = verify( trim( attribute ), int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim( attribute ) ) == 0 ) ) then

            msg_line = "Please use an integer to define a scaling factor for energy."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) scale_factor

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:12) == 'temperature ' ) then

          key_temperature = .true.

          nochar = verify( trim( attribute ), float_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use a float for temperature."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) temp

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:10) == 'potential ' ) then

          key_potential = .true.

          select case( attribute )

            case( 'lj-coul', 'bh-coul', 'none' )

              read(attribute, '(A)', iostat=ios) potential

            case default

              msg_line = "Please enter a valid potential function. Options&
                         & are 'lj-coul', 'bh-coul', and 'none'."

              call err%error('e',message="while reading INPUT file.")

              call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
                &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

              call err%error('e',tip=msg_line)

              call Display_date_time( "FINISHED AT: " )
            
              stop

          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:13) == 'write_frames ' ) then

          key_write_frames = .true.

          select case( attribute )

            case( 'XYZ', 'xyz', 'MOP', 'mop', 'none' )

              read(attribute, '(A)', iostat=ios) writeframe

            case default

              msg_line = "Please enter a valid option for frame writing.&
                         & Options are 'XYZ', 'MOP' and 'none'."

              call err%error('e',message="while reading INPUT file.")

              call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
                &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

              call err%error('e',tip=msg_line)

              call Display_date_time( "FINISHED AT: " )
            
              stop

          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:10) == 'mopac_job ' ) then 

          key_mopac_job = .true.

          read(attribute, '(A)', iostat=ios) mopac_head

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:9) == 'ref_mol1 ' ) then

          key_ref_mol1 = .true.

          nochar = verify( trim( attribute ), int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use an integer ( n > 0 ) to specify the index &
                       & of atomic site from molecule 1 for centering."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")
            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) ref1

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:9) == 'ref_mol2 ' ) then

          key_ref_mol2 = .true.

          nochar = verify( trim( attribute ), int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use an integer ( n > 0 ) to specify the index &
                       & of atomic site from molecule 1 for centering."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) ref2
              
          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:13) == 'rot_ref_mol1 ' ) then

          key_rot_ref_mol1 = .true.

          nochar = verify( trim( attribute ), int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use an integer to specify which atom from &
                       &molecule 1 will be the rotation reference."
             
            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
                 &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) vector1

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:13) == 'rot_ref_mol2 ' ) then

          key_rot_ref_mol2 = .true.

          nochar = verify( trim( attribute ), int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use an integer to specify which atom from &
                       &molecule 2 will be the rotation reference."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) vector2

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:16) == 'cutoff_distance ' ) then

          key_cutoff_distance = .true.

          nochar = verify( trim( attribute ), float_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use a float to specify the maximum site-site &
                       &distance (in Angstrom) to calculate pair energies. &
                       &"//new_line('a')   

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) cutoff

            if ( cutoff < 1.0_DP ) then
                
              cutoff_sqr = fpinf

            else

              cutoff_sqr = cutoff * cutoff

            endif

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:18) == 'shortest_distance ' ) then

          key_shortest_distance = .true.

          nochar = verify( trim( attribute ), float_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) ) then

            msg_line = "Please use a float to specify the minimum site-site &
                       &distance (in Angstrom) to consider the configuration &
                       &as valid."//new_line('a')//"    Structures with atomic &
                       &contact below such value are skipped and a strongly &
                       &repulsive energy value is considered."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)
  
            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute, *, iostat=ios) rcut

            rcut_sqr = rcut * rcut

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:10) == 'write_xtc ' ) then

          key_write_xtc = .true.

          select case( attribute )

            case( "yes", "true", "T", "no", "false", "F" )

              read(attribute, '(A)', iostat=ios) wrtxtc

            case default

              msg_line = "Please enter any valid option for this keyword: &
                         &yes/true/T or no/false/F."

              call err%error('e',message="while reading INPUT file.")

              call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
                &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

              call err%error('e',tip=msg_line)

              call Display_date_time( "FINISHED AT: " )
            
              stop

          end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer(1:18) == 'lowest_structures ' ) then

          key_lowest_structures = .true.

          nochar = verify( trim(attribute) , int_alphabet )

          if ( ( nochar > 0 ) .or. ( len( trim(attribute) ) == 0 ) .or. ( attribute == '0' ) ) then

            msg_line = "use an integer ( > 0 ) to specify the number of lowest-&
                       &energy structures to write."

            call err%error('e',message="while reading INPUT file.")

            call err%error('e',check="line "//trim(adjustl(line_number))//". Keyword '"//trim(adjustl(keyword))//"' &
              &has an invalid attribute '"//trim(adjustl(attribute))//"'.")

            call err%error('e',tip=msg_line)

            call Display_date_time( "FINISHED AT: " )
            
            stop

          else

            read(attribute,*, iostat=ios) nstruc

          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else if ( buffer == '' ) then

          !write(*,*) "blank lines are cool too"

          continue

        else 

          write(*,'(/,T3,A,A,A,i2,/)') "  ERROR: UNKNOWN KEYWORD '", trim(buffer), "' ON LINE ", line

          call Display_date_time( "FINISHED AT: " )

          stop

        endif

          continue      

      endif

    enddo

    close( file_unit )

    call check_keys

    return
  end subroutine Read_input_file

  subroutine Check_keys
    USE mod_cmd_line, only: grid_type
    implicit none

    if ( ( key_translation_factor .eqv. .false. ) .and. ( grid_type == "shell" ) ) then
        call err%error('e',message="Missing valid entry for 'translation_factor' on INPUT file!")
    else if ( ( key_translation_factor .eqv. .true. ) .and. ( grid_type /= "shell" ) ) then
        call err%error('w',message="NOTE: Ignoring unused INPUT entry 'translation_factor'")
    else if ( key_rot1_factor .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'rot1_factor' on INPUT file!")
    else if ( key_rot2_factor .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'rot2_factor' on INPUT file!")
    else if ( key_rot2_range .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'rot2_range' on INPUT file!")
    else if ( key_scale_factor .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'scaling_factor' on INPUT file!")
    else if ( key_temperature .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'temperature' on INPUT file!")
    else if ( key_potential .eqv. .false. ) then
        call err%error("Missing valid entry for 'potential' on INPUT file!")
    else if ( key_write_frames .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'write_frames' on INPUT file!")
    else if ( key_ref_mol1 .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'ref_mol1' on INPUT file!")
    else if ( key_ref_mol2 .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'ref_mol2' on INPUT file!")
    else if ( key_rot_ref_mol1 .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'rot_ref_mol1' on INPUT file!")
    else if ( key_rot_ref_mol2 .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'rot_ref_mol2' on INPUT file!")
    else if ( key_shortest_distance .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'shortest_distances' on INPUT file!")
    else if ( key_cutoff_distance .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'cutoff_distances' on INPUT file!")
    else if ( key_write_xtc .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'write_xtc' on INPUT file!")
    else if ( key_lowest_structures .eqv. .false. ) then
        call err%error('e',message="Missing valid entry for 'lowest_structures' on INPUT file!")
    endif

    if ( ( writeframe == "MOP" ) .and. ( key_mopac_job .eqv. .false. ) ) then
        call err%error('e',message="MOPAC input files will be written but no 'mopac_job' entry &
                 &was found on INPUT file!")
    else if ( ( writeframe /= "MOP" ) .and. ( key_mopac_job .eqv. .true. ) ) then
        call err%error('w',message="Ignoring unused INPUT entry 'mopac_job'")
    endif
          
    call Display_date_time( "FINISHED AT: " )

  end subroutine Check_keys

end module mod_input_read
