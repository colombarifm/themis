###################################################################################################
#                                                                                                 #
#> @file   run_path.sh                                                                            #
#> @author Felippe M. Colombari                                                                   #
#> @brief  Run multiple Themis calculations considering different scaling factors for spherical   #
#          translation grid. Thermodynamic properties are then obtained along the intermolecular  #
#          separation coordinate.                                                                 #
#> @date - Jun, 2017                                                                              #
#> - first script created                                                                         #
#> @date - Dec, 2019                                                                              #
#> - documentation and revision                                                                   #
#> @date - Jan, 2024                                                                              #
#> - new revision                                                                                 #
#> @email  bug reports to: colombarifm@hotmail.com                                                #
#> @note   usage: ./run_path.sh --min <min_dist> --max <max_dist> --step <step>                   #
#> @note   usage: ./run_path.sh -h (for help)                                                     #
#> @note   units: Angstrom                                                                        #
#                                                                                                 #
###################################################################################################

#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

min_dist=0
max_dist=0
step=0

#-------------------------------------------------------------------------------------------------#

Show_header () {

  columns=80
  printf "\n\t"
  printf "=%.0s" {1..80}

  line[1]=""
  line[2]=""
  line[3]="RUN_PATH"
  line[4]=""
  line[5]="Run multple Themis calculations considering different scaling" 
  line[6]="factors for spherical translation grid. Thermodynamic properties"
  line[7]="are then obtained along the intermolecular separation coordinate"
  line[8]=""
  line[9]="Author : Felippe Mariano Colombari"
  line[10]=""
  line[11]="Brazil"
  line[12]="colombarifm@hotmail.com"
  line[13]=""

  for i in $( seq 1 1 13 )
  do 
  
    printf "\t%*s\n" $(((${#line[$i]}+$columns)/2)) "${line[$i]}"

  done
  
  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n"

  return 1

}

#-------------------------------------------------------------------------------------------------#

Get_initial_date () {

  date_day=$( date +%m-%d-%Y )
  date_time=$( date +%T )

  printf "\tJob started at : %10s  -  %8s\n" $date_day $date_time
  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n\t"

  return 1

}

#-------------------------------------------------------------------------------------------------#

Show_usage () {

  printf "\t\t%6s\n"            "Usage:"

  line[1]="$0 --min <min_dist> --max <max_dist> --step <step>" 
  line[2]=""
  line[3]="<min_dist> is the minimum distance" 
  line[4]="<max_dist> is the maximum distance" 
  line[5]="<step> is the step increment" 
  line[6]=""
  line[7]="Note: values in Angstrom" 
  line[8]=""
  line[9]="$0 -h (shows this help)"

  for i in $( seq 1 1 9 )
  do

    printf "\t%*s\n" $(((${#line[$i]}+$columns)/2)) "${line[$i]}"

  done

  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n"

  exit -1

}

#-------------------------------------------------------------------------------------------------#

Check_arg () {

  if [ -z "$arg" ] || [[ "$arg" == -[0-9]* ]]
  then
    printf "ERROR : invalid option for %s\n\t" $flag_meaning
    printf "TIP   : real number in flag %2s must be > 0.0\n\n\t" $flag
    printf "        Enter '-h' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit
  fi

  if [ -z "$arg" ] || [[ "$arg" == -[a-z]* ]]
  then
    printf "ERROR : invalid option for %s\n\t" $flag_meaning
    printf "TIP   : no valid argument supplied for flag %2s\n\n\t" $flag
    printf "        Enter '--help' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit
  fi

  if (( $( echo "$arg <= 0.0" | bc -l ) ))
  then
    printf "ERROR : invalid option for %s\n\t" $flag_meaning
    printf "TIP   : real number in flag %2s must be > 0.0\n\n\t" $flag
    printf "        Enter '-h' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit
  fi

}

#-------------------------------------------------------------------------------------------------#

Call_cmd_line () {

  while [[ "$#" -gt 0 ]] 
  do
  
    case "$1" in
      --help) Show_usage
              exit 0
              ;;
      --min) flag=$1
             flag_meaning="min_distance"
             shift
             arg=$1
             Check_arg 
             min_dist=$arg
             ;;
      --max) flag=$1
             flag_meaning="max_distance"
             shift
             arg=$1
             Check_arg 
             max_dist=$arg
             ;;
      --step) flag=$1
              flag_meaning="step_increment"
              shift
              arg=$1
              Check_arg 
              step=$arg
              ;;
      *) Show_usage
         exit 0
         ;;
    esac
    shift
  done

  return 1

}

Check_options () {

  if (( $( echo "$min_dist >= $max_dist" | bc -l ) ))
  then

    printf "ERROR : invalid option for %s\n\t" "min_distance/max_distance"
    printf "TIP   : real number supplied for flag %2s must be smaller than one in flag %2s\n\n\t" "--min" "--max"
    printf "        Enter '--help' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n\t"
    exit
    
  fi

  printf "Calculating thermodynamic properties :\n" 
  npoints=$( echo "1 + ($max_dist - $min_dist)/$step" | bc )
  printf "\t\t > From  : %6.2f Angstrom\n" $min_dist
  printf "\t\t > To    : %6.2f Angstrom\n" $max_dist
  printf "\t\t > dr    : %6.2f Angstrom\n" $step
  printf "\t\t > Total : %6i separation distances\n" $npoints

}
  
Inquire_files () {
  
  if [ -d "files" ]; then

    printf "\n\tReading directory \"files/\"\t\t ... found!\n"
  
  else
    
    printf "\n\tMissing directory \"files/\". Aborting.\n"
    exit
  
  fi

  files_calc="conf1.xyz conf2.xyz INPUT themis parameters1 parameters2"
  
  printf "\n\tChecking Themis input files...\n"

  for i in $files_calc; do
    
    if [ -f "files/$i" ]; then
    
      printf "\n\t\t%12s%-20s%11s" "Input file " "\"files/$i\" " " ... found!"  
    
    else

      printf "\n\t\t%12s%-20s%21s" "Input file " "\"files/$i\" " " ... not found! Aborting."  
  
      Get_final_date

      exit
    
    fi
  
  done

  printf "\n\n\tOK. Starting calculations!\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n"

}

Run_themis () {

  mkdir -p r_$distance
  cd       r_$distance
  cp       ../files/conf1.xyz   .
  cp       ../files/conf2.xyz   .
  cp       ../files/INPUT       .
  cp       ../files/parameters1 .
  cp       ../files/parameters2 .

  printf "\t\t%11s\t%06.2f\t\t%3s" " Distance :" $distance " ... "
  
  ../files/themis --run --shell $i > out

  printf "%4s\n" "DONE"

  delta_A=$( grep "TOTAL OVER TRANSLATIONAL GRID" output.log | awk '{printf $6}' )
  mTdelta_S=$( grep "TOTAL OVER TRANSLATIONAL GRID" output.log | awk '{printf $7}' )
  delta_E=$( grep "TOTAL OVER TRANSLATIONAL GRID" output.log | awk '{printf $8}' )

  cd     ..

}


Get_final_date () {

  date_day=$( date +%m-%d-%Y )
  date_time=$( date +%T )

  printf "\t"
  printf "=%.0s" {1..80}
  printf "\n\n\tOK. Job finished at : %10s  -  %8s\n" $date_day $date_time
  printf "\n\t"
  printf "=%.0s" {1..80}

}

#------------------------------------------------| JOB SEQUENCE |---------------------------------------------#

Show_header
Get_initial_date
Call_cmd_line $@
Check_options
Inquire_files

printf "%60s\n"                   "############################################################" >  profile.dat
printf "%60s\n"                   "####  Thermodynamic profile along separation coodinate  ####" >> profile.dat
printf "%35s\t%18i\t%4s\n"        "####  Number of points along path: " $npoints          "####" >> profile.dat
printf "%60s\n"                   "############################################################" >> profile.dat
printf "%06s\t%16s\t%16s\t%16s\n" "# r(A)" "dA (kJ/mol)" "dE (kJ/mol)" "-TdS (kJ/mol)"           >> profile.dat
printf "%60s\n"                   "############################################################" >> profile.dat

for i in $( seq $min_dist $step $max_dist ); do

  distance=$( printf "%06.2f" $i )

  calc_number=$( printf "%04i" $npoints )

  Run_themis $i

  printf "%06.2f\t%16.5f\t%16.5f\t%16.5f\n" $distance $delta_A $delta_E $mTdelta_S               >> profile.dat

done

Get_final_date

#------------------------------------------------|   JOB DONE   |---------------------------------------------#
