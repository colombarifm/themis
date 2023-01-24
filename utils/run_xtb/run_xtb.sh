###################################################################################################
#                                                                                                 #
#> @file   run_xtb.sh                                                                             #
#> @author Felippe M. Colombari                                                                   #
#> @brief  Generate an ensemble of XYZ coordinates with Themis and loop over them performing      # 
#          single-point calculations with xTB software.                                           #
#> @date - Aug, 2019                                                                              #
#> - first version                                                                                #
#> @date - Nov, 2020                                                                              #
#> - documentation and revision                                                                   #
#> @date - Jun, 2022                                                                              #
#> - add more path verifications (XTBHOME and EXEC)                                               #
#> @date - Jul, 2022                                                                              #
#> - add automatic generation/calculation for separated molecules                                 #
#> - include internal energy for mol 2 conformational ensemble                                    #
#> @email  bug reports to: colombarifm@hotmail.com                                                #
#> @note   usage: ./run_xtb.sh --charge <charge> --uhf <uhf> (--restart)                          #
#> @note   usage: ./run_xtb.sh --help (for help)                                                  #
#                                                                                                 #
###################################################################################################

#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"

#### Set xTB environment variables according to 
#### https://xtb-docs.readthedocs.io/en/latest/setup.html#setting-up-xtb

XTBHOME=
XTBPATH=${XTBHOME}:${HOME}
MANPATH=${MANPATH}:${XTBHOME}/man
PATH=${PATH}:${XTBHOME}/bin:${XTBHOME}/python
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${XTBHOME}/lib
PYTHONPATH=${PYTHONPATH}:${XTBHOME}/python
EXEC=

export XTBHOME PATH XTBPATH MANPATH LD_LIBRARY_PATH PYTHONPATH

#### Set system settings (memory, parallelization, etc)

# To distribute the number of threads reasonable in the OMP section it is recommended to use
export OMP_NUM_THREADS=2,1

# The default linear algebra backend of xtb is the Math Kernel Library, to run it in parallel export
export MKL_NUM_THREADS=2

# You might want to deactivate nested OMP constructs by
export OMP_MAX_ACTIVE_LEVELS=1

# To calculate larger systems an appropriate OMP stacksize must be provided
export OMP_STACKSIZE=4G

# To avoid stack overflows when calculating large molecules, you should unlimit the system stack
ulimit -s unlimited

#-------------------------------------------------------------------------------------------------#

Show_header () {

  columns=80
  printf "\n\t"
  printf "=%.0s" {1..80}

  line[1]=""
  line[2]=""
  line[3]="RUN_XTB"
  line[4]=""
  line[5]="Generate an ensemble of XYZ coordinates with Themis"
  line[6]="Loop over them performing single-point calculations with xTB software"
  line[7]=""
  line[8]="Author : Felippe Mariano Colombari"
  line[9]=""
  line[10]="Brazil"
  line[11]="colombarifm@hotmail.com"
  line[12]=""

  for i in $( seq 1 1 12 )
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
  printf "\n\t"

  return 1

}

#-------------------------------------------------------------------------------------------------#

Show_usage () {

  printf "\t\t%s\n"            "Usage:"

  line[1]="$0 --charge <charge> --uhf <uhf> (--restart)" 
  line[2]=""
  line[3]="--charge <charge> is the total charge of the structures"
  line[4]="--uhf <uhf> is the number of umpaired electrons of the structure"
  line[5]="--restart is the optional flag to restart xTB calculations"
  line[6]=""
  line[7]="$0 --help (shows this help)"

  for i in $( seq 1 1 7 )
  do

    printf "\t%*s\n" $(((${#line[$i]}+$columns)/2)) "${line[$i]}"

  done

  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n"

  exit -1

}

#-------------------------------------------------------------------------------------------------#

Check_arg () {

  if [ -z "$arg" ] #|| [[ "$arg" =~ ^-?[0-9]+$ ]]
  then

    printf "\n\t%s" "$error_msg"
    printf "\n\t%s" "$help_msg"
    printf "\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit 0

  fi

  if [[ "$arg" =~ ^-?[0-9]+$ ]]
  then

    echo "" > /dev/null

  else

    printf "\n\t%s" "$error_msg"
    printf "\n\t%s" "$check_msg"
    printf "\n\t%s" "$help_msg"
    printf "\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit 0

  fi
}

#-------------------------------------------------------------------------------------------------#

Call_cmd_line () {

  help_msg="TIP   : Enter '--help' option for help"
  
  if (( $# < 1 ))
  then
    error_msg="ERROR : no arguments given"
    Check_arg
  fi
  
  while (( $# > 0 ))
  do
    if [[ "$1" != --* ]]; then
      error_msg="ERROR : Expected an option beginning with '--' but found '$1'"
      Check_arg
    elif [[ "$1" == -- ]]; then
      error_msg="ERROR : unknown command-line option '$1'"
      Check_arg
    fi
  
    case "$1" in
      --help)    Show_usage
                 exit 0
                 ;;
      --charge)  flag=$1
                 flag_meaning="total_charge"
                 found_charge=true   
                 shift
                 arg=$1
                 error_msg="ERROR : invalid option for '$flag_meaning'"
                 check_msg="CHECK : argument for flag '$flag' must be an integer"
                 Check_arg 
                 chrg=$arg
                 ;;
      --uhf)     flag=$1
                 flag_meaning="unpaired_electrons"
                 found_uhf=true
                 shift
                 arg=$1
                 error_msg="ERROR : invalid option for '$flag_meaning'"
                 check_msg="CHECK : argument for flag '$flag' must be an integer"
                 Check_arg 
                 uhf=$arg
                 ;;
      --restart) flag=$1
                 flag_meaning="restart_job"
                 restart_job=true
                 ;;
      *)         flag=$1
                 printf "\n\t%s" "ERROR : unknown option for '$flag'"
                 printf "\n\t%s"   "$help_msg"
                 printf "\n\n\t"
                 printf "=%.0s" {1..80}
                 printf "\n"
                 exit 0
                 ;;
    esac
    shift

  done

  if [ $found_charge = false ] || [ $found_uhf = false ]
  then

    printf "\n\t%s" "ERROR : some flags and arguments are missing"
    printf "\n\t%s" "$help_msg"
    printf "\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit 0

  fi

  return 1

}

Check_options () {

  if [ -z $XTBHOME ]
  then

    printf "\n\t%s" "ERROR : XTBHOME path not set" 
    printf "\n\t%s" "TIP   : edit line 21 of this script" 
    printf "\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit 0

  fi

  if [ -f $EXEC ]
  then

    printf "\n\t%s" "ERROR : EXEC path not set" 
    printf "\n\t%s" "TIP   : edit line 30 of this script" 
    printf "\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit 0

  fi

}

Inquire_files () {
  
  if [ -d "0_files" ]; then

    printf "\n\t%27s\t%20s\n" "Reading directory '0_files'" "... found!"
  
  else
    
    printf "\n\t%27s\t%20s" "Missing directory '0_files' ... ABORTING!"
    
    Get_final_date
    
    exit
  
  fi

  files_calc="conf1.xyz conf2.xyz INPUT themis grid.xyz input_xtb"
  
  printf "\n\t%s\n" "Checking Themis input files in 0_files/"

  for i in $files_calc; do
    
    if [ -f "0_files/$i" ]; then
    
      printf "\n\t\t%12s%-20s%12s" "Input file " "'0_files/$i'" " ... found!"  
    
    else

      printf "\n\t\t%12s%-20s%21s" "Input file " "'0_files/$i'" " ... not found! Aborting."  
  
      Get_final_date

      exit
    
    fi
  
  done

  printf "\n\n\t%s" "OK. Starting calculations!"
  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n"

}

Run_themis () {

  mkdir -p 1_structures/
  cd       1_structures/
  cp       ../0_files/conf1.xyz   .
  cp       ../0_files/conf2.xyz   .
  cp       ../0_files/INPUT       .
  cp       ../0_files/grid.xyz    .
  cp       ../0_files/themis      .

  printf "\n\t%s\n"       "Step 1: Ensemble generation with Themis" 
  printf "\n\t\t%s\t%42s" "Running Themis " "... "
  
  ../0_files/themis --run --user grid.xyz > ../log_themis.out

  printf "%s\n" "DONE"

  #

  mkdir -p separated_molecules/
  cd       separated_molecules/
  cp       ../conf1.xyz .
  cp       ../conf2.xyz .
  cp       ../INPUT     .
  cp       ../themis    .

  # zero here means one translation point
  sed -i 's/.*point_rot_factor :.*/\point_rot_factor : 0 /'     INPUT
  # zero here means one rotation move around ref atom of mol2
  sed -i 's/.*translation_factor :.*/\translation_factor : 0 /' INPUT
  # one here means one rotation move around axis of mol2
  sed -i 's/.*axis_rot_moves :.*/\axis_rot_moves : 1 /'         INPUT

  printf "\n\t%s\n"       "Step 1.1: Separated molecules generation with Themis" 
  printf "\n\t\t%s\t%42s" "Running Themis " "... "
  
  ../../0_files/themis --run --shell 50.0 > log_themis.out

  printf "%s\n" "DONE"

  cd ..

  #

  number_of_files=$( find -maxdepth 1 -type f -name 'point_*.xyz' | wc -l )

  printf "\t\t%s\t%10s%i" "Coordinate files written to 1_structures/ " "... " $number_of_files
  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n\t"

}

Run_continuation () {

  cd 1_structures

  printf "\n\t%s\n"   "Step 1: Ensemble generation with Themis" 
  printf "\n\t\t%s\n" "Continuation job" 
  
  number_of_files=$( find -maxdepth 1 -type f -name 'point_*.xyz' | wc -l )
  
  printf "\t\t%s\t%10s%i" "Remaining coordinate files in 1_structures/ " "... " $number_of_files

  if [ "$number_of_files" == 0 ]
  then

    printf "\n\n\t\t%s" "Nothing to do... exiting"

    Get_final_date
    
    exit
  
  fi

  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n\t"

}

Run_xtb () {

  printf "\n\t%16s%i%35s"   "Step 2: Perform " $number_of_files " single-point calculations with xTB"
  printf "\n\n\t\t%s\t%26s" "Setting up xTB calculations " "... "

  mkdir -p ../2_finished

  cp ../0_files/input_xtb .

  printf "%s"             "DONE"
  printf "\n\t\t%s\t%36s" "Starting xTB calculations " "... good luck!"
  printf "\n\n\t\t" 
  printf "*%.0s" {1..70}
  
  printf "\n\n\t\t%s" "NOTE : YOU CAN CHECK THE STATUS OF THE JOB"
  printf "\n\n\t\t%s" "* remaining coordinate files ~~~~~> 1_structures/ directory:"
  printf "\n\n\t\t%s" "  find 1_structures/ -maxdepth 1 -type f -name 'point_*.xyz' | wc -l"
  printf "\n\n\t\t%s" "* files of finished xTB calculations ~~~~~> 2_finished/ directory"
  printf "\n\n\t\t%s" "  find 2_finished/ -maxdepth 1 -type f -name 'point_*.xyz' | wc -l"
  
  printf "\n\n\t\t"
  printf "*%.0s" {1..70}
  
  ### calculations for the ensemble of interacting molecules

  for i in $( find -maxdepth 1 -type f -name 'point_*.xyz' )
  do                           
  
    file=$( basename $i .xyz ) 

    #Check_dummy

    $EXEC \
      $file.xyz \
    --uhf $uhf \
    --chrg $chrg \
    --namespace $file \
    --norestart \
    --input input_xtb > $file.log 2> /dev/null 
                       
    mv $file.xyz ../2_finished/
    mv $file.log ../2_finished/
                        
  done 

  printf "\n\n\t\t%s" " DONE! Calculating internal energies..."
  printf "\n\n\t"
  
  ### calculations for the ensemble of non-interacting molecule 2 conformations

  cd  separated_molecules/

  rm internal_energies.dat 2&> /dev/null

  for i in $( find -maxdepth 1 -type f -name 'point_*.xyz' )
  do                           
  
    file=$( basename $i .xyz ) 

    #Check_dummy
  
    $EXEC \
      $file.xyz \
    --uhf $uhf \
    --chrg $chrg \
    --namespace $file \
    --norestart \
    --input ../input_xtb > $file.log 2> /dev/null 
    
    grep "TOTAL ENERGY" $file.log | awk '{printf $4}' >> internal_energies.dat

  done
  
  cp internal_energies.dat ../2_finished/

  cd ..

  printf "\n\n\t\t%s" " DONE! Now, follow the next steps"
  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n"

  mkdir -p ../3_rerun_themis
  cp       ../0_files/conf1.xyz ../3_rerun_themis/
  cp       ../0_files/conf2.xyz ../3_rerun_themis/
  cp       ../0_files/INPUT     ../3_rerun_themis/
  cp       ../0_files/grid.xyz  ../3_rerun_themis/
  cp       ../0_files/themis    ../3_rerun_themis/

  printf "\n\t%s"   "Step 3: Extract total energies from output files and get interaction energies"
  printf "\n\n\t%s" "        * Copy xtb_extract.py to 2_finished/"
  printf "\n\n\t%s" "        * Edit xtb_extract.py properly ACCORDING TO YOUR SYSTEM"
  printf "\n\t%s"   "          ** Enter number of translation and rotation moves" 
 # printf "\n\t%s"   "          ** Enter internal energy (separated dimer energy)"
  printf "\n\n\t%s" "        * python3 xtb_extract.py"

  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n"
  
  printf "\n\t%s"   "Step 4: Calculate properties with Themis"
  printf "\n\n\t%s" "        * Copy energy.log to 3_rerun_themis/"
  printf "\n\n\t%s" "        * ./themis --rerun --user grid.xyz"

}

Check_dummy () {

  ####################### delete dummy sites (X) from XYZ files ############################
  ### use this routine if dummy sites (X) were used in conf1.xyz and/or conf2.xyz files:
  ### these atom labels are not accepted by xTB. For example, if there was 1 dummy site
  ### (center of mass) for conf1.xyz and another one (center of mass) for conf2.xyz, a
  ### total of two dummy sites will be found in the XYZ files (point_*.xyz). To run xTB, 
  ### we must:
  ###
  ### 1. read the number of atoms (with dummy sites)
  old_numat=$( head -1 $file.xyz )                  

  ### 2. read the number of dummy sites (labeled as X)
  n_dummy=$( grep -c "^X" $file.xyz )

  ### 3. subtract the number of dummy sites from the number of atoms 
  new_numat=$( echo "$old_numat - $n_dummy" | bc )  

  ### 4. change the number of atoms of input files
  sed -i "1 s/^.*$/$new_numat/" $file.xyz         

  ### 5. delete dummy sites of input files
  sed -i '/X/d' $file.xyz 
  #############################################################################################

}

Get_final_date () {

  date_day=$( date +%m-%d-%Y )
  date_time=$( date +%T )

  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n\tJob finished at : %10s  -  %8s\n\n\t" $date_day $date_time
  printf "=%.0s" {1..80}

}

#------------------------------------------------| JOB SEQUENCE |---------------------------------------------#

found_uhf=false
found_charge=false
restart_job=false

Show_header
Get_initial_date
Call_cmd_line $@
Check_options
Inquire_files
if [ "$restart_job" = true ]
then
  Run_continuation
else
  Run_themis 
fi
Run_xtb
Get_final_date

#------------------------------------------------|   JOB DONE   |---------------------------------------------#
