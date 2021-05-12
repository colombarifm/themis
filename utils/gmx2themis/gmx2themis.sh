###################################################################################################
#                                                                                                 #
#> @file   gmx2themis.sh                                                                          #
#> @author Felippe M. Colombari                                                                   #
#> @brief  Extract oplsaa charge, sigma and epsilon parameters from gromacs files and write them  #
#>         to Themis parameters file.                                                             #
#> @date - Mar, 2019                                                                              #
#> - first script created                                                                         #
#> @date - May, 2021                                                                              #
#> - documentation and revision                                                                   #
#> @email  bug reports to: colombarifm@hotmail.com                                                #
#> @note   usage: ./gmx2themis.sh --top <topology file> --ff <ffnonbonded file>                   #
#>                                --numat <number of atoms> --out <output file>                   #
#> @note   usage: ./gmx2themis.sh --help (for help)                                               #
#> @note   units: Angstrom and kJ/mol                                                             #
#> @note   units: Please check the resulting file carefully!!!                                    #
#                                                                                                 #
###################################################################################################

#!/bin/bash 

export LC_NUMERIC="en_US.UTF-8"

#-------------------------------------------------------------------------------------------------#

Show_usage () {

  columns=80

  printf "\t"
  printf "=%.0s" {1..80}
  printf "\n"
  printf "\t\t%6s"            "Usage:"

  printf "\n\n"

  line[1]="$0 --top <topology file> --ff <ffnonbonded file>" 
  line[2]="   --numat <number of atoms> --out <output file>"
  line[3]=""
  line[4]="   <topology file> is the file containing the [ atoms ] directive" 
  line[5]="<ffnonbonded file> is the file containing the [ atomtypes ] directive" 
  line[6]=" <number of atoms> is the number of atoms" 
  line[7]="     <output file> is the file in which the parameters will be written"
  line[8]=""
  line[9]="Note: units are Angstrom and kJ/mol" 
  line[10]=""
  line[11]="$0 --help (shows this help)"

  for i in $( seq 1 1 11 )
  do

    printf "\t%*s\n" $(((${#line[$i]}+$columns)/2)) "${line[$i]}"

  done

  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n"

  exit -1

}

#-------------------------------------------------------------------------------------------------#

Check_arg_num () {

  if [ -z "$arg_num" ] || [[ "$arg_num" == -[0-9]* ]] || [[ "$arg_num" == *[a-z]* ]]
  then
    printf "ERROR : invalid option for %s\n\t" $flag_meaning
    printf "TIP   : integer in flag %2s must be > 0\n\n\t" $flag
    printf "        Enter '-h' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit
  fi

}

#-------------------------------------------------------------------------------------------------#

Check_arg_str () {

  if [ -z "$arg_str" ] || [[ "$arg_str" == -[a-z]* ]]
  then
    printf "ERROR : invalid option for %s\n\t" $flag_meaning
    printf "TIP   : no valid argument supplied for flag %2s\n\n\t" $flag
    printf "        Enter '-h' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit
  fi

}

#-------------------------------------------------------------------------------------------------#

Call_cmd_line () {

  while [[ "$1" == -* ]] 
  do
  
    case "$1" in
      -h|--help) Show_usage
              exit 0
              ;;
      --top) flag=$1
          flag_meaning="topology_file"
          shift
          arg_str=$1
          Check_arg_str
          topology_file=$arg_str 
          ;;
      --ff) flag=$1
          flag_meaning="ffnonbonded_file"
          shift
          arg_str=$1
          Check_arg_str
          ffnonbonded_file=$arg_str
          ;;
      --num) flag=$1
          flag_meaning="number_of_atoms"
          shift
          arg_num=$1
          Check_arg_num
          numat=$arg_num
          ;;
      --out) flag=$1
          flag_meaning="output_file"
          shift
          arg_str=$1
          Check_arg_str
          output=$arg_str
          ;;
      --)
          shift
          break
          ;;
    esac
    shift
  done

  return 1

}


#-------------------------------------------------------------------------------------------------#

Check_topology () {

  required_string="[ atoms ]"

  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n\t"
  printf "Looking for topology file : %s" "$topology_file"

  if [ -f $topology_file ]
  then

    printf " FOUND\n\t"
    printf "Checking topology file    :" 

    if grep -Fxq "$required_string" $topology_file
    then
      printf " %s directive FOUND\n" "$required_string"
    else
      printf "\n\t"
      printf "ERROR : invalid topology file!\n\t" 
      printf "TIP   : %s directive not found in %s file \n\t" "$required_string" "$topology_file"
      printf "        Enter '-h' option for help.\n\n\t"
      printf "=%.0s" {1..80}
      printf "\n"
      exit
    fi

  else

    printf "\n\t"
    printf "ERROR : invalid topology file!\n\t" 
    printf "TIP   : %s file not found \n\t" "$topology_file"
    printf "        Enter '-h' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit

  fi

}

#-------------------------------------------------------------------------------------------------#

Read_topology () {

  word1="opls_[0-9][0-9][0-9] "    
  word2="opls_[0-9][0-9][0-9][A-Z]"  

  i=1
  while read -r line
  do

    case $line in 
            
      *${word1}*) 
  
        atom_type[$i]=$( echo $line | awk '{printf "%-9s", $2}' )
        #echo "${atom_type[$i]}" ${#atom_type[$i]} "A"
        atom_symb[$i]=$( echo $line | awk '{printf "%4s", $5}' )
        atom_chrg[$i]=$( echo $line | awk '{printf "%10.5f", $7}' )
        (( i++ ))
        ;;
      
      *${word2}*) 
  
        atom_type[$i]=$( echo $line | awk '{printf "%-9s", $2}' )
        #echo "${atom_type[$i]}" ${#atom_type[$i]} "B"
        atom_symb[$i]=$( echo $line | awk '{printf "%4s", $5}' )
        atom_chrg[$i]=$( echo $line | awk '{printf "%10.5f", $7}' )
        (( i++ ))
        ;;
      *)
        ;;

    esac
    shift
  done < $topology_file

}

#-------------------------------------------------------------------------------------------------#

Check_nonbonded () {

  required_string="[ atomtypes ]"
  
  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n\t"
  printf "Looking for forcefield file : %s " "$ffnonbonded_file"

  if [ -f $ffnonbonded_file ]
  then

    printf " FOUND\n\t"
    printf "Checking forcefield file    :" 

    if grep -Fxq "$required_string" $ffnonbonded_file
    then
      printf " %s directive FOUND\n" "$required_string"
    else
      printf "\n\t"
      printf "ERROR : invalid forcefield file!\n\t" 
      printf "TIP   : %s directive not found in %s file \n\t" "$required_string" "$ffnonbonded_file"
      printf "        Enter '-h' option for help.\n\n\t"
      printf "=%.0s" {1..80}
      printf "\n"
      exit
    fi

  else

    printf "\n\t"
    printf "ERROR : invalid forcefield file!\n\t" 
    printf "TIP   : %s file not found \n\t" "$ffnonbonded_file"
    printf "        Enter '-h' option for help.\n\n\t"
    printf "=%.0s" {1..80}
    printf "\n"
    exit

  fi

}

#-------------------------------------------------------------------------------------------------#

Read_nonbonded () {

  i=1

  printf "%4s\t%10s\t%10s\t%10s\n" "###" "chrg" "sigma (A)" "eps (kJ/mol)" > $output

  for i in $( seq 1 1 $numat )
  do

    atom_sigm[$i]=$( grep "${atom_type[$i]}" $ffnonbonded_file | awk '{printf "%10.5f", $7*10}' )
    atom_epsl[$i]=$( grep "${atom_type[$i]}" $ffnonbonded_file | awk '{printf "%10.5f", $8}' )

  done

}

#-------------------------------------------------------------------------------------------------#

Write_parameters () {

  printf "\n\t"
  printf "=%.0s" {1..80}
  printf "\n\n\t"
  printf "Writing parameters to %s file ..." $output

  for i in $( seq 1 1 $numat )
  do

    printf "%4s\t%10.5f\t%10.5f\t%10.5f\n" ${atom_symb[$i]} ${atom_chrg[$i]} ${atom_sigm[$i]} ${atom_epsl[$i]} >> $output

  done

  printf " DONE"
  printf "\n\n\t"
  printf "=%.0s" {1..80}
  printf "\n\t"

}

#-------------------------------------------------------------------------------------------------#

Call_cmd_line $@
Check_topology
Read_topology
Check_nonbonded
Read_nonbonded
Write_parameters
