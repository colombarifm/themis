###################################################################################################
#                                                                                                 #
#> @file   top2themis.py                                                                          #
#> @author Felippe M. Colombari                                                                   #
#> @brief  Extract charge, sigma and epsilon parameters from gromacs top+itp files and write them #
#>         to Themis parameters file.                                                             #
#> @date - Mar, 2025                                                                              #
#> - first python script created (better than the bash old one)                                   #
#> @email  bug reports to: colombarifm@hotmail.com                                                #
#> @note   usage: python3 top2themis.py --top <topology file> --itp <ffnonbonded file>            #
#>                                --out <output file>                                             #
#> @note   final units: Angstrom and kJ/mol                                                       #
#> @note   Please check the resulting file carefully!!!                                           #
#> @note   This version was tested for oplsaa and charmm36 forcefields!!                          #
#                                                                                                 #
###################################################################################################

import os
import re
import numpy as np
import argparse

term_size = 105
print('\n|' + '=' * term_size + '|\n')

######################################### Command line arguments parsing ##############################################

description = "Tool to convert Gromacs .top/.itp parameters to Themis format"

parser = argparse.ArgumentParser( description = description )
parser.add_argument("-t","--top", help="Gromacs .top file (must contain [ atoms ] directive)",     required=True)
parser.add_argument("-i","--itp", help="Gromacs .itp file (must contain [ atomtypes ] directive)", required=True)
parser.add_argument("-o","--out", help="Output parameter file", required=True)

args     = parser.parse_args()
top_file = args.top
itp_file = args.itp
out_file = args.out

############################################# Input files checking ####################################################

def Check_tops( file_in, directive ):
    if os.path.isfile(file_in):
        print(f"\t\tTopology file '{file_in}' found\n")
        if directive in open(file_in).read():
            print("\t\t%s directive found\n" % directive)
        else:
            print("\t\t%s directive not found! Aborting...\n" % directive)
            exit()
    else:
        print(f"\t\tTopology file '{file_in}' not found! Aborting...\n")
        exit()

######################################## Array/dictionay initialization ###############################################

label  = []
atype  = []
charge = [] 
nb_sig = []
nb_eps = []

arr_dic = dict( 
        Labels  = label,
        Atypes  = atype,
        Charges = charge,
        NBs_sig = nb_sig,
        NBs_eps = nb_eps
        )

###################################### Read atoms from molecule .top file #############################################

print("\tTOP filename read: %s\n" % top_file)
Check_tops(top_file, '[ atoms ]' )

ftop  = open(top_file)
flist = ftop.readlines()

def Get_atoms( ):
    parsing  = False
    tmp      = []
    ini_str  = "[ atoms ]"
    cmt_str  = ";"
    skp_str  = "\n"
    end_str  = "[ bonds ]"

    for line in flist:
        if line.startswith(ini_str):
            parsing = True
        elif line.startswith(end_str):
            parsing = False
        if parsing:
            if not line.startswith(ini_str):
                if not line.startswith(skp_str) and not line.startswith(cmt_str):
                    tmp  = line.strip().split() 
                    atype.append(tmp[1])
                    label.append(tmp[4])
                    charge.append(np.float(tmp[6]))
    arr_dic["Atypes"]  = atype
    arr_dic["Labels"]  = label
    arr_dic["Charges"] = charge

Get_atoms()

ftop.close()

n_atoms = len(arr_dic["Atypes"])
n_types = len(set(arr_dic["Atypes"]))

################################## Read atomtypes from forcefield .itp file ###########################################

print("\tITP filename read: %s\n" % itp_file)
Check_tops(itp_file, '[ atomtypes ]')

fitp  = open(itp_file) 
flist = fitp.readlines()

def Get_atomtypes( atom ):
    parsing  = False
    tmp      = []
    ini_str  = "[ atomtypes ]"
    cmt_str  = ";"
    skp_str  = "\n"
    end_str1 = "[ moleculetype ]"
    end_str2 = "[ pairtypes ]"

    for line in flist:
        if line.startswith(ini_str):
            parsing = True
        elif line.startswith(end_str1) or line.startswith(end_str2):
            parsing = False
        if parsing:
            if not line.startswith(ini_str):
                if not line.startswith(skp_str) and not line.startswith(cmt_str):
                    tmp = line.strip().split() 
                    nb_sig.append(10*np.float(tmp[5]))
                    nb_eps.append(np.float(tmp[6]))
    arr_dic["NBs_sig"] = nb_sig
    arr_dic["NBs_eps"] = nb_eps
   
fitp.close()

####################################### Write parameters to output file ###############################################

fout = open(out_file, "w")

fout.write("%4s\t%12s\t%12s\t%12s\n" % ( "###", "charge", "sig_(A)", "eps_(kJ/mol)"))

for i in range(0, n_atoms):

    iatom = arr_dic["Atypes"][i]

    Get_atomtypes( atom = iatom )

    fout.write("%4s\t%12.6f\t%12.6f\t%12.6f\n" % ( arr_dic["Labels"][i], (arr_dic["Charges"][i]), arr_dic["NBs_sig"][i], arr_dic["NBs_eps"][i]))

fout.close()

########################################### Write final info and exit #################################################

print("\tParameters written to file: %s\n" % out_file)
print("\t\tNumber of atoms: %i\n" % n_atoms) 
print("\t\tNumber of unique atomtypes: %i\n" % n_types) 
print("\t\tTotal charge: %12.9f\n" % np.sum(arr_dic["Charges"]))

print('|' + '=' * term_size + '|')
