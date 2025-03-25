###################################################################################################
#                                                                                                 #
#> @file   gmxdump2themis.py                                                                      #
#> @author Felippe M. Colombari                                                                   #
#> @brief  Extract charge, sigma and epsilon parameters from gromacs .tpr dumpfile and write them #
#>         to Themis parameters file.                                                             #
#> @date - Mar, 2025                                                                              #
#> - first python script created                                                                  #
#> @email  bug reports to: colombarifm@hotmail.com                                                #
#> @note   usage: python3 gmxdump2themis.py --tpr <tpr file> --gmx <gromacs executable>           #
#>                                --dump <dumpfile for .tpr> --out <final output>                 #
#> @note   final units: Angstrom and kJ/mol                                                       #
#> @note   Please check the resulting file carefully!!!                                           #
#> @note   This version was tested for charmm36 forcefield!!                                      #
#                                                                                                 #
###################################################################################################

import numpy as np
import os
import argparse
import shutil


parser = argparse.ArgumentParser(description="""Let us create a user contact.""")
parser.add_argument("-t","--tpr",  help="Gromacs TPR file",      required=True)
parser.add_argument("-g","--gmx",  help="GMX executable",        required=True)
parser.add_argument("-d","--dump", help="Dump file for TPR",     required=True)
parser.add_argument("-o","--out",  help="Output parameter file", required=True)

args = parser.parse_args()
tpr_file = args.tpr
gmx_exec = args.gmx
dmp_file = args.dump
out_file = args.out

term_size = 105
print('\n|' + '=' * term_size + '|\n')

def Check_tpr( tpr_in ):
    if os.path.isfile(tpr_in):
        print(f"\t\tTPR file '{tpr_in}' found\n")
    else:
        print(f"\t\tTPR file '{tpr_in}' not found! Aborting...\n")
        exit()

def Check_exec( exec_name ):
    path = shutil.which(exec_name) 
    if path is None:
        print(f"\t\tNo executable found for command {exec_name}! Aborting...\n")
        exit()
    else:
        print(f"\t\tExecutable {exec_name} found at {path}\n")

print("\tTPR filename read: %s\n" % tpr_file)
Check_tpr(tpr_file)

print("\tGMX executable prefix: %s\n" % gmx_exec)
Check_exec(gmx_exec)

gmxdump = gmx_exec + " dump -s " + tpr_file + " > " + dmp_file + " 2> /dev/null"

os.system( gmxdump )

print("\tTPR parameters dump to file: %s \n" % dmp_file)

label  = []
atype  = []
utype  = []
charge = [] 
nb_c6  = []
nb_c12 = []
nb_sig = []
nb_eps = []
sig_tmp = []
eps_tmp = []

arr_dic = dict( 
        Labels  = label,
        Atypes  = atype,
        Utypes  = utype,
        Charges = charge,
        NBs_c6  = nb_c6,
        NBs_c12 = nb_c12,
        NBs_sig = nb_sig,
        NBs_eps = nb_eps
    )

flist = open(dmp_file).readlines()

parsing = False
for line in flist:
    if line.startswith("topology:"):
        parsing = True
    if parsing:
        if "#atoms" in line:
            n_atoms = int(line.split("=",1)[1])
        if "atnr" in line:
            n_types = int(line.split("=",1)[1])
        if "LJ_SR" in line:
            c6_tmp  = float(line.split("c6=",1)[1].split(",",1)[0])
            c12_tmp = float(line.split("c12=",1)[1].split(",",1)[0])
            nb_c6.append(c6_tmp)
            nb_c12.append(c12_tmp)

    elif line.startswith("cmap"):
        parsing = True
    if parsing:
        if "q=" in line:
            q_tmp  = float(line.split("q=",1)[1].split(",",1)[0])
            charge.append(q_tmp)
        if "atom" in line:
            if "name" in line:
                label_tmp  = line.split("name=\"",1)[1].split("\"",1)[0] 
                label.append(label_tmp)

        if "type" in line:
            if "name" in line:
                atype_tmp  = line.split("name=\"",1)[1].split("\"",1)[0] 
                atype.append(atype_tmp)

        if "atom" in line:
            if "type" in line:
                utype_tmp  = line.split("type=",1)[1].split(",",1)[0] 
                utype.append(utype_tmp)

arr_dic["NBs_c6"] = nb_c6
arr_dic["NBs_c12"] = nb_c12
arr_dic["Charges"] = charge
arr_dic["Labels"] = label
arr_dic["Atypes"] = atype
arr_dic["Utypes"] = utype

for k in range(0, n_types * n_types):
    i = k // n_types + 1  # Calculate i (row)
    j = k % n_types + 1   # Calculate j (column)

    sig = 10 * ( arr_dic["NBs_c12"][k] / arr_dic["NBs_c6"][k] ) ** (1/6)
    eps = ( ( arr_dic["NBs_c6"][k] ** 2 ) / ( 4 * arr_dic["NBs_c12"][k] ) ) 

    if i == j:
        sig_tmp.append(sig)
        eps_tmp.append(eps)

for i in range(0,n_atoms):
    idx = int(arr_dic["Utypes"][i])

    nb_sig.append(sig_tmp[idx])
    nb_eps.append(eps_tmp[idx])

arr_dic["NBs_sig"] = nb_sig
arr_dic["NBs_eps"] = nb_eps

############################################################################################################################

f = open("parameters.txt", "w")

f.write("%4s\t%12s\t%12s\t%12s\n" % ( "###", "chrg", "sig_(A)", "eps_(kJ/mol)"))

for i in range(n_atoms):
    f.write("%4s\t%12.6f\t%12.6f\t%12.6f\n" % ( arr_dic["Labels"][i], (arr_dic["Charges"][i]), arr_dic["NBs_sig"][i], arr_dic["NBs_eps"][i]))

f.close()

print("\tParameters written to file: %s\n" % out_file)
print("\t\tNumber of atoms: %i\n" % n_atoms) 
print("\t\tNumber of unique atomtypes: %i\n" % n_types) 
print("\t\tTotal charge: %12.9f\n" % np.sum(arr_dic["Charges"]))

print('|' + '=' * term_size + '|')
