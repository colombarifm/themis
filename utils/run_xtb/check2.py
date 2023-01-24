en = []
doub = []
en_triple = []

with open("teste.dat", "r") as infile:
   readf = infile.readlines()
   for line in readf:
      energy = float(line.strip("'"))
      en.append(energy)
      doubled = energy * 2
      doub.append(doubled)
      print(f'energy is: ', energy, 'doubled is:', doubled)



print(f'energy is: ', en, 'doubled is:', doub)

