from numpy import loadtxt
import os

#lines = loadtxt("teste.dat")
#for i in lines:
#    print(i)

with open('teste.dat') as filein:
    for line in filein:
        print("here is: ", line)


text_file = open("teste.dat", "r")
lines = text_file.readlines()
print(lines)
text_file.close()

print(lines[1])

#print(lines)

#print(fval)

print("error")

data = []
with open('teste.dat') as f:
  for line in f.readlines():
      data.append(lineToData(line.split()))

print(data)

#or i in lines:
#   new_val = float(float(i) * float(2))
#   print(new_val)
#
#rint("ok")
#
#rint(new_val)

#print(new_val[1])
