import subprocess
import sys
import os

# tree = sys.argv[1]

# precision correspond to 64, 128, 256, 512 & 1024 bit precision
#precision = [28, 39, 78, 155, 309]

tree = 1
prec = sys.argv[1]
path = "/Users/serbanstan/git/TreePCG/graphs/rand_1000_u1000/"

print '*****************************************************************************************************************************'

#for tree in range(4,6):
#	for prec in precision:
#		subprocess.call(str('python worker.py ' + str(prec) + ' ' + str(tree)), shell=True)



subprocess.call(str('python worker.py ' + str(prec) + ' ' + str(tree) + ' ' + path), shell=True)

