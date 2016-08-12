import subprocess
import sys
import os

# tree = sys.argv[1]

# precision correspond to 64, 128, 256, 512 & 1024 bit precision
precision = [28, 39, 78, 155, 309]

print '*****************************************************************************************************************************'

for tree in range(1,6):
	for p in precision:
		subprocess.call(str('python worker.py ' + str(p) + ' ' + str(tree)), shell=True)

