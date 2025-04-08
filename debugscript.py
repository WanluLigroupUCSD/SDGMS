import subprocess
import sys
import time
import os
import random
OFN = "cp2kO.inp"
EFN = "cp2kE.inp"

print("PYTHON: opened")

def CP2K(OFN,EFN,typ,CN,geomXYZ):
    if typ == "O":
        file = open(str(CN) + '/opt.xyz', 'w')
        try:
            with open(str(CN) + '/input.xyz','r') as file2:
                for line in file2:
                    file.write(line.strip()+'\n')
        finally:
            file2.close()
            file.close()
    elif typ == "E":
            file = open( str(CN) + '/e.txt', 'w')
            try:
                file.write(str(random.random() - 0.5))
            finally:
                file.close()

command = sys.argv[1]
if command == "cancel":
    com = 'pkill cp2k.psmp'
    subprocess.run(com,shell=True)
else:
    CP2K(OFN,EFN,sys.argv[1],sys.argv[2],"input.xyz")

