import subprocess
import sys
import time
import os
OFN = "cp2kO.inp"
EFN = "cp2kE.inp"

print("PYTHON: opened")

def CP2K(OFN,EFN,typ,CN,geomXYZ):
    print("Python: method")
    #OFN is optimization file, EFN is energy file, geomXYZ is the geometry XYZ file name, CN is calculation number
    com = 'mkdir ' + str(CN)
    subprocess.Popen(com,shell=True)
    #print("Python: mkdir opened")
    runcom = ''
    if typ == "O":
        com = 'cp ' + OFN + ' ' + str(CN) + '/' + OFN
        subprocess.Popen(com, shell=True)
        runcom = 'cd ' + str(CN) + '\nnohup cp2k.psmp -i ' + OFN + ' -o cp2k.out &' 
    elif typ == "E":
        #print("Python: inside E block")
        com = 'cp ' + EFN + ' ' + str(CN) + '/' + EFN
        subprocess.Popen(com, shell=True)
        runcom = 'cd ' + str(CN) + '\nnohup cp2k.psmp -i ' + EFN + ' -o cp2k.out &'
    com = 'cp ' + geomXYZ + ' ' + str(CN) + '/' + geomXYZ
    subprocess.Popen(com, shell=True)
    #print("Python: escaped Popen 2")
    
    #submit the command
    subprocess.Popen(runcom,shell=True)
    print("Python escaped main popen")
    #wait for it to finish

    finished = False
    fail= False
    import time

    start_time = time.time()

    started = False
    while not started:
        if (time.time() - start_time) > 120:
            print("WARNING: cp2k did not start")
            started = True
            fail = True
        #for some reason this keeps triggering even though cp2kout is opened
        if os.path.exists(str(CN) + "/cp2k.out"):
            started = True
            fail = False
            print("PYTHON: cp2k started")
        #else:
        #    print("file on path",str(CN) + "/cp2k.out"," not yet started")
    if fail == False:
        while not finished:
            with open(str(CN) + "/cp2k.out", 'r') as file:
                for line in file:
                    if line[:41] == " The number of warnings for this run is :":
                        finished = True
                    elif "[ABORT]" in line:
                        finished = True
                        fail = True
                #    else:
                #        print("line: ",line)
                file.close()
    print("PYTHON: cp2k has finished")
    #the calculation is finished, analyze.
    with open( str(CN) + "/cp2k.out", 'r') as file:
        for line in file:
            if line[:47] == " ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):":
                energy = line[47:].strip()
            elif "[ABORT]" in line:
                fail = True
            #else:
            #    print("line: ",line)
        file.close() 
    if fail == True:
        if typ == "O":
            file = open(str(CN) + '/opt.xyz', 'w')
            try:
                file.write("-1")
            finally:
                file.close()
        elif typ == "E":
            file = open( str(CN) + '/e.txt', 'w')
            try:
                file.write("F")
            finally:
                file.close()
    else:
        if typ == "O":
            com = 'mv ' + str(CN) + '/*-pos-1.xyz ' + str(CN) + '/outxyz.xyz'
            subprocess.run(com,shell=True)

            with open( str(CN) + "/outxyz.xyz", 'r') as file:
                lines = []
                natoms = 0
                nl = False
                skipLine = False
                for line in file:
                    if natoms == 0:
                        natoms = line.strip()
                        nl = True
                    elif line.strip() == natoms:
                        lines = []
                        nl = True
                    else:
                        nl = False
                    if nl:
                        lines.append(natoms)
                        #lines.append('Frame 0 ,Energy: ' + energy + ' , time: -1')
                    elif line[:4].strip() == "i =".strip():
                        lines.append('Frame 0 ,Energy: ' + energy + ' , time: -1')
                    else:
                        lines.append(line)
                file.close()  

            file = open(str(CN) + '/opt.xyz', 'w')
            try:
                for line in lines:
                    file.write(line.strip() + '\n')
            finally:
                file.close()
        elif typ == "E":
            file = open( str(CN) + '/e.txt', 'w')
            try:
                file.write(energy)
            finally:
                file.close()

command = sys.argv[1]
if command == "cancel":
    com = 'pkill cp2k.psmp'
    subprocess.run(com,shell=True)
else:
    CP2K(OFN,EFN,sys.argv[1],sys.argv[2],"input.xyz")

