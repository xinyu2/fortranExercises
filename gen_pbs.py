import os
import time
import numpy as np
cores = [1,2,4,8,16,32,64,128,256,512]
wall = '00:10:00'
pbs = []

for c in cores:
    node = int(np.ceil(c/32.))
    filename = "milagro_c%d.pbs"%(c)
    node = str(node)
    if c < 32:
        ppn = c
    else :
        ppn = 32
    t = open('runtime.txt', 'w')
    t.write("Job_ID \t Cores \t Time_av")
    f = open(filename, 'w')
    f.write("""#!/bin/bash
### set the number of nodes
### set the number of PEs per node
#PBS -l nodes="""+node+""":ppn="""+str(ppn)+ """:xe
#PBS -q normal
### set the wallclock time
#PBS -l walltime="""+wall+"""
### set the job name
#PBS -N milagro_c"""+str(c)+"""
### set the job stdout and stderr
#PBS -e c"""+str(c)+""".${PBS_JOBID}.err
#PBS -o c"""+str(c)+""".${PBS_JOBID}.out

cd $PBS_O_WORKDIR

aprun -n """+str(c)+""" ./loop 
RUN="$PBS_JOBID"
T_STEP = $(grep -i "transport time" ${PBS_JOBID}.OU | awk '{sum += $6; count +=1} END {print sum/count}')
echo "$RUN\t"""+str(c)+"""\t$T_STEP\t" >> runtime.txt
""")

    pbs.append(filename)

for sub in pbs:
    print(sub)
    os.system("qsub %s"%sub)
    time.sleep(1)
