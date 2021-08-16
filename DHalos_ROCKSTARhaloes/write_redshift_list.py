from numpy import *
import os
import numpy as np

scalefactor=[]
redshift=[]
nsnapshot=[]
for file in os.listdir("/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs"):
    if file.endswith("list"):
        t1=os.path.join(file)
        t2=t1.split('_')[1]
        t3=t2.split('.l')[0]
        nsnapshot.append(int(t3))
        path="/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs/" + t1
        f=open(path)
        for i,line in enumerate(f):
            if i==1:
                scalefactor.append(float(line.split('= ')[1]))
                redshift.append(1/float(scalefactor[-1])-1)
            elif i>1:
                break
        f.close()
        
a1=sorted(np.array(scalefactor))
a2=sorted(np.array(redshift),reverse=True)
print(a2)
a3=sorted(np.array(nsnapshot))
print(a3)
A=np.transpose(np.matrix([a3,a1,a2]))
print(A)
header='snapshot(0) scalefactor(1) redshift(2)'
np.savetxt("redshift_list.txt",A,fmt=['%d','%.5f','%.6f'],header=header)
