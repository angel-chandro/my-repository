import os
import numpy as np
import time

dir="/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs/"
ffo='/home/chandro/output/Dhalos_ROCKSTAR_test/descendants.txt'
haloid=[]
did=[]
sh=[]
sd=[]
toth=0
for file in os.listdir(dir):
    if file.endswith("list"):
        t1=os.path.join(file)
        t2=t1.split('.')[0]
        t3=t2.split('_')[1]
        if int(t3)>28 and int(t3)<34:
            ff=open(dir+file,'r')
            print('File:',file)
            start_time = time.time()
            for line in ff:
                a=line.split()[0]
                if a.isdigit() == True:
                    haloid.append(int(line.split()[0]))
                    did.append(int(line.split()[1]))
                    sh.append(int(t3))
                    sd.append(int(t3)+1)
                    toth+=1
            ff.close()
            print('Execution time:',time.time()-start_time,'s')

A=np.transpose(np.matrix([haloid,sh,did,sd]))
ff = open(ffo,'w')
header=str(toth)+' haloes'
np.savetxt(ffo,A,fmt=['%10i','%3i','%10i','%3i'],header=header)
ff.close()
print('Output: ',ffo)
