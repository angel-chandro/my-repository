import os
import numpy as np
import time
import h5py

f=h5py.File("descendants.hdf5","w")
grp0=f.create_group("Snapshots")

dir="/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs/"
ffo='/home/chandro/output/Dhalos_ROCKSTAR_test/descendants.hdf5'
haloid=[]
did=[]
sh=[]
sd=[]
toth=0
hsnap=0
nsnap=0
grpsnap=[]
for file in os.listdir(dir):
    if file.endswith("list"):
        t1=os.path.join(file)
        t2=t1.split('.')[0]
        t3=t2.split('_')[1]
        nsnap+=1
        grpsnap.append(str(nsnap))
        if int(t3)<10:
            pathgroup="Snap_00"+t3
        elif int(t3)<100:
            pathgroup="Snap_0"+t3
        else:
            pathgroup="Snap_"+t3
        grpsnap[-1]=grp0.create_group(pathgroup)
        ff=open(dir+file,'r')
        print('File:',file)
        start_time = time.time()            
        for line in ff:
            a=line.split()[0]
            if a.isdigit() == True:
                hsnap+=1
                haloid.append(int(line.split()[0]))
                did.append(int(line.split()[1]))
                sh.append(int(t3))
                sd.append(int(t3)+1)
                toth+=1
        ff.close()
        dset=grpsnap[-1].create_dataset("ID",(hsnap,),dtype='i')
        dset[:]=haloid
        dset1=grpsnap[-1].create_dataset("Head",(hsnap,),dtype='i')
        dset1[:]=did
        dset2=grpsnap[-1].create_dataset("HeadSnap",(hsnap,),dtype='i')
        dset2[:]=sd
        grpsnap[-1].attrs['NHalos']=hsnap
        did=[]
        haloid=[]
        sh=[]
        sd=[]
        hsnap=0


grp1=f.create_group("Header")
grp1.attrs['NSnaps']=nsnap
f.close()
print('Execution time:',time.time()-start_time,'s')
print('Output: ',ffo)
