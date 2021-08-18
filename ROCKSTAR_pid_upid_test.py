from numpy import *
import os
import numpy as np

pid=[]
upid=[]
ieq=0
ineq=0
for file in os.listdir("/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs/hlists"):
    if file.endswith("list"):
        t1=os.path.join(file)
        t2=t1.split('_')[1]
        t3=t2.split('.l')[0]
        path="/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs/hlists/" + t1
        f=open(path)
        for i,line in enumerate(f):
            if i>=64:
                linemod=line.split(" ")
                linemodf=list(filter(None,linemod))
                #pid.append(int(linemodf[5]))
                #upid.append(int(linemodf[6]))
                if int(linemodf[5])==int(linemodf[6]):
                    ieq=ieq+1
                    #print(ieq)
                else:
                    ineq=ineq+1
                    #print(ineq)
        pid.append(ieq)
        upid.append(ineq)
        ieq=0
        ineq=0
        print(path)
        print(pid)
        print(upid)
        print(100*sum(upid)/(sum(pid)+sum(upid)))

A=np.transpose(np.matrix([pid,upid,100*sum(upid)/(sum(pid)+sum(upid))]))
header='eq(0) neq(1) percsum(2)'
np.savetxt("ROCKSTAR_pid_upid_test.txt",A,fmt=['%d','%d','½f'],header=header)
