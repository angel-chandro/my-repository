import os
import numpy as np
import time

# Open ROCKSTAR output halo files (in this case /data5/UNITSIM/1Gpc_2048/FixAmp_001 simulation)
for file in os.listdir("/data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs"):
    # Only consider halo files (end with word "list")
    if file.endswith("list"):
        t1=os.path.join(file)
        t2=t1.split('.')[0] # "out_{snapshot}" to name new files created
        print('File:',file)
        command='/home/gustavo/bin/find_parents /data5/UNITSIM/1Gpc_2048/FixAmp_001/ROCKSTAR/outputs/'+file+' > '+t2+'_parents.list'
        start_time = time.time()
        # Run find_parents command: /path/to/find_parents /path/to/ROCKSTAR/outfiles/out_{snapshot}.list > /path/where/to/save/output/obtained/out_{snapshot}_parents.list
        os.system(command)
        print('Execution time:',time.time()-start_time,'s')

        
