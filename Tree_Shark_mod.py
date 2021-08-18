#import numpy as np
#import math
#import random
from scipy.io import FortranFile
import os
#import numpy as np
#import math
import random
import math
#from numpy import*
#from pylab import*
#import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np



#d = np.loadtxt("/data5/UNITSIM/1Gpc_4096/fixedAmp_001/ROCKSTAR/tree/locations.dat" , comments = '#', dtype = {'names': ('TreeRootID', 'FileID', 'Offset', 'Filename '),'formats': ('i', 'f', 'f','S1')})

a_z = []
x = []
y = []
z = []
vx = []
vy = []
vz = []
ID_h = []
Mvir = []
Rvir = []


a = os.listdir("/data5/UNITSIM/1Gpc_4096/fixedAmp_001/ROCKSTAR/tree/.")
b = len(a)
#print(b)

outfile = open("Halos_tree_DOC.txt",'w')


for j in range( 0 , b , 1 ):

       u = a[j]
       #print(u)
       p = str("tree_")
       u1 = a[j][:5]
       if u1==p:
          print("yes")
          #r = r + 1
          #print(r)
          file10 = open("/data5/UNITSIM/1Gpc_4096/fixedAmp_001/ROCKSTAR/tree/tree_0_0_0.dat")
          #file10=open("/data5/UNITSIM/1Gpc_4096/fixedAmp_001/ROCKSTAR/tree/tree_1_2_7.dat")
          f10 = file10.readlines()
          f101 = len(f10)
          oki = []
          #file10 = open("tree_0_0_0.dat")
          for k in range( 0, f101 , 1):
            #aa = k
             #print(aa)
          #ac = file10.readline()
                ac = f10[k]
                if k==0:
                       print(ac[0:8])
                u2 = ac[0]
                u3 = ac[0:8]
          
                if u2 != "#" and u2 == " " and float(u3) == 1.00000:

                       q1 = " ".join(ac.split())
                       #print(q1)
                       y1 = q1.split(" ")
                       #print(y1)
                       ar = np.array(y1)
                       print(ar[0])
                       #arr = np.transpose(ar)
                       #print(len(ar))
                       #print((ar))
                       #u1 = float(ar[0])
                       #u4 = int((ar[1]))
                       #np.savetxt('test.out'+u,ar.T, fmt = '%s')

                       outfile.write('%f \n' %(float(ar[0])))
                #outfile.write('%f %d %f %f %f %f %f %f %.8e %f \n' %(float(u1),u3,float(ar[2]),float(ar[3]),float(ar[4]),float(ar[5]),float(ar[6]),float(ar[7]),float(ar[8]),float(ar[9]),float(ar[10]),float(ar[11]),float(ar[12]),float(ar[13]),float(ar[14]),float(ar[15]),float(ar[16]),float(ar[17]),float(ar[18]),float(ar[19]),float(ar[20]),float(ar[21]),float(ar[22]),float(ar[23]),float(ar[24]),float(ar[25]),float(ar[26]),int(ar[27]),int(ar[28]),int(ar[29]),int(ar[30]),float(ar[31]),float(ar[32]),float(ar[33]),float(ar[34]),float(ar[35]),float(ar[36]),float(ar[37]),float(ar[38]),float(ar[39]),float(ar[40]),float(ar[41]),float(ar[42]),float(ar[43]),float(ar[44]),float(ar[45]),float(ar[46]),float(ar[47]),float(ar[48]),float(ar[49]),float(ar[50]),float(ar[51]),float(ar[52]),float(ar[53]),float(ar[54]),float(ar[55]),float(ar[56]),float(ar[57]),float(ar[58])))
                #outfile.write((ar))i

                #lista[r-1].append(ac)

                #if (u1) == 0.4309 and u3 == ID_G_i:
                
                   #print("encontrado")
                   #lista[r-1].append(ac)
                   #a_z.append(u1)
                   #x.append(float(ar[17]))
                   #y.append(float(ar[18]))
                   #z.append(float(ar[19]))
                   #vx.append(float(ar[20]))
                   #vy.append(float(ar[21]))
                   #vz.append(float(ar[22]))
                   #ID_h.append(u3)
                   #Mvir.append(float(ar[10]))
                   #Rvir.append(float(ar[11]))
          #print("finitouno")       
#out = np.transpose(np.array([a_z,ID_h,x,y,z,vx,vy,vz,Mvir,Rvir]))
outfile.close()
#np.savetxt("Satellite_Galaxy_Halpha_All_z_1.321_halos.dat", out , header = "a_factor IDHalo x y z vx vy vz Mvir Rvir", fmt = "%f %d %f %f %f %f %f %f %f %f")
#print(lista[0])
#print(len(lista[0]))

#print(lista[1])
#print(len(lista[1]))

"""
for j in range( 0 , len(ID_G) , 1):
            ID_G_i = int(ID_G[j])
            print(ID_G_i)

            for i in range( 0 , len(lista) , 1):
                        ai = lista[i]
                        for k in range( 0 , len(a) , 1):
                                    bi = ai[k]

"""
print("fin")
#print(len(ID_h))

