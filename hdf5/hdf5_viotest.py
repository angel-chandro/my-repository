#HDF5 DOCUMENTATION FOR PYTHON
import h5py
import numpy as np

#creating file
f = h5py.File('myfile.hdf5','w')
dset = f.create_dataset("mydataset",(100,),dtype='i')

#open file for reading
f = h5py.File('myfile.hdf5','r')
#show keys file
print(list(f.keys()))
#examine dataset object
dset = f['mydataset']
print(dset.shape) #dataset shape
print(dset.dtype) #dataset data type
print(dset.name) #dataset name
print(f.name) #root group name

#writing data
dset[...] = np.arange(100)
print(dset[0])
print(dset[10])
print(dset[0:100:10])

#creating group
f = h5py.File('myfile.hdf5','a')
grp = f.create_group("subgroup")
#creating dataset in new group
dset2 = grp.create_dataset("another_dataset",(50,),dtype='f')
print(dset2.name)

#creating dataset (not intermediate group)
dset3 = f.create_dataset('subgroup2/dataset_three',(10,),dtype='i')
print(dset3.name)

print(list(f.keys())) #print objects names directly attached to the group
print(list(f.values())) #print objects contained in the gropu
print(list(f.items())) #print (name,value) for objects directly attached to the group

#VIOLETA WEBPAGE
#creating the file
thefile = 'myfilevio.hdf5'

#writing into the file
hf = h5py.File(thefile,'w')

#creating dataset
#Header
head = hf.create_dataset('header',(1,))
#adding attributes
head.attrs[u'volume']       = 1000.
head.attrs[u'units_volume'] = u'(Mpc/h)**3'
print(head.attrs[u'volume'])
print(head.attrs[u'units_volume'])

#generating data in group
#Create the edges and mid points of bins
step = 1
edges = np.array(np.arange(9.,16.,step))
mhist = edges[1:]-0.5*step
#Create some random (integer) positions
pos = np.random.randint(10,size=(3,100))
#Storing data in a group
hfdat = hf.create_group('data')

#Adding label to generated data
#creating 'mass' dataset inside 'data' group
hfdat.create_dataset('mass',data=mhist)
hfdat['mass'].dims[0].label = 'Mass (Msun/h)'
print(hfdat['mass'])
print(hfdat['mass'].dims[0])
print(hfdat['mass'].dims[0].label)
print(hfdat['mass'][:])
#creating 'pos' dataset inside 'data' group
hfdat.create_dataset('pos',data=pos)
hfdat['pos'].dims[0].label = 'x,y,z (Mpc/h)'
print(list(hfdat.items()))

#Close the output file
hf.close()

#Reading file
f = h5py.File(thefile,'r')
#accesing information in the groups
header = f['header']
data = f['data']

#reading header
print(list(header.attrs.items())) #list header attributes 

print(header.attrs['volume'])
boxsize = header.attrs['volume']**(1/3)
print(boxsize)

#reading data
mass = data['mass'][:] #or mass = f['data/mass'][:]
print(mass)
pos = data['pos'][:] #or pos = f['data/pos'][:]
print(pos)

mass5 = f['data/mass'][5:10] #or mass5 = data['mass'][5:10]

#read labels
print([dim.label for dim in data['pos'].dims])
print(len(data['pos'].dims))
print(data['pos'].dims[0].label)

#close again file
f.close()
