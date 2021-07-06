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
