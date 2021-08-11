# README

Dhalos_ROCKSTARhaloes directory: contains files related to running ROCKSTAR haloes on DHaloes (or at least try it). The different files are:

* DHalos/ROCKSTAR.test.param: parameter file for ROCKSTAR haloes that has been created.
* DHalos/Build_Trees/ROCKSTAR_reader.f90: file created to read ROCKSTAR haloes so that DHalos makes the halo trees. Still not finished, need to check the halo ID format necessary to proccess ROCKSTAR haloes (search information on DHalos papers connected with IDCorr variable)
* DHalos/Build_Trees/tree_files.f90, DHalos/Build_Trees/subhalo_data.f90, DHalos/Build_Trees/merger_trees.f90, DHalos/Build_Trees/CMakeLists.txt: files that have been modified.

* ROCKSTAR_parent_test.py: python file used to run find_parents tool over different ROCKSTAR halo list outputs (this find_parents tool adds a column with "pid": parent ID). Because the current strategy is to create the halo trees based on ROCKSTAR halo catalogues and to do that, it is necessary to know if the halo is a host halo ("pid" = -1) or is a subhalo that belongs to a host halo ("pid" = host halo ID).
The other strategy that has been left aside so far is to run DHalos over ConsistentTrees halo catalogues, what is a possibility since "pid" (different from the one above) and "upid" properties only differs in a 1% of the halos (here we have a 3 level hierarchy, not 2 level). The property to distinguish host haloes and subhaloes would be "upid" in this case.
