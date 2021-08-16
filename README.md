# README

Dhalos_ROCKSTARhaloes directory: contains files related to running ROCKSTAR haloes on DHaloes (or at least try it) that have been created or modified. The modifications and additions are briefly explained, but they can be observed more accurately in my DHalos repository here [Link](https://github.com/angel-chandro/DHalos) (they can have comments attached). The different files are:

* DHalos/ROCKSTAR.test.param: parameter file for ROCKSTAR haloes that has been created.
Modifications: subfind_format = ROCKSTAR, ROCKSTAR_path included (where you can find ROCKSTAR halo catalogues), ROCKSTAR cosmological parameters and particle mass, desc_file = none (not necessary with ROCKSTAR) and hydrorun = .false

* DHalos/Build_Trees/ROCKSTAR_reader.f90: file created to read ROCKSTAR haloes so that DHalos makes the halo trees. Copied from nifty_reader.f90, but functions names have been changed to read_ROCKSTAR and ROCKSTAR_find_main_progenitors. But only read_ROCKSTAR function has been internally varied: "nsubd", "haloid", "descid", "dsnap" variables have been removed because it doesn't need descendants file; modified halo properties just as they appeared in ROCKSTAR halo catalogues; when reading halo catalogue it is necessary to skip header but in this case header has more than 1 line (not as nifty catalogue) so modifications have been made to solve this including "nhead" counter in order to count header lines; also as ROCKSTAR halo ID can be equal to 0 some conditions have been changed.
Still not finished, need to check the halo ID format necessary to proccess ROCKSTAR haloes (search information about DHalos connected with IDCorr variable, but it seems that this IDCorr is necesary in halo IDs and descendant halo IDs because nifty_reader.f90 included it and nifty ID format is almost the same as ROCKSTAR). So in principle, IDCorr must keep as in "nifty_reader.f90" file and add it to descendants ID in this case (maybe information here [Link](http://gavo.mpa-garching.mpg.de/Millennium/pages/help/HelpSingleHTML.jsp#identifiers)).

* DHalos/Build_Trees/tree_files.f90, DHalos/Build_Trees/subhalo_data.f90: "sfactor" variable that refers to scale factor included (because redshift_list file has the format: snapshot(0), scalefactor(1), redshift (2) in 3 different columns).

* DHalos/Build_Trees/merger_trees.f90: ROCKSTAR_reader and subfind_format = ROCKSTAR included.

* DHalos/Build_Trees/CMakeLists.txt: ROCKSTAR_reader included.

* ROCKSTAR_parent_test.py: python file used to run find_parents tool over different ROCKSTAR halo list outputs (this find_parents tool adds a column with "pid": parent ID). Because the current strategy is to create the halo trees based on ROCKSTAR halo catalogues and to do that, it is necessary to know if the halo is a host halo ("pid" = -1) or is a subhalo that belongs to a host halo ("pid" = host halo ID).
The other strategy that has been left aside so far is to run DHalos over ConsistentTrees halo catalogues, what is a possibility since "pid" (different from the one above) and "upid" properties only differs in a 1% of the halos (here we have a 3 level hierarchy, not 2 level). The property to distinguish host haloes and subhaloes would be "upid" in this case.

* write_redshift_list.py: python file used to create redshift list file with 3 columns as explained above. The format is: snapshot(0), scalefactor(1), redshift (2).

Useful information about the project can be found in different links below:
* DHalos info: [Link](https://arxiv.org/pdf/1311.6649.pdf)
* ROCKSTAR info: [Link](https://arxiv.org/pdf/1110.4370.pdf), [Link](https://arxiv.org/pdf/1110.4372.pdf)
* UNIT info: [Link](https://arxiv.org/pdf/1811.02111.pdf)
* SHARK info: [Link](https://arxiv.org/pdf/1807.11180.pdf)
* SAMs info: [Link](https://arxiv.org/pdf/1412.2712.pdf), [Link](https://arxiv.org/pdf/astro-ph/0610031.pdf)
* ROCKSTAR concentration info: [Link](https://arxiv.org/pdf/2007.09012.pdf) (4.3 section)
* Calibration model info: [Link](https://arxiv.org/pdf/2103.01072.pdf)
* Approximations made on baryonic matter info: [Link](https://arxiv.org/pdf/1804.03097.pdf)
* Differences ROCKSTAR-VELOCIraptor in DHalos info: [Link](https://arxiv.org/pdf/2106.12664.pdf)