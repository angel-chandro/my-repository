cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Location of F90_Util library
include_directories(${CMAKE_BINARY_DIR}/F90_Util)
link_directories(${CMAKE_BINARY_DIR}/F90_Util)

# Location of MPI_Util library
include_directories(${CMAKE_BINARY_DIR}/MPI_Util)
link_directories(${CMAKE_BINARY_DIR}/MPI_Util)

# Location of Read_Subfind library
include_directories(${CMAKE_BINARY_DIR}/Read_Subfind)
link_directories(${CMAKE_BINARY_DIR}/Read_Subfind)

add_executable(build_trees merger_trees.f90 subhalo_data.f90
               dhalo_descendants.f90 metadata.f90 subhalo_ids.f90
               dhalo_index.f90 subhalo_interpolation.f90
               enclosing_subhalos.f90 subhalo_trees.f90
               hbt_reader.f90 nifty_reader.f90 ROCKSTAR_reader.f90
	       VELOCIraptor_reader.f90 read_subhalos.f90 sub_io.f90
               host_halos.f90 remerging.f90 main.f90
               splitting.f90 tree_files.f90 redistribute.f90 )

target_link_libraries(build_trees read_subfind f90_util mpi_util ${HDF5_WRAPPER_LIBS} ${MPI_Fortran_LIBRARIES})

install (TARGETS build_trees DESTINATION bin)
