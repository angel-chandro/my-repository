module tree_files
!
! Module to write out the merger trees
!
  use mpi
  use mpi_types
  use kind_numbers
  use hdf5_wrapper

  implicit none
  private

  ! Callable routines in this module
  public :: write_tree_files

  ! Assignment of files to processors
  integer, dimension(:), allocatable :: num_files_task
  integer, dimension(:), allocatable :: first_file
  integer, dimension(:), allocatable :: last_file

  ! Which processor each file is assigned to
  integer, dimension(:), allocatable :: file_task_id

  ! Type to store tree metadata
  type merger_tree_info_type
     sequence
     integer(kind=int8byte)   :: dhalo_index
     integer(kind=int8byte)   :: subhalo_index
     integer                  :: nnodes
     integer                  :: final_snapnum
     integer                  :: ifile
     integer                  :: itask
     integer                  :: itree_mem
     integer                  :: itree_file
  end type merger_tree_info_type

  ! Compression level (adding 10 indicates shuffle filter should be used)
  integer, parameter :: gzip_level = 18
  
contains

  integer function create_mpi_tree_type()
!
! Register an MPI type for transferring tree metadata
!
    use mpi
    implicit none
    type (merger_tree_info_type), dimension(2) :: tree
    integer(MPI_ADDRESS_KIND) :: base, offset
    integer :: extent
    integer, parameter :: countmax = 100
    integer, dimension(countmax) :: blocklen, disp, typearr
    integer :: i, ierr

    ! Figure out size of subhalo_type (may include padding)
    call MPI_GET_ADDRESS(tree(1), base,   ierr)
    call MPI_GET_ADDRESS(tree(2), offset, ierr)
    extent = int(offset-base,kind(extent))

    ! MPI type consists of a series of blocks with length, type and
    ! displacement for each one
    i = 0

    ! dhalo_index
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INT8BYTE
    call MPI_GET_ADDRESS(tree(1)%dhalo_index, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! subhalo_index
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INT8BYTE
    call MPI_GET_ADDRESS(tree(1)%subhalo_index, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! nnodes
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(tree(1)%nnodes, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! final_snapnum
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(tree(1)%final_snapnum, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! ifile
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(tree(1)%ifile, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! itask
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(tree(1)%itask, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! itree_mem
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(tree(1)%itree_mem, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! itree_file
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(tree(1)%itree_file, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Mark end of type
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_UB
    disp(i)     = extent

    ! Register the new type
    call MPI_TYPE_STRUCT(i, blocklen, disp, typearr, create_mpi_tree_type, ierr)
    call MPI_TYPE_COMMIT(create_mpi_tree_type, ierr)

    return
  end function create_mpi_tree_type


  subroutine assign_files_to_tasks(nfiles, nproc_io)
!
! Decide which MPI tasks will write which files
!
    implicit none
    integer, intent(in) :: nfiles, nproc_io
    ! MPI stuff
    integer :: myid, numprocs, ierr
    ! Loops etc
    integer :: i
    ! Which processors are allowed to do I/O
    logical, dimension(:), allocatable :: can_write 
    integer :: nstep, files_left, tasks_left

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
 
    allocate(num_files_task(0:numprocs-1))
    allocate(first_file(0:numprocs-1))
    allocate(last_file(0:numprocs-1))
    allocate(file_task_id(0:nfiles-1))
    allocate(can_write(0:numprocs-1))

    ! Decide which tasks are allowed to write
    if(nproc_io.ge.numprocs)then
       can_write = .true.
    else
       tasks_left = nproc_io
       can_write = .false.
       nstep = numprocs / nproc_io
       do i = 0, numprocs-1, nstep
          can_write(i) = .true.
          tasks_left = tasks_left - 1
          if(tasks_left.eq.0)exit
       end do
    endif

    ! Decide number of files to write on each task
    num_files_task = 0
    files_left = nfiles
    i = 0
    do while(files_left.gt.0)
       if(can_write(i))then
          num_files_task(i) = num_files_task(i) + 1
          files_left = files_left - 1
       endif
       i = i + 1
       if(i.eq.numprocs)i = 0
    end do

    if(sum(num_files_task).ne.nfiles)call panic("Total number of files is wrong!")

    ! Find first and last file for each task
    do i = 0, numprocs-1, 1
       first_file(i) = sum(num_files_task(0:i-1))
       last_file(i)  = sum(num_files_task(0:i))-1
    end do

    ! For each file, record which task will write it
    file_task_id = -1
    do i = 0, numprocs-1, 1
       file_task_id(first_file(i):last_file(i)) = i
    end do
    if(any(file_task_id.lt.0))call panic("File not assigned to a task!")

    deallocate(can_write)

    return
  end subroutine assign_files_to_tasks


  subroutine find_tree_position_in_file()
    !
    ! For each merger tree, find its position in its assigned 
    ! tree file (1=first tree in file)
    !
    ! This fills in the itree_file component in the merger_tree
    ! array.
    !
    use sort
    use multi_key_sort
    use merger_trees
    implicit none
    ! Which file we're writing
    integer :: ifile, jfile
    ! Sorting index for local trees
    ! This puts the trees in order of which file they should go in
    integer, dimension(:), allocatable :: send_idx
    ! Sorting index for received trees
    integer, dimension(:), allocatable :: recv_idx
    ! MPI stuff
    integer :: MPI_TREE_TYPE
    integer :: myid, numprocs, ierr
    ! Send/receive buffers for tree info
    type (merger_tree_info_type), dimension(:), allocatable :: sendbuf, recvbuf
    integer, dimension(:), allocatable :: nsend, nrecv, soffset, roffset
    ! Loops etc
    integer :: itree, idest, dest_ifile
    integer :: nsend_tot, nrecv_tot, ipass, i

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
    allocate(nsend(0:numprocs-1), nrecv(0:numprocs-1))
    allocate(soffset(0:numprocs-1), roffset(0:numprocs-1))

    MPI_TREE_TYPE = create_mpi_tree_type()

    ! Make sorting index (by output file) for local merger trees
    allocate(send_idx(ntrees))
    call sort_index(ntrees, merger_tree%ifile, send_idx)

    ! Loop over files to write
    do jfile = 0, maxval(num_files_task)-1, 1 

       ! Determine which file we're looking at on this iteration (may be none)
       ifile = first_file(myid) + jfile
       if(ifile.gt.last_file(myid))ifile = -1

       ! Copy info about trees to be written to this file to send buffer
       ! in order of which file (and task) they're to be written to.
       ! Also accumulate number of trees to send to each task.
       do ipass = 1, 2, 1
          nsend_tot = 0
          nsend     = 0
          do itree = 1, ntrees, 1
             ! Find out where this tree is to be sent to
             idest = file_task_id(merger_tree(send_idx(itree))%ifile)
             ! If the destination task is writing this file on this
             ! iteration, copy tree to send buffer
             dest_ifile = first_file(idest)+jfile
             if(dest_ifile.le.last_file(idest).and.dest_ifile.eq.merger_tree(send_idx(itree))%ifile)then
                nsend_tot    = nsend_tot    + 1
                nsend(idest) = nsend(idest) + 1
                if(ipass.eq.2)then
                   sendbuf(nsend_tot)%dhalo_index   = merger_tree(send_idx(itree))%dhalo_index
                   sendbuf(nsend_tot)%subhalo_index = merger_tree(send_idx(itree))%subhalo_index
                   sendbuf(nsend_tot)%nnodes        = merger_tree(send_idx(itree))%nnodes
                   sendbuf(nsend_tot)%final_snapnum = merger_tree(send_idx(itree))%final_snapnum
                   sendbuf(nsend_tot)%ifile         = merger_tree(send_idx(itree))%ifile
                   sendbuf(nsend_tot)%itask         = myid
                   sendbuf(nsend_tot)%itree_mem     = send_idx(itree)
                   sendbuf(nsend_tot)%itree_file    = -1
                endif
             endif
          end do
          if(ipass.eq.1)then
             allocate(sendbuf(max(nsend_tot,1)))
          endif
       end do

       ! Determine how many trees we're going to receive and allocate receive buffer
       call MPI_ALLTOALL(nsend, 1, MPI_INTEGER, nrecv, 1, MPI_INTEGER, &
            MPI_COMM_WORLD, ierr)
       nrecv_tot = sum(nrecv)
       allocate(recvbuf(nrecv_tot))

       ! Exchange the merger tree metadata
       soffset(0) = 0
       roffset(0) = 0
       do i = 1, numprocs-1, 1
          soffset(i) = soffset(i-1) + nsend(i-1)
          roffset(i) = roffset(i-1) + nrecv(i-1)
       end do
       call MPI_ALLTOALLV(sendbuf, nsend, soffset, MPI_TREE_TYPE, &
            recvbuf, nrecv, roffset, MPI_TREE_TYPE, MPI_COMM_WORLD, ierr)
       !
       ! Recvbuf now contains the metadata for the trees this processor is going to
       ! write to this file.
       !
       ! Determine the position in the file where each tree will be written
       ! Trees are sorted by final dhalo index, then within a dhalo by subhalo index.
       allocate(recv_idx(nrecv_tot))
       call two_key_sort_index(nrecv_tot, recvbuf%dhalo_index, recvbuf%subhalo_index, recv_idx)
       do i = 1, nrecv_tot, 1
          recvbuf(recv_idx(i))%itree_file = i
       end do

       ! Return updated information to the originating processor
       call MPI_ALLTOALLV(recvbuf, nrecv, roffset, MPI_TREE_TYPE, &
            sendbuf, nsend, soffset, MPI_TREE_TYPE, MPI_COMM_WORLD, ierr)
          
       ! Update merger tree array
       do i = 1, nsend_tot, 1
          merger_tree(sendbuf(i)%itree_mem)%itree_file = sendbuf(i)%itree_file
       end do

       deallocate(sendbuf, recvbuf, recv_idx)

       ! Next file
    end do

    ! Done with tree metadata type
    call MPI_TYPE_FREE(MPI_TREE_TYPE, ierr)

    deallocate(send_idx, nsend, nrecv, soffset, roffset)

    ! All trees should have been assigned a position
    if(any(merger_tree%itree_file.lt.0))&
         call panic("find_tree_position_in_file() - Entry in merger_tree%itree_file not set!")

    return
  end subroutine find_tree_position_in_file
  

  character(len=500) function tree_file_name(treedir, tree_basename, ifile)
!
! Generate the name of a tree file
!
    use subhalo_data, only : node_at_snapshot
    implicit none
    character(len=*), intent(in) :: treedir, tree_basename
    integer                      :: ifile
    character(len=10)            :: filenum, snapnum
    integer                      :: isnap
    
    isnap = ubound(node_at_snapshot,1)

    write(filenum, '(1i10)')ifile
    if(isnap.gt.999)then
       write(snapnum,'(i4.4)')isnap
    else
       write(snapnum,'(i3.3)')isnap
    endif
    tree_file_name = &
         trim(treedir)//"/treedir_"// &
         trim(snapnum)//"/"// &
         trim(tree_basename)//"_"// &
         trim(snapnum)//"."// &
         trim(adjustl(filenum))//".hdf5"
    
    return
  end function tree_file_name


  function create_tree_file(treedir, tree_basename, nfiles, ifile, &
       mpart, lbox, remerge) result(ihdf5)
!
! Create a new tree file
!
    use metadata
    use subhalo_data, only : zred, sfactor ! include
    implicit none
    ! Where to put the output
    character(len=*), intent(in) :: treedir, tree_basename
    ! Which file to create
    integer, intent(in) :: ifile, nfiles
    ! Box size and particle mass
    real, intent(in) :: lbox, mpart
    logical, intent(in) :: remerge
    ! Output file handle
    integer :: ihdf5
    ! Internal
    character(len=500) :: fname, dirname
    integer, dimension(:), allocatable :: snapnum
    integer :: ifirst, ilast, i

    ifirst = lbound(zred,1)
    ilast  = ubound(zred,1)
    allocate(snapnum(ifirst:ilast))
    do i = ifirst, ilast, 1
       snapnum(i) = i
    end do
    fname = tree_file_name(treedir, tree_basename, ifile)

    ! Make sure the directory exists
    !if(ifile.eq.0)then
    i = scan(fname, "/", .true.)
    if(i.gt.1)then
       dirname = fname(:i-1)
       call make_directory(dirname)
    endif
    !endif

    call hdf5_create_file(ihdf5, fname)
    call hdf5_create_group(ihdf5, "haloTrees")
    call hdf5_create_group(ihdf5, "treeIndex")
    call hdf5_create_group(ihdf5, "hostIndex")
    call hdf5_create_group(ihdf5, "treeConstruction")
    call hdf5_create_group(ihdf5, "outputTimes")
    call hdf5_write_data(ihdf5, "outputTimes/snapshotNumber", snapnum)
    call hdf5_write_data(ihdf5, "outputTimes/redshift",       zred)
    call hdf5_write_data(ihdf5, "outputTimes/expansion",      sfactor) ! include

    ! Record which file this is in the set
    call hdf5_create_group(ihdf5, "fileInfo")
    call hdf5_write_attribute(ihdf5, "fileInfo/numberOfFiles", nfiles)
    call hdf5_write_attribute(ihdf5, "fileInfo/thisFile",      ifile)

    deallocate(snapnum)

    call write_metadata(ihdf5)
    call hdf5_write_attribute(ihdf5, "simulation/boxSize",      lbox)
    call hdf5_write_attribute(ihdf5, "simulation/particleMass", mpart)

    if(remerge)then
       call hdf5_write_attribute(ihdf5, "/haloTrees/treesAreSelfContained", 1)
    else
       call hdf5_write_attribute(ihdf5, "/haloTrees/treesAreSelfContained", 0)
    endif
    call hdf5_write_attribute(ihdf5, "/haloTrees/treesHaveSubhalos", 1)
    call hdf5_write_attribute(ihdf5, "/haloTrees/velocitiesIncludeHubbleFlow", 0)

    return
  end function create_tree_file


  subroutine write_tree_files(nfiles, treedir, tree_basename, &
       read_wdm, mpart, lbox, remerge, nproc_io, use_total_mass)
    !
    ! Write out the merger tree files
    !
    use sort
    use multi_key_sort
    use merger_trees
    implicit none
    ! How many files to write
    integer, intent(in) :: nfiles, nproc_io
    ! Where to put the output
    character(len=*), intent(in) :: treedir, tree_basename
    logical, intent(in) :: read_wdm
    real,    intent(in) :: mpart, lbox
    logical, intent(in) :: remerge
    logical, intent(in) :: use_total_mass
    ! Which file we're writing
    integer :: ifile, jfile
    ! Sorting index for local trees
    ! This puts the trees in order of which file they should go in
    integer, dimension(:), allocatable :: send_idx
    ! Sorting index for received trees
    integer, dimension(:), allocatable :: recv_idx
    ! MPI stuff
    integer :: MPI_TREE_TYPE
    integer :: myid, numprocs, ierr
    ! Send/receive buffers for tree info
    type (merger_tree_info_type), dimension(:), allocatable :: sendbuf, recvbuf
    integer, dimension(:), allocatable :: nsend, nrecv, soffset, roffset
    ! Loops etc
    integer :: itree, idest, dest_ifile
    integer :: nsend_tot, nrecv_tot, ipass, i
    ! HDF5 file handle
    integer :: ihdf5
    ! Output buffere for integer quantities
    integer, dimension(:), allocatable :: ioutbuf
    ! Number of final hosts in this file
    integer :: nhost
    integer(kind=int8byte) :: current_dhalo
    integer, dimension(:), allocatable :: host_ntrees
    ! Number of files to write on this iteration
    integer :: files_to_write

    ! Decide which tasks write which files
    call assign_files_to_tasks(nfiles, nproc_io)

    ! Decide where to write each tree in each file
    call find_tree_position_in_file()    

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
    allocate(nsend(0:numprocs-1), nrecv(0:numprocs-1))
    allocate(soffset(0:numprocs-1), roffset(0:numprocs-1))

    MPI_TREE_TYPE = create_mpi_tree_type()

    ! Make sorting index (by output file then position in file) for local merger trees
    allocate(send_idx(ntrees))
    call two_key_sort_index(ntrees, merger_tree%ifile, merger_tree%itree_file, send_idx)

    ! Loop until all files have been written
    do jfile = 0, maxval(num_files_task)-1, 1 

       ! Count how many files to write this time
       if(myid.eq.0)then
          files_to_write = 0
          do i = 0, numprocs-1, 1
             if(num_files_task(i).gt.jfile)then
                !write(*,*)"    Rank ", i, " writing"
                files_to_write = files_to_write + 1
             endif
          end do
          write(*,'(a,i6,a,i6,a,i6,a)')" Output iteration ", jfile+1, " of ", maxval(num_files_task), &
               ": writing ", files_to_write, " tree files"
       endif

       ! Determine which file we're writing on this iteration (may be none)
       ifile = first_file(myid) + jfile
       if(ifile.gt.last_file(myid))ifile = -1

       !
       ! Calculate indexing information for this file
       !
       ! Quantities we need:
       ! -------------------
       !
       ! /haloTrees/:
       !   * numberOfTrees: number of "trees" (may actually be branches) in this file
       ! /treeIndex/:
       !   * numberOfNodes: number of nodes in each tree in this file
       !   * firstNode:     offset to the first node in each tree in this file
       !   * finalSnapshot: final snapshot number for each tree
       ! /hostIndex/:
       !   * firstTree:     first tree for each final host
       !   * numberOfTrees: number of trees for each final host
       !
       ! Gather tree metadata for this file
       ! ----------------------------------
       !
       ! Copy info about trees to be written out to send buffer
       ! in order of which file and position in file.
       ! Also accumulate number of trees to send to each task.
       do ipass = 1, 2, 1
          nsend_tot = 0
          nsend     = 0
          do itree = 1, ntrees, 1
             ! Find out where this tree is to be sent to
             idest = file_task_id(merger_tree(send_idx(itree))%ifile)
             ! If the destination task is writing this file on this
             ! iteration, copy tree to send buffer
             dest_ifile = first_file(idest)+jfile
             if(dest_ifile.le.last_file(idest).and.dest_ifile.eq.merger_tree(send_idx(itree))%ifile)then
                nsend_tot    = nsend_tot    + 1
                nsend(idest) = nsend(idest) + 1
                if(ipass.eq.2)then
                   sendbuf(nsend_tot)%dhalo_index   = merger_tree(send_idx(itree))%dhalo_index
                   sendbuf(nsend_tot)%subhalo_index = merger_tree(send_idx(itree))%subhalo_index
                   sendbuf(nsend_tot)%nnodes        = merger_tree(send_idx(itree))%nnodes
                   sendbuf(nsend_tot)%final_snapnum = merger_tree(send_idx(itree))%final_snapnum
                   sendbuf(nsend_tot)%ifile         = merger_tree(send_idx(itree))%ifile
                   sendbuf(nsend_tot)%itask         = myid
                   sendbuf(nsend_tot)%itree_mem     = send_idx(itree)
                   sendbuf(nsend_tot)%itree_file    = merger_tree(send_idx(itree))%itree_file
                endif
             endif
          end do
          if(ipass.eq.1)allocate(sendbuf(max(nsend_tot,1)))
       end do

       ! Determine how many trees we're going to receive and allocate receive buffer
       call MPI_ALLTOALL(nsend, 1, MPI_INTEGER, nrecv, 1, MPI_INTEGER, &
            MPI_COMM_WORLD, ierr)
       nrecv_tot = sum(nrecv)
       allocate(recvbuf(nrecv_tot))

       ! Exchange the merger tree metadata
       soffset(0) = 0
       roffset(0) = 0
       do i = 1, numprocs-1, 1
          soffset(i) = soffset(i-1) + nsend(i-1)
          roffset(i) = roffset(i-1) + nrecv(i-1)
       end do
       call MPI_ALLTOALLV(sendbuf, nsend, soffset, MPI_TREE_TYPE, &
            recvbuf, nrecv, roffset, MPI_TREE_TYPE, MPI_COMM_WORLD, ierr)

       !
       ! Recvbuf now contains the metadata for the trees this processor is going to
       ! write to this file.
       !
       ! Get index to sort trees into output order
       allocate(recv_idx(nrecv_tot))
       call sort_index(nrecv_tot, recvbuf%itree_file, recv_idx)

       if(ifile.ge.0)then
          !
          ! Create the merger tree file
          !
          ihdf5 = create_tree_file(treedir, tree_basename, nfiles, ifile, &
               mpart, lbox, remerge)
          !
          ! Write treeIndex data
          !
          ! treeIndex/numberOfNodes
          allocate(ioutbuf(nrecv_tot))
          do i = 1, nrecv_tot, 1
             ioutbuf(i) = recvbuf(recv_idx(i))%nnodes
          end do
          call hdf5_write_data(ihdf5, "treeIndex/numberOfNodes", ioutbuf, &
               gzip=gzip_level)
          ! treeIndex/firstNode
          if(nrecv_tot.gt.0)ioutbuf(1) = 0
          do i = 2, nrecv_tot, 1
             ioutbuf(i) = ioutbuf(i-1) + recvbuf(recv_idx(i-1))%nnodes
          end do
          call hdf5_write_data(ihdf5, "treeIndex/firstNode", ioutbuf, &
               gzip=gzip_level)
          ! treeIndex/finalSnapshot
          do i = 1, nrecv_tot, 1
             ioutbuf(i) = recvbuf(recv_idx(i))%final_snapnum
          end do
          call hdf5_write_data(ihdf5, "treeIndex/finalSnapshot", ioutbuf, &
               gzip=gzip_level)
          deallocate(ioutbuf)
          !
          ! Write hostIndex data
          !
          ! First count number of unique final dhalos in this file
          nhost = 0
          current_dhalo = -1
          do i = 1, nrecv_tot, 1
             if(recvbuf(recv_idx(i))%dhalo_index.ne.current_dhalo.or.nhost.eq.0)then
                nhost = nhost + 1
                current_dhalo = recvbuf(recv_idx(i))%dhalo_index
             endif
          end do
          ! Count number of trees per final host
          allocate(host_ntrees(nhost))
          host_ntrees = 0
          nhost = 0
          current_dhalo = -1
          do i = 1, nrecv_tot, 1
             if(recvbuf(recv_idx(i))%dhalo_index.ne.current_dhalo.or.nhost.eq.0)then
                nhost = nhost + 1
                current_dhalo = recvbuf(recv_idx(i))%dhalo_index
             endif
             host_ntrees(nhost) = host_ntrees(nhost) + 1
          end do
          ! hostIndex/numberOfTrees
          call hdf5_write_data(ihdf5, "hostIndex/numberOfTrees", host_ntrees, &
               gzip=gzip_level)
          ! hostIndex/firstTree
          allocate(ioutbuf(nhost))
          if(nhost.gt.0)ioutbuf(1) = 0
          do i = 2, nhost, 1
             ioutbuf(i) = ioutbuf(i-1) + host_ntrees(i-1)
          end do
          call hdf5_write_data(ihdf5, "hostIndex/firstTree", ioutbuf, &
               gzip=gzip_level)
          deallocate(ioutbuf, host_ntrees)
            
          ! Write trees per file
          call hdf5_write_attribute(ihdf5, "haloTrees/numberOfTrees", nrecv_tot)

       endif

       ! Write the merger tree data
       call write_tree_data(nsend_tot, sendbuf, nrecv_tot, recvbuf, recv_idx, &
            ihdf5, read_wdm, mpart, use_total_mass)

       ! Close the tree file (if we wrote one)
       if(ifile.ge.0)then
          call add_comments(ihdf5)
          call hdf5_close_file(ihdf5)
       endif
          
       deallocate(sendbuf, recvbuf, recv_idx)

       ! Next file
    end do

    ! Done with tree metadata type
    call MPI_TYPE_FREE(MPI_TREE_TYPE, ierr)

    deallocate(send_idx, nsend, nrecv, soffset, roffset)

    return
  end subroutine write_tree_files
     

  subroutine write_tree_data(nsend_tot, trees_send, &
       nrecv_tot, trees_recv, recv_idx, &
       ihdf5, read_wdm, mpart, use_total_mass)
    !
    ! Write out the merger trees for the current set of files
    ! (zero or one files per processor, ifile=-1 on tasks which
    ! are not writing a file on this iteration).
    !
    ! Input:
    !
    ! nsend_tot: total number of local trees we need to send 
    !            for output on this iteration
    !
    ! trees_send: array of merger_tree_info_type sorted by ordering
    !             in the output (i.e. by file then position in file).
    !             Describes the merger trees stored locally which are
    !             to be written on this iteration.
    !
    ! nrecv_tot: total number of merger trees we need to receive and
    !            write to file ifile
    !
    ! trees_recv: array of merger_tree_info_type (NOT sorted).
    !             Describes the merger trees to be written to file
    !             ifile, most of which will not be stored locally.
    !
    ! recv_idx: sorting index to sort trees_recv into output order
    !
    ! ihdf5: hdf5 wrapper file handle for the output file
    !
    ! ifile: index of the file to be written by this task (-1 for none)
    !
    use merger_trees
    use subhalo_data

    implicit none
    ! Input parameters
    integer, intent(in) :: nsend_tot, nrecv_tot
    type (merger_tree_info_type), dimension(:), intent(in) :: trees_recv
    type (merger_tree_info_type), dimension(:), intent(in) :: trees_send
    integer, dimension(:), intent(in) :: recv_idx
    integer, intent(in) :: ihdf5
    logical, intent(in) :: read_wdm
    real,    intent(in) :: mpart
    logical, intent(in) :: use_total_mass
    ! Range of trees we need to send to each task
    type task_info
       integer :: first_tree
       integer :: last_tree
       integer :: itree
       integer :: isub
       type(subhalo_tree_node), pointer :: next_node
    end type task_info
    type(task_info), dimension(:), allocatable :: task
    ! Next tree and subhalo to import
    integer :: itree_import, jtree_import, ktree_import
    integer :: isub_import, jsub_import, ksub_import
    ! MPI stuff
    integer :: myid, numprocs, ierr, status(MPI_STATUS_SIZE)
    integer :: ptask, ngrp
    ! Maximum number of subhalos to write per iteration
    integer, parameter :: nsubmax = 100000
    ! Loops etc
    integer :: i, idest, itask
    integer :: nleft, nleft_max, nwrite, nwrite_left, nnodes, nodes_left
    integer :: n_to_buffer, offset, nsub
    ! Number of subhalos to receive
    integer, dimension(:), allocatable :: nsub_recv, nsub_send
    integer, dimension(:), allocatable :: soffset, roffset
    ! Send buffer for subhalos
    type(subhalo_type), dimension(:), allocatable :: sendbuf, recvbuf, iobuf
    ! Pointer to the base of a tree
    type(subhalo_tree_node), pointer :: base
    integer :: file_offset

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

    allocate(nsub_recv(0:numprocs-1), nsub_send(0:numprocs-1))
    allocate(soffset(0:numprocs-1), roffset(0:numprocs-1))

    ! Initialise task info:
    ! - Range of trees to send to each task
    ! - Index of current tree and pointer to (and index of) next subhalo to write out
    allocate(task(0:numprocs-1))
    task%first_tree = 0
    task%last_tree  = -1
    do i = 1, nsend_tot, 1
       idest = file_task_id(trees_send(i)%ifile)
       if(task(idest)%first_tree.le.0)task(idest)%first_tree = i
       task(idest)%last_tree = i
    end do
    do i = 0, numprocs-1, 1
       task(i)%itree = task(i)%first_tree
       task(i)%isub  = 1       
       nullify(task(i)%next_node)
    end do

    !
    ! Loop over blocks of subhalos to write to the output file
    !
    nleft = sum(trees_recv(1:nrecv_tot)%nnodes) ! Total subhalos to write to current file
    itree_import = 1 ! Next tree to import 
    isub_import  = 1 ! Next subhalo in the tree to import
    file_offset  = 0 ! Position in file to write next set of subhalos
    do while(.true.)
       
       ! Decide how many subhalos to write on this iteration
       nwrite = min(nleft, nsubmax)
       
       ! Calculate how many subhalos we need to import from each
       ! other task
       nsub_recv    = 0
       nwrite_left  = nwrite
       jtree_import = itree_import
       jsub_import  = isub_import
       do while(nwrite_left.gt.0)
          ! Check if we're writing all remaining subhalos of the next tree
          nnodes = trees_recv(recv_idx(jtree_import))%nnodes
          itask  = trees_recv(recv_idx(jtree_import))%itask
          nodes_left = nnodes - jsub_import + 1
          if(nwrite_left.ge.nodes_left)then
             ! Writing the rest of this tree
             nsub_recv(itask) = nsub_recv(itask) + nodes_left
             nwrite_left  = nwrite_left - nodes_left
             ! Go to first subhalo in next tree
             jtree_import = jtree_import + 1
             jsub_import  = 1 
          else
             ! Writing part of this tree
             nsub_recv(itask) = nsub_recv(itask) + nwrite_left
             ! Advance to remaining part of this tree
             jsub_import = jsub_import + nwrite_left
             nwrite_left = 0
          endif   
       end do

       ! Find out how many subhalos we need to send to each task
       call MPI_ALLTOALL(nsub_recv, 1, MPI_INTEGER, nsub_send, 1, MPI_INTEGER, &
            MPI_COMM_WORLD, ierr)
       
       !
       ! Exchange subhalo data
       !
       ! This effectively an MPI_ALLTOALLV, but here we implement it using
       ! point to point communication so that we don't have to allocate send 
       ! buffers for all destinations simultaneously.
       !
       ptask = 0
       do while(2**ptask.lt.numprocs)
          ptask = ptask + 1
       end do
       allocate(recvbuf(sum(nsub_recv)))
       roffset(0) = 0
       do i = 1, numprocs-1, 1
          roffset(i) = roffset(i-1) + nsub_recv(i-1)
       end do
       allocate(sendbuf(nsubmax))
       do ngrp = 0, (2**ptask)-1, 1
          idest = ieor(myid, ngrp)
          if(idest.lt.numprocs)then
             ! Populate send buffer for this destination
             n_to_buffer = nsub_send(idest)
             offset      = 1
             do while(n_to_buffer.gt.0)
                ! Get pointer to base node in next tree to send data from
                base => merger_tree(trees_send(task(idest)%itree)%itree_mem)%ptr
                ! Retrieve the subhalo merger tree data
                ! Returns the rest of the current tree or n_to_buffer subhalos,
                ! whichever is smaller. 
                call get_tree_data(base, task(idest)%next_node, n_to_buffer, sendbuf(offset:), nsub)
                ! Advance to next tree if we read all remaining subhalos in this one
                if(.not.associated(task(idest)%next_node))then
                   task(idest)%itree = task(idest)%itree + 1
                endif
                ! Advance to next part of send buffer
                offset = offset + nsub
                n_to_buffer = n_to_buffer - nsub
             end do
             ! Exchange subhalos
             if(nsub_send(idest).gt.0.and.nsub_recv(idest).gt.0)then
                call MPI_SENDRECV(sendbuf(1:nsub_send(idest)), nsub_send(idest), &
                     SUBHALO_MPI_TYPE, idest, 0, &
                     recvbuf(roffset(idest)+1:roffset(idest)+nsub_recv(idest)), &
                     nsub_recv(idest), SUBHALO_MPI_TYPE, idest, 0, &
                     MPI_COMM_WORLD, status, ierr)
             else if(nsub_send(idest).gt.0)then
                call MPI_SEND( &
                     sendbuf(1:nsub_send(idest)), &
                     nsub_send(idest), SUBHALO_MPI_TYPE, idest, 0, &
                     MPI_COMM_WORLD, ierr)
             else if(nsub_recv(idest).gt.0)then
                call MPI_RECV( &
                     recvbuf(roffset(idest)+1:roffset(idest)+nsub_recv(idest)), &
                     nsub_recv(idest), SUBHALO_MPI_TYPE, idest, 0, &
                     MPI_COMM_WORLD, status, ierr)
             endif
          endif
       end do
       deallocate(sendbuf)

       ! Now need to assemble array of subhalos to write in the
       ! appropriate order
       allocate(iobuf(sum(nsub_recv)))
       nsub_recv    = 0
       nwrite_left  = nwrite
       ktree_import = itree_import
       ksub_import  = isub_import
       offset       = 0
       do while(nwrite_left.gt.0)
          ! Check if we're writing all remaining nodes of the next tree
          nnodes = trees_recv(recv_idx(ktree_import))%nnodes
          itask  = trees_recv(recv_idx(ktree_import))%itask
          nodes_left = nnodes - ksub_import + 1
          if(nwrite_left.ge.nodes_left)then
             ! Writing the rest of this tree
             iobuf(offset+1:offset+nodes_left) = &
                  recvbuf(roffset(itask)+nsub_recv(itask)+1:roffset(itask)+nsub_recv(itask)+nodes_left)
             nsub_recv(itask) = nsub_recv(itask) + nodes_left
             nwrite_left  = nwrite_left - nodes_left
             offset = offset + nodes_left
             ! Go to first subhalo in next tree
             ktree_import = ktree_import + 1
             ksub_import  = 1 
          else
             ! Writing part of this tree
             iobuf(offset+1:offset+nwrite_left) = &
                  recvbuf(roffset(itask)+nsub_recv(itask)+1:roffset(itask)+nsub_recv(itask)+nwrite_left)
             nsub_recv(itask) = nsub_recv(itask) + nwrite_left
             offset = offset + nwrite_left
             ! Advance to remaining part of this tree
             ksub_import = ksub_import + nwrite_left
             nwrite_left = 0
          endif   
       end do

       ! Write out the subhalos
       if(sum(nsub_recv).gt.0)then
          call write_subhalos(ihdf5, file_offset, iobuf, read_wdm, mpart, use_total_mass)
          file_offset = file_offset + sum(nsub_recv)
       endif

       deallocate(iobuf, recvbuf)
       
       ! Advance to next tree to import
       itree_import = jtree_import
       isub_import  = jsub_import

       ! Update number of subhalos left to write
       nleft = nleft - nwrite

       ! Check if we're done
       call MPI_ALLREDUCE(nleft, nleft_max, 1, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       if(myid.eq.0)write(*,*)"  Maximum subhalos left to write: ", nleft_max
       if(nleft_max.eq.0)exit
    end do

    deallocate(task, nsub_recv, nsub_send, soffset, roffset)

    return
  end subroutine write_tree_data




  subroutine write_subhalos(ihdf5, file_offset, subhalo, read_wdm, mpart, use_total_mass)
    !
    ! Write an array of subhalos to the output file at the specified offset
    !
    use subhalo_data

    implicit none
    integer, intent(in) :: ihdf5, file_offset
    type (subhalo_type), dimension(:), intent(in) :: subhalo
    logical, intent(in) :: read_wdm, use_total_mass
    real,    intent(in) :: mpart
    integer, dimension(7) :: start, count, start2d, count2d
    integer :: nwrite
    integer :: i
    ! Output buffer for vectors
    real,    dimension(:,:), allocatable :: buf
    integer, dimension(:),   allocatable :: flag
    real,    dimension(:),   allocatable :: redshift
    integer(kind=int8byte), dimension(:), allocatable :: mbid

    nwrite = size(subhalo)
    start(1)      = file_offset + 1
    count(1)      = size(subhalo)
    start2d(1:2)  = (/1, file_offset + 1 /)
    count2d(1:2)  = (/3, size(subhalo)   /)

    ! Quantities defining the tree structure
    call hdf5_write_data(ihdf5, "/haloTrees/nodeIndex", &
         subhalo%id, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/descendantIndex", &
         subhalo%descendant_id, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/hostIndex", &
         subhalo%dhalo_index, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/snapshotNumber", &
         subhalo%snapnum, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/descendantSnapshot", &
         subhalo%descendant_snapnum, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/enclosingIndex", &
         subhalo%parent_id, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/fofIndex", &
         subhalo%fof_index, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/mbpsContributed", &
         subhalo%mbps_contributed, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/descendantHost", &
         subhalo%next_dhalo_index, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    call hdf5_write_data(ihdf5, "/haloTrees/mainProgenitorIndex", &
         subhalo%main_prog_id, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Number of particles in the subfind group
    call hdf5_write_data(ihdf5, "/haloTrees/particleNumber", &
         subhalo%np, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    ! Mass of the subfind group (Note: ignores contamination)
    if(use_total_mass)then
       ! Use total mass from subfind
       call hdf5_write_data(ihdf5, "/haloTrees/nodeMass", &
            subhalo%mass, &
            start=start, count=count, extensible=.true., gzip=gzip_level)
    else
       ! Use mass = np*mpart
       call hdf5_write_data(ihdf5, "/haloTrees/nodeMass", &
            subhalo%np * mpart, &
            start=start, count=count, extensible=.true., gzip=gzip_level)
    endif

    ! Number of particles in the group when it was isolated
    call hdf5_write_data(ihdf5, "/haloTrees/originalParticleNumber", &
         subhalo%np_orig, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! MostboundID
    allocate(mbid(nwrite))
    mbid = subhalo%mostboundid
    call hdf5_write_data(ihdf5, "/haloTrees/mostBoundID", &
         mbid, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    deallocate(mbid)
    allocate(buf(1:3,nwrite))
    ! Position
    do i = 1, nwrite, 1
       buf(1:3, i) = subhalo(i)%pos(1:3)
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/position", &
         buf(1:3, 1:nwrite), start=start2d, count=count2d, &
         extensible=.true., gzip=gzip_level)
    ! Centre of mass
    do i = 1, nwrite, 1
       buf(1:3, i) = subhalo(i)%cofm(1:3)
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/centreOfMass", &
         buf(1:3, 1:nwrite), start=start2d, count=count2d, &
         extensible=.true., gzip=gzip_level)
    ! Velocity
    do i = 1, nwrite, 1
       buf(1:3, i) = subhalo(i)%vel(1:3)
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/velocity", &
         buf(1:3, 1:nwrite), start=start2d, count=count2d, &
         extensible=.true., gzip=gzip_level)
    ! Spin (angular momentum as output by subfind, not the spin parameter)
    do i = 1, nwrite, 1
       buf(1:3, i) = subhalo(i)%spin(1:3)
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/angularMomentum", &
         buf(1:3, 1:nwrite), start=start2d, count=count2d, &
         extensible=.true., gzip=gzip_level)
    deallocate(buf)
    ! Vmax
    call hdf5_write_data(ihdf5, "/haloTrees/maximumCircularVelocity", &
         subhalo%vmax, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    ! RVmax
    call hdf5_write_data(ihdf5, "/haloTrees/maximumCircularVelocityRadius", &
         subhalo%rvmax, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    ! Velocity dispersion
    call hdf5_write_data(ihdf5, "/haloTrees/velocityDispersion", &
         subhalo%veldisp, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Concentration
    call hdf5_write_data(ihdf5, "/haloTrees/cnfw", &
         subhalo%cnfw, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Lambda
    call hdf5_write_data(ihdf5, "/haloTrees/lambda", &
         subhalo%lambda, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Vvir
    call hdf5_write_data(ihdf5, "/haloTrees/Vvir", &
         subhalo%Vvir, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Mgas
    call hdf5_write_data(ihdf5, "/haloTrees/Mgas", &
         subhalo%mgas, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Rhalf
    call hdf5_write_data(ihdf5, "/haloTrees/halfMassRadius", &
         subhalo%rhalf, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Main branch information:
    ! Maximum past vmax on main branch
    call hdf5_write_data(ihdf5, "/haloTrees/mainBranchMaximumVmax", &
         subhalo%max_vmax, &
         start=start, count=count, extensible=.true., gzip=gzip_level)
    ! Maximum past mass on main branch
    call hdf5_write_data(ihdf5, "/haloTrees/mainBranchMaximumMass", &
         subhalo%max_mass, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    ! Main progenitor flag
    allocate(flag(nwrite))
    do i = 1, nwrite, 1
       if(iand(subhalo(i)%flags, MAIN_PROGENITOR_FLAG).ne.0)then
          flag(i) = 1
       else
          flag(i) = 0
       endif
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/isMainProgenitor", &
         flag, start=start, count=count, extensible=.true., gzip=gzip_level)
    ! FoF centre flag
    do i = 1, nwrite, 1
       if(iand(subhalo(i)%flags, FOF_CENTRE_FLAG).ne.0)then
          flag(i) = 1
       else
          flag(i) = 0
       endif
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/isFoFCentre", &
         flag, start=start, count=count, extensible=.true., gzip=gzip_level)
    ! DHalo centre flag
    do i = 1, nwrite, 1
       if(iand(subhalo(i)%flags, DHALO_CENTRE_FLAG).ne.0)then
          flag(i) = 1
       else
          flag(i) = 0
       endif
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/isDHaloCentre", &
         flag, start=start, count=count, extensible=.true., gzip=gzip_level)
    ! Interpolated flag
    do i = 1, nwrite, 1
       if(iand(subhalo(i)%flags, INTERPOLATED_FLAG).ne.0)then
          flag(i) = 1
       else
          flag(i) = 0
       endif
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/isInterpolated", &
         flag, start=start, count=count, extensible=.true., gzip=gzip_level)
    ! Remerged flag
    do i = 1, nwrite, 1
       if(iand(subhalo(i)%flags, REMERGED_FLAG).ne.0)then
          flag(i) = 1
       else
          flag(i) = 0
       endif
    end do
    call hdf5_write_data(ihdf5, "/haloTrees/isRemerged", &
         flag, start=start, count=count, extensible=.true., gzip=gzip_level)
    if(read_wdm)then
       ! WDM flag
       do i = 1, nwrite, 1
          if(iand(subhalo(i)%flags, WDM_SPURIOUS_FLAG).ne.0)then
             flag(i) = 1
          else
             flag(i) = 0
          endif
       end do
       call hdf5_write_data(ihdf5, "/haloTrees/WDMisSpuriousFlag", &
            flag, start=start, count=count, extensible=.true., gzip=gzip_level)
       ! Index of subhalo which set the WDM flag
       call hdf5_write_data(ihdf5, "/haloTrees/halfMaxMassIndex", &
            subhalo%half_max_mass_id, &
            start=start, count=count, extensible=.true., gzip=gzip_level)
    endif
    deallocate(flag)
    ! Redshift
    allocate(redshift(nwrite))
    do i = 1, nwrite, 1
       redshift(i) = zred(subhalo(i)%snapnum)
     end do
    call hdf5_write_data(ihdf5, "/haloTrees/redshift", &
         redshift, start=start, count=count, extensible=.true., gzip=gzip_level)
    deallocate(redshift)

    ! Position in subfind output
    call hdf5_write_data(ihdf5, "/haloTrees/positionInCatalogue", &
         subhalo%index, &
         start=start, count=count, extensible=.true., gzip=gzip_level)

    return
  end subroutine write_subhalos


  subroutine add_comments(ihdf5)

    use subhalo_reader
    implicit none
    integer, intent(in) :: ihdf5
    ! Date and time
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer   :: values(8)
    character(len=500) :: date_str

    ! Add comments
    ! haloTrees group
    call hdf5_write_attribute(ihdf5, "/haloTrees/nodeIndex/comment", &
         "Identifier for each node in the tree, unique within the simulation")
    call hdf5_write_attribute(ihdf5, "/haloTrees/descendantIndex/comment", &
         "nodeIndex of this node's descendant")
    call hdf5_write_attribute(ihdf5, "/haloTrees/hostIndex/comment", &
         "nodeIndex of the top level parent halo")
    call hdf5_write_attribute(ihdf5, "/haloTrees/snapshotNumber/comment", &
         "Snapshot at which this node exists")
    call hdf5_write_attribute(ihdf5, "/haloTrees/descendantSnapshot/comment", &
         "Snapshot at which this node's descendant exists")
    call hdf5_write_attribute(ihdf5, "/haloTrees/enclosingIndex/comment", &
         "nodeIndex of the node which contains this node")
    call hdf5_write_attribute(ihdf5, "/haloTrees/fofIndex/comment", &
         "Identifier for the FoF group which contains this node")
    call hdf5_write_attribute(ihdf5, "/haloTrees/mbpsContributed/comment", &
         "How many of the descendant node's most bound particles came from this node")
    call hdf5_write_attribute(ihdf5, "/haloTrees/particleNumber/comment", &
         "Number of particles in the SubFind group corresponding to this node")
    call hdf5_write_attribute(ihdf5, "/haloTrees/originalParticleNumber/comment", &
         "Number of particles in the node when it was last an independent halo")
    call hdf5_write_attribute(ihdf5, "/haloTrees/mostBoundID/comment", &
         "Most bound particle ID for the SubFind group corresponding to this node")
    call hdf5_write_attribute(ihdf5, "/haloTrees/position/comment", &
         "Position of potential minimum of the node as calculated by SubFind (comoving)")
    call hdf5_write_attribute(ihdf5, "/haloTrees/velocity/comment", &
         "Velocity of this node as calculated by SubFind (peculiar)")
    call hdf5_write_attribute(ihdf5, "/haloTrees/angularMomentum/comment", &
         "Angular momentum of this node as calculated by SubFind")
    call hdf5_write_attribute(ihdf5, "/haloTrees/maximumCircularVelocity/comment", &
         "Maximum circular velocity of this node as calculated by SubFind "//&
         "(physical units inc. hubble flow)")
    call hdf5_write_attribute(ihdf5, "/haloTrees/velocityDispersion/comment", &
         "Velocity dispersion of this node as calculated by SubFind "//&
         "(physical units inc. hubble flow)")
    call hdf5_write_attribute(ihdf5, "/haloTrees/halfMassRadius/comment", &
         "Half mass radius of the subfind group corresponding to this node (comoving)")
    call hdf5_write_attribute(ihdf5, "/haloTrees/isMainProgenitor/comment", &
         "1 if this node is the main progenitor of its descendant, 0 otherwise")
    call hdf5_write_attribute(ihdf5, "/haloTrees/isFoFCentre/comment", &
         "1 if this node is the top level subhalo in its FoF group, 0 otherwise")
    call hdf5_write_attribute(ihdf5, "/haloTrees/isDHaloCentre/comment", &
         "1 if this node is the most massive in its DHalo, 0 otherwise")
    call hdf5_write_attribute(ihdf5, "/haloTrees/isInterpolated/comment", &
         "1 if this node was added by the merger tree code, 0 otherwise")
    call hdf5_write_attribute(ihdf5, "/haloTrees/isRemerged/comment", &
         "1 if this node was merged back onto its host due to a de-merger, 0 otherwise")
    call hdf5_write_attribute(ihdf5, "/haloTrees/redshift/comment", &
         "redshift at which this node exists")
    call hdf5_write_attribute(ihdf5, "/haloTrees/descendantHost/comment", &
         "Which host the descendant of this node exists at")
    call write_linking_parameters(ihdf5, "treeConstruction")

    ! Add date and time
    call date_and_time(date, time, zone, values)
    date_str = date//" "//time//" "//zone
    call hdf5_write_attribute(ihdf5, "fileInfo/creationTime", date_str)

    return
  end subroutine add_comments

end module tree_files
