program build_trees
!
! Construct subhalo merger trees
!
  use mpi
  use get_argumentsmod
  use read_parameters
  use read_subfind
  use subhalo_data
  use subhalo_reader
  use subhalo_trees
  use redistribute
  use enclosing_subhalos
  use splitting
  use dhalo_index_module
  use tree_files
  use merger_trees
  use remerging
  use subhalo_interpolation
  use metadata
  use host_halos
  use hbt_reader
  use nifty_reader
  use ROCKSTAR_reader ! include ROCKSTAR_reader
  use relabel_snapshots
  use VELOCIraptor_reader

  implicit none

  ! Input parameters
  character(LEN=500) :: outdir, basedir, paramfile, subfind_format
  character(LEN=500) :: subtab_path, subids_path, zfile, tree_basename
  character(LEN=500) :: hbt_path, nifty_path, ROCKSTAR_path,snapnum_file,VELOCIraptor_path,desc_file
  ! include ROCKSTAR_path
  character(LEN=500), dimension(10) :: gadget_param
  integer :: nproc_io
  integer :: nparam
  integer :: ifirst, ilast
  integer :: id_size, float_size, nsteptrace, ntreefile
  real    :: lbox, rfactor, mfrac, mpart, dum, t1, t2, h0
  logical :: remerge, interpolate
  integer :: nsub_tot, nsum
  logical :: separate_branches, hydrorun
  logical :: use_total_mass
  logical :: discard_spurious = .false.
  real    :: sphericity_cut = 0.0, mass_cut = 0.0
  character(len=500) :: wdm_file, descfil
  ! MPI startup
  integer :: myid, numprocs, ierr
  ! Subhalo data
  type (subhalo_type), dimension(:), pointer :: subhalo
  ! Loops etc
  integer :: isnap, isub, i
  type (subhalo_tree_node), pointer :: node
  integer :: new_trees

  ! Descendants.txt information holder, used at east by VELOCIraptor
  integer(kind=int8byte), dimension(:), allocatable :: halo_ids, desc_ids
  integer, dimension(:), allocatable :: desc_snaps
  integer :: n_halos
!
! Start MPI and find out which processor this is
!
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  if(myid.eq.0) &
       write(*,*)'build_trees starting on ', numprocs, ' processors'
!
! Read in parameters
!
  if(myid.eq.0)then
     call set_usage("build_trees.exe <parameterfile>")
     call get_argument(1,paramfile)
     call read_parameterfile(paramfile,err_unit=-1,verbose=.true.)
     call get_parameter("basedir",basedir)
     call get_parameter("subids_path",subids_path)
     call get_parameter("subtab_path",subtab_path)
     call get_parameter("hbt_path",hbt_path)
     call get_parameter("nifty_path",nifty_path)
     call get_parameter("ROCKSTAR_path",ROCKSTAR_path) ! include
     call get_parameter("VELOCIraptor_path",VELOCIraptor_path)
     call get_parameter("desc_file",desc_file)
     call get_parameter("hydrorun",hydrorun)
     call get_parameter("snapnum_file",snapnum_file)
     call get_parameter("id_size",id_size)
     call get_parameter("float_size",float_size)
     call get_parameter("nproc_io", nproc_io)
     call get_parameter("ifirst",ifirst)
     call get_parameter("ilast",ilast)
     call get_parameter("treedir",outdir)
     call get_parameter("tree_basename",tree_basename)
     call get_parameter("subfind_format",subfind_format)
     call get_parameter("nsteptrace",nsteptrace)
     call get_parameter("lbox",lbox)
     call get_parameter("mpart",mpart)
     call get_parameter("h0",h0)
     call get_parameter("zfile",zfile)
     call get_parameter("rfactor",rfactor)
     call get_parameter("mfrac",mfrac)
     call get_parameter("remerge", remerge)
     call get_parameter("interpolate", interpolate)
     call get_parameter("ntreefile", ntreefile)
     call get_parameter("gadget_parameters", gadget_param, &
          nelements=nparam, default=(/"none"/))
     call get_parameter("separate_branches", separate_branches)
     call get_parameter("wdm_file", wdm_file)
     call get_parameter("discard_spurious", discard_spurious, required=wdm_file.ne."none")
     call get_parameter("sphericity_cut", sphericity_cut, required=wdm_file.ne."none")
     call get_parameter("mass_cut",       mass_cut,       required=wdm_file.ne."none")
     call get_parameter("use_total_mass", use_total_mass)
     call finalize_parameters(abort_on_unused=.false., abort_on_default=.false.)
     if(remerge.and.(.not.interpolate)) &
          call panic("Remerging without interpolation is not supported!")
     if(discard_spurious.and.(.not.separate_branches)) &
          call panic("Discarding spurious groups only works with separate_branches=.true.!")        
     if(.not.use_total_mass.and.mpart.le.0.0)then
        call panic("mpart must be specified if use_total_mass=.false.")
     endif

  end if
!
! Broadcast parameters to all processors
!
  call MPI_BCAST(basedir,len(basedir),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(subids_path,len(subids_path),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(subtab_path,len(subtab_path),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(hbt_path,len(hbt_path),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nifty_path,len(nifty_path),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ROCKSTAR_path,len(ROCKSTAR_path),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr) ! include
  call MPI_BCAST(VELOCIraptor_path,len(VELOCIraptor_path),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(desc_file,len(desc_file),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(snapnum_file,len(snapnum_file),MPI_CHARACTER,0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ifirst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ilast,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nproc_io,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(outdir,len(outdir),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tree_basename,len(tree_basename),MPI_CHARACTER,0,&
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(subfind_format,len(subfind_format),MPI_CHARACTER,0,&
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(lbox,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(mpart,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rfactor,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(id_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(float_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(mfrac,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(remerge, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(hydrorun, 1,MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(interpolate, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ntreefile, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nsteptrace, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(separate_branches, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(wdm_file, len(wdm_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(discard_spurious, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(sphericity_cut, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mass_cut, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nparam, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(gadget_param, len(gadget_param)*size(gadget_param), &
       MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(use_total_mass, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
 
!
! Read file with snapshot renumbering info
!
  call read_snapnum_file(snapnum_file)
  if(trim(adjustl(snapnum_file)).ne."none") &
       call get_relabelled_range(ifirst, ilast)
!
! Read information about the simulation and linking algorithm so it
! can be added to the output files
!
  if(gadget_param(1).eq."none")nparam = 0
  call get_metadata(gadget_param(1:nparam))
  if(subfind_format.ne."HBT" .and. subfind_format.ne."nifty" .and. subfind_format.ne."VELOCIraptor" &
       .and. subfind_format.ne."ROCKSTAR")then ! include
     call read_linking_parameters(outdir, ifirst, &
          rfactor, mfrac, remerge, interpolate)
  endif
!
! Get redshifts on processor 0 and broadcast
!
  allocate(zred(ifirst:ilast))
  allocate(sfactor(ifirst:ilast)) ! include
  if(myid.eq.0)then
     ! Read from ascii list
     open(unit=1,file=zfile,status='old',form='formatted')
     if(subfind_format .eq. "nifty")then
          read(1,*)
          do i=ifirst,ilast,1
             isnap=-1
             do while(isnap.lt.file_snapnum(i))
                read(1,*)isnap,dum,zred(i)
             end do
             if(isnap.ne.file_snapnum(i))call panic('Redshift/snapshot indices do not match!')
             if(i.gt.ifirst)then
                if(zred(i).ge.zred(i-1)) &
                     call panic("Redshifts must be in descending order!")
             end if
            end do
      else if(subfind_format .eq. "VELOCIraptor")  then 
             read(1,*)
             do i=ifirst,ilast,1
                isnap=-1
                do while(isnap.lt.file_snapnum(i))
                   read(1,*)isnap,dum,zred(i)
                end do
                if(isnap.ne.file_snapnum(i))call panic('Redshift/snapshot indices do not match!')
                if(i.gt.ifirst)then
                   if(zred(i).ge.zred(i-1)) &
                        call panic("Redshifts must be in descending order!")
                end if
             end do
      ! include
      else if(subfind_format .eq. "ROCKSTAR")  then
             read(1,*)
             do i=ifirst,ilast,1
                isnap=-1
                do while(isnap.lt.file_snapnum(i))
                   read(1,*)isnap,sfactor(i),zred(i)
                   ! redshift file with 3 columns: snapshot (0), scale factor (1), redshift (2)
                end do
                if(isnap.ne.file_snapnum(i))call panic('Redshift/snapshot indices do not match!')
                if(i.gt.ifirst)then
                   if(zred(i).ge.zred(i-1)) &
                        call panic("Redshifts must be in descending order!")
                   if(sfactor(i).le.sfactor(i-1)) &
                        call panic("Scale factors must be in ascending order!")
                end if
             end do
       else
          do i=ifirst,ilast,1
             isnap=-1
             do while(isnap.lt.file_snapnum(i))
                read(1,*)isnap,zred(i)
             end do
             if(isnap.ne.file_snapnum(i))call panic('Redshift/snapshot indices do not match!')
             if(i.gt.ifirst)then
                if(zred(i).gt.zred(i-1)) &
                     call panic("Redshifts must be in descending order!")
             end if
          end do
     endif
     close(1)
  end if
  call MPI_BCAST(zred,(ilast-ifirst+1),MPI_REAL,0,MPI_COMM_WORLD,ierr)
!
! Initialise module for reading subfind outputs
!
  if(subfind_format.ne."HBT".and.subfind_format.ne."nifty".and.subfind_format.ne."VELOCIraptor" &
       .and.subfind_format.ne."ROCKSTAR")then ! include
     call set_subfind_format(subfind_format,id_size,float_size)
     call set_subfind_path(basedir,subtab_path,subids_path)
  endif
!
! Initialise subhalo data module
!
  call subhalo_data_init(ifirst, ilast)

!
! Read the descendants.txt whole for VELOCIraptor,
! which consumes more memory but is far more efficient
!
  ! ROCKSTAR does not need read descendants because they are in the same file
  if( subfind_format == "VELOCIraptor" .or. subfind_format == "nifty" ) then
    descfil = trim(desc_file)
    print '(A)', "Reading tree file whole"
    call read_descendants(descfil, n_halos, halo_ids, desc_ids, desc_snaps)
  endif
!
! Loop over earlier snapshots and add subhalos to the merger trees
!
  do isnap = ilast, ifirst, -1
     !print*,'isnapmain2:',isnap
     if(myid.eq.0) &
          write(*,*)"Adding subhalos for snapshot ", isnap
     ! Read subhalos for this snapshot
     if(subfind_format.ne."HBT".and.subfind_format.ne."nifty".and.subfind_format.ne."VELOCIraptor" &
          .and.subfind_format.ne."ROCKSTAR")then ! include
        ! Use SubFind + Find_Descendants output
        call read_subhalos(outdir, isnap, ilast, subhalo, nproc_io, wdm_file, sphericity_cut)
     else if (subfind_format.eq."HBT") then
        ! Use Jiaxin's HBT output which already contains descendant info
        if(wdm_file.ne."none")call panic("WDM mode not implemented for HBT groups!")
        ! Need to add code to distribute catalogue to run >1 MPI process
        if(numprocs.gt.1)call panic("Can only run on one processor with HBT groups!")
        call read_hbt(hbt_path, isnap, ilast, subhalo)
     else if (subfind_format.eq."nifty") then
        ! Use nIFTy/AFH output which already contains descendant info
        if(wdm_file.ne."none")call panic("WDM mode not implemented for nifty groups!")
        ! Need to add code to distribute catalogue to run >1 MPI process
        if(numprocs.gt.1)call panic("Can only run on one processor with nifty groups!")

        call cpu_time ( t1 )
        call read_nifty(nifty_path, isnap, ilast, subhalo, mpart, lbox, &
                              &n_halos, halo_ids, desc_ids, desc_snaps)
        call cpu_time ( t2 )
        write ( *, * ) 'Elapsed CPU time in read_nifty = ', t2 - t1

     ! include
     else if (subfind_format.eq."ROCKSTAR") then
	! Use ROCKSTAR output which already contains descendant info
        if(wdm_file.ne."none")call panic("WDM mode not implemented for ROCKSTAR groups!")
        ! Need to add code to distribute catalogue to run >1 MPI process
        if(numprocs.gt.1)call panic("Can only run on one processor with ROCKSTAR groups!")
        call cpu_time ( t1 )
        call read_ROCKSTAR(ROCKSTAR_path, isnap, ilast, subhalo, mpart, lbox)
        call cpu_time ( t2 )
        write ( *, * ) 'Elapsed CPU time in read_ROCKSTAR = ', t2 - t1

     else if (subfind_format.eq."VELOCIraptor") then
        ! Use nIFTy/AFH output which already contains descendant info
        if(wdm_file.ne."none")call panic("WDM mode not implemented for VELOCIraptor groups!")
        ! Need to add code to distribute catalogue to run >1 MPI process
        if(numprocs.gt.1)call panic("Can only run on one processor with VELOCIraptor groups!")
        call cpu_time ( t1 )
        call read_VELOCIraptor(VELOCIraptor_path, isnap, ilast, subhalo, mpart, lbox, &
                              &n_halos, halo_ids, desc_ids, desc_snaps, h0, hydrorun)
        call cpu_time ( t2 )
        write ( *, * ) 'Elapsed CPU time in read_VELOCIraptor = ', t2 - t1
     endif
     call MPI_REDUCE(size(subhalo), nsub_tot, 1, MPI_INTEGER, MPI_SUM, &
          0, MPI_COMM_WORLD, ierr)
     if(myid.eq.0) then
        write(*,*)"  Read in ", nsub_tot, " subhalos"
     endif
     ! Redistribute final subhalos among mpi tasks to improve memory load balancing
     if(isnap.eq.ilast)then
        if(myid.eq.0)write(*,*)"  Redistributing subhalos between tasks"
        call redistribute_subhalos(subhalo)
     endif
     ! Add subhalos to the merger trees
     if(isnap.lt.ilast) &
          call add_to_tree(isnap, subhalo, nsteptrace, ilast,&
          find_all=(subfind_format.ne."HBT".and.subfind_format.ne."nifty".and.subfind_format.ne."VELOCIraptor" &
          .and.subfind_format.ne."ROCKSTAR"))
     ! Any remaining subhalos become new merger trees rooted at
     ! this snapshot
     new_trees = 0
     do isub = 1, size(subhalo), 1
        if(iand(subhalo(isub)%flags, IN_TREE_FLAG).eq.0)then
           ! These subhalos should not have descendants
           if(subhalo(isub)%descendant_id.ge.0)then
              if(subfind_format.ne."HBT".and.subfind_format.ne."nifty".and.subfind_format.ne."VELOCIraptor" &
                   .and.subfind_format.ne."ROCKSTAR")then
                 call panic("Subhalo with descendant not assigned to a tree!")
              else
                 subhalo(isub)%descendant_id = -1
              endif
           endif
           ! Make a new node with no descendant
           node => new_tree_node(subhalo(isub), isnap)
           new_trees = new_trees + 1
        endif
     end do
     call MPI_REDUCE(new_trees, nsum, 1, MPI_INTEGER, MPI_SUM, 0, &
          MPI_COMM_WORLD, ierr)
     if(myid.eq.0) then          
        write(*,*)"  There are ", nsum, " new merger trees at this redshift"
     endif
     ! Next snapshot
     deallocate(subhalo)
     call write_memory_usage()

  end do

  ! Free up a bit of memory
  if( subfind_format == "VELOCIraptor" .or. subfind_format == "nifty" ) then
    deallocate(halo_ids, desc_ids, desc_snaps)
  endif
  ! ROCKSTAR doesn't need descendants file
!
! Identify and flag main progenitor of each subhalo
!
  if(myid.eq.0) &
       write(*,*)"Finding main progenitor for each subhalo"
  call set_tree_pointers()
  !call find_main_progenitors()
  if (subfind_format.eq."nifty") then
     call nifty_find_main_progenitors()
  else if (subfind_format.eq."ROCKSTAR") then ! include
     call ROCKSTAR_find_main_progenitors()
  else if (subfind_format.eq."VELOCIraptor") then 
     call VELOCIraptor_find_main_progenitors()
  else
     call find_main_progenitors()
  endif
!
! Discard WDM spurious subhalos if necessary
! If discard_spurious=.false. this just makes sure that
! WDM flags are consistent along branches.
!
  if(wdm_file.ne."none")then
     if(myid.eq.0) &
          write(*,*)"Checking for spurious WDM tree branches"
     call discard_spurious_groups(discard_spurious, mass_cut, mpart)
     ! Removing subhalos invalidates progenitor linked lists,
     ! so recalculate them
     call set_tree_pointers()
     if (subfind_format.eq."nifty") then
        call nifty_find_main_progenitors()
     else if (subfind_format.eq."ROCKSTAR") then ! include
        call ROCKSTAR_find_main_progenitors()
     else if (subfind_format.eq."VELOCIraptor") then
        call VELOCIraptor_find_main_progenitors()
     else
        call find_main_progenitors()
     endif
  endif
!
! Find out which subhalos "enclose" each other
!
  if(myid.eq.0) &
       write(*,*)"Finding enclosing subhalos"
  call group_by_radius(lbox,rfactor,use_total_mass)
!
! Split off subhalos which have retained enough mass
!
  if(myid.eq.0) &
       write(*,*)"Splitting off subhalos which have retained mass"
  call split_by_mass(mfrac,use_total_mass)
!
! For each subhalo, find max np and vmax along main branch
!
  call find_max_mass_and_vmax(mpart, use_total_mass)
!
! Add interpolated subhalos
!
  if(interpolate)then
     if(myid.eq.0) &
          write(*,*)"Adding interpolated subhalos"
     call interpolate_subhalos(ifirst, ilast)
  endif
  call write_memory_usage()
!
! Do remerging if necessary
!
  if(remerge)then
     ! DHalo indexes are calculated as a side effect of the remerging
     ! procedure, so there's no need to do it separately in this case.
     ! This currently only works if interpolation is enabled.
     if(myid.eq.0)write(*,*)"Remerging subhalos"
     call remerge_subhalos(ifirst, ilast)
  else
     ! Just calculate DHalo indexes
     if(myid.eq.0)write(*,*)"Assigning subhalos to dhalos"
     call find_dhalo_index(ifirst, ilast)  
  endif
  call write_memory_usage()
!
! Identify the base of each merger tree and find final
! dhalo, subhalo and snapnum
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myid.eq.0) &
       write(*,*)"Identifying merger trees"
  call find_merger_trees(separate_branches,subfind_format)
  call write_memory_usage()
!
! Find the host halo ID for each subhalo merger tree. Since trees don't
! always survive to the final time this is not just the DHalo ID of the
! final subhalo in the tree. The result is stored in the dhalo_index
! for each subhalo tree.
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myid.eq.0) &
       write(*,*)"Finding final host halos"
  call find_host_halos()
!
! Decide which DHalo trees go in which files
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myid.eq.0) &
       write(*,*)"Assigning merger trees to files"
  call assign_final_dhalos_to_files(ntreefile)
  call write_memory_usage()
!
! Output the merger trees
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myid.eq.0) &
       write(*,*)"Writing out merger trees"
  call write_tree_files(ntreefile, outdir, trim(tree_basename), &
       (wdm_file.ne."none"), mpart, lbox, remerge, nproc_io, use_total_mass)
!
! Deallocate subhalo merger tree data
!
  call subhalo_data_deallocate()
!
! Write file to indicate completion and shut down
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(myid.eq.0)then
     open(unit=1,file=trim(outdir)//"/build_trees_done", &
          status="unknown", form="formatted")
     write(1,*)"Ran to completion"
     close(1)
  endif
  deallocate(zred)
  deallocate(sfactor) ! include
  call MPI_FINALIZE(ierr)
  
  if(myid.eq.0) &
       write(*,*)"Finished"

contains

  subroutine write_memory_usage()
!
! Write maximum memory usage to stdout
!
    implicit none
    real :: vmem, rss
    real :: max_vmem, min_vmem, max_rss, min_rss
    integer :: ierr

    call getmemuse(vmem, rss)
    call MPI_REDUCE(vmem, max_vmem, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(rss,  max_rss,  1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(vmem, min_vmem, 1, MPI_REAL, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(rss,  min_rss,  1, MPI_REAL, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

    if(myid.eq.0)then
       if(max_vmem.gt.0.or.max_rss.gt.0)write(*,*)
       if(max_vmem.gt.0)&
            write(*,*)"  Virtual memory usage:  min = ", min_vmem, ", max = ", max_vmem, " Mb"
       if(max_rss.gt.0)&
            write(*,*)"  Resident set size   :  min = ", min_rss,  ", max = ", max_rss,  " Mb"
       if(max_vmem.gt.0.or.max_rss.gt.0)write(*,*)
    endif

    return
  end subroutine write_memory_usage

end program build_trees



