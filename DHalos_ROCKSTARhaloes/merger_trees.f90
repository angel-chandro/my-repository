module merger_trees

  use mpi
  use subhalo_data
  use subhalo_trees
  use sort
  use multi_key_sort
  use mpi_parallel_sort
  use assignment

  implicit none
  private
  save

  public :: walk_tree
  public :: find_merger_trees
  public :: get_tree_data
  public :: count_tree_nodes
  public :: discard_spurious_groups
  public :: assign_final_dhalos_to_files

  ! Information about subhalo merger trees
  ! (or branches if we have separate_branches=true)
  integer, public :: ntrees
  type merger_tree_type
     type(subhalo_tree_node), pointer :: ptr
     integer                  :: nnodes
     integer(kind=int8byte)   :: dhalo_index
     integer(kind=int8byte)   :: subhalo_index
     integer                  :: final_snapnum
     integer                  :: ifile
     integer                  :: itree_file
  end type merger_tree_type
  type(merger_tree_type), dimension(:), allocatable, public :: merger_tree

contains


  function walk_tree(node)
!
! Visit next node in the tree, returns NULL if there are no more
!
    type (subhalo_tree_node), pointer :: node, tmp, walk_tree
    
    tmp => node
    if(associated(node%progenitor))then
       walk_tree => node%progenitor
    else
       do while(.not.associated(tmp%next_progenitor))
          tmp => tmp%descendant
          if(.not.associated(tmp))then
             nullify(walk_tree)
             return
          endif
       end do
       walk_tree => tmp%next_progenitor
    endif
    
    return
  end function walk_tree


  subroutine find_merger_trees(separate_branches,subfind_format)
!
! Find all subhalos with no descendant - each one is the base of
! a merger tree
!
    use nifty_reader
    use VELOCIraptor_reader
    use ROCKSTAR_reader ! include

    implicit none
    logical :: separate_branches
    integer :: isnap
    type(subhalo_tree_node), pointer :: node
    integer :: i
    character(LEN=500) :: subfind_format

    ! Sanity check
    do isnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
       ! Loop over subhalos at this snapshot
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
          if(.not.associated(node%descendant))then
             if(node%subhalo%descendant_id.ne.-1)then
                call panic("Failed to find descendant for node!")
             endif
          else
             if(node%subhalo%descendant_id.eq.-1)then
                call panic("Found descendant when descendant ID=-1!")
             endif
             if(node%subhalo%descendant_id.ne.node%descendant%subhalo%id)then
                call panic("Incorrect descendant found for node!")
             endif
          endif
          node => node%next 
       end do
    end do

    ! Chop merger trees into branches if necessary
    if(separate_branches)then
       do isnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
          ! Loop over subhalos at this snapshot
          node => node_at_snapshot(isnap)%first
          do while(associated(node))
             if(associated(node%descendant))then
                if(iand(MAIN_PROGENITOR_FLAG, node%subhalo%flags).eq.0)then
                   nullify(node%descendant)
                   !node%subhalo%descendant_id = -1
                endif
             endif
             node => node%next 
          end do
       end do
       ! Recalculate progenitor pointers
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

    ! Count trees
    ntrees = 0
    do isnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
       ! Loop over subhalos at this snapshot
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
          if(.not.associated(node%descendant))ntrees = ntrees + 1
          node => node%next 
       end do
    end do

    allocate(merger_tree(ntrees))
    ntrees = 0
    do isnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
          if(.not.associated(node%descendant))then
             ! Found a new tree
             ntrees = ntrees + 1
             ! Store a pointer to the tree
             merger_tree(ntrees)%ptr => node
          endif
          node => node%next 
       end do
       ! Next snapshot
    end do

    ! Find final DHalo, subhalo and snapnum for each tree
    do i = 1, ntrees, 1
       merger_tree(i)%dhalo_index   = merger_tree(i)%ptr%subhalo%dhalo_index
       merger_tree(i)%subhalo_index = merger_tree(i)%ptr%subhalo%id
       merger_tree(i)%final_snapnum = merger_tree(i)%ptr%subhalo%snapnum
       merger_tree(i)%itree_file    = -1 ! Position in file, not known yet
    end do

    ! Store number of nodes in each tree
    call count_tree_nodes()

    return
  end subroutine find_merger_trees
  

  subroutine get_tree_data(base, next_node, max_subhalos, subarr, nsub)
!
! Return an array of subhalos from the specified tree. If there are more
! than max_subhalos then return just this many and a pointer to the next
! node to output.
!
    implicit none
    ! Pointer to the base of the tree
    type(subhalo_tree_node), pointer :: base
    ! Pointer to the next node to return. Will return NULL if there
    ! are no more nodes
    type(subhalo_tree_node), pointer :: next_node
    ! Maximum no. of subhalos to return
    integer, intent(in) :: max_subhalos
    ! Data array to return. Should be allocated 1:max_subhalos
    type (subhalo_type), dimension(:) :: subarr
    ! Number of subhalos returned
    integer, intent(out) :: nsub
    
    nsub = 0
    ! Start from the root if next_node is null
    if(.not.associated(next_node))next_node => base

    do while(associated(next_node).and.nsub.lt.max_subhalos)
       ! Copy this subhalo to the output array
       nsub = nsub + 1
       subarr(nsub) = next_node%subhalo
       ! Advance to the next subhalo
       next_node => walk_tree(next_node)
    end do

    return
  end subroutine get_tree_data


  subroutine count_tree_nodes()

    implicit none
    type(subhalo_tree_node), pointer :: node
    integer :: i, nbase

    ! Count nodes in each tree
    do i = 1, ntrees, 1
       nbase = 0
       merger_tree(i)%nnodes = 0
       node => merger_tree(i)%ptr
       do while(associated(node))
          if(.not.associated(node%descendant))nbase = nbase + 1
          merger_tree(i)%nnodes = merger_tree(i)%nnodes + 1
          node => walk_tree(node)
          ! Next node
       end do
       if(nbase.ne.1)call panic("Tree with <>1 node with no descendant!")
       ! Next tree
    end do

  end subroutine count_tree_nodes


  
  subroutine discard_spurious_groups(discard_spurious, mass_cut, mpart)
    !
    ! Discard tree branches which are tagged  as fake at
    ! the point where they have half of their maximum mass.
    !
    ! WDM flags from input files: 0=real, 1=not real
    !
    implicit none
    logical :: discard_spurious
    real, intent(in) :: mass_cut, mpart
    type(subhalo_tree_node), pointer :: node, prev_node, next_node
    type(subhalo_tree_node), pointer :: prog, desc
    integer :: isnap
    integer :: max_np
    logical :: is_spurious
    integer :: nbranch, nbranch_all
    integer :: nremove, nremove_all
    integer :: myid, numprocs, ierr
    integer*8 :: half_max_mass_id
    logical :: new_branch

    ! Find which processor this is
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

    nbranch = 0
    nremove = 0

    ! Loop over nodes starting at early times
    do isnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
        
          ! If node has no main progenitor, its the start of a branch
          new_branch = .true.
          prog => node%progenitor
          do while(associated(prog))
             if(iand(prog%subhalo%flags, MAIN_PROGENITOR_FLAG).ne.0)&
                  new_branch = .false.
             prog => prog%next_progenitor
          end do

          if(new_branch)then
             nbranch = nbranch + 1
             ! Find maximum mass along this branch
             desc => node
             max_np = 0
             do
                ! Record maximum mass
                max_np = max(max_np, desc%subhalo%np)
                ! Go to next node on this branch, if there is one
                if(iand(desc%subhalo%flags, MAIN_PROGENITOR_FLAG).eq.0)then
                   ! End of branch
                   exit
                else
                   ! Go to next node on branch
                   if(.not.associated(desc%descendant))then
                      write(0,*)"Main progenitor has no descendant?!"
                      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                   endif
                   desc=>desc%descendant
                endif
             end do
             ! Find flag at half maximum mass
             half_max_mass_id = -1
             desc => node
             do
                if(desc%subhalo%np.ge.0.5*max_np)then
                   is_spurious = (iand(desc%subhalo%flags, WDM_SPURIOUS_FLAG).ne.0)
                   half_max_mass_id = desc%subhalo%id
                   exit
                endif
                desc=>desc%descendant                                
             end do
             ! We also want to flag branches where the maximum mass is less than
             ! mass_cut
             if(max_np*mpart.lt.mass_cut)is_spurious=.true.
             ! Set flags for whole branch
             if(is_spurious)nremove = nremove + 1
             desc => node
             do
                if(is_spurious)then
                   ! Set flag
                   desc%subhalo%flags = ior(desc%subhalo%flags, WDM_SPURIOUS_FLAG)
                else
                   ! Clear flag
                   desc%subhalo%flags = iand(desc%subhalo%flags, not(WDM_SPURIOUS_FLAG))
                endif
                ! Record ID of subhalo that set the flag
                desc%subhalo%half_max_mass_id = half_max_mass_id
                ! Go to next node on this branch, if there is one
                if(iand(desc%subhalo%flags, MAIN_PROGENITOR_FLAG).eq.0)then
                   ! End of branch
                   exit
                else
                   ! Go to next node on branch
                   desc=>desc%descendant
                endif
             end do
          endif
          ! Next starting node
          node => node%next 
       end do
    end do

    ! Just tag branches if this is not set
    if(.not.discard_spurious)return

    ! Remove spurious subhalos from linked list and deallocate
    do isnap = ubound(node_at_snapshot,1), lbound(node_at_snapshot,1), -1
       ! Loop over nodes at this snapshot
       node => node_at_snapshot(isnap)%first
       nullify(prev_node)
       do while(associated(node))
          ! Store pointer to next node
          next_node => node%next
          ! Check if we keep this one
          if(iand(node%subhalo%flags, WDM_SPURIOUS_FLAG).ne.0)then
             ! Need to discard this node
             if(.not.associated(prev_node))then
                ! Haven't found any real nodes on this snapshot yet,
                ! so need to update pointer to first node on snapshot
                ! (this will eventually nullify %first if there are no
                ! real nodes at this snapshot)
                node_at_snapshot(isnap)%first => next_node
             else
                ! Make next pointer on last real node point to
                ! next node
                prev_node%next => next_node
             endif
             ! Nullify descendant pointers of any progenitors of this node
             prog => node%progenitor
             do while(associated(prog))
                nullify(prog%descendant)
                prog%subhalo%descendant_id = -1
                prog => prog%next_progenitor
             end do
             ! Can't deallocate nodes if we're allocating them in blocks,
             ! but if they're removed from the linked list they wont be
             ! used in the output merger trees and they'll get deallocated
             ! when we free all the blocks at the end.
             !deallocate(node)
          else
             ! This is a real one. Keep a pointer to last real
             ! node we found.
             prev_node => node
          end if
          ! Advance to next node
          node => next_node
       end do
       ! Set pointer to last real node
       node_at_snapshot(isnap)%last => prev_node
       ! Next snapshot
    end do    

    ! Report number of branches removed
    call MPI_REDUCE(nremove, nremove_all, 1, MPI_INTEGER, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(nbranch, nbranch_all, 1, MPI_INTEGER, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    if(myid.eq.0)then
       write(*,*)"Removed ", nremove_all, " branches of ", nbranch_all
    endif

    return
  end subroutine discard_spurious_groups


  subroutine assign_final_dhalos_to_files(nfiles)
!
! Randomly assign each final DHalo to a file
!
    implicit none
    ! Parameters
    integer, intent(in) :: nfiles
    ! Sorting
    integer(kind=int8byte), dimension(:), allocatable :: dhalo_index, idx
    integer,                dimension(:), allocatable :: ifile, ifile_send
    ! MPI stuff
    integer :: myid, numprocs, ierr, status(MPI_STATUS_SIZE), iproc
    integer :: nrecv
    ! Current file and dhalo
    integer(kind=int8byte) :: current_dhalo
    integer                :: current_file
    ! Loops etc
    integer :: i
    real(kind=real8byte) :: r8

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

    ! Do parallel sort on final DHalo indexes
    allocate(dhalo_index(ntrees), idx(ntrees), ifile(ntrees))
    dhalo_index = merger_tree%dhalo_index
    call parallel_sort(ntrees, dhalo_index, idx)
    
    if(myid.gt.0)then
       ! Processors other than zero send their dhalo indexes to zero
       call MPI_SEND(dhalo_index, ntrees, MPI_INTEGER8, 0, 0, MPI_COMM_WORLD, ierr)
       ! And then receive the file index in return
       call MPI_RECV(ifile, ntrees, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, ierr)
    else
       ! Processor 0 draws a random file number for each unique dhalo index
       current_file  = -1
       current_dhalo = -1
       do iproc = 0, numprocs-1, 1
          ! Get dhalo indexes from next processor
          if(iproc.gt.0)then
             deallocate(dhalo_index)
             call MPI_PROBE(iproc, 0, MPI_COMM_WORLD, status, ierr)
             call MPI_GET_COUNT(status, MPI_INTEGER8, nrecv, ierr)
             allocate(dhalo_index(nrecv))
             call MPI_RECV(dhalo_index, nrecv, MPI_INTEGER8, iproc, 0, &
                  MPI_COMM_WORLD, status, ierr)
          else
             nrecv = ntrees
          endif
          ! Assign file indexes
          allocate(ifile_send(nrecv))
          do i = 1, nrecv, 1
             if(current_dhalo.ne.dhalo_index(i))then
                ! New final Dhalo - pick a file to put it in
                call random_number(r8)
                current_file  = floor(r8*nfiles)
                if(current_file.lt.0)current_file = 0
                if(current_file.ge.nfiles)current_file = nfiles-1
                current_dhalo = dhalo_index(i)
             endif
             ifile_send(i) = current_file
          end do
          ! Send file index back
          if(iproc.gt.0)then
             call MPI_SEND(ifile_send, nrecv, MPI_INTEGER, iproc, 0, MPI_COMM_WORLD, ierr)
          else
             ifile = ifile_send
          endif
          deallocate(ifile_send)
          ! Next processor
       end do
    endif
    deallocate(dhalo_index)

    ! Now need to return ifile to original order
    call assign_elements(merger_tree%ifile, ifile, idx)
    
    deallocate(idx, ifile)

    return
  end subroutine assign_final_dhalos_to_files

end module merger_trees
