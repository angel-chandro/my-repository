module subhalo_data

  use kind_numbers
  use sort

  implicit none
  public
  save

  type subhalo_type
     ! If this is modified must make corresponding changes to MPI type
     ! definition in subhalo_data_init() below.
     sequence
     integer(kind=int8byte) :: mostboundid      ! INPUT: ID of most bound particle in this subhalo
     integer(kind=int8byte) :: id               ! INPUT: ID of this subhalo
     integer(kind=int8byte) :: descendant_id    ! INPUT: ID of the descendant subhalo
     integer(kind=int8byte) :: parent_id        ! OUTPUT: used to find which subhalos "enclose" which others
     integer(kind=int8byte) :: dhalo_index      ! OUTPUT: which DHalo this subhalo belongs to
     integer(kind=int8byte) :: next_dhalo_index ! OUTPUT: which DHalo descendant of this subhalo belongs to
     integer(kind=int8byte) :: half_max_mass_id ! OUTPUT: ID of subhalo on main prog branch with half of max. mass
     integer(kind=int8byte) :: main_prog_id     ! OUTPUT: ID of main progenitor
     integer            :: snapnum              ! OUTPUT: snapshot at which this subhalo exists
     integer            :: descendant_snapnum   ! INPUT: snapshot at which descendant exists
     integer            :: fof_index            ! INPUT: which FoF group this subhalo belongs to
     integer            :: mbps_contributed     ! OUTPUT: number of descendants most bound particles from this subhalo 
     real, dimension(3) :: pos                  ! INPUT: position of subhalo (potential minimum, Mpc/h)
     real, dimension(3) :: vel                  ! INPUT: velocity of subhalo (peculiar, km/sec)
     real, dimension(3) :: spin                 ! INPUT (not used yet): angular momentum of subhalo
     real, dimension(3) :: cofm                 ! INPUT: centre of mass of subhalo (Mpc/h)
     real               :: mvir                 ! INPUT : virial mass of the subhalo (Msun/h)
     real               :: mgas                 ! INPUT : gas mass of the subhalo (Msun/h)
     real               :: veldisp              ! INPUT (not used yet): velocity dispersion of subhalo (km/sec)
     real               :: cnfw                 ! INPUT: concentration NFW profile.
     real               :: Vvir                 ! INPUT: virial velocity (km/s)
     real               :: lambda               ! INPUT: lambda Bullock.
     real               :: vmax                 ! INPUT (not used yet): max circ.vel of subhalo (km/sec)
     real               :: rvmax                ! INPUT (not used yet): radius of max circ.vel of subhalo (Mpc/h, comoving)
     real               :: rhalf                ! INPUT: half mass radius of subhalo (Mpc/h, comoving)
     integer            :: np                   ! INPUT: number of particles in subhalo
     integer            :: np_orig              ! OUTPUT: number of particles last time subhalo was isolated
     integer            :: flags                ! INPUT and OUTPUT: various flags stored as single bits
                                                ! Read routine should set this to 2 for central subhalo in FoF
                                                ! group, 0 otherwise.
     integer            :: index                ! INPUT (not used): position of subhalo in input catalogue
     real               :: max_mass             ! OUTPUT: maximum mass along past main branch
     real               :: max_vmax             ! OUTPUT: maximum vmax along past main branch (km/sec)
     real               :: mass                 ! INPUT: total subhalo mass
     real               :: mass_orig            ! OUTPUT: mass last time subhalo was isolated
     real               :: dummy
  end type subhalo_type

  type subhalo_tree_node
     type (subhalo_type) :: subhalo
     type (subhalo_tree_node), pointer :: descendant
     type (subhalo_tree_node), pointer :: progenitor
     type (subhalo_tree_node), pointer :: next_progenitor
     type (subhalo_tree_node), pointer :: next
     integer :: itree ! Which subhalo tree this subhalo belongs to
  end type subhalo_tree_node

  ! Type to store a linked list of blocks of subhalo tree nodes
  integer, parameter :: blocksize = 10000
  type node_block_type
     integer                                         :: nnodes
     type(node_block_type), pointer                  :: next
     type (subhalo_tree_node), dimension(:), pointer :: node
  end type node_block_type

  ! MPI data type handle
  integer :: SUBHALO_MPI_TYPE

  ! Bit masks for flags
  integer, parameter :: MAIN_PROGENITOR_FLAG = ishft(1, 0)
  integer, parameter :: FOF_CENTRE_FLAG      = ishft(1, 1)
  integer, parameter :: DHALO_CENTRE_FLAG    = ishft(1, 2)
  integer, parameter :: IN_TREE_FLAG         = ishft(1, 3)
  integer, parameter :: INTERPOLATED_FLAG    = ishft(1, 4)
  integer, parameter :: REMERGED_FLAG        = ishft(1, 5)
  integer, parameter :: WDM_SPURIOUS_FLAG    = ishft(1, 6)
  integer, parameter :: CHECKED_FLAG         = ishft(1, 7)

  ! Pointers to subhalo tree nodes at the start and end of the linked
  ! list for each snapshot
  type linked_list_start
     type (subhalo_tree_node), pointer :: first
     type (subhalo_tree_node), pointer :: last
  end type linked_list_start

  ! Pointers to first subhalo at each snapshot
  type (linked_list_start), dimension(:), allocatable :: node_at_snapshot

  ! This is used to make an array of pointers to nodes
  type tree_node_ptr
     type (subhalo_tree_node), pointer :: ptr
  end type tree_node_ptr

  ! Redshift list
  real, dimension(:), allocatable :: zred
  real, dimension(:), allocatable :: sfactor ! include

  ! 
  ! Internal types / data
  !

  ! This allows us to make an array of arrays of pointers (!)
  type tree_node_ptr_array
     type (tree_node_ptr), dimension(:), pointer :: node
  end type tree_node_ptr_array
  private :: tree_node_ptr_array

  ! Pointers to possible descendants
  type (tree_node_ptr_array), dimension(:), allocatable :: desc_at_snapshot
  private :: desc_at_snapshot

  ! Pointer to start of linked list of node blocks and current block
  ! Current node block
  type(node_block_type), pointer :: current_block
  type(node_block_type), pointer :: first_block

contains

  subroutine subhalo_data_init(ifirst, ilast)
!
! Register an MPI type to store a subhalo and allocate pointers
! to the first subhalo at each snapshot
!
    use mpi
    implicit none
    integer, intent(in) :: ifirst, ilast
    type (subhalo_type), dimension(2) :: sub
    integer(MPI_ADDRESS_KIND) :: base, offset
    integer :: extent
    integer, parameter :: countmax = 100
    integer, dimension(countmax) :: blocklen, disp, typearr
    integer :: i, ierr, isnap

    ! Figure out size of subhalo_type (may include padding)
    call MPI_GET_ADDRESS(sub(1), base,   ierr)
    call MPI_GET_ADDRESS(sub(2), offset, ierr)
    extent = int(offset-base,kind(extent))

    ! MPI type consists of a series of blocks with length, type and
    ! displacement for each one
    i = 0

    ! MostboundID
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%mostboundid, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

   ! ID
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%id, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! DescendantID
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%descendant_id, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! ParentID
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%parent_id, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! DHaloID
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%dhalo_index, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Next DHalo index
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%next_dhalo_index, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Half max mass index
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%half_max_mass_id, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Main progenitor ID
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER8
    call MPI_GET_ADDRESS(sub(1)%main_prog_id, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Snapnum
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%snapnum, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Descendant snapnum
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%descendant_snapnum, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! FoF index
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%fof_index, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! No. of MBPs contributed to descendant
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%mbps_contributed, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Position
    i = i + 1
    blocklen(i) = 3
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%pos, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Velocity
    i = i + 1
    blocklen(i) = 3
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%vel, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Spin
    i = i + 1
    blocklen(i) = 3
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%spin, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! CofM
    i = i + 1
    blocklen(i) = 3
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%cofm, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Velocity dispersion
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%veldisp, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Concentration
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%cnfw, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Gas mass
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%mgas, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! lambda
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%lambda, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Vvir
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%Vvir, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Vmax
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%vmax, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! RVmax
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%rvmax, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Rhalf
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%rhalf, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Np
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%np, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Np_orig
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%np_orig, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Flags
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%flags, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Index
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_INTEGER
    call MPI_GET_ADDRESS(sub(1)%index, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Max mass
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%max_mass, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Max vmax
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%max_vmax, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Mass
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%mass, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Mass_orig
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_REAL
    call MPI_GET_ADDRESS(sub(1)%mass_orig, offset, ierr)
    disp(i)     = int(offset-base,kind(disp))

    ! Mark end of type
    i = i + 1
    blocklen(i) = 1
    typearr(i)  = MPI_UB
    disp(i)     = extent

    ! Register the new type
    call MPI_TYPE_STRUCT(i, blocklen, disp, typearr, SUBHALO_MPI_TYPE, ierr)
    call MPI_TYPE_COMMIT(SUBHALO_MPI_TYPE, ierr)

    !
    ! Allocate pointers to tree nodes
    !
    allocate(node_at_snapshot(ifirst:ilast))
    do isnap = ifirst, ilast, 1
       nullify(node_at_snapshot(isnap)%first)
       nullify(node_at_snapshot(isnap)%last)
    end do

    ! Allocate array of arrays of pointers to possible descendants
    allocate(desc_at_snapshot(ifirst:ilast))
    do isnap = ifirst, ilast, 1
       nullify(desc_at_snapshot(isnap)%node)
    end do

    ! Initialise list of blocks
    nullify(first_block)
    nullify(current_block)

    return
  end subroutine subhalo_data_init


  function new_tree_node(subhalo, isnap)
!
! Allocate a new tree node and return a pointer to it
!
    implicit none
    integer,                  intent(in) :: isnap
    type (subhalo_type),      intent(in) :: subhalo
    type (subhalo_tree_node), pointer    :: new_tree_node

    ! Allocate first block of nodes if necessary
    if(.not.associated(first_block))then
       allocate(first_block)
       allocate(first_block%node(blocksize))
       current_block => first_block
       current_block%nnodes = 0
       nullify(current_block%next)
    endif

    ! Make sure there's at least one node left in the current
    ! block - may need to allocate a new one
    if(current_block%nnodes.ge.blocksize)then
       ! Need to make a new block
       allocate(current_block%next)
       allocate(current_block%next%node(blocksize))
       nullify(current_block%next%next)
       current_block%next%nnodes = 0
       current_block => current_block%next
    endif

    ! Use the next available node in the current block,
    current_block%nnodes = current_block%nnodes + 1
    new_tree_node => current_block%node(current_block%nnodes)

    ! Initialise this node
    new_tree_node%subhalo = subhalo
    nullify(new_tree_node%descendant)
    nullify(new_tree_node%progenitor)
    nullify(new_tree_node%next_progenitor)
    nullify(new_tree_node%next)

    ! Add this node to the linked list for this snapshot
    if(associated(node_at_snapshot(isnap)%first))then
       node_at_snapshot(isnap)%last%next => new_tree_node
    else
       node_at_snapshot(isnap)%first => new_tree_node
    endif
    node_at_snapshot(isnap)%last => new_tree_node

    return
  end function new_tree_node


  subroutine subhalo_data_deallocate()
!
! Dealloacte the subhalo merger tree data
!
    implicit none
    integer :: ierr, jsnap
    type(node_block_type),    pointer :: this_block, prev_block

    ! Deallocate pointers to descendants
    do jsnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
       if(associated(desc_at_snapshot(jsnap)%node))then
          deallocate(desc_at_snapshot(jsnap)%node)
          nullify(desc_at_snapshot(jsnap)%node)
       endif
    end do
    deallocate(desc_at_snapshot)

    ! Deallocate pointers to start/end of linked lists
    deallocate(node_at_snapshot)

    ! Deallocate blocks of nodes
    this_block => first_block
    do while(associated(this_block))
       prev_block => this_block
       this_block => this_block%next
       deallocate(prev_block)
    end do

    ! Deallocate MPI type we defined
    call MPI_TYPE_FREE(SUBHALO_MPI_TYPE, ierr)

    return
  end subroutine subhalo_data_deallocate


  subroutine add_to_tree(isnap, subhalo, nsteptrace, ilast,&
       find_all)
!
! Add subhalos for snapshot isnap to the merger trees. The IN_TREE flag
! gets set for any subhalos which are added to a tree. This is so that the
! remaining subhalos can be made into the "roots" of new trees later.
!
    use mpi
    use sort
    implicit none
    ! Parameters
    integer, intent(in) :: isnap, nsteptrace, ilast
    type (subhalo_type), dimension(:), intent(inout) :: subhalo
    logical :: find_all
    ! Subhalo send/receive buffers
    type (subhalo_type), dimension(:), allocatable :: sub_send
    type (subhalo_type), dimension(:), allocatable :: sub_recv
    ! Number of subhalos
    integer :: nsub, nsub_max, nsub_recv
    ! MPI stuff
    integer :: myid, numprocs, ierr
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: isrc, idest, iproc
    ! Sorting
    integer, dimension(:), allocatable :: idx, prog_idx
    integer(kind=int8byte), dimension(:), allocatable :: id
    ! Pointers to tree nodes
    type (subhalo_tree_node), pointer :: node
    type (tree_node_ptr), dimension(:), pointer :: node_array
    integer :: ndesc
    ! Loops etc
    integer :: isub, jsnap, i, j
    integer :: last_snapshot, last_snapshot_local
    ! New tree node pointer
    type (subhalo_tree_node), pointer    :: new_node
    integer :: nmatched, n_with_desc
    integer :: nmatched_tot, n_with_desc_tot, nsub_tot

    ! Find which processor this is
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

    ! Deallocate any descendant pointer arrays we don't need any more -
    ! any which are more than nsteptrace snapshots later than snapshot isnap.
    do jsnap = isnap + nsteptrace + 1, ilast, 1
       if(associated(desc_at_snapshot(jsnap)%node))then
          deallocate(desc_at_snapshot(jsnap)%node)
          nullify(desc_at_snapshot(jsnap)%node)
       endif
    end do
    
    ! Get maximum number of subhalos across all processors and allocate buffers
    nsub = size(subhalo)
    call MPI_ALLREDUCE(nsub, nsub_max, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierr)
    allocate(sub_send(nsub_max))
    allocate(sub_recv(nsub_max))

    ! Populate send buffer with progenitor subhalos ordered by descendant ID.
    ! Keep the sort index so we can copy updated flags back into the subhalo
    ! array later.
    n_with_desc = 0
    allocate(prog_idx(nsub))
    call sort_index(nsub, subhalo%descendant_id, prog_idx)
    do isub = 1, nsub, 1
       sub_send(isub) = subhalo(prog_idx(isub))
       if(sub_send(isub)%descendant_id.ge.0)n_with_desc = n_with_desc + 1
    end do

    ! Get index of processor to send to
    idest = myid + 1
    if(idest.ge.numprocs)idest = idest - numprocs
    ! Get index of processor to receive from
    isrc = myid - 1
    if(isrc.lt.0)isrc = isrc + numprocs

    ! Find the last snapshot at which there is a descendant
    last_snapshot_local = maxval(sub_send(1:nsub)%descendant_snapnum)
    call MPI_ALLREDUCE(last_snapshot_local, last_snapshot, &
         1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Make arrays of pointers to subhalos at possible descendant snapshots.
    ! These need to be sorted by subhalo ID.
    ! Loop over descendant snapshots
    do jsnap = isnap+1, last_snapshot, 1
       ! May have already set up array for this snapshot on a previous
       ! iteration
       if(associated(desc_at_snapshot(jsnap)%node))cycle 
       ! Count nodes at this snapshot
       ndesc = 0
       node => node_at_snapshot(jsnap)%first
       do while(associated(node))
          ndesc = ndesc + 1
          node => node%next
       end do
       ! Make array of pointers to descendants at this snapshot
       allocate(desc_at_snapshot(jsnap)%node(ndesc))
       ndesc = 0
       node => node_at_snapshot(jsnap)%first
       do while(associated(node))
          ndesc = ndesc + 1
          desc_at_snapshot(jsnap)%node(ndesc)%ptr => node
          node => node%next
       end do
       ! Sort pointers so we can access subhalos in order of ID
       allocate(node_array(ndesc), idx(ndesc), id(ndesc))
       do i = 1, ndesc, 1
          node_array(i)%ptr => desc_at_snapshot(jsnap)%node(i)%ptr
          id(i)             =  desc_at_snapshot(jsnap)%node(i)%ptr%subhalo%id
       end do
       call sort_index(ndesc, id, idx)
       do i = 1, ndesc, 1
          desc_at_snapshot(jsnap)%node(i)%ptr => node_array(idx(i))%ptr
       end do
       deallocate(node_array, idx, id)
    end do

    ! Pass the subhalos around all processors in sequence
    nmatched = 0
    do iproc = 0, numprocs-1, 1
       ! Pass local data to the next processor and receive from the
       ! previous one
       call MPI_SENDRECV(&
            sub_send, nsub,     SUBHALO_MPI_TYPE, idest, 0, &
            sub_recv, nsub_max, SUBHALO_MPI_TYPE, isrc,  0, &
            MPI_COMM_WORLD, status, ierr)
       call MPI_GET_COUNT(status, SUBHALO_MPI_TYPE, nsub_recv, ierr)
       ! Check if these imported subhalos can be matched up to descendants
       ! on this processor
       do jsnap = isnap+1, last_snapshot, 1
          ! May have no possible descendants on this processor
          if(size(desc_at_snapshot(jsnap)%node).lt.1)cycle
          j = 1 ! Descendant subhalo index
          ! Loop over imported (progenitor) subhalos
          do i = 1, nsub_recv, 1
             if(sub_recv(i)%descendant_snapnum.eq.jsnap)then
                do while(desc_at_snapshot(jsnap)%node(j)%ptr%subhalo%id.lt.&
                     sub_recv(i)%descendant_id.and. &
                     j.lt.size(desc_at_snapshot(jsnap)%node))
                   j = j + 1
                end do
                if(desc_at_snapshot(jsnap)%node(j)%ptr%subhalo%id.eq.&
                     sub_recv(i)%descendant_id)then
                   ! Found a match - add this subhalo to the merger tree
                   new_node => new_tree_node(sub_recv(i), isnap)
                   new_node%descendant => desc_at_snapshot(jsnap)%node(j)%ptr
                   ! Check we didn't already assign this one to a tree
                   if(iand(sub_recv(i)%flags, IN_TREE_FLAG).ne.0) &
                        call panic("Attempt to assign a subhalo to >1 tree!")
                   ! Flag it as matched
                   sub_recv(i)%flags = &
                        ior(sub_recv(i)%flags, IN_TREE_FLAG)
                   nmatched = nmatched + 1
                endif
             end if
             ! Next progenitor subhalo
          end do
          ! Next descendant snapshot
       end do

       ! Received subhalos will be sent on next iteration
       nsub = nsub_recv
       sub_send(1:nsub) = sub_recv(1:nsub)
    end do

    call MPI_ALLREDUCE(nmatched, nmatched_tot, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(nsub, nsub_tot, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(n_with_desc, n_with_desc_tot, 1, MPI_INTEGER, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(myid.eq.0)then
       write(*,*)"  Found descendants for ", nmatched_tot, " of ", nsub_tot, &
            "subhalos"
    endif

    if(find_all)then
       ! SubFind input - should be able to find all descendants
       if(nmatched_tot.ne.n_with_desc_tot) &
            call panic("No. subhalos matched not equal to no. with descendants")
    endif

    ! Should have original set of subhalos back on this processor now.
    ! Copy updated flags indicating which subhalos have been added to a tree.
    do isub = 1, nsub, 1
       subhalo(prog_idx(isub))%flags = sub_send(isub)%flags
    end do
    deallocate(prog_idx)

    ! Deallocate communication buffers
    deallocate(sub_send, sub_recv)

    return
  end subroutine add_to_tree


  subroutine sort_subhalos(isnap, sort_nodes, next_dhalo)
!
! Return an array of pointers to access subhalos sorted by DHaloID.
! Allocates sort_nodes.
!
    implicit none
    ! Snapshot to do
    integer, intent(in) :: isnap
    logical :: next_dhalo
    ! Pointers to subhalos
    type (tree_node_ptr), dimension(:), allocatable :: sort_nodes, tmp_nodes
    ! Count of local subhalos and max over all processors
    integer :: nsub
    ! Current/next dhalo index
    integer(kind=int8byte), dimension(:), allocatable :: cur_dhalo_index
    ! Sorting index
    integer, dimension(:), allocatable :: idx
    ! Pointer to a node
    type (subhalo_tree_node), pointer :: node
    ! Loops etc
    integer :: i

    ! Generate array of pointers to access local subhalos in order of
    ! current dhalo index. First count subhalos.
    nsub = 0
    node => node_at_snapshot(isnap)%first
    do while(associated(node))
       nsub = nsub + 1
       if(associated(node%descendant))then
          ! Check interpolation has been done
          if(node%subhalo%descendant_snapnum.ne.node%subhalo%snapnum+1) &
               call panic( &
               "Descendant of subhalo is not at the next snapshot!")
       endif
       node => node%next
    end do

    ! Then store dhalo indexes
    allocate(cur_dhalo_index(nsub), idx(nsub))
    nsub = 0
    node => node_at_snapshot(isnap)%first
    do while(associated(node))
       nsub = nsub + 1
       if(next_dhalo)then
          cur_dhalo_index(nsub) = node%subhalo%next_dhalo_index
       else
          cur_dhalo_index(nsub) = node%subhalo%dhalo_index
       endif
       node => node%next
    end do
       
    ! And make sorting index by current dhalo
    call sort_index(nsub, cur_dhalo_index, idx)
       
    ! Store pointers to subhalos
    allocate(tmp_nodes(nsub))
    nsub = 0
    node => node_at_snapshot(isnap)%first
    do while(associated(node))
       nsub = nsub + 1
       tmp_nodes(nsub)%ptr => node
       node => node%next
    end do

    ! Sort pointers by current dhalo index
    allocate(sort_nodes(nsub))
    do i = 1, nsub, 1
       sort_nodes(i)%ptr => tmp_nodes(idx(i))%ptr
    end do
    deallocate(tmp_nodes, idx, cur_dhalo_index)

    return
  end subroutine sort_subhalos

end module subhalo_data
