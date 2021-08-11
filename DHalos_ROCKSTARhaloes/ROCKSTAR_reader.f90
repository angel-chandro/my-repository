! For reading the Cluster information
module ROCKSTAR_reader

  use mpi
  use sub_io
  use subhalo_properties
  use subhalo_data
  use read_subfind
  use hdf5_wrapper
  use kind_numbers
  use subhalo_ids
  !use code_version
  use relabel_snapshots
  use VELOCIraptor_reader

  implicit none
  private

  public :: read_ROCKSTAR, ROCKSTAR_find_main_progenitors

  integer, parameter :: MAXFILES = 4096

contains


  subroutine read_ROCKSTAR(ROCKSTAR_path, isnap, ilast, subhalo, mpart, lbox)
    !
    ! Read in subhalos for the specified snapshot
    !
    ! Input units (as described in the AFH documentation)
    ! -----------------------------------------------
    !
    ! Positions              - comoving kpc/h
    ! Sizes (e.g. Rhalf)     - comoving kpc/h
    ! Velocities of subhalos - km/sec peculiar (a \dot{comoving position})
    ! Velocity dispersion    - km/sec
    ! Max. circular velocity - km/sec
    ! Spin                   - ?
    ! Masses                 - Msun/h     
    !
    ! Output units:
    ! -------------
    !
    ! All positions and lengths  - comoving Mpc/h
    ! Velocities of subhalos     - km/sec peculiar
    ! Veldisp/Vmax               - km/sec
    ! Spin                       - ?
    ! Masses                 - Msun/h     
    !
    implicit none
    ! Input parameters
    character(len=*), intent(in) :: ROCKSTAR_path
    integer,          intent(in) :: isnap, ilast 
    real,             intent(in) :: mpart,lbox
    ! Output subhalo catalogue
    type (subhalo_type), dimension(:), pointer :: subhalo
    ! Groups to read on this processor
    !integer, intent(in) :: nsubd
    integer :: nhead, ihead
    integer :: isub, nsub
    !integer, intent(in) :: dsnap(:)
    !integer(kind=int8byte), intent(in) :: haloid(:), descid(:)
    ! Loops etc
    integer :: i1, i2, i, nsnap, ios, idesc, innpart
    !
    character(Len=1024):: subfile,descfil,snum,sf,zz,toth
    integer(kind=int8byte) :: inhalo,hhalo,id,id_inhalo,id_hhalo,idhalo
    integer(kind=int8byte), parameter :: IDCorr = 1000000000000_8

    real :: t1,t2, t3
    real :: inbmass,inx,iny,inz,invx,invy,invz, &
         rvir,rmax,r2,mbp_offset,com_offset,invmax,v_esc,vdispersion,lambda,&
         lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Eax500c,Eay500c,Eaz500c,&
         b500c,c500c,T_U_ratio,M_200Mean,M_200Crit,M500c,M2500c,&
         sam_Mvir,MpeBehroozi,MpeDIemer
    
 
    ! Get the snapshot input file
    if(isnap.le.ilast) then  
       nsnap = int(file_snapnum(isnap))
       write(snum,'(I3)') nsnap
     else
       return
    endif

    !One specific file to test
    subfile=trim(ROCKSTAR_path)//'/out_'//trim(adjustl(snum))//'_parents.list'
    !subfile='/home/chandro/ROCKSTAR_parents/out_30_parents.list'
    print*,trim(subfile)


    !----------------------------------------
    !Read subhalo information
    open(unit=20,file=trim(subfile),status='old')
    isub = 0 ; idesc = 0
    ios = 1 ; nhead = 0
    ! Count the subhaloes at a given redshift
    do while (ios.ne.0)
       read(20,*,iostat=ios) inhalo
       nhead = nhead + 1
    enddo
    close(20)
    ! nhead: number of lines of header in haloes file
    open(unit=20,file=trim(subfile),status='old')
    ! skip header
    do ihead = 1,nhead-1
       read(20,*,iostat=ios)
    enddo
    ! read haloes
    do while (ios.eq.0)
       read(20,*,iostat=ios) inhalo,idhalo
       if(ios.ne.0)exit
       if(inhalo.ge.0)then
          isub = isub +1
          if(idhalo.ge.0)then
             ! because if idhalo=-1 it means that halo hasn't descendants
             ! Count number of haloes with descendants
             idesc = idesc + 1 
          endif

       !else if(isfactor.gt.sfactor(nsnap))then ! Descendant's file ordered by increasing snapshot
       !   exit
       !endif
       !isub = isub + 1
       
       else !Sanity check
          print*,'STOP! (B) Halo id<0: ',isub,inhalo  ; stop
       endif
    enddo
    close(20)
    
    ! Number of subhaloes here should be equal to that in the descendants file
    nsub = isub ;
    print*,nsub,' subhaloes in this snapshot, ',&
         idesc,' with descendants.' 
    !if(nsub.ne.nsubd)then
    !   print*,'STOP: ',nsub,' is not equal to ',nsubd ; stop
    !endif
    

    ! Allocate output array
    allocate(subhalo(nsub))
    subhalo%flags       = 0

    !Get the information from the file
    open(unit=20,file=trim(subfile),status='old')
    ios = 0
    isub = 0 ; idesc = 0 ; t3 = 0.0  
    ! skip header
    do ihead = 1,nhead-1
       read(20,*,iostat=ios)
    enddo

    do while (ios.eq.0)

       !read(20,*,iostat=ios) inhalo,hhalo,nss,inbmass,innpart, &
       !     inx,iny,inz,invx,invy,invz,rvir,rmax,r2, &
       !     mbp_offset,com_offset,invmax,v_esc,vdispersion,&
       !     lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,Ecy,Ecz,&
       !     ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW,FoFMass,M_200Mean,&
       !     M_200Crit,M_TopHat,R_200Mean,R_200Crit
       read(20,*,iostat=ios) inhalo,idhalo,inbmass, &
            invmax,vdispersion,rvir,rmax,innpart, &
            inx,iny,inz,invx,invy,invz,Lx,Ly,Lz, &
            lambda,r2,sam_Mvir,M_200Mean,M_200Crit,M500c,M2500c, &
            mbp_offset,com_offset,lambdaE,b,c,Eax,Eay,Eaz, &
            b500c,c500c,Eax500c,Eay500c,Eaz500c,T_U_ratio, &
            MpeBehroozi,MpeDIemer,hhalo
       if(ios.ne.0)exit

       isub = isub + 1 
       !if(mod(isub,10000).eq.0)print*,isub,' subhaloes read.'

       !Subhalo mass
       subhalo(isub)%mvir    = inbmass

       !Positions
       subhalo(isub)%pos(1)  = inx/1000.
       subhalo(isub)%pos(2)  = iny/1000.
       subhalo(isub)%pos(3)  = inz/1000.

       !Peculiar velocities
       subhalo(isub)%vel(1)  = invx
       subhalo(isub)%vel(2)  = invy
       subhalo(isub)%vel(3)  = invz

       !Velocity dispersion
       subhalo(isub)%veldisp = vdispersion
       subhalo(isub)%cnfw    = -1 
       subhalo(isub)%lambda  = lambda 
       subhalo(isub)%Vvir    = -1
       subhalo(isub)%mgas    = -1
       subhalo(isub)%vmax    = invmax

       !Number of particles
       subhalo(isub)%np      = innpart !innpart. Rockstar usually write
                                                  !the number of particles that are uniquely assigned to a halo. However, after
                                                  !identifying this number, there is an unbinding procedure and Mvir is calculated
                                                  !after that. Thus from these inputs Mvir/innpart is different to mpart. The right
                                                  !thing here is therefore to redefine the number of particles that are attached to
                                                  !that halo.
       
       !Spin
       subhalo(isub)%spin(1)  = Lx
       subhalo(isub)%spin(2)  = Ly
       subhalo(isub)%spin(3)  = Lz
       
       ! Rhalf with an rfactor=1
       ! Comoving coordinates (Mpc/h) to be passed
       subhalo(isub)%rhalf = rvir/2000.0

       !No information on most bounded particle
       subhalo(isub)%mostboundid = 0
       subhalo(isub)%mbps_contributed = 0

       ! Calculate the full ID of this subhalo if it the one
       ! read from the file doesn't look like a full ID
       if(inhalo.lt.1e9) then
          id = inhalo + nsnap*IDCorr
       else
          id = inhalo
       endif
       ! if DHalos needs a special ID format, IDCorr will have to be used
       ! (search in ROCKSTAR info)
       
       ! Set the halo index
       subhalo(isub)%id = inhalo
        
       ! Set the host halo index
       if (hhalo.lt.0)then
          hhalo = inhalo      
          subhalo(isub)%flags = ior(subhalo(isub)%flags, FOF_CENTRE_FLAG)
          print*,'hhalo,inhalo,subhalo(isub)%flags',hhalo,inhalo,subhalo(isub)%flags
       endif

       ! Remove the snapshot information of the ID of this subhalo if the
       ! one
       ! read from the file looks like a full ID
       if(inhalo.gt.1e9) then
          id_inhalo = inhalo - nsnap*IDCorr
          id_hhalo  = hhalo  - nsnap*IDCorr
          !print*,'También pasa por aquí:',id_inhalo,id_hhalo
       else
          id_inhalo = inhalo
          id_hhalo  = hhalo
       end if
       ! again if DHalos needs a special ID format, IDCorr will have to be used 
       ! (search in ROCKSTAR info)
       
       subhalo(isub)%fof_index   = int(hhalo,4)

       ! Position of subhalo in input catalogue
       subhalo(isub)%index       = int(inhalo,4)

       ! print*,'subhalo(isub)%index,id',subhalo(isub)%index,id
       ! Set various indexes
       subhalo(isub)%dhalo_index = 0
       subhalo(isub)%snapnum     = isnap

       !Read in descendant info
       !call cpu_time ( t1 )
       !call find_in_descendants(subhalo(isub), nsubd, idesc, haloid, descid, dsnap)
       !call cpu_time ( t2 ) ; t3 = t3 + t2 - t1
       ! ROCKSTAR does not need descendants file, because descendants are included in the proper file
       
       subhalo(isub)%descendant_id = idhalo
       subhalo(isub)%descendant_snapnum = isnap+1

       idesc = idesc+1
       
       if(subhalo(isub)%id.lt.0) then !Sanity check
          print*,'STOP! (C) Halo id<=0: ',isub,subhalo(isub)%id ; stop
       endif
    enddo
    print*,'Elapsed CPU time finding descendants = ', t3
    print*,maxval(subhalo(1:nsubd)%descendant_snapnum)
    close(20)

    if (idesc.ne.isub)then !Sanity check
       print*,'STOP! The number of descendants does not match the number of subhaloes' ; stop
    endif
!
    return
  end subroutine read_ROCKSTAR

  subroutine ROCKSTAR_find_main_progenitors()
!
! Sort progenitor linked lists into descending order of number of
! most massive
!
    use sort
    use multi_key_sort

    implicit none

    type (subhalo_tree_node), pointer :: node, this_prog
    integer :: isnap, i
    integer :: nprog
    type (tree_node_ptr), dimension(:), allocatable :: progenitor
    integer(int8byte), dimension(:), allocatable :: sort_key, subid
    integer(int8byte) :: max_np
    integer, dimension(:), allocatable :: idx

    ! Loop over snapshots
    do isnap = ubound(node_at_snapshot,1), lbound(node_at_snapshot,1)+1, -1
       ! Loop over subhalos at this snapshot
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
          if(associated(node%progenitor))then
             ! Check if node has >1 progenitor
             if(associated(node%progenitor%next_progenitor))then
                ! Count progenitors of this node
                nprog = 0
                this_prog => node%progenitor
                do while(associated(this_prog))
                   nprog = nprog + 1
                   this_prog => this_prog%next_progenitor
                end do
                ! Store pointers to progenitors
                allocate(progenitor(nprog))
                nprog = 0
                this_prog => node%progenitor
                do while(associated(this_prog))
                   nprog = nprog + 1
                   progenitor(nprog)%ptr => this_prog
                   this_prog => this_prog%next_progenitor
                end do
                ! Sort progenitors by Mvir, and then
                ! ID
                allocate(sort_key(nprog), idx(nprog), subid(nprog))
                max_np = -1
                do i = 1, nprog, 1
                   max_np = max(max_np, progenitor(i)%ptr%subhalo%np)
                end do
                do i = 1, nprog, 1
                   sort_key(i) = -1*int(progenitor(i)%ptr%subhalo%np, int8byte) ! Descending order
                   subid(i) = progenitor(i)%ptr%subhalo%id
                end do
                call two_key_sort_index(nprog, sort_key, subid, idx)
                ! Rebuild linked list in sorted order
                node%progenitor => progenitor(idx(1))%ptr
                this_prog => node%progenitor
                do i = 2, nprog, 1
                   this_prog%next_progenitor => progenitor(idx(i))%ptr
                   this_prog => this_prog%next_progenitor
                end do
                ! Mark end of linked list
                nullify(this_prog%next_progenitor)
                deallocate(progenitor, sort_key, idx, subid)
             endif
          endif
          ! Next node at this snapshot
          node => node%next
       end do
       ! Next snapshot
    end do

    ! Set subhalo main progenitor flags
    ! First set them all to zero
    do isnap = lbound(node_at_snapshot,1), ubound(node_at_snapshot,1), 1
       ! Loop over subhalos at this snapshot
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
          node%subhalo%flags = iand(node%subhalo%flags, not(MAIN_PROGENITOR_FLAG))
          ! Next node at this snapshot
          node => node%next
       end do
       ! Next snapshot
    end do

    ! Loop over snapshots
    do isnap = ubound(node_at_snapshot,1), lbound(node_at_snapshot,1)+1, -1
       ! Loop over subhalos at this snapshot
       node => node_at_snapshot(isnap)%first
       do while(associated(node))
          ! Check if subhalo has a progenitor
          if(associated(node%progenitor))then
             if(node%progenitor%subhalo%np.gt.0)then
                ! First progenitor is main if it contributed some mass
                node%progenitor%subhalo%flags = &
                     ior(node%progenitor%subhalo%flags, &
                     MAIN_PROGENITOR_FLAG)                
             endif
          endif
          ! Next node at this snapshot
          node => node%next
       end do
       ! Next snapshot
    end do

    return
  end subroutine ROCKSTAR_find_main_progenitors

end module ROCKSTAR_reader
