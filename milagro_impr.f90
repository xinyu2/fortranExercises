!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
program milagroimpr
  implicit none
  INCLUDE 'mpif.h'
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer,parameter :: impi0=0 !the master rank
  integer,parameter :: maxpars=8*2!1024*1024
  integer,parameter :: BUFFSIZE=8,MAXSTEP=1
  integer,parameter :: dimx=8,dimy=2,dimz=1
  integer,parameter :: rmpbtag=5,rpctag=10,rftag=15,sndtag=20,rbftag=20
  integer,parameter :: dmpi=1 !debug this mpi rank

  integer,dimension(:),allocatable :: maxbfsize  ! max particle buffer size on each ranks
  integer,dimension(:),allocatable :: pcomplete  ! array for number of particle complete on master
  logical :: globalFinish  ! global finish tag on slaves
  integer :: pcmplt        ! number of particle complete on slave(cross boundary)

  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,it,istat(MPI_STATUS_SIZE),i
  !********************************
  !* reqs:
  !* reqmbs: max buffer size
  !* reqrb:  receive buffer
  !* reqpc:  particle complete
  !* reqgf:  global finish
  !* rbflag: receive buffer flag
  !********************************
  integer,dimension(:),allocatable :: reqrb,reqpc,cmpind
  integer :: reqgf,outreq
  logical,dimension(:),allocatable :: rbflag,rpcflag
  logical :: rgfflag
  !*********
  !* domain
  !*********
  integer :: dsize !domain size(total cells=dimx*dimy*dimz)
  integer :: rsize !rank size(cells per rank=dimx*dimy*dimz/nmpi)
  !integer :: parsub !subset of pars, particle per rank
  type cell
     integer :: gid  ! global id
     integer :: rid  ! rank id
  end type cell
  type(cell),dimension(:),allocatable,target :: dd ! reordered cell index
  integer,dimension(:,:),allocatable :: nbrs !neighbor list for each rank
  integer :: flag ! adjacency flag between neighbors
  !***********
  !* particle
  !***********
  type par
     integer :: x, y, z   ! coordinate
     integer :: r         ! mpi rank
     integer :: e         ! positive or negative energy
     integer :: step      ! step count
     integer :: pid       ! particle id
     logical :: lalive    ! alive
  end type par
  !***********************************************************************
  !*pars:       init total particles
  !*scattpars : reorder particles to be scattered to all ranks
  !*ps:         subset of particles on each rank
  !*sndbuff:    buffer of particles that need transfer(hit rank boundary)
  !*rcvbuff:    buffer of particles that received from neibors
  !*sndidx:     index of the last particle in the send buffer
  !*rcvidx:     index of the last particle in the send buffer
  !***********************************************************************
  type(par),allocatable,target :: pars(:),scattpars(:),ps(:),&
       &sndbuff(:,:),rcvbuff(:,:)
  integer,allocatable,target :: sndidx(:)

  integer,dimension(:),allocatable :: counts,displs
  integer ::ttlPars   ! number of total particles, reduced to master
  integer ::subPars   ! number of particles scattered to each slave
  integer ::ttlEnergy ! total energy, reduced to master
  integer ::subEnergy ! subtotal of energy scattered to each slave
  logical ::existLocalPar ! exist local active particles

  integer particletype, celltype, oldtypes(0:1)   ! required variables
  integer blockcounts(0:1), offsets(0:1), extent

  !************************
  !* random generator
  !************************
  integer :: values(1:8), k
  integer, dimension(:), allocatable :: seed
  !real(8) :: r
  !***********
  !* timing
  !***********
  real*8 :: t0,t1,t2
  !***********
  !* tree
  !***********
  type trnode
     integer :: parent
     integer :: lchild,rchild
  end type trnode
  type(trnode),target::tn

  !********************
  !* mpi setups
  !********************
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  globalFinish=.false.
  pcmplt=0
  it=0
  outreq=0
  !**************************************
  !* mpi structure for passing particles
  !**************************************
  ! setup description of the 7 MPI_INTEGER fields x, y, z, r, e, step, pid
  offsets(0) = 0
  oldtypes(0) = MPI_INTEGER
  blockcounts(0) = 7

  ! setup description of the  MPI_REAL fields n, type
  ! need to first figure offset by getting size of MPI_REAL
  call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
  offsets(1) = 7 * extent
  oldtypes(1) = MPI_LOGICAL !MPI_REAL8
  blockcounts(1) = 1

  ! define structured type and commit it
  call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
       particletype, ierr)
  call MPI_TYPE_COMMIT(particletype, ierr)

  !**************************************************
  !* mpi structure for passing domain decompositions
  !**************************************************
  ! setup description of the 1 MPI_INTEGER fields gid:globalid
  offsets(0) = 0
  oldtypes(0) = MPI_INTEGER
  blockcounts(0) = 1

  ! setup description of the  MPI_INTEGER fields rid:rankid
  call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
  offsets(1) = 1 * extent
  oldtypes(1) = MPI_INTEGER
  blockcounts(1) = 1

  ! define structured type and commit it
  call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
       celltype, ierr)
  call MPI_TYPE_COMMIT(celltype, ierr)
  !**************************
  !* allocate working arrays
  !**************************
  call setArrays
  !*****************************
  !* construct communicate tree
  !*****************************
  tn%parent=getParent(impi)
  tn%lchild=getLChild(impi,nmpi)
  tn%rchild=getRChild(impi,nmpi)
  write(6,'(A1,I4,A8,I4,A6,I4,A6,I4)') '@',impi,&
       &' parent=',tn%parent,' lchd=',tn%lchild,' rchd=',tn%rchild
!
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
!
  t0=mpi_wtime()
  if(lmpi0) then
     !call initParticles(maxpars,pars)
     call initParPerCell(maxpars,pars)
     call domainDecompose(nmpi,dd,dsize,rsize)
     call scatterParticles(nmpi,maxpars,pars,scattpars,&
          &counts,displs)
     !do i=1,nmpi
     !   write(6,*) '0>>',counts(i),displs(i)
     !enddo
  endif

  call bcastDD
  call getNeighbors(nmpi,nbrs) ! everyone get its neighbor lists

  if(impi/=impi0) then
     allocate(scattpars(maxpars))
  endif
  allocate(ps(maxpars))
  !allocate(sndbuff(nmpi,BUFFSIZE))
  allocate(sndbuff(BUFFSIZE,nmpi))
  !allocate(rcvbuff(nmpi,BUFFSIZE))
  allocate(rcvbuff(BUFFSIZE,nmpi))
  allocate(sndidx(nmpi))
  call initRank(impi,ps,pars,sndbuff,sndidx,rcvbuff)
  !call mpi_scatterv(scattpars,counts,displs,particletype,&
  !     &        ps,maxpars,particletype,&
  !     &        impi0,MPI_COMM_WORLD,ierr)
  !do while(.not.globalFinish)
  call recvMaxParBuff(impi,nbrs)
  if(tn%parent>0) then
     call recvChdParComplete
     call recvParentFinish
  endif
  t1=mpi_wtime()
  if(impi==0) then
     write(6,*) 'setup time=',t1-t0
  endif
  do while(.not.globalFinish)
     it=it+1
     !if (impi==0) then
     !   write(6,*) '===> step ',it,globalFinish
     !   call printPars(impi,ps)
     !endif
     existLocalPar=getExistLocalPar(ps)
     if (existLocalPar) then
        call movePar(ps,pcmplt)
     endif
     call recvBuffer(impi,ps,nbrs)
     call tallyPcomplete(tn,pcmplt)
     existLocalPar=getExistLocalPar(ps)
     if (.not.existLocalPar) then
        do i=1,nmpi
           if(sndidx(i)>0) then
              write(6,*) '@no local',impi,sndidx
              call sendBuffer(impi,sndbuff,sndidx,i-1)
           endif
        enddo
        if(impi==0) then
           if(pcmplt==maxpars) then
              globalFinish=.true.
              if(tn%lchild>0) then
                 call mpi_send(globalFinish,1,MPI_LOGICAL,tn%lchild,&
                      &rftag,MPI_COMM_WORLD,ierr)
              endif
              if(tn%rchild>0) then
                 call mpi_send(globalFinish,1,MPI_LOGICAL,tn%rchild,&
                      &rftag,MPI_COMM_WORLD,ierr)
              endif
           endif
        else  ! end master
           call mpi_send(pcmplt,1,MPI_integer,tn%parent,&
                &rpctag,MPI_COMM_WORLD,ierr)
           pcmplt=0
           call mpi_test(reqgf,rgfflag,istat,ierr)
           !write(6,*) '@testgf',rgfflag
           if(rgfflag) then
              !write(6,'(A5,I2,A20)') 'rank(',impi,') recv global finish'
              if(tn%lchild>0) then
                 call mpi_send(globalFinish,1,MPI_LOGICAL,tn%lchild,&
                      &rftag,MPI_COMM_WORLD,ierr)
              endif
              if(tn%rchild>0) then
                 call mpi_send(globalFinish,1,MPI_LOGICAL,tn%rchild,&
                      &rftag,MPI_COMM_WORLD,ierr)
              endif
           endif
        endif ! end slaves
     endif
     if(globalFinish) then
        write(6,*) '@globalfinish',globalFinish
        do i=1,nmpi
           flag=nbrs(impi+1,i)
           if(flag==1)then
              call mpi_cancel(reqrb(i),ierr)
           endif
        enddo
        if(lmpi0) then
           do i = 1,nmpi-1
              call mpi_cancel(reqpc(i),ierr)
           enddo
        endif
     endif
  enddo !tsp_it
  t2=mpi_wtime()
  write(6,'(A17,I5,A4,F6.2)') 'transport time, @',impi,' is ',t2-t1
  subPars=getSubPars(ps)
  subEnergy=getSubEnergy(ps)
  !write(6,'(A5,I2,A9,I5,A11,I2)') '@rank',impi,' subPars=',subPars,' subEnergy=',subEnergy

  call reduceTotalPars ! check total number of praticles on master
  if (impi==impi0) then
     write(6,*) 'after reduce>',ttlPars,ttlEnergy
  endif
  !if(impi==nmpi-1) then
  !call printPars(impi,ps)
  !endif
!
!
!--
!-- FINISH UP:
!=============
! call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output
  !call mpi_barrier(MPI_COMM_WORLD,ierr)
  if(lmpi0) then
     write(6,*)
     write(6,*) 'milagro finished'
  endif

  call MPI_TYPE_FREE(particletype, ierr)
  call MPI_TYPE_FREE(celltype, ierr)

  !****************************
  !* deallocate working arrays
  !****************************
  call mpi_finalize(ierr) !MPI
  call freeArrays


contains
  subroutine initParticles(size,pars)
    implicit none
    integer,intent(in)::size
    type(par),dimension(:),allocatable,intent(out)::pars
    integer :: i
    integer  :: temp,tot ! use tot to check sum of energy is zero
    allocate(pars(size))
    tot=0
    call random_seed()
    write(6,*) 'init pars:   x           y           z        energy'
    do i=1,size
       pars(i)%lalive=.true.
       pars(i)%x=getCoor(dimx)
       !pars(i)%x=1 ! test
       pars(i)%y=getCoor(dimy)
       pars(i)%y=1 ! test
       pars(i)%z=getCoor(dimz)
       if(i<=size/2) then
          temp=getEnergy()
          pars(i)%e=temp
          pars(size-i+1)%e=-1*temp
       endif
       if (mod(size,2)==0) then
       else
          if(i==size/2+1) then
             pars(i)%e=0
          endif
       endif
       pars(i)%pid=i ! test energy as label
       tot=tot+pars(i)%e
       if(i==size) then
          write(6,'(A2,I4,I4,I4,I3,I10)') '->',&
               &pars(i)%x,pars(i)%y,pars(i)%z,pars(i)%e,pars(i)%pid
       endif
    enddo
    write(6,'(A50,I10)') ' ================================ total energy: ',tot
  end subroutine initParticles

  subroutine initParPerCell(size,pars)
    implicit none
    integer,intent(in)::size
    type(par),dimension(:),allocatable,intent(out)::pars
    integer :: i,j,l
    integer  :: temp,tot ! use tot to check sum of energy is zero
    allocate(pars(size))
    tot=0
    call random_seed()
    write(6,*) '  init pars:   x           y           z        energy'
    do i=1,dimy
       do j=1,dimx
          l=(i-1)*dimx+j
          pars(l)%lalive=.true.
          pars(l)%x=j
          pars(l)%y=i ! test
          pars(l)%z=getCoor(dimz)
          if(l<=size/2) then
             temp=getEnergy()
             pars(l)%e=temp
             pars(size-l+1)%e=-1*temp
          endif
          if (mod(size,2)==0) then
          else
             if(l==size/2+1) then
                pars(l)%e=0
             endif
          endif
          pars(l)%step=MAXSTEP
          pars(l)%pid=l ! test energy as label
          tot=tot+pars(l)%e
          if(l==size) then
          write(6,'(A2,12X,I5,7X,I5,5X,I5,10X,I3,7X,I5,I10)') '->',&
               &pars(l)%x,pars(l)%y,pars(l)%z,pars(l)%e,&
               &pars(l)%step,pars(l)%pid
          endif
       enddo
    enddo
    write(6,'(A50,I10)') ' ================================ total energy: ',tot
  end subroutine initParPerCell

  subroutine domainDecompose(mpi_size,dd,v,c)
    implicit none
    integer,intent(in)::mpi_size
    type(cell),dimension(:),intent(out),allocatable :: dd
    integer,intent(out) :: v ! volume, domain size, number of cells
    integer,intent(out) :: c ! cells per rank
    integer :: i!,j,k,idx
    v = dimx*dimy*dimz
    c = v/mpi_size
    allocate(dd(v))
    !**********************
    ! strip decomposition
    !**********************
    do i=1,v
       dd(i)%gid=i
       dd(i)%rid=(i-1)/c
    enddo
  end subroutine domainDecompose

  subroutine scatterParticles(mpi_size,psize,pars,scattpars,&
       & counts,displs)
    implicit none
    integer,intent(in)::mpi_size
    integer,intent(in)::psize ! particle size
    type(par),dimension(:),intent(inout)::pars
    integer,dimension(:),allocatable,intent(out) :: counts,displs
    integer,dimension(:),allocatable :: offset
    type(par),dimension(:),allocatable,intent(out)::scattpars
    integer :: i,rankid,pointer,stage
    allocate(counts(mpi_size),displs(mpi_size),offset(mpi_size))
    allocate(scattpars(psize))
    stage=0
    do i=1,mpi_size
       counts(i)=0
       displs(i)=0
       offset(i)=0
    enddo
    do i=1,psize
       rankid=getRankId(pars(i)%x,pars(i)%y,pars(i)%z,pars(i)%pid)
       pars(i)%r=rankid
       counts(rankid+1)=counts(rankid+1)+1
    enddo
    displs(1)=0
    do i=2,mpi_size
       displs(i)=counts(i-1)+displs(i-1)
    enddo
    do i=1,psize
       rankid=pars(i)%r
       offset(rankid+1)=offset(rankid+1)+1
       pointer=displs(rankid+1)+offset(rankid+1)
       scattpars(pointer)=pars(i)
    enddo
    if(allocated(offset))then
       deallocate(offset)
    endif
  end subroutine scatterParticles

  subroutine printPars(myrank,ps)
    implicit none
    integer,intent(in)::myrank
    type(par),dimension(:),intent(in)::ps
    integer::i
    if(impi==myrank)then
       write(6,*) '>print pars'
       do i=1,maxpars
          !i=maxpars ! only print one particles
          write(6,'(1X,A4,I3,A3,I3,A1,I3,A1,I2,A6,I2,A8,I2,A6,I5,A5,I10,A7,L1)') &
               &'rnk(',impi,') (',ps(i)%x,',',ps(i)%y,',',ps(i)%z,&
               &') rnk=',ps(i)%r,' energy=',ps(i)%e,' step=',ps(i)%step,&
               &' pid=',ps(i)%pid,' alive=',ps(i)%lalive
       enddo
       write(6,*)
    endif
  end subroutine printPars


  subroutine bcastDD
    implicit none
    integer :: v
    v=dimx*dimy*dimz
    if(impi/=impi0) then
       allocate(dd(v))
       allocate(pars(maxpars))
    endif
    call mpi_bcast(dd,v,celltype,impi0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(pars,maxpars,particletype,impi0,MPI_COMM_WORLD,ierr)
  end subroutine bcastDD

  subroutine printDD(myrank,dd)
    implicit none
    integer,intent(in)::myrank
    type(cell),dimension(:),intent(in)::dd
    integer::i,v
    v=dimx*dimy*dimz
    if(impi==myrank)then
       do i=1,v
          write(6,*) 'dd>>',impi,dd(i)%gid,dd(i)%rid
       enddo
       write(6,*)
    endif
  end subroutine printDD

  subroutine initRank(myrank,ps,pars,sndbuff,sndidx,rcvbuff)
    implicit none
    integer, intent(in)::myrank
    type(par),dimension(:),  intent(inout)::ps
    type(par),dimension(:),  intent(in)::pars
    type(par),dimension(:,:),intent(inout)::sndbuff,rcvbuff
    integer,  dimension(:),intent(inout)::sndidx
    integer::i,j
    do i=1,maxpars
       ps(i)%x=0
       ps(i)%y=0
       ps(i)%z=0
       ps(i)%r=-1
       ps(i)%e=0
       ps(i)%step=MAXSTEP
       ps(i)%pid=0
       ps(i)%lalive=.false.
       if(pars(i)%r==myrank) then
          ps(i)=pars(i)
       endif
    enddo
    do i=1,nmpi
       do j=1,BUFFSIZE
          sndbuff(j,i)%x=0
          sndbuff(j,i)%y=0
          sndbuff(j,i)%z=0
          sndbuff(j,i)%r=-1
          sndbuff(j,i)%e=0
          sndbuff(j,i)%step=MAXSTEP
          sndbuff(j,i)%pid=0
          sndbuff(j,i)%lalive=.false.
          rcvbuff(j,i)%x=0
          rcvbuff(j,i)%y=0
          rcvbuff(j,i)%z=0
          rcvbuff(j,i)%r=-1
          rcvbuff(j,i)%e=0
          rcvbuff(j,i)%step=MAXSTEP
          rcvbuff(j,i)%pid=0
          rcvbuff(j,i)%lalive=.false.
       enddo !j
       sndidx(i)=0
    enddo !i
  end subroutine initRank

  subroutine getNeighbors(nmpi,nbrs)
    implicit none
    integer,intent(in)::nmpi
    integer,dimension(:,:),intent(out),allocatable::nbrs
    integer::i,j
    allocate(nbrs(nmpi,nmpi))
    !******************
    !* initialize nbrs
    !******************
    do i=1,nmpi
       do j=1,nmpi
          nbrs(i,j)=0
       enddo
    enddo
    !********************
    !* nbrs for strip-dd
    !********************
    do i=1,nmpi
       if(i>1) nbrs(i,i-1)=1
       if(i<nmpi) nbrs(i,i+1)=1
    enddo
    nbrs(1,nmpi)=1 ! loop
    nbrs(nmpi,1)=1 ! loop
  end subroutine getNeighbors

  subroutine recvMaxParBuff(myrank,nbrs)
    implicit none
    integer,intent(in)::myrank
    integer,dimension(:,:),intent(in)::nbrs
    integer::i,flag,source
    do i=1,nmpi
       source=i-1
       flag=nbrs(myrank+1,i)
       if(flag==1)then
          call mpi_irecv(rcvbuff(1,i),maxbfsize(i),particletype,source,rbftag,&
               & MPI_COMM_WORLD,reqrb(i),ierr)
       endif
    enddo
  end subroutine recvMaxParBuff

  subroutine recvParComplete()
    implicit none
    integer::i,source
    do i=1,nmpi-1
       source=i
       call mpi_irecv(pcomplete(i),1,MPI_INTEGER,source,rpctag,&
            & MPI_COMM_WORLD,reqpc(i),ierr)
    enddo
  end subroutine recvParComplete

  subroutine recvChdParComplete()
    implicit none
    if (tn%lchild>0) then
       call mpi_irecv(pcomplete(1),1,MPI_INTEGER,tn%lchild,rpctag,&
            & MPI_COMM_WORLD,reqpc(1),ierr)
    endif
    if (tn%rchild>0) then
       call mpi_irecv(pcomplete(2),1,MPI_INTEGER,tn%rchild,rpctag,&
            & MPI_COMM_WORLD,reqpc(2),ierr)
    endif
  end subroutine recvChdParComplete

  subroutine recvFinish()
    implicit none
    call mpi_irecv(globalFinish,1,MPI_LOGICAL,0,rftag,&
         & MPI_COMM_WORLD,reqgf,ierr)
  end subroutine recvFinish

  subroutine recvParentFinish()
    implicit none
    call mpi_irecv(globalFinish,1,MPI_LOGICAL,tn%parent,rftag,&
         & MPI_COMM_WORLD,reqgf,ierr)
  end subroutine recvParentFinish

  subroutine reduceTotalPars()
    implicit none
    call mpi_reduce(subPars,ttlPars,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(subEnergy,ttlEnergy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  end subroutine reduceTotalPars

  subroutine setArrays()
    implicit none
    integer::i
    allocate(pcomplete(nmpi),maxbfsize(nmpi))
    allocate(reqrb(nmpi),reqpc(2),cmpind(2),&
         &rbflag(nmpi),rpcflag(nmpi-1))
    !do i=1,nmpi
    pcomplete=0
    maxbfsize=BUFFSIZE ! test should init to 0 and eval to maxpar/nmpi
    reqrb=0
    reqpc=0
    cmpind=0
    rbflag=.false.
    rpcflag=.false.
    !enddo
    allocate(counts(nmpi))
    allocate(displs(nmpi))
    do i=1,nmpi
       counts(i)=0
       displs(i)=0
    enddo
  end subroutine setArrays

  subroutine freeArrays()
    implicit none
    if(allocated(pcomplete)) deallocate(pcomplete)
    if(allocated(reqrb)) deallocate(reqrb)
    if(allocated(reqpc)) deallocate(reqpc)
    if(allocated(cmpind)) deallocate(cmpind)
    if(allocated(rbflag)) deallocate(rbflag)
    if(allocated(rpcflag)) deallocate(rpcflag)
    if(allocated(maxbfsize)) deallocate(maxbfsize)
    if(allocated(counts)) deallocate(counts)
    if(allocated(displs)) deallocate(displs)
    if(allocated(ps)) deallocate(ps)
    if(allocated(sndbuff)) deallocate(sndbuff)
    if(allocated(rcvbuff)) deallocate(rcvbuff)
    if(allocated(sndidx)) deallocate(sndidx)
    if(allocated(dd)) deallocate(dd)
  end subroutine freeArrays

  subroutine movePar(ps,pcmplt)
    implicit none
    type(par), dimension(:),intent(inout)::ps
    integer,intent(out)::pcmplt
    !type(cell),dimension(:),intent(in)   ::dd
    integer :: dir  ! x:0, y:1, z:2
    integer :: step ! -1, 0 , +1
    integer :: oldrank,newrank,stage,idle
    !********************************
    !* test random number generator
    !********************************
    call date_and_time(values=values)
    call random_seed(size=k)
    if(.not.allocated(seed)) allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
    do i=1,maxpars
       if(ps(i)%lalive) then
          do idle=1,10000
          enddo ! pretending computation
          dir=getDirection()
          step=getStep()
          oldrank=ps(i)%r
          !write(6,*) '@m1>',dir,step,ps(i)%x,ps(i)%y,ps(i)%z
          select case(dir)
          case(0) !x direction
             ps(i)%x=getNewCoor(ps(i)%x,step,dimx)
          case(1) !y direction
             ps(i)%y=getNewCoor(ps(i)%y,step,dimy)
          case(2) !z direction
             ps(i)%z=getNewCoor(ps(i)%z,step,dimz)
          case default
             stop 'invalid direction'
          end select
          stage=1
          newrank=getRankId(ps(i)%x,ps(i)%y,ps(i)%z,ps(i)%pid)
          ps(i)%r=newrank
          ps(i)%step=ps(i)%step-1
          if(ps(i)%step==0) then
             pcmplt=pcmplt+1
             ps(i)%lalive=.false.
             !write(6,'(A6,I2,A5,I5,A9)') '@rank(',impi,') pr.',i,' finished'
          endif
          !if(impi==0) then
          !write(6,'(5X,A6,I2,A6,I2,A1,I2,I2,I2,I10)') &
          !     &'move>(',oldrank,') to (',newrank,')',&
          !     &ps(i)%x,ps(i)%y,ps(i)%z,ps(i)%pid
          !endif
          if(newrank/=oldrank) then
             !write(6,'(1X,A6,I2,A6,I2,A6,I5,A5,I5,A1)') 'from (',oldrank,&
             !     &') to (',newrank,') top(',sndidx(newrank+1),') mx(',&
             !     &maxbfsize(newrank+1),')'
             if(sndidx(newrank+1)<maxbfsize(newrank+1)) then
                call writeBuffer(sndbuff,sndidx,newrank,ps(i))
             else
                !call printBuffer(impi,sndbuff,sndidx)
                call sendBuffer(impi,sndbuff,sndidx,newrank)
                call writeBuffer(sndbuff,sndidx,newrank,ps(i))
             endif
             call delParticle(i,ps) ! remove it from local ps
          endif
       endif
    enddo
    !****************************
    !* print leftover snd buffer
    !****************************
    !if (impi==dmpi-1) then
    !   call printBuffer(impi,sndbuff,sndidx)
    !endif
  end subroutine movePar

  subroutine writeBuffer(sndbuff,sndidx,newrank,p)
    implicit none
    type(par),dimension(:,:),intent(inout)::sndbuff
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::newrank
    type(par),intent(in)::p
    integer::j
    sndidx(newrank+1)=sndidx(newrank+1)+1 ! increment number of buffered par
    j=sndidx(newrank+1)
    sndbuff(j,newrank+1)=p
  end subroutine writeBuffer

  subroutine printBuffer(myrank,sndbuff,sndidx)
    implicit none
    integer,intent(in)::myrank
    type(par),dimension(:,:),intent(in)::sndbuff
    integer,dimension(:),intent(in)::sndidx
    integer::i,j,top
    write(6,'(A12,I2,A2)') 'buffer@rank(',myrank,')'
    do i=1,nmpi
       top=sndidx(i)
       do j=1,top
          if(sndbuff(j,i)%lalive) then
             write(6,'(A9,I2,A5,I2,A3,I2,A1,I2,A1,I2,A6,I2,A8,I2,A6,I5,A5,I10)') &
                  &' nbr rank(',i-1,')par(',j,&
                  &')=(',sndbuff(j,i)%x,',',sndbuff(j,i)%y,',',&
                  &sndbuff(j,i)%z,') rnk=',sndbuff(j,i)%r,&
                  &' energy=',sndbuff(j,i)%e,' step=',sndbuff(j,i)%step,&
                  &' pid=',sndbuff(j,i)%pid
          endif
       enddo
    enddo
  end subroutine printBuffer

  subroutine sendBuffer(myrank,sndbuff,sndidx,newrank)
    implicit none
    integer,intent(in)::myrank
    type(par),dimension(:,:),intent(inout)::sndbuff
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::newrank
    integer::count
    count=sndidx(newrank+1) ! increment number of buffered par
    if(myrank<0) then ! don't print
       write(6,'(A12,I2,A2,I5,A8,I2,A2)') 'send@rank(',myrank,') ',&
            &count,'prs to (',newrank,')'
    endif
    call mpi_send(sndbuff(1,newrank+1),BUFFSIZE,particletype,newrank,&
         &sndtag,MPI_COMM_WORLD,ierr)
    call cleanBuffer(sndbuff,sndidx,newrank)
  end subroutine sendBuffer

  subroutine cleanBuffer(sndbuff,sndidx,newrank)
    implicit none
    !integer,intent(in)::myrank
    type(par),dimension(:,:),intent(inout)::sndbuff
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::newrank
    integer :: j
    !write(6,'(A12,I2,A5,I2,A2)') 'clean@rank(',myrank,').bf(',newrank,')'
    do j=1,sndidx(newrank+1)
       sndbuff(j,newrank+1)%x=0
       sndbuff(j,newrank+1)%y=0
       sndbuff(j,newrank+1)%z=0
       sndbuff(j,newrank+1)%r=-1
       sndbuff(j,newrank+1)%e=0
       sndbuff(j,newrank+1)%step=MAXSTEP
       sndbuff(j,newrank+1)%pid=0
       sndbuff(j,newrank+1)%lalive=.false.
    enddo
    sndidx(newrank+1)=0 ! clean buffer top
  end subroutine cleanBuffer

  subroutine delParticle(idx,ps)
    implicit none
    integer,intent(in)::idx
    type(par),dimension(:),intent(inout)::ps
    ps(idx)%x=0
    ps(idx)%y=0
    ps(idx)%z=0
    ps(idx)%e=0
    ps(idx)%r=-1
    ps(idx)%step=MAXSTEP
    ps(idx)%pid=0
    ps(idx)%lalive=.false.
  end subroutine delParticle

  subroutine recvBuffer(myrank,ps,nbrs)
    implicit none
    integer,intent(in)::myrank
    type(par),dimension(:),intent(inout)::ps
    integer,dimension(:,:),intent(in)::nbrs
    integer::i,flag,source,rcvidx,psid
    do i=1,nmpi
       source=i-1
       flag=nbrs(myrank+1,i)
       if(flag==1)then
          call mpi_test(reqrb(i),rbflag(i),istat,ierr)
          if(rbflag(i)) then
!-- add incoming particles to main particle array on this rank
             do rcvidx = 1, maxbfsize(i)
                !write(6,'(I2,I2,I2,I2,I3,I10)') impi,rcvbuff(rcvidx,i)%x,rcvbuff(rcvidx,i)%y,&
                !     &rcvbuff(rcvidx,i)%z,rcvbuff(rcvidx,i)%e,rcvbuff(rcvidx,i)%pid
                if(rcvbuff(rcvidx,i)%x/=0) then
                   psid=rcvbuff(rcvidx,i)%pid
                   ps(psid)=rcvbuff(rcvidx,i)
                end if
                rcvbuff(rcvidx,i)%x=0
                rcvbuff(rcvidx,i)%y=0
                rcvbuff(rcvidx,i)%z=0
                rcvbuff(rcvidx,i)%r=-1
                rcvbuff(rcvidx,i)%e=0
                rcvbuff(rcvidx,i)%step=MAXSTEP
                rcvbuff(rcvidx,i)%pid=0
                rcvbuff(rcvidx,i)%lalive=.false.
             enddo
!-- repost the nonblocking receive
             call mpi_irecv(rcvbuff(1,i),maxbfsize(i),particletype, &
                  source,rbftag,MPI_COMM_WORLD,reqrb(i),ierr)
          endif
       endif
    enddo
  end subroutine recvBuffer

  subroutine tallyPcomplete(tn,pcmplt)
    implicit none
    type(trnode),intent(in)::tn
    integer,intent(inout)::pcmplt
    integer::i,idx,chd
    call mpi_testsome(2,reqpc,outreq,cmpind,istat,ierr)
    if(outreq>0) then
       do i=1,outreq
          idx=cmpind(i)
          pcmplt=pcmplt+pcomplete(i)
          if (idx==1) then
             chd=tn%lchild
          else
             chd=tn%rchild
          endif
          call mpi_irecv(pcomplete(i),1,MPI_INTEGER,chd,rpctag,&
               & MPI_COMM_WORLD,reqpc(i),ierr)
       enddo
    endif
  end subroutine tallyPcomplete

  integer function getCoor(dim)
    implicit none
    integer,intent(in)::dim
    integer :: c
    real*8  :: temp
    call random_number(temp)
    c=int(temp*dim)+1
    getCoor=c
  end function getCoor

  integer function getEnergy()
    implicit none
    real*8   :: temp
    integer  :: e
    call random_number(temp)
    e=INT((temp*2-1)*10)
    getEnergy=e
  end function getEnergy

  integer function getRankId(x,y,z,pid)
    implicit none
    integer,intent(in) ::x,y,z,pid
    integer::grdidx,r
    grdidx=x+(y-1)*dimx+(z-1)*dimx*dimy
    if(grdidx<0) then
       write(6,*) '@getrankid',x,y,z,pid
    endif
    r=dd(grdidx)%rid
    !if(r>nmpi-1) then
    !write(6,'(A9,I4,A1,I4,A1,I4,A7,I8,A4,I10)') &
    !     &'@getrid>(',x,',',y,',',z,' grdid=',grdidx,&
    !     &'pid=',pid
    !endif
    getRankId=r
  end function getRankId

  integer function getSubPars(ps)
    implicit none
    type(par),dimension(:),intent(in)::ps
    integer::i,s
    s=0
    do i=1,maxpars
       if(ps(i)%x>0) then
          s=s+1
       endif
    enddo
    getSubPars=s
  end function getSubPars

  integer function getSubEnergy(ps)
    implicit none
    type(par),dimension(:),intent(in)::ps
    integer::i,e
    e=0
    do i=1,maxpars
       if(ps(i)%lalive) then
          e=e+ps(i)%e
       endif
    enddo
    getSubEnergy=e
  end function getSubEnergy

  logical function getExistLocalPar(ps)
    implicit none
    type(par),dimension(:),intent(in)::ps
    integer::i,s
    logical::l
    s=0
    l=.FALSE.
    do i=1,maxpars
       if(ps(i)%lalive) then
          s=s+1
       endif
    enddo
    if(s/=0) then
       l=.TRUE.
    endif
    getExistLocalPar=l
  end function getExistLocalPar

  integer function getDirection()
    implicit none
    real*8   :: temp
    integer  :: d
    call random_number(temp)
    !d=FLOOR(temp*3)
    !d=FLOOR(temp*2) ! limit to 0:x,1:y for test
    d=1 ! test, move right
    getDirection=d
    !write(6,*) 'dir>',d
  end function getDirection

  integer function getStep()
    implicit none
    real*8   :: temp
    integer  :: s
    call random_number(temp)
    !s=FLOOR(temp*3)-1
    s=1 ! test, move 1
    getStep=s
    !write(6,*) 'step>',s
  end function getStep

  integer function getNewCoor(oldcoor,step,boundary)
    implicit none
    integer,intent(in) :: oldcoor, step, boundary
    integer :: newcoor
    newcoor=oldcoor+step
    if(newcoor<1)then
       newcoor=boundary
    elseif(newcoor>boundary) then
       newcoor=1
    endif
    getNewCoor=newcoor
  end function getNewCoor

  integer function getParent(myrank)
    implicit none
    integer,intent(in)::myrank
    integer::p
    if(myrank==0) then
       p=-1
    else
       if(isOdd(myrank)) then
          p=(myrank-1)/2
       else
          p=myrank/2-1
       endif
    endif
    getParent=p
  end function getParent

  integer function getLChild(myrank,maxrank)
    implicit none
    integer,intent(in)::myrank,maxrank
    integer::lchd
    lchd=myrank*2+1
    if (lchd>=maxrank) then
       lchd=-1
    endif
    getLChild=lchd
  end function getLChild

  integer function getRChild(myrank,maxrank)
    implicit none
    integer,intent(in)::myrank,maxrank
    integer::rchd
    rchd=(myrank+1)*2
    if (rchd>=maxrank) then
       rchd=-1
    endif
    getRChild=rchd
  end function getRChild

  logical function isOdd(myrank)
    implicit none
    integer,intent(in)::myrank
    logical::l
    l=.false.
    if(mod(myrank,2)==0) then
    else
       l=.true.
    endif
    isOdd=l
  end function isOdd

end program milagroimpr
! vim: fdm=marker
