!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
program minimilagro
  implicit none
  INCLUDE 'mpif.h'
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer,parameter :: impi0=0 !the master rank
  integer,parameter :: maxpars=17
  integer,parameter :: dimx=4,dimy=8,dimz=1
  integer,parameter :: rmpbtag=5,rpctag=10,rftag=15,sndtag=20,rbftag=20
  integer,parameter :: dmpi=2 !debug this mpi rank

  integer,dimension(:),allocatable :: maxbfsize  ! max particle buffer size on each ranks
  integer,dimension(:),allocatable :: pcomplete  ! array for number of particle complete on master
  logical :: globalFinish  ! global finish tag on slaves
  integer :: pcmplt        ! number of particle complete on each slave

  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,it,istat(MPI_STATUS_SIZE,4),i
  integer :: tsp_start,tsp_end !time step
  !********************************
  !* reqs:
  !* reqmbs: max buffer size
  !* reqrb:  receive buffer
  !* reqpc:  particle complete
  !* reqgf:  global finish
  !* rbflag: receive buffer flag
  !********************************
  integer,dimension(:),allocatable :: req,reqmbs,reqrb,reqpc
  integer :: reqgf
  logical,dimension(:),allocatable :: rbflag
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
  !***********
  !* particle
  !***********
  type par
     integer :: x, y, z  ! coordinate
     !integer :: dx,dy,dz ! delta coordinate
     integer :: r        ! mpi rank
     integer  :: e        ! positive or negative energy
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
  integer,allocatable,target :: sndidx(:),rbuffidx(:)

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

  !********************
  !* mpi setups
  !********************
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  tsp_start = 1
  tsp_end   = 1
  !**************************************
  !* mpi structure for passing particles
  !**************************************
  ! setup description of the 7 MPI_INTEGER fields x, y, z, r
  offsets(0) = 0
  oldtypes(0) = MPI_INTEGER
  blockcounts(0) = 4

  ! setup description of the  MPI_REAL fields n, type
  ! need to first figure offset by getting size of MPI_REAL
  call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
  offsets(1) = 4 * extent
  oldtypes(1) = MPI_INTEGER !MPI_REAL8
  blockcounts(1) = 1

  ! define structured type and commit it
  call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
       particletype, ierr)
  call MPI_TYPE_COMMIT(particletype, ierr)

  !**************************************************
  !* mpi structure for passing domain decompositions
  !**************************************************
  ! setup description of the 7 MPI_INTEGER fields x, y, z, r
  offsets(0) = 0
  oldtypes(0) = MPI_INTEGER
  blockcounts(0) = 1

  ! setup description of the  MPI_REAL fields n, type
  ! need to first figure offset by getting size of MPI_REAL
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
!
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
!
  if(lmpi0) then
     call initParticles(maxpars,pars)
     call domainDecompose(nmpi,dd,dsize,rsize)
     call scatterParticles(nmpi,maxpars,pars,scattpars,&
          &counts,displs)
     !call printPars(impi0,scattpars)
     do i=1,nmpi
        write(6,*) '0>>',counts(i),displs(i)
     enddo
  endif

  call bcastDD
  !call printDD(1,dd)
  call getNeighbors(nmpi,nbrs) ! everyone get its neighbor lists

  if(impi/=impi0) allocate(scattpars(maxpars))
  allocate(ps(maxpars))
  allocate(sndbuff(nmpi,int(maxpars/nmpi)))
  allocate(rcvbuff(nmpi,int(maxpars/nmpi)))
  allocate(sndidx(nmpi),rbuffidx(nmpi))
  call initRank(ps,sndbuff,sndidx)
  call mpi_scatterv(scattpars,counts,displs,particletype,&
       &        ps,maxpars,particletype,&
       &        impi0,MPI_COMM_WORLD,ierr)

  do it=tsp_start,tsp_end
     !if(impi==dmpi) then
     !   call printPars(impi,ps)
     !endif
     subPars=getSubPars(ps)
     subEnergy=getSubEnergy(ps)
     !write(6,*) 'subPars=',impi,subPars
     call postRecvMaxParBuff(impi,nbrs)

     if(impi==impi0) then
        call postRecvParComplete
     else
        call postRecvFinish
     endif
     call reduceTotalPars ! check total number of praticles on master
     if (impi==impi0) then
        write(6,*) 'after reduce>',ttlPars,ttlEnergy
     endif
     if (impi==impi0) then
        write(6,*) '===> step ',it
     endif
     existLocalPar=getExistLocalPar(ps)
     if (existLocalPar.eqv..TRUE.) then
        call movePar(ps)
     endif
     call mpi_waitall(nmpi,reqrb,istat,ierr)
     call recvBuffer(impi,ps,nbrs)
     if(impi==dmpi) then
        call printPars(impi,ps)
     endif
  enddo !tsp_it
!
!
!--
!-- FINISH UP:
!=============
! call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output
  call mpi_barrier(MPI_COMM_WORLD,ierr)
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
      pars(i)%x=getCoor(dimx)
      pars(i)%y=getCoor(dimy)
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
      pars(i)%e=i ! test energy as label
      tot=tot+pars(i)%e
      write(6,*) '->',pars(i)%x,pars(i)%y,pars(i)%z,pars(i)%e
    enddo
    write(6,'(A50,I1)') ' ================================ total energy: ',tot
  end subroutine initParticles

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
    integer :: i,rankid,pointer
    allocate(counts(mpi_size),displs(mpi_size),offset(mpi_size))
    allocate(scattpars(psize))
    do i=1,mpi_size
       counts(i)=0
       displs(i)=0
       offset(i)=0
    enddo
    do i=1,psize
       rankid=getRankId(pars(i)%x,pars(i)%y,pars(i)%z)
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
          write(6,'(1X,A4,I2,A3,I2,A1,I2,A1,I2,A6,I2,A8,I2)') &
               &'rnk(',impi,') (',ps(i)%x,',',ps(i)%y,',',ps(i)%z,&
               &') rnk=',ps(i)%r,' energy=',ps(i)%e
       enddo
       write(6,*)
    endif
  end subroutine printPars


  subroutine bcastDD
    implicit none
    integer :: v
    v=dimx*dimy*dimz
    if(impi/=impi0) allocate(dd(v))
    call mpi_bcast(dd,v,celltype,impi0,MPI_COMM_WORLD,ierr)
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

  subroutine initRank(ps,sndbuff,sndidx)
    implicit none
    type(par),dimension(:),  intent(inout)::ps
    type(par),dimension(:,:),intent(inout)::sndbuff
    integer,  dimension(:),intent(inout)::sndidx
    integer::i,j,bfsize
    bfsize=int(maxpars/nmpi)
    do i=1,maxpars
       ps(i)%x=0
       ps(i)%y=0
       ps(i)%z=0
       ps(i)%r=-1
       ps(i)%e=0
    enddo
    do i=1,nmpi
       do j=1,bfsize
          sndbuff(i,j)%x=0
          sndbuff(i,j)%y=0
          sndbuff(i,j)%z=0
          sndbuff(i,j)%r=-1
          sndbuff(i,j)%e=0
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
    nbrs(1,2)=1
    nbrs(2,1)=1
    nbrs(2,3)=1
    nbrs(3,2)=1
    nbrs(3,4)=1
    nbrs(4,3)=1
  end subroutine getNeighbors

  subroutine postRecvMaxParBuff(myrank,nbrs)
    implicit none
    integer,intent(in)::myrank
    integer,dimension(:,:),intent(in)::nbrs
    integer::i,flag,source
    do i=1,nmpi
       source=i-1
       flag=nbrs(myrank+1,i)
       if(flag==1)then
          call mpi_irecv(maxbfsize(i),1,MPI_INTEGER,source,rmpbtag,&
            & MPI_COMM_WORLD,reqmbs(i),ierr)
          !write(6,'(A17,I2,A9,I2,A2,I1,A7,I5,A1)') ' recvMaxParBuff@(',&
          !     &myrank,')<-nbrnk(',source,'):',flag,' maxbf(',maxbfsize(i),')'
          call mpi_irecv(rcvbuff(i,1),maxbfsize(i),particletype,source,rbftag,&
            & MPI_COMM_WORLD,reqrb(i),ierr)
       endif
    enddo
  end subroutine postRecvMaxParBuff

  subroutine postRecvParComplete()
    implicit none
    integer::i,source
    do i=1,nmpi-1
       source=i
       call mpi_irecv(pcomplete(i),1,MPI_INTEGER,source,rpctag,&
            & MPI_COMM_WORLD,reqpc(i),ierr)
       !write(6,*) 'par complete size>>',pcomplete(i)
    enddo
  end subroutine postRecvParComplete

  subroutine postRecvFinish()
    implicit none
    !use the last req to receive global finish flag
    !do i=1,nmpi-1
    call mpi_irecv(globalFinish,1,MPI_LOGICAL,0,rftag,&
    !     & MPI_COMM_WORLD,req(nmpi),ierr)
         & MPI_COMM_WORLD,reqgf,ierr)
    !write(6,*) 'global finish>>',globalFinish
    !enddo
  end subroutine postRecvFinish

  subroutine reduceTotalPars()
    implicit none
    call mpi_reduce(subPars,ttlPars,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(subEnergy,ttlEnergy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  end subroutine reduceTotalPars

  subroutine setArrays()
    implicit none
    integer::i
    allocate(pcomplete(nmpi),maxbfsize(nmpi))
    allocate(req(nmpi),reqmbs(nmpi),reqrb(nmpi),reqpc(nmpi),rbflag(nmpi))
    do i=1,nmpi
       pcomplete(i)=0
       maxbfsize(i)=3 ! test should init to 0 and eval to maxpar/nmpi
    enddo
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
    if(allocated(req)) deallocate(req)
    if(allocated(reqmbs)) deallocate(reqmbs)
    if(allocated(reqrb)) deallocate(reqrb)
    if(allocated(reqpc)) deallocate(reqpc)
    if(allocated(rbflag)) deallocate(rbflag)
    if(allocated(maxbfsize)) deallocate(maxbfsize)
    if(allocated(counts)) deallocate(counts)
    if(allocated(displs)) deallocate(displs)
    if(allocated(ps)) deallocate(ps)
    if(allocated(sndbuff)) deallocate(sndbuff)
    if(allocated(rcvbuff)) deallocate(rcvbuff)
    if(allocated(sndidx)) deallocate(sndidx)
    if(allocated(rbuffidx)) deallocate(rbuffidx)
    if(allocated(dd)) deallocate(dd)
  end subroutine freeArrays

  subroutine movePar(ps)
    implicit none
    type(par), dimension(:),intent(inout)::ps
    !type(cell),dimension(:),intent(in)   ::dd
    integer :: dir  ! x:0, y:1, z:2
    integer :: step ! -1, 0 , +1
    integer :: oldrank,newrank
    !********************************
    !* test random number generator
    !********************************
     call date_and_time(values=values)
     call random_seed(size=k)
     if(.not.allocated(seed)) allocate(seed(1:k))
     seed(:) = values(8)
     call random_seed(put=seed)
     do i=1,maxpars
        if(ps(i)%r/=-1) then
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
           newrank=getRankId(ps(i)%x,ps(i)%y,ps(i)%z)
           ps(i)%r=newrank
           !write(6,'(5X,A6,I2,A6,I2,A1)') 'move>(',oldrank,') to (',newrank,')'
           if(newrank/=oldrank) then
              !write(6,'(1X,A6,I2,A6,I2,A6,I5,A5,I5,A1)') 'from (',oldrank,&
              !     &') to (',newrank,') top(',sndidx(newrank+1),') mx(',&
              !     &maxbfsize(newrank+1),')'
              if(sndidx(newrank+1)<maxbfsize(newrank+1)) then
                 call writeBuffer(sndbuff,sndidx,newrank,ps(i))
              else
                 call sendBuffer(impi,sndbuff,sndidx,newrank)
                 call cleanBuffer(sndbuff,sndidx,newrank)
                 call writeBuffer(sndbuff,sndidx,newrank,ps(i))
              endif
              call removeParticle(i,ps) ! remove it from local ps
           endif
           pcmplt=pcmplt+1
        endif
     enddo
     !****************************
     !* print leftover snd buffer
     !****************************
     if (impi==dmpi-1) then
        call printBuffer(impi,sndbuff,sndidx)
     endif
  end subroutine movePar

  subroutine writeBuffer(sndbuff,sndidx,newrank,p)
    implicit none
    type(par),dimension(:,:),intent(inout)::sndbuff
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::newrank
    type(par),intent(in)::p
    integer::i
    sndidx(newrank+1)=sndidx(newrank+1)+1 ! increment number of buffered par
    i=sndidx(newrank+1)
    sndbuff(newrank+1,i)=p
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
          if(sndbuff(i,j)%r/=-1) then
             write(6,'(A9,I2,A5,I2,A3,I2,A1,I2,A1,I2,A6,I2,A8,I2)') &
                  &' nbr rank(',i-1,')par(',j,&
                  &')=(',sndbuff(i,j)%x,',',sndbuff(i,j)%y,',',&
                  &sndbuff(i,j)%z,') rnk=',sndbuff(i,j)%r,&
                  &' energy=',sndbuff(i,j)%e
          endif
       enddo
    enddo
  end subroutine printBuffer

  subroutine sendBuffer(myrank,sndbuff,sndidx,newrank)
    implicit none
    integer,intent(in)::myrank
    type(par),dimension(:,:),intent(in)::sndbuff
    integer,dimension(:),intent(in)::sndidx
    integer,intent(in)::newrank
    integer::count
    count=sndidx(newrank+1) ! increment number of buffered par
    write(6,'(A12,I2,A2,I5,A5,I2,A2)') 'send@rank(',myrank,') ',&
         &count,' to (',newrank,')'
    call mpi_send(sndbuff(newrank+1,1),count,particletype,newrank,&
         &sndtag,MPI_COMM_WORLD,ierr)
  end subroutine sendBuffer

  subroutine cleanBuffer(sndbuff,sndidx,newrank)
    implicit none
    !integer,intent(in)::myrank
    type(par),dimension(:,:),intent(inout)::sndbuff
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::newrank
    integer :: i
    !write(6,'(A12,I2,A5,I2,A2)') 'clean@rank(',myrank,').bf(',newrank,')'
    do i=1,sndidx(newrank+1)
       sndbuff(newrank+1,i)%x=0
       sndbuff(newrank+1,i)%y=0
       sndbuff(newrank+1,i)%z=0
       sndbuff(newrank+1,i)%r=-1
       sndbuff(newrank+1,i)%e=0
    enddo
    sndidx(newrank+1)=0 ! clean buffer top
  end subroutine cleanBuffer

  subroutine removeParticle(idx,ps)
    implicit none
    integer,intent(in)::idx
    type(par),dimension(:),intent(inout)::ps
    ps(idx)%x=0
    ps(idx)%y=0
    ps(idx)%z=0
    ps(idx)%r=-1
    ps(idx)%e=0
  end subroutine removeParticle

  subroutine recvBuffer(myrank,ps,nbrs)
    implicit none
    integer,intent(in)::myrank
    type(par),dimension(:),intent(inout)::ps
    integer,dimension(:,:),intent(in)::nbrs
    integer::i,flag,source,rcvidx,pstop
    do i=1,nmpi
       source=i-1
       flag=nbrs(myrank+1,i)
       if(flag==1)then
          call mpi_test(reqrb(i),rbflag(i),istat,ierr)
          if(rbflag(i)) then
!-- add incoming particles to main particle array on this rank
             pstop=getSubPars(ps)
             do rcvidx = 1, maxbfsize(i)
                if(rcvbuff(i,rcvidx)%r/=-1) then
                   ps(pstop+rcvidx)=rcvbuff(i,rcvidx)
                end if
             enddo
!-- repost the nonblocking receive
             call mpi_irecv(rcvbuff(i,1),maxbfsize(i),particletype, &
                  source,rbftag,MPI_COMM_WORLD,reqrb(i),ierr)
          endif
       endif
    enddo
  end subroutine recvBuffer

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

  integer function getRankId(x,y,z)
    implicit none
    integer,intent(in) ::x,y,z
    integer::paridx,r
    paridx=x+(y-1)*dimx+(z-1)*dimx*dimy
    r=dd(paridx)%rid
    getRankId=r
  end function getRankId

  integer function getSubPars(ps)
    implicit none
    type(par),dimension(:),intent(in)::ps
    integer::i,s
    s=0
    do i=1,maxpars
       if(ps(i)%r/=-1) then
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
       if(ps(i)%r/=-1) then
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
       if(ps(i)%r/=-1) then
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
       newcoor=1
    elseif(newcoor>boundary) then
       newcoor=boundary
    endif
    getNewCoor=newcoor
  end function getNewCoor

end program minimilagro
! vim: fdm=marker
