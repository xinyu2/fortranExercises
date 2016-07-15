!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
program minimilagro
  implicit none
  INCLUDE 'mpif.h'
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer,parameter :: impi0=0 !the master rank
  integer,parameter :: maxpars=17,totalValue=0
  integer,parameter :: dimx=4,dimy=8,dimz=1
  integer,parameter :: rmpbtag=5,rpctag=10,rftag=15

  integer :: maxbsize   ! max particle buffer size
  integer :: parcmplt   ! number of particle complete on slaves
  logical :: glbfinsh   ! global finish tag

  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,it,i
  integer :: tsp_start,tsp_end !time step
  !*********
  !* domain
  !*********
  integer :: dsize !domain size(total cells=dimx*dimy*dimz)
  integer :: rsize !rank size(cells per rank=dimx*dimy*dimz/nmpi)
  integer :: parsub !subset of pars, particle per rank
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
     integer :: dx,dy,dz ! delta coordinate
     integer :: r        ! mpi rank
     integer  :: e        ! positive or negative energy
  end type par
  type(par),allocatable,target :: pars(:),scattpars(:),ps(:)
  integer,dimension(:),allocatable :: counts,displs

  integer particletype, oldtypes(0:1)   ! required variables
  integer blockcounts(0:1), offsets(0:1), extent

!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  tsp_start = 1
  tsp_end   = 1

  ! setup description of the 7 MPI_INTEGER fields x, y, z, dx, dy, dz, r
  offsets(0) = 0
  oldtypes(0) = MPI_INTEGER
  blockcounts(0) = 7

  ! setup description of the  MPI_REAL fields n, type
  ! need to first figure offset by getting size of MPI_REAL
  call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr)
  offsets(1) = 7 * extent
  oldtypes(1) = MPI_INTEGER !MPI_REAL8
  blockcounts(1) = 1

  ! define structured type and commit it
  call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
       particletype, ierr)
  call MPI_TYPE_COMMIT(particletype, ierr)

!
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
!
  if(lmpi0) then
     call initParticles(nmpi,maxpars,pars)
     call domainDecompose(nmpi,dd,dsize,rsize)
     !call getNeighbors(nmpi,nbrs)
     call scatterParticles(nmpi,maxpars,pars,dd,scattpars,&
          &counts,displs)
     !call printPars(impi0,scattpars)
     !do i=1,nmpi
     !   write(6,*) '>>',counts(i),displs(i)
     !enddo
  endif

  call getNeighbors(nmpi,nbrs) ! everyone get its neighbor lists

  if(impi/=impi0) allocate(scattpars(maxpars))

  allocate(ps(maxpars))
  call initPs(ps)

  call mpi_scatterv(scattpars,counts,displs,particletype,&
       &        ps,maxpars,particletype,&
       &        impi0,MPI_COMM_WORLD,ierr)

  !call printPars(impi,ps)
  do it=tsp_start,tsp_end
     write(6,*) '===> step ',it

     call postRecvMaxParBuff(impi,nbrs)
     if(impi==impi0) then
        call postRecvParComplt
     else
        call postRecvFinish
     endif


     if(impi==impi0)then
        call scatterParticles(nmpi,maxpars,pars,dd,scattpars,&
             &counts,displs)
     endif
     call mpi_scatterv(scattpars,counts,displs,particletype,&
       &        ps,maxpars,particletype,&
       &        impi0,MPI_COMM_WORLD,ierr)
     !call printPars(impi0,ps)
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
  call mpi_finalize(ierr) !MPI

contains
  subroutine initParticles(mpi_size,size,pars)
    implicit none
    integer,intent(in)::mpi_size,size
    type(par),dimension(:),allocatable,intent(out)::pars
    integer :: i
    integer  :: temp,tot ! use tot to check sum of energy is zero
    allocate(pars(size))
    tot=0
    call random_seed()
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
      tot=tot+pars(i)%e
      write(6,*) '->',pars(i)%x,pars(i)%y,pars(i)%z,pars(i)%e
    enddo
    write(6,*) '>>',tot
  end subroutine initParticles

  subroutine domainDecompose(mpi_size,dd,v,c)
    implicit none
    integer,intent(in)::mpi_size
    type(cell),dimension(:),intent(out),allocatable :: dd
    integer,intent(out) :: v ! volume, domain size, number of cells
    integer,intent(out) :: c ! cells per rank
    integer :: i,j,k,idx
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

  subroutine scatterParticles(mpi_size,psize,pars,dd,scattpars,&
       & counts,displs)
    implicit none
    integer,intent(in)::mpi_size
    integer,intent(in)::psize ! particle size
    type(par),dimension(:),intent(inout)::pars
    type(cell),dimension(:),intent(in) :: dd
    integer,dimension(:),allocatable,intent(out) :: counts,displs
    integer,dimension(:),allocatable :: offset
    type(par),dimension(:),allocatable,intent(out)::scattpars
    integer :: i,paridx,rankid,pointer
    allocate(counts(mpi_size),displs(mpi_size),offset(mpi_size))
    allocate(scattpars(psize))
    do i=1,mpi_size
       counts(i)=0
       displs(i)=0
       offset(i)=0
    enddo
    do i=1,psize
       paridx=pars(i)%x+(pars(i)%y-1)*dimx
       pars(i)%r=dd(paridx)%rid
       counts(dd(paridx)%rid+1)=counts(dd(paridx)%rid+1)+1
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
       do i=1,maxpars
          write(6,*) '>>>',impi,'x',ps(i)%x,'y',ps(i)%y,'e',ps(i)%e
       enddo
       write(6,*)
    endif
  end subroutine printPars

  subroutine initPs(ps)
    implicit none
    type(par),dimension(:),intent(inout)::ps
    integer::i
    do i=1,maxpars
       ps(i)%x=0
       ps(i)%y=0
       ps(i)%z=0
       ps(i)%dx=0
       ps(i)%dy=0
       ps(i)%dz=0
       ps(i)%r=-1
       ps(i)%e=0
    enddo
  end subroutine initPs

  subroutine getNeighbors(nmpi,nbrs)
    implicit none
    integer,intent(in)::nmpi
    integer,dimension(:,:),intent(out),allocatable::nbrs
    integer::i,j
    allocate(nbrs(nmpi,nmpi-1))
    !******************
    !* initialize nbrs
    !******************
    do i=1,nmpi
       do j=1,nmpi-1
          nbrs(i,j)=-1
       enddo
    enddo
    !********************
    !* nbrs for strip-dd
    !********************
    nbrs(1,1)=1
    nbrs(2,1)=0
    nbrs(2,2)=2
    nbrs(3,1)=1
    nbrs(3,2)=3
    nbrs(4,1)=2
  end subroutine getNeighbors

  subroutine postRecvMaxParBuff(myrank,nbrs)
    implicit none
    integer,intent(in)::myrank
    integer,dimension(:,:),intent(in)::nbrs
    integer::i,source

    do i=1,nmpi-1
       source=nbrs(myrank+1,i)
       if(source>-1)then
          call mpi_irecv(maxbsize,1,MPI_INTEGER,source,rmpbtag,&
            & MPI_COMM_WORLD,ierr)
          !write(*,*) 'prmpb>>',myrank,source
       endif

    enddo
  end subroutine postRecvMaxParBuff

  subroutine postRecvParComplt()
    implicit none
    integer::i,source
    do i=1,nmpi-1
       source=i
       call mpi_irecv(parcmplt,1,MPI_INTEGER,source,rpctag,&
            & MPI_COMM_WORLD,ierr)
       !write(6,*) 'prpc>>',source
    enddo
  end subroutine postRecvParComplt

  subroutine postRecvFinish()
    implicit none
    !call mpi_irecv(rf,1,MPI_LOGICAL,0,rftag,&
    !     & MPI_COMM_WORLD,ierr)
    !write(6,*) 'prf>>',0
  end subroutine postRecvFinish

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
    e=(temp*2-1)*10
    getEnergy=e
  end function getEnergy

end program minimilagro
! vim: fdm=marker
