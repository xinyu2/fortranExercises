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
  integer,parameter :: rmpbtag=5,rpctag=10,rftag=15

  integer,dimension(:),allocatable :: maxbfsize  ! max particle buffer size on each ranks
  integer,dimension(:),allocatable :: pcomplete  ! array for number of particle complete on master
  logical :: globalFinish  ! global finish tag on slaves
  integer :: pcmplt        ! number of particle complete on each slave

  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,it,i
  integer,dimension(:),allocatable :: req
  integer :: tsp_start,tsp_end !time step
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
  !*pbuff:      subset of particles that need transfer(hit rank boundary)
  !***********************************************************************
  type(par),allocatable,target :: pars(:),scattpars(:),ps(:),pbuff(:)

  integer,dimension(:),allocatable :: counts,displs
  integer ::ttlPars   ! number of total particles, reduced to master
  integer ::subPars   ! number of particles scattered to each slave
  integer ::ttlEnergy ! total energy, reduced to master
  integer ::subEnergy ! subtotal of energy scattered to each slave
  logical ::existLocalPar ! exist local active particles

  integer particletype, oldtypes(0:1)   ! required variables
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
     call scatterParticles(nmpi,maxpars,pars,dd,scattpars,&
          &counts,displs)
     !call printPars(impi0,scattpars)
     do i=1,nmpi
        write(6,*) '0>>',counts(i),displs(i)
     enddo

  endif

  call getNeighbors(nmpi,nbrs) ! everyone get its neighbor lists

  if(impi/=impi0) allocate(scattpars(maxpars))
  allocate(ps(maxpars))
  allocate(pbuff(maxpars))
  call initPs(ps)
  call mpi_scatterv(scattpars,counts,displs,particletype,&
       &        ps,maxpars,particletype,&
       &        impi0,MPI_COMM_WORLD,ierr)

  !call printPars(impi,ps)
  do it=tsp_start,tsp_end
     write(6,*) '===> step ',it
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
     existLocalPar=getExistLocalPar(ps)
     if (existLocalPar.eqv..TRUE.) then
        call movePar(ps)
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

  !****************************
  !* deallocate working arrays
  !****************************
  call freeArrays

  call mpi_finalize(ierr) !MPI



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
       !paridx=pars(i)%x+(pars(i)%y-1)*dimx
       !pars(i)%r=dd(paridx)%rid
       pars(i)%r=getRankId(pars(i)x,pars(i)%y,pars(i)%z)
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
          call mpi_irecv(maxbfsize(i),1,MPI_INTEGER,source,rmpbtag,&
            & MPI_COMM_WORLD,req(i),ierr)
          !write(*,*) 'RecvMaxParBuff>>',myrank,source
       endif
       !write(*,*) '              >>',myrank,i,maxbfsize(i)
    enddo
  end subroutine postRecvMaxParBuff

  subroutine postRecvParComplete()
    implicit none
    integer::i,source
    do i=1,nmpi-1
       source=i
       call mpi_irecv(pcomplete(i),1,MPI_INTEGER,source,rpctag,&
            & MPI_COMM_WORLD,req(i),ierr)
       !write(6,*) 'par complete size>>',pcomplete(i)
    enddo
  end subroutine postRecvParComplete

  subroutine postRecvFinish()
    implicit none
    integer::i
    do i=1,nmpi-1
       call mpi_irecv(globalFinish,1,MPI_LOGICAL,0,rftag,&
            & MPI_COMM_WORLD,req(i),ierr)
       !write(6,*) 'global finish>>',globalFinish
    enddo
  end subroutine postRecvFinish

  subroutine reduceTotalPars()
    implicit none
    call mpi_reduce(subPars,ttlPars,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call mpi_reduce(subEnergy,ttlEnergy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  end subroutine reduceTotalPars

  subroutine setArrays()
    implicit none
    integer::i
    allocate(pcomplete(nmpi-1))
    allocate(maxbfsize(nmpi-1))
    allocate(req(nmpi-1))
    do i=1,nmpi-1
       pcomplete(i)=0
       maxbfsize(i)=0
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
    if(allocated(maxbfsize)) deallocate(maxbfsize)
    if(allocated(counts)) deallocate(counts)
    if(allocated(displs)) deallocate(displs)
    if(allocated(ps)) deallocate(ps)
    if(allocated(pbuff)) deallocate(pbuff)
  end subroutine freeArrays

  subroutine movePar(ps)
    implicit none
    type(par),dimension(:),intent(inout)::ps
    integer :: dir  ! x:0, y:1, z:2
    integer :: step ! -1, 0 , +1
    integer :: oldrank,newrank
    !********************************
    !* test random number generator
    !********************************
     call date_and_time(values=values)
     call random_seed(size=k)
     allocate(seed(1:k))
     seed(:) = values(8)
     call random_seed(put=seed)
     do i=1,maxpars
        if(ps(i)%r/=-1) then
           dir=getDirection()
           step=getStep()
           oldrank=ps(i)%r
           select case(dir)
           case(0)
              ps(i)%x=ps(i)%x+step
           case(1)
              ps(i)%y=ps(i)%y+step
           case(2)
              ps(i)%z=ps(i)%z+step
           case default
              stop 'invalid direction'
           end select
           newrank=getRankId(ps(i)%x,ps(i)%y,ps(i)%z)



        pcmplt=pcmplt+1
        endif
     enddo

  end subroutine movePar

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
    d=FLOOR(temp*3)
    getDirection=d
    !write(6,*) 'dir>',d
  end function getDirection

  integer function getStep()
    implicit none
    real*8   :: temp
    integer  :: s
    call random_number(temp)
    s=FLOOR(temp*3)-1
    getStep=s
    !write(6,*) 'step>',s
  end function getStep

end program minimilagro
! vim: fdm=marker
