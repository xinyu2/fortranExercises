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
  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,it
  integer :: tsp_start,tsp_end !time step
  integer :: dsize !domain size(total cells=dimx*dimy*dimz)
  integer :: rsize !rank size(cells per rank=dimx*dimy*dimz/nmpi)
  type par
     integer :: x, y, z  ! coordinate
     integer :: dx,dy,dz ! delta coordinate
     integer :: r        ! mpi rank
     real*8  :: e        ! positive or negative energy
  end type par
  type(par),allocatable,target :: pars(:)
  type cell
     integer :: gid  ! global id
     integer :: rid  ! rank id
  end type cell
  type(cell),dimension(:),allocatable,target :: dd ! reordered cell index
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  tsp_start = 1
  tsp_end   = 4
!
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
!
  if(lmpi0) then
     call initParticles(nmpi,maxpars,pars)
     call domainDecompose(nmpi,dd,dsize,rsize)
     call scatterParticles(nmpi,maxpars,pars,dd)
     do it=tsp_start,tsp_end

     enddo !tsp_it
  endif
!
!
!--
!-- FINISH UP:
!=============
! call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output
  if(lmpi0) then
     write(6,*)
     write(6,*) 'milagro finished'
  endif

  call mpi_finalize(ierr) !MPI

contains
  subroutine initParticles(mpi_size,size,pars)
    implicit none
    integer,intent(in)::mpi_size,size
    type(par),dimension(:),allocatable,intent(out)::pars
    integer :: i
    real*8  :: temp,tot ! use tot to check sum of energy is zero
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

  subroutine scatterParticles(mpi_size,psize,pars,dd)
    implicit none
    integer,intent(in)::mpi_size
    integer,intent(in)::psize ! particle size
    type(par),dimension(:),intent(inout)::pars
    type(cell),dimension(:),intent(in) :: dd
    integer,dimension(:),allocatable :: counts,displ,offset
    type(par),dimension(:),allocatable::scatter
    integer :: i,paridx,rankid,pointer
    allocate(counts(mpi_size),displ(mpi_size),offset(mpi_size))
    allocate(scatter(psize))
    do i=1,mpi_size
       counts(i)=0
       displ(i)=0
       offset(i)=0
    enddo
    do i=1,psize
       paridx=pars(i)%x+(pars(i)%y-1)*dimx
       pars(i)%r=dd(paridx)%rid
       counts(dd(paridx)%rid+1)=counts(dd(paridx)%rid+1)+1
    enddo
    displ(1)=1
    do i=2,mpi_size
       displ(i)=counts(i-1)+displ(i-1)
    enddo
    do i=1,psize
       rankid=pars(i)%r
       pointer=displ(rankid+1)+offset(rankid+1)
       offset(rankid+1)=offset(rankid+1)+1
       scatter(pointer)=pars(i)
    enddo
    do i=1,psize
       write(6,*) '>',i,scatter(i)%r
    enddo

  end subroutine scatterParticles

  integer function getCoor(dim)
    implicit none
    integer,intent(in)::dim
    integer :: c
    real*8  :: temp
    call random_number(temp)
    c=int(temp*dim)+1
    getCoor=c
  end function getCoor

  real*8 function getEnergy()
    implicit none
    real*8  :: temp, e
    call random_number(temp)
    e=temp*2-1
    getEnergy=e
  end function getEnergy

end program minimilagro
! vim: fdm=marker
