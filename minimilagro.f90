!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
program minimilagro
  implicit none
  INCLUDE 'mpif.h'
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer,parameter :: impi0=0 !the master rank
  integer,parameter :: maxpars=640,totalValue=0
  integer,parameter :: dimx=4,dimy=8,dimz=1
  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,it
  integer :: tsp_start,tsp_end
  type par
     integer :: x, y, z  ! coordinate
     integer :: dx,dy,dz ! delta coordinate
     real*8  :: energy   ! positive or negative
     logical :: cmplt    ! flag for complete
  end type par
  type(par),allocatable,target :: pars(:)
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  tsp_start = 1
  tsp_end   = 40
!
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
!
  if(lmpi0) then
     call initParticles(nmpi,maxpars,pars)
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
  subroutine initParticles(mpi_size,par_size,out_pars)
    implicit none
    integer,intent(in)::mpi_size,par_size
    type(par),dimension(:),allocatable,intent(out)::out_pars
    integer :: i
    real*8  :: temp,tot ! use tot to check sum of energy is zero
    allocate(out_pars(par_size))
    tot=0
    call random_seed()
    do i=1,par_size
      call random_number(temp)
      out_pars(i)%x=int(temp*dimx)+1
      call random_number(temp)
      out_pars(i)%y=int(temp*dimy)+1
      call random_number(temp)
      out_pars(i)%z=int(temp*dimz)+1
      if (i<par_size) then
         call random_number(temp)
         out_pars(i)%energy=temp*2-1
         tot=tot+out_pars(i)%energy
      else
         out_pars(i)%energy=0-tot
         tot=tot+out_pars(i)%energy
         write(6,*) '>',i,'energy=',out_pars(i)%energy
      endif
    enddo
    write(6,*) '>>',tot
  end subroutine initParticles

end program minimilagro
! vim: fdm=marker
