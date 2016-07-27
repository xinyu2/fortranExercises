program milagro_strip

  implicit none
  INCLUDE 'mpif.h'
!--------------------
!- A milagro mini-app
!--------------------
!-- mpi variables
  integer,parameter :: impi0=0 !the master rank
  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr, istat(MPI_STATUS_SIZE)
  integer :: irequestall,irequest(2)!3
  logical :: iflag(3)
  integer,allocatable :: ireqslav(:)
  logical,allocatable :: iflagslav(:)
!-- non-mpi variables
  integer :: ixoffset,iyoffset
  integer :: iimpi
  integer :: i,j, ii,jj
  integer :: ipart, ibuff
  integer :: src_nspr ! number of source particles per MPI rank
  integer :: itotiter
!-- grid
  integer,parameter :: grd_nx = 12
  integer,parameter :: grd_ny = 10
  integer :: grd_xarr(grd_nx)
  integer :: grd_yarr(grd_ny)
  integer :: grd_icell(grd_nx,grd_ny)
!-- decomposition
  integer :: gas_nx,gas_ny
  integer :: gas_icell(grd_nx/4,grd_ny)
  integer :: gas_neighbors(2)
  integer :: gas_ighostx(grd_ny)
  integer :: gas_ighosty(grd_nx/4)
!-- particles
  integer,parameter :: prt_nmax = 1 !max array size per MPI rank
  integer,parameter :: src_ns = 1   !total number of source particles
  type packet
     logical :: lalive
     integer :: ix, iy
  end type packet
  type(packet) :: prt_particles(prt_nmax) ! particle array (per rank)
!-- particle buffers
  integer,parameter :: prt_nbuf = 2
  integer*8 :: nbuffsize
  integer :: nbuf
  type(packet) :: prt_rcvbuffer(prt_nbuf,2)
  type(packet) :: prt_sndbuffer(prt_nbuf,2)
!-- number of particles completed per rank (stored on master)
  logical :: prt_ldone, prt_lalldone
  integer :: prt_ndeadall,prt_ndead
  integer :: prt_ndonex,prt_ndoney,idonex,idoney
  integer, allocatable :: prt_ncompleted(:)

!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
!-- sanity check
  if(nmpi /= 4) stop 'nmpi /= 4'

!
!-- INITIALIZE:

!-- grid
  forall(i = 1:grd_nx) grd_xarr(i)=i
  forall(j = 1:grd_ny) grd_yarr(j)=j
!-- global grid mapping
  do j = 1, grd_ny
     do i = 1, grd_nx
        grd_icell(i,j) = i+grd_nx*(j-1)
     enddo
  enddo

!-- decomposition
!-----------------------!
!- r -!- r -!- r -!- r -!
!- a -!- a -!- a -!- a -!
!- n -!- n -!- n -!- n -!
!- k -!- k -!- k -!- k -!
!- 0 -!- 1 -!- 2 -!- 3 -!
!-----------------------!
  gas_ny = grd_ny
  gas_nx = grd_nx/4
  ixoffset = gas_nx*impi
  iyoffset = gas_ny*0
  do jj = 1, gas_ny
     do ii = 1, gas_nx
!-- convert local 2d index to global
        i = ii+ixoffset
        j = jj+iyoffset
        gas_icell(ii,jj) = grd_icell(i,j)
     enddo
  enddo
!-- sanity checks
  if(lmpi0 .and. any(gas_icell/=grd_icell(:grd_nx/4,:grd_ny))) &
       stop 'rank 0 does not have correct sub-domain'
  if(impi==1 .and. any(gas_icell/=grd_icell(grd_nx/4+1:grd_nx/2,:grd_ny))) &
       stop 'rank 1 does not have correct sub-domain'
  if(impi==2 .and. any(gas_icell/=grd_icell(grd_nx/2+1:grd_nx*3/4,:grd_ny))) &
       stop 'rank 2 does not have correct sub-domain'
  if(impi==3 .and. any(gas_icell/=grd_icell(grd_nx*3/4+1:grd_nx,:grd_ny))) &
       stop 'rank 3 does not have correct sub-domain'
!-- calculate the neighbors of this rank
!-- (1st index in x-right, 2nd index in x-left) loop 0-1-2-3-0
  if(impi==0) gas_neighbors = [1,3]
  if(impi==1) gas_neighbors = [2,0]
  if(impi==2) gas_neighbors = [3,1]
  if(impi==3) gas_neighbors = [0,2]
!-- get global indices of ghost cells from neighbors
!-- x-boundary
  ixoffset = gas_nx*(impi+1)
  iyoffset = gas_ny*0
  if(impi<nmpi-1) ixoffset=ixoffset+1
  do jj=1,gas_ny
     i = ixoffset
     j = jj+iyoffset
     gas_ighostx(jj) = grd_icell(i,j)
  enddo
!-- y-boundary
  ixoffset = gas_nx*(impi)
  iyoffset = gas_ny
  if(impi<impi) iyoffset=iyoffset+1 ! no y ghostcells
  do ii=1,gas_nx
     i = ii+ixoffset
     j = iyoffset
     gas_ighosty(ii) = grd_icell(i,j)
  enddo

!-- distribute particle numbers
  src_nspr = 0
!-- inactive particles have negative index
  prt_particles%lalive = .false.
  prt_particles%ix = -1
  prt_particles%iy = -1
  do ipart = 1, src_ns
!-- leftmost cells (using global index)
     i = 1
     j = ipart
!-- check serialized global index against sub-domain map
     if(any(gas_icell==grd_icell(i,j))) then
!-- increment set size
        src_nspr=src_nspr + 1
!-- instantiate a particle
        prt_particles(src_nspr)%lalive = .true.
        prt_particles(src_nspr)%ix=i
        prt_particles(src_nspr)%iy=j
     endif
  enddo
!-- print particle set sizes
  write(0,*) 'impi =', impi,': living particles =',count(prt_particles%lalive)

!--==== MILAGRO ===============================================================

!-- master rank will send the completion logical to all slave ranks if true
  prt_lalldone = .false.

!-- initialize buffers
  prt_sndbuffer%lalive = .false.
  prt_sndbuffer%ix = -1
  prt_sndbuffer%iy = -1
  prt_rcvbuffer%lalive = .false.
  prt_rcvbuffer%ix = -1
  prt_rcvbuffer%iy = -1

!-- post nonblocking receive for each neighbor
  nbuffsize = sizeof(prt_rcvbuffer(:,1))
  call mpi_irecv(prt_rcvbuffer(:,1),nbuffsize, &
       MPI_BYTE,gas_neighbors(1),1,MPI_COMM_WORLD,irequest(1),ierr)
  call mpi_irecv(prt_rcvbuffer(:,2),nbuffsize, &
       MPI_BYTE,gas_neighbors(2),2,MPI_COMM_WORLD,irequest(2),ierr)

!-- set-up completion condition
  if(lmpi0) then
     allocate(prt_ncompleted(nmpi-1))
     allocate(ireqslav(nmpi-1))
     allocate(iflagslav(nmpi-1))
!-- post nonblocking receives for number of completed particles
     do iimpi = 1,nmpi-1
        call mpi_irecv(prt_ncompleted(iimpi),1,MPI_INTEGER, &
             iimpi,3,MPI_COMM_WORLD,ireqslav(iimpi),ierr)
     enddo
  else
     call mpi_irecv(prt_lalldone,1,MPI_LOGICAL,impi0,4, &
          !MPI_COMM_WORLD,irequest(3),ierr)
          MPI_COMM_WORLD,irequestall,ierr)
  endif

!-- TRANSPORT:

  prt_ndeadall = 0
  itotiter = 0
  do while (.not. prt_lalldone .and. itotiter<10)
     if (lmpi0) then
        write(6,'(A12,I9,I9)') '> iter=',itotiter,ipart
     endif
     prt_ndonex = 0
     prt_ndoney = 0
     prt_ndead = 0

!-- buffer indices
     idonex = 0
     idoney = 0

!-- simply move particles to the right
     do ipart = 1, prt_nmax
!-- ignore dead and outbound particles
        if(.not.prt_particles(ipart)%lalive) cycle
        prt_ldone = .false.
        do while(.not.prt_ldone)
           if(prt_particles(ipart)%ix==grd_nx) then
!-- made it to the right of the global domain, particle is finished
              prt_ndeadall = prt_ndeadall+1
              prt_ndead = prt_ndead+1
              prt_particles(ipart)%lalive = .false.
              prt_ldone = .true.
           else
!-- go right
              prt_particles(ipart)%ix=prt_particles(ipart)%ix+1
!-- short-cuts
              i = prt_particles(ipart)%ix
              j = prt_particles(ipart)%iy
!-- check if particle is still in sub-domain
              if(all(gas_icell/=grd_icell(i,j))) then
                 !-- add particle to a send buffer and de-activate it locally
                 prt_ldone = .true.
                 if(any(gas_ighostx==grd_icell(i,j))) then

                    prt_ndonex = prt_ndonex + 1
                    idonex = idonex + 1
                    prt_sndbuffer(idonex,1)=prt_particles(ipart)
                    if(idonex == prt_nbuf) then
!-- send the buffer to the x-neighbor
                       call mpi_send(prt_sndbuffer(:,1), &
                            sizeof(prt_sndbuffer(:,1)), MPI_BYTE, &
                            gas_neighbors(1),1,MPI_COMM_WORLD,ierr)
!-- reset x-buffer counter
                       idonex = 0
                       prt_sndbuffer(:,1)%lalive = .false.
                       prt_sndbuffer(:,1)%ix = -1
                       prt_sndbuffer(:,1)%iy = -1
                    endif
                 else
                    prt_ndoney = prt_ndoney + 1
                    idoney = idoney + 1
                    prt_sndbuffer(idoney,2)=prt_particles(ipart)
                    if(idoney == prt_nbuf) then
!-- send the buffer to the y-neighbor
                       call mpi_send(prt_sndbuffer(:,2), &
                            sizeof(prt_sndbuffer(:,2)), MPI_BYTE, &
                            gas_neighbors(2),2,MPI_COMM_WORLD,ierr)
!-- reset y-buffer counter
                       idoney = 0
                       prt_sndbuffer(:,2)%lalive = .false.
                       prt_sndbuffer(:,2)%ix = -1
                       prt_sndbuffer(:,2)%iy = -1
                    endif
                 endif
!-- deactivate local particle if buffered for sending
                 prt_particles(ipart)%lalive = .false.
                 prt_particles(ipart)%ix = -1
                 prt_particles(ipart)%iy = -1
              endif
           endif
        enddo
     enddo

!-- correct local particle count for losses
     src_nspr = src_nspr-prt_ndonex-prt_ndoney-prt_ndead
     !call mpi_waitall(1,irequest,istat,ierr)
!-- test for incoming particle buffer

!-- x-neighbor
     call mpi_test(irequest(1),iflag(1),istat,ierr)
     if(iflag(1)) then
        !if (impi==1) then
        write(6,*) ' here@@'
        !endif
        nbuf = count(prt_rcvbuffer(:,1)%lalive)
        src_nspr = src_nspr+nbuf
!-- add incoming particles to main particle array on this rank
        do ibuff = 1, nbuf
           do ipart = 1, prt_nmax
              if(.not.prt_particles(ipart)%lalive .and. &
                   prt_particles(ipart)%ix==-1) exit
           enddo
           prt_particles(ipart)=prt_rcvbuffer(ibuff,1)
           if(.not.prt_particles(ipart)%lalive) &
                stop 'added dead particle from x-bound'
        enddo
!-- repost the nonblocking receive
        call mpi_irecv(prt_rcvbuffer(:,1),nbuffsize,MPI_BYTE, &
             gas_neighbors(1),1,MPI_COMM_WORLD,irequest(1),ierr)
     endif

!-- y-neighbor
     call mpi_test(irequest(2),iflag(2),istat,ierr)
     if(iflag(2)) then
        nbuf = count(prt_rcvbuffer(:,2)%lalive)
        src_nspr = src_nspr+nbuf
!-- add incoming particles to main particle array on this rank
        do ibuff = 1, nbuf
           do ipart = 1, prt_nmax
              if(.not.prt_particles(ipart)%lalive .and. &
                   prt_particles(ipart)%ix==-1) exit
           enddo
           prt_particles(ipart)=prt_rcvbuffer(ibuff,2)
           if(.not.prt_particles(ipart)%lalive) &
                stop 'added dead particle from y-bound'
        enddo
!-- repost the nonblocking receive
        call mpi_irecv(prt_rcvbuffer(:,2),nbuffsize,MPI_BYTE, &
             gas_neighbors(2),2,MPI_COMM_WORLD,irequest(2),ierr)
     endif

!-- check if there are no more active local particles
     if(src_nspr==0) then
!-- send any partially full buffers
!-- over the x-boundary

        if(any(prt_sndbuffer(:,1)%lalive)) then
           if (impi==0) then
              do i=1,count(prt_sndbuffer(:,1)%lalive)
                 write(6,'(A12,I9,I9)') '> =',prt_sndbuffer(i,1)%ix,&
                      &prt_sndbuffer(i,1)%iy
              enddo
           endif
           call mpi_send(prt_sndbuffer(:,1), nbuffsize, &
                MPI_BYTE, gas_neighbors(1),1,MPI_COMM_WORLD,ierr)
           prt_sndbuffer(:,1)%lalive = .false.
           prt_sndbuffer(:,1)%ix = -1
           prt_sndbuffer(:,1)%iy = -1
        endif
!-- over the y-boundary
        if(any(prt_sndbuffer(:,2)%lalive)) then
           call mpi_send(prt_sndbuffer(:,2),nbuffsize, &
                MPI_BYTE,gas_neighbors(2),2,MPI_COMM_WORLD,ierr)
           prt_sndbuffer(:,2)%lalive = .false.
           prt_sndbuffer(:,2)%ix = -1
           prt_sndbuffer(:,2)%iy = -1
        endif

!-- collect counts of dead particles to master rank
        if(lmpi0) then
!-- count dead particles
           !prt_ndeadall = prt_ndeadall + prt_ndeadall
           do iimpi = 1,nmpi-1
              call mpi_test(ireqslav(iimpi),iflagslav(iimpi),istat,ierr)
              if(iflagslav(iimpi)) then
                 prt_ndeadall=prt_ndeadall+prt_ncompleted(iimpi)
!-- respost nonblocking receive for particle count
                 call mpi_irecv(prt_ncompleted(iimpi),1,MPI_INTEGER, &
                      iimpi,3,MPI_COMM_WORLD,ireqslav(iimpi),ierr)
              endif
           enddo
!-- check if total dead count equals the original number
           if(prt_ndeadall==src_ns) then
              prt_lalldone = .true.
              do iimpi = 1,nmpi-1
                 call mpi_send(prt_lalldone,1,MPI_LOGICAL,iimpi,4, &
                      MPI_COMM_WORLD,ierr)
              enddo
           endif
        else
!-- send the local particle completed count
           call mpi_send(prt_ndeadall,1,MPI_INTEGER,impi0,3, &
                MPI_COMM_WORLD,ierr)
           prt_ndeadall = 0
        endif

!-- cancel open nonblocking receive requests if finished
        if(prt_lalldone) then
!-- the transport is finished, cancel all non-blocking receives
           call mpi_cancel(irequest(1),ierr)
           call mpi_cancel(irequest(2),ierr)
           if(lmpi0) then
              do iimpi = 1,nmpi-1
                 call mpi_cancel(ireqslav(iimpi),ierr)
              enddo
           else
              !call mpi_cancel(irequest(3),ierr)
              call mpi_cancel(irequestall,ierr)
           endif
        endif

     endif
!-- increment outer iteration counter
     itotiter = itotiter + 1
  enddo

!-- now reprint particle set sizes to check transport
  write(0,*) 'impi =', impi,': finished particles =',count(prt_particles%ix/=-1)

!-- simulation finished
  call mpi_finalize(ierr)



end program milagro_strip
