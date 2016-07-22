program testmpi
  implicit none
  include 'mpif.h'
!--------------------
!- A milagro mini-app
!--------------------
!-- mpi variables
  integer,parameter :: impi0=0 !the master rank
  logical :: lmpi0 !true for the master rank
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,istat(MPI_STATUS_SIZE,4),i,nb
  integer,allocatable,dimension(:)::ireq
  logical,allocatable,dimension(:)::iflag
  integer,allocatable,dimension(:)::rnks
  logical :: alldone
  !integer :: left,right
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  allocate(rnks(nmpi),ireq(nmpi),iflag(nmpi))
  alldone=.false.
  nb=0
  rnks=-1
  do i=0,3
     !if(i/=impi) then
     call mpi_irecv(rnks(i+1),4, &
          MPI_BYTE,i,1,MPI_COMM_WORLD,ireq(i+1),ierr)
     !endif
  enddo
  do i=0,3
     !if(i/=impi) then
     call mpi_send(impi*10,4, &
          MPI_BYTE,i,1,MPI_COMM_WORLD,ierr)
     !endif
  enddo
  call mpi_waitall(4,ireq,istat,ierr)


  do i=0,3
     !if(i/=impi) then
     call mpi_test(ireq(i+1),iflag(i+1),istat,ierr)
     if(iflag(i+1)) then
        nb=nb+1
        write(6,'(A5,I2,A7,I2,I2)') &
             &'rank(',impi,') recv:',rnks(i+1),nb
     endif
     !endif
  enddo

!-- simulation finished
  call mpi_finalize(ierr)
  deallocate(rnks,ireq,iflag)
end program testmpi
