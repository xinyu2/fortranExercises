program testmpi
  implicit none
  include 'mpif.h'
!--------------------
!- A milagro mini-app
!--------------------
!-- mpi variables
  integer,parameter :: impi0=0 !the master rank
  logical :: lmpi0 !true for the master rank
  logical :: lalldone,ldone
  integer :: impi !mpi rank
  integer :: nmpi !number of mpi tasks
  integer :: ierr,istat(MPI_STATUS_SIZE,4),i,nb,sumnb
  integer,allocatable,dimension(:)::ireq
  logical,allocatable,dimension(:)::iflag
  logical,allocatable,dimension(:)::lflag
  integer,allocatable,dimension(:)::rnks
  integer,parameter :: alldonetag=10
  integer :: ireqalldone,istatalldone(MPI_STATUS_SIZE)
  logical :: iflagalldone

!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
  allocate(rnks(nmpi),ireq(nmpi),iflag(nmpi),lflag(nmpi))
  lalldone=.false.
  ldone=.false.
  nb=0
  rnks=-1
  lflag=.false.
  do while(.not.lalldone )
     if(impi/=0) then
        call mpi_irecv(lalldone,4, &
             MPI_BYTE,0,alldonetag,MPI_COMM_WORLD,ireqalldone,ierr)
        call mpi_test(ireqalldone,iflagalldone,istatalldone,ierr)
        if(iflagalldone) then
           write(6,'(A5,I2,A20)') 'rank(',impi,') recv alldone'
        endif
     endif
     if (.not.ldone) then
        do i=0,3
           call mpi_irecv(rnks(i+1),4, &
                MPI_BYTE,i,1,MPI_COMM_WORLD,ireq(i+1),ierr)
        enddo
     endif
     if(.not.ldone) then
        do i=0,3
           call mpi_send(impi*10,4, &
                MPI_BYTE,i,1,MPI_COMM_WORLD,ierr)
        enddo
     endif
     !call mpi_waitall(4,ireq,istat,ierr)
     if(.not.ldone) then
        do i=0,3
           call mpi_test(ireq(i+1),iflag(i+1),istat,ierr)
           !write(6,*) '>test',impi,iflag
           if(iflag(i+1)) then
              lflag(i+1)=iflag(i+1)
              nb=count(lflag .eqv. .true.)
              if(nb==4) ldone=.true.
              !if(nb<4) then
              !   nb=nb+1
              !endif
              write(6,'(A5,I2,A7,I5,I5)') &
                   &'rank(',impi,') recv:',rnks(i+1),nb
           endif
        enddo
     endif
     call mpi_reduce(nb,sumnb,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

     if(impi==0) then
        if (sumnb==16) lalldone=.true.
        write(6,*) 'sum nb=',sumnb
        do i = 1,nmpi-1
           call mpi_send(lalldone,1,MPI_LOGICAL,i,alldonetag, &
                MPI_COMM_WORLD,ierr)
        enddo
     endif
  enddo
!-- simulation finished
  call mpi_finalize(ierr)
  deallocate(rnks,ireq,iflag)
end program testmpi
