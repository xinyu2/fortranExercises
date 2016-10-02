!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
module milagromod
  use particlemod
!     -------------------
  implicit none
  INCLUDE 'mpif.h'
!
  integer,parameter :: BUFFSIZE=3
  integer,private :: ierr,it,istat(MPI_STATUS_SIZE),i,ii
  integer,parameter :: rmpbtag=5,rpctag=10,rftag=15,&
       &sndtag=20,snd2tag=22,rbftag=20,rbf2tag=22
  !********************************
  !* reqs:
  !* reqmbs: max buffer size
  !* reqrb:  receive buffer
  !* reqpc:  particle complete
  !* reqgf:  global finish
  !* rbflag: receive buffer flag
  !********************************
  integer,dimension(:),allocatable :: pcomplete  ! array for number of particle complete on master
  logical :: globalFinish    ! global finish tag on slaves
  integer :: ipcmplt         ! number of particle complete on slave(buff, dead, flux)
  integer,dimension(:),allocatable :: reqrb,reqrb2,reqpc,cmpind
  integer :: reqgf,outreq
  logical,dimension(:),allocatable :: rbflag,rb2flag,rpcflag
  logical :: rgfflag
  logical,dimension(:),allocatable :: lhasBuff
  !**************
  !* dictionary
  !**************
  type cell
     integer :: gid            ! global id
     integer :: rid            ! rank id
  end type cell
  type(cell),dimension(:),allocatable,target :: dd ! reordered cell index
  type(cell), allocatable, dimension(:),target :: myGhosts

  !***************************
  !* neighbor adjacent matrix
  !***************************
  integer,dimension(:),allocatable :: nbrs !neighbor list for each rank
  integer :: totnbr ! total number of neighbors for each rank
  integer :: flag           ! adjacency flag between neighbors
  !***************************
  !* send/receive buffer
  !***************************
  type(packet),allocatable,target  :: sndbuff(:,:), rcvbuff(:,:)
  type(packet2),allocatable,target :: sndbuff2(:,:),rcvbuff2(:,:)
  integer,dimension(:),allocatable :: db2r ! (totnbr): buffIdx->rankNbr
  integer,dimension(:),allocatable :: sndidx ! (totnbr): capacity pointer
  integer pckttyp, pckttyp2, pbtyp, celltype, oldtypes(0:2) ! required variables
  integer blockcounts(0:2), offsets(0:2), extent
  type pb
     type(packet) :: pb1
     type(packet2):: pb2
  end type pb
  type(pb),allocatable,target :: sndbf(:,:),rcvbf(:,:)
  !***********
  !* tree
  !***********
  type trnode
     integer :: parent
     integer :: lchild,rchild
  end type trnode
  type(trnode),target::tn

contains

! -----------------------------------------------------
! subroutines
! -----------------------------------------------------
  !**********************
  ! ghost & dictionary  *
  !**********************
  subroutine check_neighbor2(myghosts, dd, a_dum, counter, x_d, y_d, z_d, &
       & side, outer, inner,   box, x,y,z)
    implicit none
    integer, intent(in) :: x_d, y_d, z_d, side
    integer, intent(in) :: outer, inner,  x, y, z
    integer, intent(in) :: box(x,y,z)
    type(cell), intent(in),    dimension(:) :: dd
    integer, intent(in) :: a_dum(x_d,y_d,z_d)
    type(cell), intent(inout), dimension(:) :: myghosts
    integer, intent(inout) :: counter
    integer :: escaped, look_i, look_n
    integer :: i, j, rank, one, two
    do j=1, outer
       do i=1, inner
          escaped = -1

          !check if the cells are on the domain boundary (side)
          !which means the ghost cells are out of the domain (-1)
          select case(side)
          case(1)
             look_i = box(i,j,1)
             look_n = look_i-x_d*y_d
             yolo1: do one = 1, y_d
                do two = 1, x_d
                   if(look_i==a_dum(two,one,1)) then
                      escaped = 1
                      exit yolo1
                   endif
                enddo
             enddo yolo1
          case(2)
             look_i = box(i,j,z)
             look_n = look_i+x_d*y_d
             yolo2: do one = 1, y_d
                do two = 1, x_d
                   if(look_i==a_dum(two,one,z_d)) then
                      escaped = 1
                      exit yolo2
                   endif
                enddo
             enddo yolo2
          case(3)
             look_i = box(i,1,j)
             look_n = look_i-y_d+x_d !+x_d chenx
             yolo3: do one = 1, z_d
                do two = 1, x_d
                   if(look_i==a_dum(two,1,one)) then
                      escaped = 1
                      exit yolo3
                   endif
                enddo
             enddo yolo3
          case(4)
             look_i = box(i,y,j)
             look_n = look_i+y_d-x_d !-x_d chenx
             yolo4: do one = 1, z_d
                do two = 1, x_d
                   if(look_i==a_dum(two,y_d,one)) then
                      escaped = 1
                      exit yolo4
                   endif
                enddo
             enddo yolo4
          case(5)
             look_i = box(1,i,j)
             look_n = look_i-1
             yolo5: do one = 1, z_d
                do two = 1, y_d
                   if(look_i==a_dum(1,two,one)) then
                      escaped = 1
                      exit yolo5
                   endif
                enddo
             enddo yolo5
          case(6)
             look_i = box(x,i,j)
             look_n = look_i+1
             yolo6: do one = 1, z_d
                do two = 1, y_d
                   if(look_i==a_dum(x_d,two,one)) then
                      escaped = 1
                      exit yolo6
                   endif
                enddo
             enddo yolo6
          case default
             stop 'invalid side'
          end select
          !put ghost cell index into dictionary

          !if particles escaped, assign rank=-1, hence out of boundary
          if(escaped == 1) then
             !else, find the rank, which has the ghost cell
          else
             myghosts(counter)%gid = look_n
             rank=dd(look_n)%rid
             !put rank #, which contains the ghost cell index into dictionary
             myghosts(counter)%rid = rank
             counter = counter + 1
          endif
       enddo
    enddo
  end subroutine check_neighbor2

  subroutine getGhost(dd_a, dd, x_d, y_d, z_d, counts, displs, nmpi, impi, x,y,z, myghosts)
    use gridmod
    implicit none
    integer,intent(in) :: x_d, y_d, z_d
    integer,intent(in) :: nmpi,impi
    integer,intent(in),dimension(:) :: dd_a
    type(cell),intent(in),dimension(:) :: dd
    integer,intent(in),dimension(nmpi) :: counts, displs
    integer,intent(in) :: x,y,z

    type(cell), intent(out),dimension(:),allocatable :: myghosts

    integer :: i, n, num_ghosts, counter
    integer, dimension(x_d*y_d*z_d) :: a
    integer :: a_dum(x_d,y_d,z_d)
    integer, allocatable :: box(:,:,:)
    !create a 3D version of the domain
    a(:) = (/((i), i=1,size(a))/)
    a_dum = reshape(a(:),(/x_d,y_d,z_d/))

    n=impi+1
    !iterate through all of the ranks


    allocate(box(x,y,z))
    box = reshape(dd_a(displs(n):displs(n)+counts(n)-1),(/x,y,z/))

    select case(grd_igeom)
    case(1,11)
       num_ghosts=4
    case(2)
       num_ghosts=4*(x+y)
    case(3)
       num_ghosts=4*(x*y+y*z+x*z)
    case default
       stop 'invalid grd-igeom'
    end select

    allocate(myghosts(num_ghosts))
    myghosts%gid = -1
    myghosts%rid = -1
    counter = 1

    !get ghosts of all 6 sides
    call check_neighbor2(myghosts, dd, a_dum, counter, &
         & x_d, y_d, z_d, 1, y, x,   box, x,y,z)

    call check_neighbor2(myghosts, dd, a_dum, counter, &
         & x_d, y_d, z_d, 2, y, x,   box, x,y,z)

    call check_neighbor2(myghosts, dd, a_dum, counter, &
         & x_d, y_d, z_d, 3, z, x,   box, x,y,z)

    call check_neighbor2(myghosts, dd, a_dum, counter, &
         & x_d, y_d, z_d, 4, z, x,   box, x,y,z)

    call check_neighbor2(myghosts, dd, a_dum, counter, &
         & x_d, y_d, z_d, 5, z, y,   box, x,y,z)

    call check_neighbor2(myghosts, dd, a_dum, counter, &
         & x_d, y_d, z_d, 6, z, y,   box, x,y,z)
    deallocate(box)
  end subroutine getGhost

  subroutine myDimensions(dd_a, x, y, z, counts, displs, nmpi, impi, dim_x, dim_y,dim_z)
    implicit none
    integer, intent(in)  :: x,y,z
    integer, intent(in)  :: nmpi,impi
    integer, intent(in), dimension(z*y*x) :: dd_a
    integer, intent(in), dimension(nmpi) :: counts, displs
    integer, intent(out) ::dim_x, dim_y, dim_z
    integer :: i, j, k

    i=impi+1

    dim_x = 0
    dim_y = 0

    !calculate dimension 'x' of the rank
    do j=displs(i), displs(i)+counts(i)-1
       if(dd_a(j)==x*y*z) then
          dim_x = dim_x + 1
          exit
       endif

       if (dd_a(j)<dd_a(displs(i))+x) then
          dim_x = dim_x + 1
       else
          exit
       endif
       if (dd_a(j+1)-dd_a(j)/=1) exit
    enddo

    !get absolute 'z' boundary of the rank to calculate 'y'
    do k=1, z
       if(x*y*(k-1)<dd_a(displs(i)) .and. dd_a(displs(i))<=x*y*k) exit
    enddo

    !calculate dimension 'y' of the rank
    do j=displs(i), displs(i)+counts(i)-1, dim_x
       if(dd_a(j)>k*x*y) exit
       dim_y = dim_y + 1
    enddo

    !calculate dimension 'z' of the rank
    dim_z = counts(i)/(dim_x*dim_y)
  end subroutine myDimensions

  subroutine getNeighbors(nmpi,myghosts,nbrs,totnbr)
    implicit none
    integer,intent(in)::nmpi
    type(cell),intent(in),dimension(:) :: myghosts
    integer,dimension(:),intent(out),allocatable::nbrs
    integer,intent(out) :: totnbr
    integer::i,r
    allocate(nbrs(nmpi))
    !******************
    !* initialize nbrs
    !******************
    nbrs=0
    totnbr=0
    !********************
    !* nbrs for strip-dd
    !********************
    do i=1,size(pack(myghosts%rid,myghosts%rid>-1))
       r=myghosts(i)%rid
       nbrs(r+1)=1
    enddo
    totnbr=size(pack(nbrs,nbrs==1))
  end subroutine getNeighbors

  subroutine nbrDict(nmpi,totnbr,nbrs,db2r)
    implicit none
    integer,intent(in)::nmpi
    integer,intent(in) :: totnbr
    integer,dimension(:),intent(inout),allocatable::nbrs
    integer,dimension(:),intent(out),allocatable::db2r
    integer::i,temptot
    allocate(db2r(totnbr))
    temptot=totnbr
    do i=nmpi,1,-1
       if(nbrs(i)==1)then
          nbrs(i)=temptot
          db2r(temptot)=i-1
          temptot=temptot-1
       endif
    enddo
    if(temptot/=0)then
       stop 'build neighbor dictionary error'
    endif
  end subroutine nbrDict

  subroutine buildDictionary(dd_indexes,counts,displs,dd)
    implicit none
    integer, intent(in),  dimension(:) :: dd_indexes
    integer, intent(in),  dimension(:) :: counts, displs
    type(cell), intent(out), dimension(:) :: dd
    integer :: i,dict_idx,rnk,cnt
    do rnk=0,size(counts)-1
       do cnt = 1, counts(rnk+1)
          i = displs(rnk+1)+cnt - 1
          dict_idx=dd_indexes(i)
          dd(dict_idx)%gid=dict_idx
          dd(dict_idx)%rid=rnk
       enddo
    enddo
  end subroutine buildDictionary
  !************
  ! mpi APIs  *
  !************
  subroutine define_mpi_str
    implicit none
    !------------------
    !*********************************************
    !* mpi structure for passing particle-packet1
    !*********************************************
    ! setup description of the 7 MPI_REAL fields x, y, z,
    ! mu, om, t, e, e0, wl
    offsets(0) = 0
    oldtypes(0) = MPI_REAL
    blockcounts(0) = 9

    ! define structured type and commit it
    call MPI_TYPE_STRUCT(1, blockcounts, offsets, oldtypes, &
         pckttyp, ierr)
    call MPI_TYPE_COMMIT(pckttyp, ierr)

    !*********************************************
    !* mpi structure for passing particle-packet2
    !*********************************************
    ! setup description of the 4 MPI_REAL fields:
    ! mux, muy, muz, dist
    offsets(0) = 0
    oldtypes(0) = MPI_REAL
    blockcounts(0) = 4

    ! setup description of the  MPI_CHAR fields: stat
    call MPI_TYPE_EXTENT(MPI_REAL, extent, ierr)
    offsets(1) =  4 * extent
    oldtypes(1) = MPI_CHAR
    blockcounts(1) = 4

    ! setup description of the  MPI_INTEGER fields:
    ! ix, iy, iz, ic, ig, itype, ipart, istep, idist
    call MPI_TYPE_EXTENT(MPI_CHAR, extent, ierr)
    offsets(2) =  offsets(1) + 4 * extent
    oldtypes(2) = MPI_INTEGER
    blockcounts(2) = 9

    ! define structured type and commit it
    call MPI_TYPE_STRUCT(3, blockcounts, offsets, oldtypes, &
         pckttyp2, ierr)
    call MPI_TYPE_COMMIT(pckttyp2, ierr)
    !*********************************************
    !* mpi structure for passing pb
    !*********************************************
    ! setup description of the 4 MPI_REAL fields:
    ! mux, muy, muz, dist
    offsets(0) = 0
    oldtypes(0) = MPI_REAL
    blockcounts(0) = 13

    ! setup description of the  MPI_CHAR fields: stat
    call MPI_TYPE_EXTENT(MPI_REAL, extent, ierr)
    offsets(1) =  4 * extent
    oldtypes(1) = MPI_CHAR
    blockcounts(1) = 4

    ! setup description of the  MPI_INTEGER fields:
    ! ix, iy, iz, ic, ig, itype, ipart, istep, idist
    call MPI_TYPE_EXTENT(MPI_CHAR, extent, ierr)
    offsets(2) =  offsets(1) + 4 * extent
    oldtypes(2) = MPI_INTEGER
    blockcounts(2) = 9

    ! define structured type and commit it
    call MPI_TYPE_STRUCT(3, blockcounts, offsets, oldtypes, &
         pbtyp, ierr)
    call MPI_TYPE_COMMIT(pbtyp, ierr)
  end subroutine define_mpi_str

  subroutine recvMaxParBuff(myrank,db2r)
    implicit none
    !------------------
    integer,intent(in)::myrank
    integer,dimension(:),intent(in)::db2r
    integer::i,source
    do i=1,size(db2r)
       source=db2r(i)
       !call mpi_irecv(rcvbuff(1,i),BUFFSIZE,pckttyp,&
       !     &source,rbftag,MPI_COMM_WORLD,reqrb(i),ierr)
       !call mpi_irecv(rcvbuff2(1,i),BUFFSIZE,pckttyp2,&
       !     &source,rbf2tag,MPI_COMM_WORLD,reqrb2(i),ierr)
       call mpi_irecv(rcvbf(1,i),BUFFSIZE,pbtyp,&
            &source,rbftag,MPI_COMM_WORLD,reqrb(i),ierr)
       !write(6,*) 'rcv pckt2-> buff2',source
    enddo
  end subroutine recvMaxParBuff

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

  subroutine recvParentFinish()
    implicit none
    call mpi_irecv(globalFinish,1,MPI_LOGICAL,tn%parent,rftag,&
         & MPI_COMM_WORLD,reqgf,ierr)
  end subroutine recvParentFinish

  subroutine checkBuffer(myrank,db2r,lhas)
    implicit none
    integer,intent(in)::myrank
    integer,dimension(:),intent(in)::db2r
    logical,dimension(:),intent(out)::lhas
    integer ::i
    lhas =.true.
    do i=1,size(db2r)
       call mpi_test(reqrb(i),rbflag(i),istat,ierr)
       lhas(i)=lhas(i) .and. rbflag(i)
       !call mpi_test(reqrb2(i),rb2flag(i),istat,ierr)
       !lhas(i)=lhas(i) .and. rb2flag(i)
    enddo
  end subroutine checkBuffer

  !subroutine writeBuffer(myrank,mythrd,sndbuff,sndbuff2,&
  subroutine writeBuffer(myrank,mythrd,sndbf,&
       &sndidx,r,p,p2)
    implicit none
    integer,intent(in)::myrank,mythrd
    !type(packet),dimension(:,:),intent(inout)::sndbuff
    !type(packet2),dimension(:,:),intent(inout)::sndbuff2
    type(pb),dimension(:,:),intent(inout)::sndbf
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::r
    type(packet),intent(in)::p
    type(packet2),intent(in)::p2
    integer::j
    sndidx(r)=sndidx(r)+1 ! increment number of buffered par
    j=sndidx(r)
    !if(myrank==0) then
    !   write(6,*) '@wrtbuff',mythrd,j,r
    !endif
    !sndbuff(j,r)=p
    !sndbuff2(j,r)=p2
    sndbf(j,r)%pb1=p
    sndbf(j,r)%pb2=p2
  end subroutine writeBuffer

  !subroutine sendBuffer(myrank,mythrd,sndbuff,sndbuff2,&
  subroutine sendBuffer(myrank,mythrd,sndbf,&
       &sndidx,r,newrank)
    implicit none
    integer,intent(in)::myrank,mythrd
    !type(packet),dimension(:,:),intent(inout)::sndbuff
    !type(packet2),dimension(:,:),intent(inout)::sndbuff2
    type(pb),dimension(:,:),intent(inout)::sndbf
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::r,newrank
    integer::count
    count=sndidx(r) ! increment number of buffered par
    !if(myrank==0) then ! don't print
    !   write(6,'(A15,I2,A2,I5,A8,I2,A3,I3)') '>>>send@rank(',myrank,&
    !        &') ',count,'prs to (',newrank,') @',mythrd
    !   write(6,*) '@send',mythrd,newrank,r,'>bf',sndbuff(:,r)
    !endif
    !call mpi_send(sndbuff(1,r),BUFFSIZE,pckttyp,newrank,&
    !     &sndtag,MPI_COMM_WORLD,ierr)
    !call mpi_send(sndbuff2(1,r),BUFFSIZE,pckttyp2,newrank,&
    !     &snd2tag,MPI_COMM_WORLD,ierr)
    call mpi_send(sndbf(1,r),BUFFSIZE,pbtyp,newrank,&
         &sndtag,MPI_COMM_WORLD,ierr)
    !call cleanBuffer(sndbuff,sndbuff2,sndidx,r)
    call cleanBuffer(sndbf,sndidx,r)
  end subroutine sendBuffer

  !subroutine cleanBuffer(sndbuff,sndbuff2,sndidx,r)
  subroutine cleanBuffer(sndbf,sndidx,r)
    implicit none
    !integer,intent(in)::myrank
    !type(packet),dimension(:,:),intent(inout)::sndbuff
    !type(packet2),dimension(:,:),intent(inout)::sndbuff2
    type(pb),dimension(:,:),intent(inout)::sndbf
    integer,dimension(:),intent(inout)::sndidx
    integer,intent(in)::r
    integer :: j
    !write(6,'(A12,I2,A5,I2,A2)') 'clean@rank(',myrank,').bf(',newrank,')'
    !do j=1,sndidx(r)
    !sndbuff(:,r) =packet(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    !sndbuff2(:,r)=packet2(0.0,0.0,0.0,0.0,' ',0,0,0,0,0,0,0,0,0)
    sndbf(:,r)%pb1 =packet(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    sndbf(:,r)%pb2=packet2(0.0,0.0,0.0,0.0,' ',0,0,0,0,0,0,0,0,0)
    !enddo
    sndidx(r)=0 ! clean buffer top
  end subroutine cleanBuffer

  subroutine tallyPcomplete(myrank,tn,ncmplt)
    implicit none
    !------------------
    integer,intent(in)::myrank
    type(trnode),intent(in)::tn
    integer,intent(inout)::ncmplt
    integer::i,idx,chd
    if((tn%lchild>0).or.(tn%rchild>0)) then
       call mpi_testsome(2,reqpc,outreq,cmpind,istat,ierr)
       if(outreq>0) then
          do i=1,outreq
             idx=cmpind(i)
             ncmplt=ncmplt+pcomplete(i)
             if (idx==1) then
                chd=tn%lchild
             else
                chd=tn%rchild
             endif
             call mpi_irecv(pcomplete(idx),1,MPI_INTEGER,chd,rpctag,&
                  & MPI_COMM_WORLD,reqpc(idx),ierr)
             if(myrank<0) write(6,*) '@tally',myrank,ncmplt
          enddo
       endif
    endif
  end subroutine tallyPcomplete
  !*****************
  ! buffer manage  *
  !*****************
  subroutine initBuffer(impi,totnbr,nomp)
    implicit none
    integer,intent(in) :: impi,totnbr,nomp
    if(.not.allocated(rcvbuff)) then
       allocate(rcvbuff(BUFFSIZE,totnbr),rcvbuff2(BUFFSIZE,totnbr))
    endif
    if(.not.allocated(sndbuff)) then
       allocate(sndbuff(BUFFSIZE,totnbr),sndbuff2(BUFFSIZE,totnbr))
    endif
    if(.not.allocated(sndbf)) then
       allocate(sndbf(BUFFSIZE,totnbr))
    endif
    if(.not.allocated(reqrb)) then
       allocate(reqrb(totnbr),reqrb2(totnbr))
    endif
    if(.not.allocated(pcomplete)) then
       allocate(pcomplete(2))
    endif
    if(.not.allocated(reqpc)) then
       allocate(reqpc(2))
    endif
    if(.not.allocated(rbflag)) then
       allocate(rbflag(totnbr),rb2flag(totnbr))
    endif
    if(.not.allocated(lhasBuff)) then
       allocate(lhasBuff(totnbr))
    endif
    if(.not.allocated(sndidx)) then
       allocate(sndidx(totnbr))
    endif
    reqrb=0
    reqrb2=0
    pcomplete=0
    reqpc=0
    rbflag=.false.
    rb2flag=.false.
    lhasBuff=.false.
    rcvbuff=packet(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    rcvbuff2=packet2(0.0,0.0,0.0,0.0,' ',0,0,0,0,0,0,0,0,0)
    sndbuff=packet(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    sndbuff2=packet2(0.0,0.0,0.0,0.0,' ',0,0,0,0,0,0,0,0,0)
    sndbf%pb1=packet(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    sndbf%pb2=packet2(0.0,0.0,0.0,0.0,' ',0,0,0,0,0,0,0,0,0)
    sndidx=0
  end subroutine initBuffer
! -----------------------------------------------------
! functions
! -----------------------------------------------------
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

  integer function getRankId(myghosts,pic)
    implicit none
    type(cell),intent(in),dimension(:) :: myghosts
    integer,intent(in) ::pic
    integer::i,r
    do i=1,size(myghosts)
       if(myghosts(i)%gid==pic) then
          r=myghosts(i)%rid
          exit
       endif
    enddo
    getRankId=r
  end function getRankId
end module milagromod
! vim: fdm=marker
