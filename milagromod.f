*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module milagromod
      use particlemod
c     -------------------
      implicit none
      INCLUDE 'mpif.h'
c
      integer,parameter :: BUFFSIZE=128
      integer,private :: ierr,it,istat(MPI_STATUS_SIZE),i,ii
      integer,parameter :: rmpbtag=5,rpctag=10,rftag=15,sndtag=20,
     &rbftag=20
c     !********************************
c     !* reqs:
c     !* reqmbs: max buffer size
c     !* reqrb:  receive buffer
c     !* reqpc:  particle complete
c     !* reqgf:  global finish
c     !* rbflag: receive buffer flag
c     !********************************
      integer,dimension(:),allocatable :: pcomplete  ! array for number of particle complete on master
      logical :: globalFinish   ! global finish tag on slaves
      integer :: pcmplt         ! number of particle complete on slave(cross boundary)
      integer,dimension(:),allocatable :: reqrb,reqpc,cmpind
      integer :: reqgf,outreq
      logical,dimension(:),allocatable :: rbflag,rpcflag
      logical :: rgfflag
c
c     !**************
c     !* dictionary
c     !**************
      type cell
      integer :: gid            ! global id
      integer :: rid            ! rank id
      end type cell
      type(cell),dimension(:),allocatable,target :: dd ! reordered cell index
c
c     !***************************
c     !* neighbor adjacent matrix
c     !***************************
      integer,dimension(:),allocatable :: nbrs !neighbor list for each rank
      integer :: totnbr ! total number of neighbors for each rank
      integer :: flag           ! adjacency flag between neighbors
c     !***************************
c     !* send/receive buffer
c     !***************************
      type(packet),allocatable,target  :: sndbuff(:,:), rcvbuff(:,:)
      type(packet2),allocatable,target :: sndbuff2(:,:),rcvbuff2(:,:)
c
c     !***********
c     !* tree
c     !***********
      type trnode
      integer :: parent
      integer :: lchild,rchild
      end type trnode
      type(trnode),target::tn
c
c-- explicit interfaces
c     interface
c
      contains
c
c
c
c
      subroutine tallyPcomplete(myrank,tn,ncmplt)
      implicit none
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
               call mpi_irecv(pcomplete(idx),1,MPI_INTEGER,chd,rpctag,
     & MPI_COMM_WORLD,reqpc(idx),ierr)
               if(myrank<0) write(6,*) '@tally',myrank,ncmplt
            enddo
         endif
      endif
      end subroutine tallyPcomplete
c
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
c
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
      end module milagromod
c vim: fdm=marker
