*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module milagromod
c     -------------------
      implicit none
      INCLUDE 'mpif.h'
c
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
      end module milagromod
c vim: fdm=marker
