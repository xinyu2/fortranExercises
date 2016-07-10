*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module mpimod
c     -------------
      implicit none
      INCLUDE 'mpif.h'
c
      integer,parameter :: impi0=0 !the master rank
      logical :: lmpi0 !true for the master rank
      integer :: impi !mpi rank
      integer :: nmpi !number of mpi tasks
      integer,private :: ierr
c
      integer,private,allocatable :: counts(:)
      integer,private,allocatable :: displs(:)
c
      save
c
      contains
c
c
c
      subroutine bcast_permanent
c     --------------------------!{{{
      use inputparmod
      use inputstrmod
      use sourcemod
      use ionsmod
      use ffxsmod
      use bfxsmod
      use bbxsmod
      use gridmod
      use gasmod
      use groupmod
      use fluxmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
************************************************************************
      integer :: i,n
      integer :: il,ii,ir,ic
      integer :: nx,ny,nz
      logical,allocatable :: lsndvec(:)
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
      character(4),allocatable :: csndvec(:)
c
c-- inputparmod variables
c========================
c-- create pointer arrays
      call inputpar_create_pointers(il,ii,ir,ic)
c-- broadcast logicals
      n = il
      allocate(lsndvec(n))
      forall(i=1:n) lsndvec(i) = in_l(i)%p
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_l(i)%p = lsndvec(i)
      deallocate(lsndvec)
c-- broadcast integers
      n = ii
      allocate(isndvec(n))
      forall(i=1:n) isndvec(i) = in_i(i)%p
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_i(i)%p = isndvec(i)
      deallocate(isndvec)
c-- broadcast real*8
      n = ir
      allocate(sndvec(n))
      forall(i=1:n) sndvec(i) = in_r(i)%p
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_r(i)%p = sndvec(i)
      deallocate(sndvec)
c-- broadcast characters
      n = ic
      allocate(csndvec(n))
      forall(i=1:n) csndvec(i) = in_c(i)%p
      call mpi_bcast(csndvec,n*4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
!     forall(i=1:n) in_c(i)%p = csndvec(i) !avoid gfortran 4.6.3 compiler bug
      do i=1,n
       in_c(i)%p = csndvec(i)
      enddo
      deallocate(csndvec)
c
c-- set number of threads, i.e. non-automatic
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c
c-- everything else
c==================
c-- broadcast constants
c-- logical
      n = 3
      allocate(lsndvec(n))
      if(lmpi0) lsndvec = [str_lvoid,str_ltemp,str_lye]
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      str_lvoid = lsndvec(1)
      str_ltemp = lsndvec(2)
      str_lye = lsndvec(3)
      deallocate(lsndvec)
c
c-- integer
      n = 5
      allocate(isndvec(n))
      if(lmpi0) isndvec = [ion_nion,ion_iionmax,bb_nline,
     &  str_nc,str_nabund]
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      ion_nion     = isndvec(1)
      ion_iionmax  = isndvec(2)
      bb_nline     = isndvec(3)
      str_nc       = isndvec(4)
      str_nabund   = isndvec(5)
      deallocate(isndvec)
cc
cc-- real*8
c      n = 1
c      allocate(sndvec(n))
c      if(lmpi0) sndvec = [dummy]
c      call mpi_bcast(sndvec,n,MPI_REAL8,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc-- copy back
c      dummy        = sndvec(1)
c      deallocate(sndvec)
c
c-- dimenstions
      nx = in_ndim(1)
      ny = in_ndim(2)
      nz = in_ndim(3)
c
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       if(bb_nline>0) allocate(bb_xs(bb_nline))
       allocate(str_xleft(nx+1))
       allocate(str_yleft(ny+1))
       allocate(str_zleft(nz+1))
       allocate(str_idcell(str_nc))
       if(str_nabund>0) allocate(str_iabund(str_nabund))
      endif
c
c-- inputstr
      call mpi_bcast(str_xleft,nx+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_yleft,ny+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_zleft,nz+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_idcell,str_nc,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      if(str_nabund>0) then
       call mpi_bcast(str_iabund,str_nabund,MPI_INTEGER,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c
c-- broadcast data
c-- bound-bound
      if(bb_nline>0) then
       call mpi_bcast(bb_xs,sizeof(bb_xs),MPI_BYTE,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c-- bound-free
      call mpi_bcast(bf_ph1,6*7*30*30,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(bf_ph2,7*30*30,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- free-free
      call mpi_bcast(ff_gff,ff_nu*ff_ngg,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call bcast_ions
c
c
      contains
c
      subroutine bcast_ions
c     ---------------------!{{{
      implicit none
************************************************************************
* broadcast the ions data structure
************************************************************************
      integer :: ii,iz,iion,n
      real*8 :: vec(1000)
      integer :: nion(gas_nelem)
      integer :: nlev(ion_nion)
      real*8 :: e(ion_nion)
c
c-- evaluate shape info
      if(lmpi0) then
       iion = 0
       do iz=1,gas_nelem
        nion(iz) = ion_el(iz)%ni
        do ii=1,ion_el(iz)%ni
         iion = iion + 1
         nlev(iion) = ion_el(iz)%i(ii)%nlev
         e(iion) = ion_el(iz)%i(ii)%e
        enddo !ii
       enddo !iz
c-- sanity check
       if(iion/=ion_nion) stop "bcast_perm: ion_nion problem"
      endif
c
c-- bcast shape info and allocate
      call mpi_bcast(nion,gas_nelem,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nlev,ion_nion,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- allocate structure
      if(impi/=impi0) call ions_alloc_el(gas_nelem,nion,ion_nion,nlev)
c
c-- fill structure
      call mpi_bcast(e,ion_nion,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      iion = 0
      do iz=1,gas_nelem
       do ii=1,ion_el(iz)%ni
        iion = iion + 1
        n = nlev(iion)
c-- eion
        if(impi/=impi0) ion_el(iz)%i(ii)%e = e(iion)
c-- elev
        if(lmpi0) vec(:n) = ion_el(iz)%i(ii)%elev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%elev = vec(:n)
c-- glev
        if(lmpi0) vec(:n) = ion_el(iz)%i(ii)%glev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%glev = vec(:n)
       enddo !ii
      enddo !iz
c-- sanity check
      if(iion/=ion_nion) stop "bcast_perm: ion_nion problem"
c!}}}
      end subroutine bcast_ions
c!}}}
      end subroutine bcast_permanent
c
c
c
      subroutine scatter_inputstruct(ndim,icell1,ncell,igeom)
c     -------------------------------------------------!{{{
      use inputstrmod
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3),igeom
      integer,intent(out) :: icell1,ncell
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
************************************************************************
      !-- milagro dd vars
      integer :: i,n,nc,nx,ny,nz
      integer :: edgelen,nboxes !,reorderidx
      integer, dimension(:), allocatable :: corners
!     milagro: reorder the origin arrays: str_mass,str_temp,etc
      integer, dimension(:), allocatable :: arr_reorder
!     milagro: store everything to be scattered in reordered index
      real*8,allocatable :: arr_scatter(:) !(nc)
      real*8,allocatable :: reorder_massfrdc(:,:) !(nabund,nc)
c
      nx = ndim(1)
      ny = ndim(2)
      nz = ndim(3)

      write(6,*) '>> nx=',nx, 'ny=',ny,'nz=',nz

      if(impi==impi0) then
         edgelen=getedge(nmpi,ndim,igeom) ! cube edge's length
      endif

c
c-- calculate offsets
      allocate(counts(nmpi),displs(nmpi))
      displs(1) = 0
      nc = str_nc

!     milagro dd: put arrays to arr_reorder
!     then scatter them by mpiscatterv
!     this is done by master
      if(impi==impi0) then
         select case(igeom)
         case(1,11) !arr_reorder=raw order(xyz)
         case(2)
!           call getCorners(nx,ny,edgelen,nboxes,corners)
!           call getReorder(nx,edgelen,nboxes,corners,arr_reorder)
            call getCorners3d(igeom,nx,ny,nz,edgelen,nboxes,corners)
            call getReorder3d(igeom,nx,ny,nz,edgelen,nboxes,corners,
     &                        arr_reorder)
         case(3)
            call getCorners3d(igeom,nx,ny,nz,edgelen,nboxes,corners)
            call getReorder3d(igeom,nx,ny,nz,edgelen,nboxes,corners,
     &                        arr_reorder)
         case default
            stop 'igeom invalid'
         endselect
      endif

      do i=1,nmpi
       n = ceiling(nc/(nmpi-i+1d0))
       counts(i) = n
       if(i-1==impi) ncell = n
       if(i-1==impi) icell1 = displs(i) + 1
       if(i<nmpi) displs(i+1) = displs(i) + n
       nc = nc - n
      enddo

      if(sum(counts)/=str_nc) stop 'scatter_inputstr: counts/=str_nc'
      if(nc/=0) stop 'scatter_inputstr: nc/=0'

c
c-- allocate domain decomposed and domain compressed
!     if(impi/=impi0) allocate(str_massdc(str_nc))
!     milagro dd: use arr_scatter replace str_massdc
      if(impi==impi0) call getScatter(igeom,str_nc,str_massdc,
     &     arr_reorder,arr_scatter)
      if(impi/=impi0) allocate(arr_scatter(str_nc))
      allocate(str_massdd(ncell))

!     call mpi_scatterv(str_massdc,counts,displs,MPI_REAL8,
!     milagro dd: use arr_scatter replace str_massdc
      call mpi_scatterv(arr_scatter,counts,displs,MPI_REAL8,
     &        str_massdd,ncell,MPI_REAL8,
     &        impi0,MPI_COMM_WORLD,ierr)
c
c-- mass fractions if available
      if(str_nabund>0) then
!        if(impi/=impi0) allocate(str_massfrdc(str_nabund,str_nc))
         if(impi==impi0) call getMassfrdc(igeom,str_nabund,
     &     str_nc,str_massfrdc,arr_reorder,reorder_massfrdc)
         if(impi/=impi0) allocate(reorder_massfrdc(str_nabund,str_nc))
         allocate(str_massfrdd(str_nabund,ncell))
         n = str_nabund
!        call mpi_scatterv(str_massfrdc,n*counts,n*displs,MPI_REAL8,
         call mpi_scatterv(reorder_massfrdc,n*counts,n*displs,MPI_REAL8,
     &        str_massfrdd,n*ncell,MPI_REAL8,
     &        impi0,MPI_COMM_WORLD,ierr)
      endif
c
c-- gas temperature structure if available
      if(str_ltemp) then
        if(impi==impi0) call getScatter(igeom,str_nc,str_tempdc,
     &        arr_reorder,arr_scatter)
!      if(impi/=impi0) allocate(str_tempdc(str_nc))
!      milagro dd: use arr_scatter replace str_tempdc
       if(impi/=impi0) allocate(arr_scatter(str_nc))
       allocate(str_tempdd(ncell))
!      call mpi_scatterv(str_tempdc,counts,displs,MPI_REAL8,
!      milagro dd: use arr_scatter replace str_tempdc
       call mpi_scatterv(arr_scatter,counts,displs,MPI_REAL8,
     &  str_tempdd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      endif
c-- ye structure if available
      if(str_lye) then
         if(impi==impi0) call getScatter(igeom,str_nc,str_yedc,
     &        arr_reorder,arr_scatter)
!        if(impi/=impi0) allocate(str_yedc(str_nc))
!        milagro dd: use arr_scatter replace str_yedc
         if(impi/=impi0) allocate(arr_scatter(str_nc))
         allocate(str_yedd(ncell))
!        call mpi_scatterv(str_yedc,counts,displs,MPI_REAL8,
         call mpi_scatterv(arr_scatter,counts,displs,MPI_REAL8,
     &        str_yedd,ncell,MPI_REAL8,
     &        impi0,MPI_COMM_WORLD,ierr)
      endif
!}}}

      deallocate(arr_scatter)   !deallocate arr_scatter


!     milagro dd function: get the length of edge
!     2D, the length is sqrt of cells per rank
!     3D, the length is cbrt of cells per rank
      contains
      integer function getedge(nmpi,ndim,igeom)
      implicit none
      integer,intent(in) :: nmpi,ndim(3),igeom
      integer::l,nx,ny,nz,volume
      real*8::tasks,temp

      nx=ndim(1)
      ny=ndim(2)
      nz=ndim(3)

      volume=nx*ny*nz
      tasks=nmpi+0.0
      temp=volume/tasks
      write(6,*) '>> temp',temp
      if(igeom==2) then
         l=ceiling(temp**(1.0/2.0))
         if (l*l*nmpi>volume) then
            l=l-1               ! 2d l round up by 1 due to precision
         endif
      else if(igeom==3) then
         l=ceiling(temp**(1.0/3.0))
         if (l*l*l*nmpi>volume) then
            l=l-1               ! 3d l round up by 1 due to precision
         endif
      else
         l=ceiling(volume/tasks)
         if (l*nmpi>volume) then
            l=l-1               ! 1d l round up by 1 due to precision
         endif
      endif
      getedge=l
      write(6,*) '>> getEdge',l
      end function getedge

      subroutine getCorners(nx,ny,len,nboxes,pos)
      implicit none
      integer, intent(in)::nx,ny,len
      integer, intent(out)::nboxes
      integer,dimension(:),allocatable,intent(out) :: pos
      integer::i,j,cc,p

      nboxes=0
      cc=(nx*ny)/(len*len)
      allocate(pos(cc))
      do j=1,ny,len
         do i=1,nx,len
            p=(j-1)*nx+i
            nboxes=nboxes+1
            pos(nboxes)=p
         enddo                  !end x
      enddo                     !end y
      end subroutine getCorners


      subroutine getReorder(nx,len,nboxes,corners,arr_reorder)
      implicit none
      integer,intent(in)::nx,len,nboxes
      integer,dimension(:),intent(in)::corners
      integer,dimension(:),allocatable,intent(out)::arr_reorder
      integer::arr_size,i,j,k,l,idx,box_offset
      arr_size=len*len*nboxes
      allocate(arr_reorder(arr_size))
      box_offset=0
      l=0                       !l should == arr_size
      do k=1,nboxes
         box_offset=corners(k)-corners(1)
         do j=1,len
            do i=1,len
               l=l+1
               idx=(j-1)*nx+i+box_offset
               arr_reorder(l)=idx
            enddo
         enddo
      enddo
      end subroutine getReorder

      subroutine getCorners3d(igeom,nx,ny,nz,len,nboxes,pos)
      implicit none
      integer, intent(in)::igeom,nx,ny,nz,len
      integer, intent(out)::nboxes
      integer,dimension(:),allocatable,intent(out) :: pos
      integer::i,j,k,cc,p,nodevol

      nboxes=0
      nodevol=getVolume(igeom,len)
      cc=(nx*ny*nz)/nodevol ! number of processors
      write(6,*) '>> corner > nv',nodevol,'igeom',igeom
      allocate(pos(cc))
      do k=1,nz,len
         do j=1,ny,len
            do i=1,nx,len
               p=(k-1)*nx*ny+(j-1)*nx+i
               nboxes=nboxes+1
               pos(nboxes)=p
            enddo               !end x
         enddo                  !end y
      enddo                     !end z
      end subroutine getCorners3d

      subroutine getReorder3d(igeom,nx,ny,nz,len,nboxes,corners,
     &                        arr_reorder)
      implicit none
      integer,intent(in)::igeom,nx,ny,nz,len,nboxes
      integer,dimension(:),intent(in)::corners
      integer,dimension(:),allocatable,intent(out)::arr_reorder
      integer::arr_size,i,j,k,l,m,idx,box_offset,nodevol
      nodevol=getVolume(igeom,len)
      arr_size=nodevol*nboxes ! total cells in domain

      allocate(arr_reorder(arr_size))
      box_offset=0
      l=0                       !l should == arr_size
      do m=1,nboxes
         box_offset=corners(m)-corners(1)
         do k=1,min(len,nz)
            do j=1,len
               do i=1,len
                  l=l+1
c--               write(6,*) '>> reorder m=',m,'k=',k,'j=',j,
c-- &                       'i=',i,'l=',l
                  idx=(k-1)*nx*ny+(j-1)*nx+i+box_offset
                  arr_reorder(l)=idx
               enddo            !end x
            enddo               !end y
         enddo                  !end z
      enddo
      end subroutine getReorder3d

      subroutine getScatter(igeom,arr_size,in_arr,in_order,out_arr)
      implicit none
      integer,intent(in)::igeom,arr_size
      real*8,dimension(:),intent(in)::in_arr
      integer,dimension(:),intent(in)::in_order
      real*8,dimension(:),allocatable,intent(out)::out_arr
      integer::i,idx
      allocate(out_arr(arr_size)) !allocate array to be scattered
      if((igeom==1).or.(igeom==11)) then
         do i=1,arr_size
            out_arr(i)=in_arr(i)
         enddo
      else
         do i=1,arr_size
            idx=in_order(i)
            out_arr(i)=in_arr(idx)
         enddo
      endif
      end subroutine getScatter

!     getMassfrdc(igeom,str_nabund,str_nc,str_massfrdc,
!     &           arr_reorder,reorder_massfrdc)
      subroutine getMassfrdc(igeom,nabund,arr_size,in_arr,
     &                       in_order,out_arr)
      implicit none
      integer,intent(in)::igeom,nabund,arr_size
      real*8,dimension(:,:),intent(in)::in_arr
      integer,dimension(:),intent(in)::in_order
      real*8,dimension(:,:),allocatable,intent(out)::out_arr
      integer::i,idx
      allocate(out_arr(nabund,arr_size)) !allocate array to be scattered
      if((igeom==1).or.(igeom==11)) then
         do i=1,arr_size
            out_arr(:,i)=in_arr(:,i)
         enddo
      else
         do i=1,arr_size
            idx=in_order(i)
            out_arr(:,i)=in_arr(:,idx)
         enddo
      endif
      end subroutine getMassfrdc

      integer function getVolume(igeom,len)
      implicit none
      integer,intent(in)::igeom,len
      integer::v
      if (igeom==2) then
         v=len*len
      else if(igeom==3) then
         v=len*len*len
      endif
      getVolume=v
      end function getVolume
      end subroutine scatter_inputstruct
c
c
c
      subroutine allgather_initialrad
c     -------------------------------
      use gridmod
      use gasmod
      implicit none
************************************************************************
* gather initial gas_eraddens to grd_evolinit
************************************************************************
      call mpi_allgatherv(gas_eraddens,gas_ncell,MPI_REAL8,
     &  grd_evolinit,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      end subroutine allgather_initialrad
c
c
      subroutine allgather_gammacap
c     -----------------------------!{{{
      use gridmod
      use gasmod
      implicit none
************************************************************************
* gather gas_capgam to grd_capgam
************************************************************************
      call mpi_allgatherv(gas_emitex,gas_ncell,MPI_REAL8,
     &  grd_emitex,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_capgam,gas_ncell,MPI_REAL8,
     &  grd_capgam,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)!}}}
      end subroutine allgather_gammacap
c
c
      subroutine allreduce_gammaenergy
c     ------------------------------------!{{{
      use gridmod,nx=>grd_nx,ny=>grd_ny,nz=>grd_nz
      use timingmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
************************************************************************
      real*8 :: snd(2,grd_ncell)
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
      snd = grd_tally
      call mpi_allreduce(snd,grd_tally,2*grd_ncell,MPI_REAL8,MPI_SUM,
     &  MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpimisc, t1-t0)
c!}}}
      end subroutine allreduce_gammaenergy
c
c
      subroutine bcast_nonpermanent
c     -----------------------------!{{{
      use gridmod
      use gasmod
      use groupmod
      use sourcemod
      use totalsmod
      use particlemod
      use timingmod
      use transportmod, only:trn_noampfact
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
************************************************************************
      real*8 :: t0,t1
      real*8 :: sndgas(gas_ncell)
      real*8 :: sndgrd(grd_ncell)
      integer :: nvacant
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
c-- gather
      sndgas = 1d0/gas_temp
      call mpi_allgatherv(sndgas,gas_ncell,MPI_REAL8,
     &  grd_tempinv,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_fcoef,gas_ncell,MPI_REAL8,
     &  grd_fcoef,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_capgrey,gas_ncell,MPI_REAL8,
     &  grd_capgrey,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_emit,gas_ncell,MPI_REAL8,
     &  grd_emit,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_emitex,gas_ncell,MPI_REAL8,
     &  grd_emitex,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgatherv(gas_sig,gas_ncell,MPI_REAL8,
     &  grd_sig,counts,displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(gas_cap,grp_ng*gas_ncell,MPI_REAL,
     &  grd_cap,grp_ng*counts,grp_ng*displs,MPI_REAL,
     &  MPI_COMM_WORLD,ierr)
c
c-- broadcast
      call mpi_bcast(tot_esurf,1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- allreduce
      if(.not.trn_noampfact) then
       sndgrd = grd_eamp
       call mpi_allreduce(sndgrd,grd_eamp,grd_ncell,MPI_REAL8,MPI_SUM,
     &   MPI_COMM_WORLD,ierr)
      endif
c
c-- allgather
      nvacant = count(prt_isvacant) !count returns the same type as prt_isvacant
      call mpi_allgather(nvacant,1,MPI_INTEGER,
     &  src_nvacantall,1,MPI_INTEGER8,
     &  MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpibcast, t1-t0)
c!}}}
      end subroutine bcast_nonpermanent
c
c
      subroutine allgather_leakage
c     ------------------------------------------!{{{
      use gridmod
      use timingmod
      implicit none
************************************************************************
      integer :: n
      real*8,allocatable :: snd(:,:)
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
      allocate(snd(9,grd_ndd))
      snd = grd_opaclump(:,grd_idd1:grd_idd1+grd_ndd-1)
      call mpi_allgatherv(snd,9*grd_ndd,MPI_REAL8,
     &  grd_opaclump,9*counts,9*displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
      n = grd_nep
      allocate(snd(n,grd_ndd))
      snd = grd_emitprob(:,grd_idd1:grd_idd1+grd_ndd-1)
      call mpi_allgatherv(snd,n*grd_ndd,MPI_REAL8,
     &  grd_emitprob,n*counts,n*displs,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
      t1 = t_time()
      call timereg(t_mpimisc, t1-t0)
c!}}}
      end subroutine allgather_leakage
c
c
c
      subroutine reduce_gridtally
c     -----------------------!{{{
      use gridmod,nx=>grd_nx,ny=>grd_ny,nz=>grd_nz
      use gasmod
      use timingmod
      use fluxmod
      use sourcemod
      implicit none
************************************************************************
* Reduce the results from particle_advance that are needed for the
* temperature correction.
************************************************************************
      integer :: n
      integer :: isnd2f(flx_nmu,flx_nom)
      real*8 :: snd2f(flx_nmu,flx_nom)
      integer,allocatable :: isnd(:)
      real*8,allocatable :: snd(:)
      real*8,allocatable :: snd2(:,:)
      real*8 :: help
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
c-- flux dim==2
      n = flx_nmu*flx_nom
      isnd2f = flx_gamlumnum
      call mpi_reduce(isnd2f,flx_gamlumnum,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd2f = flx_gamluminos
      call mpi_reduce(snd2f,flx_gamluminos,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd2f = flx_gamlumdev
      call mpi_reduce(snd2f,flx_gamlumdev,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd2f = flx_gamlumtime
      call mpi_reduce(snd2f,flx_gamlumtime,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- dim==3
      allocate(isnd(grd_ncell))
      n = grd_ncell
      isnd = grd_numcensimc
      call mpi_reduce(isnd,grd_numcensimc,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      isnd = grd_numcensddmc
      call mpi_reduce(isnd,grd_numcensddmc,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      isnd = grd_methodswap
      call mpi_reduce(isnd,grd_methodswap,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(isnd)
c
      allocate(snd2(2,grd_ncell))
      snd2 = grd_tally
      call mpi_reduce(snd2,grd_tally,2*n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(snd2)
c
c-- bcast
      call mpi_bcast(src_nflux,1,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- scatter
      allocate(snd(grd_ncell))
      snd = grd_tally(1,:)
      call mpi_scatterv(snd,counts,displs,MPI_REAL8,
     &  gas_edep,gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd = grd_tally(2,:)
      call mpi_scatterv(snd,counts,displs,MPI_REAL8,
     &  gas_eraddens,gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(snd)
c
c-- timing statistics
      help = t_pckt_stat(1)
      call mpi_reduce(help,t_pckt_stat(1),1,MPI_REAL8,MPI_MIN,
     &  impi0,MPI_COMM_WORLD,ierr)
      help = t_pckt_stat(2)/nmpi
      call mpi_reduce(help,t_pckt_stat(2),1,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      help = t_pckt_stat(3)
      call mpi_reduce(help,t_pckt_stat(3),1,MPI_REAL8,MPI_MAX,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpireduc, t1-t0)
c!}}}
      end subroutine reduce_gridtally
c
c
c
      subroutine reduce_fluxtally
c     ---------------------------
      use fluxmod
      implicit none
************************************************************************
* Reduce flux arrays for output
************************************************************************
      integer :: n
      integer :: isnd3f(flx_ng,flx_nmu,flx_nom)
      real*8 :: snd3f(flx_ng,flx_nmu,flx_nom)
c
c-- flux dim==3
      n = flx_ng*flx_nmu*flx_nom
      isnd3f = flx_lumnum
      call mpi_reduce(isnd3f,flx_lumnum,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_luminos
      call mpi_reduce(snd3f,flx_luminos,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_lumdev
      call mpi_reduce(snd3f,flx_lumdev,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      end subroutine reduce_fluxtally
c
c
c
      subroutine reduce_gastemp
c     -------------------------------------!{{{
      use gridmod
      use totalsmod
      use gasmod
      use timingmod
      implicit none
************************************************************************
* for output
************************************************************************
      integer :: n
      real*8,allocatable :: sndvec(:),rcvvec(:)
      real*8 :: sndgas(gas_ncell)
      real*8 :: t0,t1
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      t0 = t_time()
c
      sndgas = 1d0/gas_temp
      call mpi_gatherv(sndgas,gas_ncell,MPI_REAL8,
     &   grd_tempinv,counts,displs,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c
c-- dim==0
      n = 12
      allocate(sndvec(n))
      allocate(rcvvec(n))
      sndvec = [
     &  tot_erad,tot_eout,tot_eext,tot_evelo,tot_emat,tot_eext0,
     &  tot_sthermal,tot_smanufac,tot_sdecaygamma,tot_sdecaybeta,
     &  tot_sfluxgamma,tot_sflux]
c
      call mpi_reduce(sndvec,rcvvec,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      if(lmpi0) then
         tot_erad = rcvvec(1)
         tot_eout = rcvvec(2)
         tot_eext = rcvvec(3)
         tot_evelo = rcvvec(4)
         tot_emat = rcvvec(5)
         tot_eext0 = rcvvec(6)
         tot_sthermal = rcvvec(7)
         tot_smanufac = rcvvec(8)
         tot_sdecaygamma = rcvvec(9)
         tot_sdecaybeta = rcvvec(10)
         tot_sfluxgamma = rcvvec(11)
         tot_sflux = rcvvec(12)
      else
c-- zero out cumulative values on all other ranks to avoid double counting.
         tot_eout = 0d0
         tot_eext = 0d0
         tot_evelo = 0d0
      endif
      deallocate(sndvec)
      deallocate(rcvvec)
c
      t1 = t_time()
      call timereg(t_mpimisc, t1-t0)
c!}}}
      end subroutine reduce_gastemp
c
c
c
      subroutine mpimod_dealloc
      deallocate(counts,displs)
      end subroutine mpimod_dealloc
c
      end module mpimod
c vim: fdm=marker
