!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
subroutine particle_advance_gamgrey(impi,nmpi)

!$ use omp_lib
  use miscmod
  use randommod
  use sourcemod
  use particlemod
  use transportmod
  use gridmod
  use physconstmod
  use inputparmod
  use timingmod
  use timestepmod
  use totalsmod
  use fluxmod
  use gasmod
  use milagromod
  implicit none
  integer,intent(in) :: impi,nmpi
!###### ############################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or be censused for continued
  !propagation in the next time step.  Currently DDMC and IMC particle events
  !are being handled in separate subroutines but this may be changed to reduce
  !total subroutine calls in program.
  !##################################################
  integer :: istat(MPI_STATUS_SIZE)
  integer :: ierr
  !integer :: nhere, nemit, ndmy -- chenx
  real*8 :: r1, edep, help
  integer :: i,j,k,l, ii, iimpi !,gas_idx -- chenx
  integer :: imu, iom, icold
  integer :: r,newrank,temprb,temprb2,nbf,nbf2 ! chenx
  integer,pointer :: ic
  integer,pointer :: ix, iy, iz
  real*8,pointer :: x,y,z,mu,om,e,e0
  real*8 :: eta, xi
  real*8 :: t0,t1  !timing
  real*8 :: labfact, cmffact, mu1, mu2, gm
  real*8 :: etot,pwr,sdtot !subdomain tot energy
  real*8 :: om0, mu0, x0, y0, z0
!
  integer :: nvol(grd_ncell)
!
  type(rnd_t) :: rndstate
  integer,save :: iomp=0
!
  real*8,parameter :: basefrac=.1d0
  real*8 :: base,edone,einv,invn,en
  integer :: n, ndone, mpart, npart, ipart
  integer*8 :: nstot,nsavail,nsbase
!
  integer,allocatable :: ipospart(:,:) !(3,npart)
!
  type(packet),target :: ptcl
  type(packet2),target :: ptcl2
! chenx
  integer,parameter::highint=1000000000
  integer::rnkhigh
!
!-- start clock
  t0 = t_time()

  rnkhigh = 0 !chenx
  grd_tally = 0d0
  flx_gamlumtime = 0d0
  flx_gamluminos = 0d0
  flx_gamlumdev = 0d0
  flx_gamlumnum = 0

!-- initializing volume numbers
  nvol = 0

!-- shortcut
  pwr = in_srcepwr

!-- total energy (converted by pwr)
  etot = sum(grd_emitex**pwr)
  sdtot = sum(gas_emitex**pwr)
  if(etot/=etot) stop 'particle_advance_gamgrey: etot nan'
!-- total particle number
  !nstot = nmpi*int(src_ns,8)
  nstot = nint(nmpi*src_ns*(sdtot/etot),8)
  !write(6,*) '@impi,',impi,'@nstot2',nstot


!-- base (flat,constant) particle number per cell over ALL RANKS
  !n = count(grd_emitex>0d0)  !number of cells that emit
  n = count(gas_emitex>0d0)  !number of cells that emit
  base = dble(nstot)/n  !uniform distribution
  base = basefrac*base

!-- number of particles available for proportional distribution
  nsbase = int(n*base,8)  !total number of base particles
  nsavail = nstot - nsbase

  !write(6,*) '@rank',impi,'nsbase',nsbase,'nsavail',nsavail
  !mpart = nint(10 + 1.1d0*src_ns) !big enough
  mpart = nint(10 + 1.1d0*nstot) !big enough
  allocate(ipospart(3,mpart))

  iimpi = 0
  npart = 0
!-- total particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/nstot
  !einv = 1d0/etot
  einv = 1d0/sdtot
  !do k=1,grd_nz
  !do j=1,grd_ny
  !do i=1,grd_nx
  do l=1,size(gas_idcell)
     !l = grd_icell(i,j,k)
     !en = grd_emitex(l)**pwr
     en = gas_emitex(l)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = int(en*nsavail*einv + base)  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     nvol(gas_idcell(l)) = n
     edone = edone + en
     ndone = ndone + n
!-- particle count on this rank
     !call sourcenumbers_roundrobin(iimpi,grd_emitex(l), &
     !     0d0,nvol(l),nemit,nhere,ndmy)
     call reverseMapping(gas_idcell(l),i,j,k)
     do ii=1,nvol(gas_idcell(l))!nhere
        npart = npart + 1
        if(npart>mpart) stop 'particle_advance_gamgrey: npart>npartmax'
        ipospart(:,npart) = [i,j,k]
     enddo !ii
  enddo
  !enddo
  !enddo
  !enddo
  call mpi_reduce(npart,maxpars,1,MPI_INTEGER,MPI_SUM,&
       &  0,MPI_COMM_WORLD,ierr) !chenx
  call recvMaxParBuff(impi,db2r)
  call recvChdParComplete(impi)
  if(tn%parent>=0) then
     call recvParentFinish(impi)
  endif

!$omp parallel &
!$omp shared(nvol) &
!$omp private(ptcl,ptcl2,x0,y0,z0,mu0,om0,cmffact,gm,mu1,mu2,eta,xi,labfact,iom,imu, &
!$omp    rndstate,edep,ierr, iomp, &
!$omp    x,y,z,mu,om,e,e0,ix,iy,iz,ic,icold,r1, &
!$omp    i,j,k, newrank, r, sndbuff, sndbuff2, sndidx) &
!$omp reduction(+:grd_tally,flx_gamluminos,flx_gamlumnum, &
!$omp    flx_gamlumdev,flx_gamlumtime,ipcmplt)

!-- thread id
!$ iomp = omp_get_thread_num()
  rndstate = rnd_states(iomp+1)
  sndidx=0                                                ! chenx
  sndbuff=packet(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)     ! chenx
  sndbuff2=packet2(0.0,0.0,0.0,0.0,' ',0,0,0,0,0,0,0,0,0) ! chenx
  ipcmplt = 0                                             ! chenx
  nbf=0
!
!-- primary particle properties
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
!-- secondary particle properties
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic

!-- default values in bounds
  y = grd_yarr(1)
  z = grd_zarr(1)

!=========================================
! omp loop for local gamma
! 1)initialize
! 2)move while 'live'
! 3)tally
!=========================================
!$omp do schedule(static,1) !round-robin
!  do ipart=1,npart
  do ipart=npart,1,-1
     i = ipospart(1,ipart)
     j = ipospart(2,ipart)
     k = ipospart(3,ipart)
!-- adopt position (get rid of this copy after merged with wlT branch)
     ix = i
     iy = j
     iz = k
     ic = grd_icell(ix,iy,iz)

!-- x position
     select case(in_igeom)
     case(1,11)
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
!-- must be inside cell
        x = min(x,grd_xarr(i+1))
        x = max(x,grd_xarr(i))
     case(2)
        call rnd_r(r1,rndstate)
        x = sqrt(r1*grd_xarr(i+1)**2 + (1d0-r1)*grd_xarr(i)**2)
!-- must be inside cell
        x = min(x,grd_xarr(i+1))
        x = max(x,grd_xarr(i))
     case(3)
        call rnd_r(r1,rndstate)
        x = r1*grd_xarr(i+1) + (1d0-r1) * grd_xarr(i)
     endselect

!-- y,z position
     call rnd_r(r1,rndstate)
     y = r1*grd_yarr(j+1) + (1d0-r1) * grd_yarr(j)
     call rnd_r(r1,rndstate)
     z = r1*grd_zarr(k+1) + (1d0-r1) * grd_zarr(k)

!-- direction cosine (comoving)
     call rnd_r(r1,rndstate)
     mu0 = 1d0-2d0*r1
     call rnd_r(r1,rndstate)
     om0 = pc_pi2*r1

!-- transform direction
     if(.not.grd_isvelocity) then
         mu = mu0
         om = om0
     else
        select case(in_igeom)!{{{
        case(1,11)
           x0 = x
           cmffact = 1d0+mu0*x0/pc_c !-- 1+dir*v/c
           mu = (mu0+x0/pc_c)/cmffact
           om = om0
        case(2)
           x0 = x
           y0 = y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           gm = 1d0/sqrt(1d0-(x**2+y**2)/pc_c**2)
!-- om
           om = atan2(sqrt(1d0-mu0**2)*sin(om0) , &
                sqrt(1d0-mu0**2)*cos(om0)+(gm*x/pc_c) * &
                (1d0+gm*(cmffact-1d0)/(gm+1d0)))
           if(om<0d0) om = om+pc_pi2
!-- mu
           mu = (mu0+(gm*y/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
                (gm*cmffact)
        case(3)
           x0 = x
           y0 = y
           z0 = z
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           mu = (mu0+z0/pc_c)/cmffact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- om
           om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(om<0d0) om = om+pc_pi2
        endselect!}}}
     endif

!-- velocity components in cartesian basis
     if(grd_igeom==1) then
!-- spherical projections
        eta = sqrt(1d0-mu**2)*cos(om)
        xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
        ptcl2%mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
        ptcl2%muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
        ptcl2%muz = mu*y-eta*sqrt(1d0-y**2)
     endif
!-- update invariant direction quantities
     if(grd_igeom==2) then
        ptcl2%mux = x*sin(om)/sin(z+om)  !-- intercept
        ptcl2%muy = x*sin(z)/sin(z+om)  !-- distance to intercept
        ptcl2%muz = pc_pi-(z+om)  !-- direction angle
        if(ptcl2%muz<0d0) ptcl2%muz = ptcl2%muz+pc_pi2
        if(ptcl2%muz<0d0) ptcl2%muz = ptcl2%muz+pc_pi2
     endif
!
!-- emission energy per particle
     e = grd_emitex(ic)/nvol(ic)
     if(grd_isvelocity) e = e*cmffact
     e0 = e

!-----------------------------------------------------------------------
!-- Advancing particle until census, absorption, or escape from domain
     if(rnkhigh==0) then
        rnkhigh=(impi+1)*(highint/(10**ceiling(log10(nmpi+0.0))))
     endif !chenx
     ptcl2%ipart = ipart + rnkhigh
     ptcl2%istep = 0
     ptcl2%idist = 0

     ptcl2%stat = 'live'
!
     do while (ptcl2%stat=='live')
        ptcl2%istep = ptcl2%istep + 1
        icold = ptcl2%ic
        call transport_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
        !if(ptcl2%stat=='buff') then
        !   write(6,*) '@> p',ptcl2%ipart,'ic',ptcl2%ic,'icold',icold,'stat',ptcl2%stat
        !endif
!-- tally
        grd_tally(1,icold) = grd_tally(1,icold) + edep

!-- Russian roulette for termination of exhausted particles
        if(ptcl%e<1d-6*e0 .and. ptcl2%stat=='live' .and. grd_capgam(ptcl2%ic)>0d0) then
           call rnd_r(r1,rndstate)!{{{
           if(r1<0.5d0) then
!-- transformation factor
              if(grd_isvelocity) then
                 select case(grd_igeom)
                 case(1,11)
                    labfact = 1.0d0 - ptcl%mu*x/pc_c
                 case(2)
                    labfact = 1d0-(ptcl%mu*ptcl%y+sqrt(1d0-ptcl%mu**2) * &
                         cos(ptcl%om)*ptcl%x)/pc_c
                 case(3)
                    labfact = 1d0-(ptcl%mu*ptcl%z+sqrt(1d0-ptcl%mu**2) * &
                         (cos(ptcl%om)*ptcl%x+sin(ptcl%om)*ptcl%y))/pc_c
                 endselect
              else
                 labfact = 1d0
              endif
!
              ptcl2%stat = 'dead'
              grd_tally(1,ptcl2%ic) = grd_tally(1,ptcl2%ic) + ptcl%e*labfact
           else
              ptcl%e = 2d0*ptcl%e
              ptcl%e0 = 2d0*ptcl%e0
           endif!}}}
        endif

!-- verify position
        if(ptcl2%stat=='live') then
           if(ptcl%x>grd_xarr(ptcl2%ix+1) .or. ptcl%x<grd_xarr(ptcl2%ix) .or. ptcl%x/=ptcl%x) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv_ggrey: x not in cell', &
                   ptcl2%ix,ptcl%x,grd_xarr(ptcl2%ix),grd_xarr(ptcl2%ix+1)
           endif
           if(ptcl%y>grd_yarr(ptcl2%iy+1) .or. ptcl%y<grd_yarr(ptcl2%iy) .or. ptcl%y/=ptcl%y) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv_ggrey: y not in cell', &
                   ptcl2%iy,ptcl%y,grd_yarr(ptcl2%iy),grd_yarr(ptcl2%iy+1)
           endif
           if(ptcl%z>grd_zarr(ptcl2%iz+1) .or. ptcl%z<grd_zarr(ptcl2%iz) .or. ptcl%z/=ptcl%z) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv_ggrey: z not in cell', &
                   ptcl2%iz,ptcl%z,grd_zarr(ptcl2%iz),grd_zarr(ptcl2%iz+1)
           endif
        endif
!-- check for errors
        if(ierr/=0 .or. ptcl2%istep>1000) then
           write(0,*) 'pagg: ierr,ipart,istep,idist:',ierr,ptcl2%ipart,ptcl2%istep,ptcl2%idist
           write(0,*) 'dist:',ptcl2%dist
           write(0,*) 't:',ptcl%t
           write(0,*) 'ix,iy,iz,ic,ig:',ptcl2%ix,ptcl2%iy,ptcl2%iz,ptcl2%ic,ptcl2%ig
           write(0,*) 'x,y,z:',ptcl%x,ptcl%y,ptcl%z
           write(0,*) 'mu,om:',ptcl%mu,ptcl%om
           write(0,*) 'mux,muy,muz:',ptcl2%mux,ptcl2%muy,ptcl2%muz
           write(0,*)
           if(ierr>0) then
              if(trn_errorfatal) stop 'particle_advance_gg: fatal transport error'
              ptcl2%stat = 'dead'
              exit
           endif
        endif
     enddo
!-- write to sndbuff chenx ---------------------------------------------!
     if(ptcl2%stat=='buff') then
        nbf=nbf+1
        newrank=getRankId(myGhosts,ptcl2%ic)
        r=nbrs(newrank+1)
        !if(impi==0) then ! don't print
        !   write(6,'(A12,I2,A2,I2,A2,I5,A8,I3,A2)') 'buff@rank(',impi,&
        !     &') ',sndidx(r),'@  thrd(',iomp,')'
        !endif
        !write(6,*) 'wrt',impi,newrank,ptcl2%ipart,nbf
        if(sndidx(r)<BUFFSIZE) then
           call writeBuffer(impi,iomp,sndbuff,sndbuff2,&
                &sndidx,r,ptcl,ptcl2)
        else
           call sendBuffer(impi,iomp,sndbuff,sndbuff2, &
                &sndidx,r,newrank)
           call writeBuffer(impi,iomp,sndbuff,sndbuff2,&
                &sndidx,r,ptcl,ptcl2)
        endif
     endif
!-- end write to sndbuff -----------------------------------------------!
!
!-- outbound luminosity tally
     if(ptcl2%stat=='flux') then
        ipcmplt = ipcmplt + 1
!-- lab frame flux group, polar, azimuthal bin
        imu = binsrch(mu,flx_mu,flx_nmu+1,.false.)
        iom = binsrch(om,flx_om,flx_nom+1,.false.)
!-- observer corrected time
        help=tsp_t+.5d0*tsp_dt
        select case(grd_igeom)
        case(1,11)
           labfact = mu*x/pc_c
        case(2)
           labfact = (mu*y+sqrt(1d0-mu**2) * &
                cos(om)*x)/pc_c
        case(3)
           labfact = (mu*z+sqrt(1d0-mu**2) * &
                (cos(om)*x+sin(om)*y))/pc_c
        endselect
        if(grd_isvelocity) labfact=labfact*tsp_t
        help=help-labfact
        !-- tally outbound luminosity
        flx_gamlumtime(imu,iom) = flx_gamlumtime(imu,iom)+help
        flx_gamluminos(imu,iom) = flx_gamluminos(imu,iom)+e
        flx_gamlumdev(imu,iom) = flx_gamlumdev(imu,iom)+e**2
        flx_gamlumnum(imu,iom) = flx_gamlumnum(imu,iom)+1
     endif
     if(ptcl2%stat=='dead') then
        ipcmplt = ipcmplt + 1
     endif
     if(ipart==1)then ! last local particle
        do r=1,size(db2r)
           if(sndidx(r)>0) then
              call sendBuffer(impi,iomp,sndbuff,sndbuff2, &
                   &sndidx,r,db2r(r))
           endif
        enddo
     endif
  enddo !ipart
!$omp end do
!
!-- save state
  rnd_states(iomp+1) = rndstate
!$omp end parallel
  !write(6,*) 'after local',impi,ipcmplt,maxpars,'wrtbf',nbf,globalFinish
! milagro transport particles in buffer chenx
  i=0
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
!-- secondary particle properties
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  nbf=0
  nbf2=0
!  if(tn%lchild>0)then
!     call mpi_irecv(pcomplete(1),1,MPI_INTEGER,tn%lchild,rpctag,&
!          & MPI_COMM_WORLD,reqpc(1),ierr)
!  endif
!  if(tn%rchild>0)then
!     call mpi_irecv(pcomplete(2),1,MPI_INTEGER,tn%rchild,rpctag,&
!          & MPI_COMM_WORLD,reqpc(2),ierr)
!  endif
  do while((i<60000).and.(.not.globalFinish))
  !do while(.not.globalFinish)
     i = i + 1
     call tallyPcomplete(impi,tn,ipcmplt)
     do j=1, size(db2r)
        temprb=reqrb(j)
        temprb2=reqrb2(j)
        call mpi_test(reqrb(j),rbflag(j),istat,ierr)
        call mpi_test(reqrb2(j),rb2flag(j),istat,ierr)
        lhas= rbflag(j).and. rb2flag(j)
        if(lhas) then
           do k=1,BUFFSIZE
              ptcl=rcvbuff(k,j)
              ptcl2=rcvbuff2(k,j)
              if(ptcl2%ipart/=0) then
                 ptcl2%stat='live'
                 nbf=nbf+1
                 !write(6,*) 'rcv',impi,'rb',temprb,temprb2,'ipart',ptcl2%ipart
              endif
!----- transport buffer
              do while (ptcl2%stat=='live')
                 ptcl2%istep = ptcl2%istep + 1
                 icold = ptcl2%ic
                 call transport_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
!-- tally
                 grd_tally(1,icold) = grd_tally(1,icold) + edep

!-- Russian roulette for termination of exhausted particles
                 if(ptcl%e<1d-6*e0 .and. ptcl2%stat=='live' .and. grd_capgam(ptcl2%ic)>0d0) then
                    call rnd_r(r1,rndstate)!{{{
                    if(r1<0.5d0) then
!-- transformation factor
                       if(grd_isvelocity) then
                          select case(grd_igeom)
                          case(1,11)
                             labfact = 1.0d0 - ptcl%mu*x/pc_c
                          case(2)
                             labfact = 1d0-(ptcl%mu*ptcl%y+sqrt(1d0-ptcl%mu**2) * &
                                  cos(ptcl%om)*ptcl%x)/pc_c
                          case(3)
                             labfact = 1d0-(ptcl%mu*ptcl%z+sqrt(1d0-ptcl%mu**2) * &
                                  (cos(ptcl%om)*ptcl%x+sin(ptcl%om)*ptcl%y))/pc_c
                          endselect
                       else
                          labfact = 1d0
                       endif
!
                       ptcl2%stat = 'dead'
                       grd_tally(1,ptcl2%ic) = grd_tally(1,ptcl2%ic) + ptcl%e*labfact
                    else
                       ptcl%e = 2d0*ptcl%e
                       ptcl%e0 = 2d0*ptcl%e0
                    endif!}}}
                 endif

!-- verify position
                 if(ptcl2%stat=='live') then
                    if(ptcl%x>grd_xarr(ptcl2%ix+1) .or. ptcl%x<grd_xarr(ptcl2%ix) .or. ptcl%x/=ptcl%x) then
                       if(ierr==0) ierr = -99
                       write(0,*) 'prt_adv_ggrey: x not in cell', &
                            ptcl2%ix,ptcl%x,grd_xarr(ptcl2%ix),grd_xarr(ptcl2%ix+1)
                    endif
                    if(ptcl%y>grd_yarr(ptcl2%iy+1) .or. ptcl%y<grd_yarr(ptcl2%iy) .or. ptcl%y/=ptcl%y) then
                       if(ierr==0) ierr = -99
                       write(0,*) 'prt_adv_ggrey: y not in cell', &
                            ptcl2%iy,ptcl%y,grd_yarr(ptcl2%iy),grd_yarr(ptcl2%iy+1)
                    endif
                    if(ptcl%z>grd_zarr(ptcl2%iz+1) .or. ptcl%z<grd_zarr(ptcl2%iz) .or. ptcl%z/=ptcl%z) then
                       if(ierr==0) ierr = -99
                       write(0,*) 'prt_adv_ggrey: z not in cell', &
                            ptcl2%iz,ptcl%z,grd_zarr(ptcl2%iz),grd_zarr(ptcl2%iz+1)
                    endif
                 endif
!-- check for errors
                 if(ierr/=0 .or. ptcl2%istep>1000) then
                    write(0,*) 'pagg: ierr,ipart,istep,idist:',ierr,ptcl2%ipart,ptcl2%istep,ptcl2%idist
                    write(0,*) 'dist:',ptcl2%dist
                    write(0,*) 't:',ptcl%t
                    write(0,*) 'ix,iy,iz,ic,ig:',ptcl2%ix,ptcl2%iy,ptcl2%iz,ptcl2%ic,ptcl2%ig
                    write(0,*) 'x,y,z:',ptcl%x,ptcl%y,ptcl%z
                    write(0,*) 'mu,om:',ptcl%mu,ptcl%om
                    write(0,*) 'mux,muy,muz:',ptcl2%mux,ptcl2%muy,ptcl2%muz
                    write(0,*)
                    if(ierr>0) then
                       if(trn_errorfatal) stop 'particle_advance_gg: fatal transport error'
                       ptcl2%stat = 'dead'
                       exit
                    endif
                 endif
                 !write(6,*) 'rcv buff@rnk',impi,'@',nbf,ptcl2%ipart,ptcl2%stat
              enddo
!-- write to sndbuff chenx ---------------------------------------------!
              if(ptcl2%stat=='buff') then
                 nbf2=nbf2+1
                 newrank=getRankId(myGhosts,ptcl2%ic)
                 r=nbrs(newrank+1)
                 !write(6,*) 'wrt buff2@rnk',impi,'@',nbf2,ptcl2%ipart,ptcl2%stat,sndidx
                 if(sndidx(r)<BUFFSIZE) then
                    call writeBuffer(impi,iomp,sndbuff,sndbuff2,&
                         &sndidx,r,ptcl,ptcl2)
                 else
                    call sendBuffer(impi,iomp,sndbuff,sndbuff2, &
                         &sndidx,r,newrank)
                    call writeBuffer(impi,iomp,sndbuff,sndbuff2,&
                         &sndidx,r,ptcl,ptcl2)
                 endif
              endif
!-- end write to sndbuff -----------------------------------------------!
!
!-- outbound luminosity tally
              if(ptcl2%stat=='flux') then
                 ipcmplt = ipcmplt + 1
!-- lab frame flux group, polar, azimuthal bin
                 imu = binsrch(mu,flx_mu,flx_nmu+1,.false.)
                 iom = binsrch(om,flx_om,flx_nom+1,.false.)
!-- observer corrected time
                 help=tsp_t+.5d0*tsp_dt
                 select case(grd_igeom)
                 case(1,11)
                    labfact = mu*x/pc_c
                 case(2)
                    labfact = (mu*y+sqrt(1d0-mu**2) * &
                         cos(om)*x)/pc_c
                 case(3)
                    labfact = (mu*z+sqrt(1d0-mu**2) * &
                         (cos(om)*x+sin(om)*y))/pc_c
                 endselect
                 if(grd_isvelocity) labfact=labfact*tsp_t
                 help=help-labfact
        !-- tally outbound luminosity
                 flx_gamlumtime(imu,iom) = flx_gamlumtime(imu,iom)+help
                 flx_gamluminos(imu,iom) = flx_gamluminos(imu,iom)+e
                 flx_gamlumdev(imu,iom) = flx_gamlumdev(imu,iom)+e**2
                 flx_gamlumnum(imu,iom) = flx_gamlumnum(imu,iom)+1
              endif
              if(ptcl2%stat=='dead') then
                 ipcmplt = ipcmplt + 1
              endif
!----- end transport buffer
           enddo ! receive and transport incoming buffer
           !write(6,*) 'rcvd',impi,nbf,'2ndbuff',nbf2
           call mpi_irecv(rcvbuff(1,j),BUFFSIZE,pckttyp,&
                &db2r(j),rbftag,MPI_COMM_WORLD,reqrb(j),ierr)
           call mpi_irecv(rcvbuff2(1,j),BUFFSIZE,pckttyp2,&
                &db2r(j),rbf2tag,MPI_COMM_WORLD,reqrb2(j),ierr)
        endif
        if(sndidx(j)>0) then
           !write(6,*) '> send incoming > last par',impi,sndidx
           call sendBuffer(impi,iomp,sndbuff,sndbuff2, &
                &sndidx,j,db2r(j))
        endif
     enddo
     !write(6,*) '@rnk', impi,ipcmplt,maxpars,'tr',tn
     if(impi==0) then ! master rank
        if(ipcmplt==maxpars) then
           globalFinish=.true.
           !write(6,*) '@0 global finish',globalFinish
           if(tn%lchild>0) then
              call mpi_send(globalFinish,1,MPI_LOGICAL,tn%lchild,&
                   &rftag,MPI_COMM_WORLD,ierr)
           endif
           if(tn%rchild>0) then
              call mpi_send(globalFinish,1,MPI_LOGICAL,tn%rchild,&
                   &rftag,MPI_COMM_WORLD,ierr)
           endif
        endif
     else             ! worker rank
        if(ipcmplt>0) then
           call mpi_send(ipcmplt,1,MPI_integer,tn%parent,&
                &rpctag,MPI_COMM_WORLD,ierr)
           ipcmplt=0
        endif
        call mpi_test(reqgf,rgfflag,istat,ierr)
        if(rgfflag) then
           !write(6,'(A5,I2,A20)') 'rank(',impi,') recv global finish'
           if(tn%lchild>0) then
              call mpi_send(globalFinish,1,MPI_LOGICAL,tn%lchild,&
                   &rftag,MPI_COMM_WORLD,ierr)
           endif
           if(tn%rchild>0) then
              call mpi_send(globalFinish,1,MPI_LOGICAL,tn%rchild,&
                   &rftag,MPI_COMM_WORLD,ierr)
           endif
        endif
     endif            ! end check globalfinish
     if(globalFinish) then
        do j=1,size(db2r)
           if(reqrb(j)>0) then
              call mpi_cancel(reqrb(j),ierr)
           endif
           if(reqrb2(j)>0) then
              call mpi_cancel(reqrb2(j),ierr)
           endif
        enddo
        do j = 1,2
           if(reqpc(j)>0) then
              call mpi_cancel(reqpc(j),ierr)
           endif
        enddo
     endif
  enddo
  write(6,*) 'finish',impi,i,globalFinish
  tot_sfluxgamma = -sum(flx_gamluminos)

!-- convert to flux per second
  help = 1d0/tsp_dt
  flx_gamluminos = flx_gamluminos*help
  flx_gamlumdev = flx_gamlumdev*help**2

  deallocate(ipospart)

  t1 = t_time()
  call timereg(t_pcktgam, t1-t0)

contains
  subroutine reverseMapping(idx,ix,iy,iz)
    implicit none
    integer,intent(in)::idx
    integer,intent(out)::ix,iy,iz
    iz=ceiling(idx/real(grd_nx*grd_ny))
    iy=ceiling((idx-(iz-1)*grd_nx*grd_ny)/real(grd_nx))
    ix=idx-(iz-1)*grd_nx*grd_ny-(iy-1)*grd_nx
  end subroutine reverseMapping

end subroutine particle_advance_gamgrey
! vim: fdm=marker
