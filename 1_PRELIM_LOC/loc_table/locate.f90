real trmin(500),tobmn(500),almn(500)
integer ngdmn(500)
integer ist222(500),ips222(500)
real tob222(500),dtold(500)
real tall(20),hall(20),aall(20)
integer kmin(2000),nall(2000),kodes(10,20000)
real dkode(20000)
real zstart(10),zst_it(10,10), w_qual(10)

character*8 ar,md,line

real fst(9000),tst(9000)

common/stations/ xst(9000),yst(9000),zst(9000)
common/krat/nkrat,istkr(500),tobkr(500),ipskr(500),qualkr(500),trfkr(500),ngood(500),alkr(500),diskr(500)
common/stat_level/stat_level
common/loc_param/wgs,res_loc1,res_loc2,dist_limit,n_pwr_dist,ncyc_av,w_P_S_diff

common/general/k_re1_syn2,VPSX_key,koe,kref,key_ft1_xy2,key_true1,key_flat1
common/loc_table/krat_min,dist_max,wgs0,dist_limit0,n_pwr_dist0,ncyc_av0,bad_max,res_1_km,sss_max,ifreq0, &
	niter_loc,dx_it(10),dy_it(10),dz_it(10),res_it1(10),res_it2(10),wps_it(10)
common/center/fi0,tet0
common/ref_table/stat_lev,dmin,depmax,distmax,nlay,zcrce_lev(20),dzsrce_lev(20),zztmax


one=1.d0
pi=asin(one)*2.d0
per=pi/180.d0
iprint=0

! Read codes of the model

koe=0
open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
221 close(1)

call read_param(ar,md)
call read_z_lim(ar,md)
call read_vref(ar,md)
call read_topo(ar)
call prepare_ref(ar,md)


dist_limit=dist_limit0; n_pwr_dist=n_pwr_dist0; ncyc_av=ncyc_av0; stat_level=stat_lev; wgs=wgs0
ifreq=ifreq0


write(*,*)' area=',ar,' model=',md,' koe=',koe,' k_re1_syn2=',k_re1_syn2
i=system('mkdir -p ../../../TMP_files/rays')


! Read the coordinates of the stations
open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
open(12,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
nst=0
33	read(1,*,end=44)fi,tet,zstat
	call SFDEC(fi,tet,0.,X,Y,Z,fi0,tet0)
	nst=nst+1
	fst(nst)=fi
	tst(nst)=tet
	xst(nst)=x
	yst(nst)=y
	zst(nst)=zstat
	write(12,*)xst(nst),yst(nst),zst(nst)
	!write(*,*)xst(nst),yst(nst),zst(nst)
	goto 33
44	close(12)
close(1)
!write(*,*)' nst=',nst

if(k_re1_syn2.eq.1) then
	open(1,file='../../../DATA/'//ar//'/inidata/rays.dat')
	write(*,*)' ******************************************'
	write(*,*)' REAL DATA INVERSION'
else if(k_re1_syn2.eq.2) then
	write(*,*)' ******************************************'
	write(*,*)' SYNTHETIC DATA INVERSION'
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays_syn.dat')
else
	write(*,*)' k_re1_syn2=',k_re1_syn2
	pause
end if

open(11,file='../../../DATA/'//ar//'/'//md//'/data/rays0.dat',form='unformatted')
open(12,file='../../../DATA/'//ar//'/'//md//'/data/ztr0.dat')
open(14,file='../../../DATA/'//ar//'/'//md//'/data/ztr_old.dat')
open(21,file='../../../TMP_files/rays/dt_dist_p0.dat')
open(22,file='../../../TMP_files/rays/dt_dist_s0.dat')
open(31,file='../../../TMP_files/rays/rays_p0.bln')
open(32,file='../../../TMP_files/rays/rays_s0.bln')
open(33,file='../../../TMP_files/rays/srces0.dat')

izt=0
nray=0
nrp=0
nrs=0
nztgood=0
dtot=0
ntot=0

!ifreq=1

! Read the sources:
992	continue
	read(1,*,end=991)fini,tini,zini,nkrat
    if(nkrat.eq.0) goto 992

    zold=zini
	zlim=z_lim(fini,tini)
	if(zold.ge.zlim) zold=zlim-1

	zlim_up=h_lim(fini,tini)
	if (zold.le.zlim_up) zold=zlim_up


	!write(*,*)fini,tini,zold,nkrat
	call SFDEC(fini,tini,0.,xold,yold,Z,fi0,tet0)
	if(zold.le.0.) zold=0
	!write(*,*)xold,yold,zold,nkrat 
	!read(1)imt,idy,ihr,imn		
	!write(*,*)imt,idy,ihr,imn		
	!xbl=xbl-100.5
	!ybl=ybl+80.7

	izt=izt+1
	!write(*,*)' nkrat=',nkrat

! Read all the records:

	do i=1,nkrat
		read(1,*)ips,ist,tobs
		!write(*,*)ips,ist,tobs
		iqua=1
		!read(1,end=991)ist,ips,iqua,tobs
		istkr(i)=ist	!ist: code of station, 
		ipskr(i)=ips
		tobkr(i)=tobs	! tobs: observerd arrival time
		qualkr(i)=w_qual(iqua)
	end do
	!write(*,*)fini,tini,zold 
!	if(izt.lt.554) goto 992
!    ifreq=1
!	write(*,*)izt,fini,tini,zold,nkrat
	if(nkrat.lt.krat_min) goto 992

		res_loc1=res_it1(1)
		res_loc2=res_it2(1)
		w_P_S_diff=wps_it(1)
		dx_loc=dx_it(1)
		dy_loc=dy_it(1)
		dz_loc=dz_it(1)


	call goal_new(xold,yold,zold, disp,aver,nk,aold)

	!write(*,*)' xold=',xold,' yold=',yold,' zold=',zold,' aold=',aold

    !do i=1,nkrat
	!	write(*,'(2i4,4f10.3)')ipskr(i),istkr(i),tobkr(i),trfkr(i),tobkr(i)-trfkr(i)-aold,diskr(i)
	!end do
!    pause

	if(koe.eq.1.and.mod(izt,2).eq.0) goto 992
	if(koe.eq.2.and.mod(izt,2).eq.1) goto 992
	!if(zold.gt.650) zold=500
	!write(*,*)izt,xold,yold,zold,nkrat 
	!do i=1,nkrat
		!write(*,*)istkr(i),ipskr(i),tobkr(i)
	!end do
	!pause


	dismin=9999999
	do i=1,nst
		hordist=sqrt((xst(i)-xold)*(xst(i)-xold)+(yst(i)-yold)*(yst(i)-yold))
		if(hordist.lt.dismin) dismin=hordist
	end do
	!write(*,*)' dismin 1111=',dismin
	!if(dismin.gt.dist_max) goto 992

	xmin=xold	!+300
	ymin=yold	!-300
	zmin=zold	!+500

	call decsf(xmin,ymin,0.,fi0,tet0,fff,ttt,h)
	zlim=z_lim(fff,ttt)
	if(zmin.ge.zlim) zmin=zlim-1

	zlim_up=h_lim(fff,ttt)
	if (zmin.le.zlim_up) zmin=zlim_up


	do iter=1,niter_loc

		res_loc1=res_it1(iter)
		res_loc2=res_it2(iter)
		w_P_S_diff=wps_it(iter)
		dx_loc=dx_it(iter)
		dy_loc=dy_it(iter)
		dz_loc=dz_it(iter)

		!write(*,*)res_loc1,res_loc2,w_P_S_diff,dx_loc,dy_loc,dz_loc

		!write(*,*)xmin,ymin,zmin,amax
		call goal_new(xmin,ymin,zmin, disp,aver,nk,amax)

		!write(*,*)' amax=',amax

		!write(*,*)res_loc1,res_loc2,w_P_S_diff


		nkode=1
		kodes(1,nkode)=0
		kodes(2,nkode)=0
		kodes(3,nkode)=0
		dkode(nkode)=amax
		ixmin1=0
		iymin1=0
		izmin1=0

	282 continue
		index=0
		do iix=1,5
			ix=ixmin1+iix-3
			dx=dx_loc*ix
			do iiy=1,5
				iy=iymin1+iiy-3
				dy=dy_loc*iy
				do iiz=1,5
					iz=izmin1+iiz-3
					do ik=1,nkode
						if(kodes(1,ik).eq.ix.and.kodes(2,ik).eq.iy.and.kodes(3,ik).eq.iz) goto 281
					end do
					dz=dz_loc*iz
					zzz=zmin+dz

					call decsf(xmin+dx,ymin+dx,0.,fi0,tet0,fff,ttt,h)
					zlim=z_lim(fff,ttt)
					if(zzz.gt.zlim) cycle

					zlim_up=h_lim(fff,ttt)
					if(zzz.lt.zlim_up) cycle

!write(*,*)xmin+dx,ymin+dy,zzz
					!write(*,'(3i3,3f6.1,f7.3)')ix,iy,iz,xmin+dx,ymin+dy,zzz
                    !write(*,*)' zzz=',zzz,' zlim=',zlim,' zlim_up=',zlim_up
					call goal_new(xmin+dx,ymin+dy,zzz, disp,aver,nk,ank)
					!write(*,*)xmin+dx,ymin+dy,zzz,' ank=',ank

					nkode=nkode+1
					kodes(1,nkode)=ix
					kodes(2,nkode)=iy
					kodes(3,nkode)=iz
					dkode(nkode)=ank
					if(ank.le.amax) cycle
					index=1
					ixmin=ix
					iymin=iy
					izmin=iz
					amax=ank
					!write(*,*)ix,iy,iz,' ank=',ank
	281				continue
				end do
			end do
		end do
!		write(*,*)ixmin,iymin,izmin,amax
        !stop
		if(index.eq.1) then
			ixmin1=ixmin
			iymin1=iymin
			izmin1=izmin
			goto 282
		end if

		xmin=xmin+dx_loc*(ixmin1)
		ymin=ymin+dy_loc*(iymin1)
		zmin=zmin+dz_loc*(izmin1)
		!write(*,*)' after iteration:',iter
		!write(*,*)xmin,ymin,zmin,amax
	end do

	call goal_new(xold,yold,zold, disp,aver,nk,aold)

!write(*,*)' xold=',xold,' yold=',yold,' zold=',zold,' aold=',aold
!	do i=1,nkrat
!		write(*,'(2i4,4f8.3)')ipskr(i),istkr(i),tobkr(i),trfkr(i),tobkr(i)-trfkr(i)-aver,diskr(i)
!	end do

	call goal_new(xmin,ymin,zmin, disp,aver,nk,ank)
!	write(*,*)' z=',zmin,' ank=',ank,' d=',disp,' nk=',nk
!write(*,*)' xmin=',xmin,' ymin=',ymin,' zmin=',zmin,' ank=',ank
!	do i=1,nkrat
!		write(*,'(2i4,3f8.3)')ipskr(i),istkr(i),tobkr(i)-trfkr(i)-aver,diskr(i)
!	end do
!    stop


	dismin=9999999
	do i=1,nst
		hordist=sqrt((xst(i)-xmin)*(xst(i)-xmin)+(yst(i)-ymin)*(yst(i)-ymin))
		if(hordist.lt.dismin) dismin=hordist
	end do
	!write(*,*)' dismin 2222=',dismin
	!pause
	if(dismin.gt.dist_max) goto 992


	do i=1,nkrat
		tobkr(i)=tobkr(i)-aver
		dt=tobkr(i)-trfkr(i)
		!write(*,347)istkr(i),ipskr(i),diskr(i),dtold(i),dt
347	format(' ist=',i4,' ips=',i3,' dist=',f7.1,' dt_old=',f9.3,' dt_new=',f9.3)
	end do

! Determine the GOOD and BAD residuals:

	nbad=0
	ngood=1
	do i=1,nkrat
		ist=istkr(i)
		xs=xst(ist)
		ys=yst(ist)
		dhor=sqrt((xs-xmin)*(xs-xmin)+(ys-ymin)*(ys-ymin))
		dist=sqrt(dhor*dhor+zmin*zmin)
		if(dist.gt.sss_max) dist=sss_max
		res_limit=dist*res_1_km
		if(ipskr(i).eq.2)res_limit=res_limit*wgs
		dt=tobkr(i)-trfkr(i)
        if(ipskr(i).eq.1)write(21,*)dhor,dt
        if(ipskr(i).eq.2)write(22,*)dhor,dt
		!write(*,*)ipskr(i),istkr(i),dt,res_limit
		if(abs(dt).lt.res_limit) cycle
		ngood(i)=0
		nbad=nbad+1
	end do
    !pause

	nk=nkrat-nbad
	abad=nbad
	akrat=nkrat
	ratio_bad=(abad/akrat)

	!write(*,*)' nbad=',nbad,' ngood=',nk
	if(ratio_bad*100.gt.bad_max) then
		!write(*,*)' BAD event!!!'
		!write(*,488)izt,xold,yold,zold,aold
		goto 992
	end if
	if(nk.lt.krat_min) goto 992

	!write(*,*)' disp=',disp,' nk=',nk,' ank=',ank
!	if(nk.lt.krat_min) goto 992


1515 continue
    call decsf(xold,yold,0.,fi0,tet0,fold,told,h)
    call decsf(xmin,ymin,0.,fi0,tet0,fnew,tnew,h)
    write(11)xmin,ymin,zmin,nk
    !write(*,*)xmin,ymin,zmin,nk
    write(12,*)fnew,tnew,zmin
    write(33,*)fnew,tnew,zmin
    write(14,*)fini,tini,zini
    !write(22,*)2
    !write(22,*)fold,told
    !write(22,*)fnew,tnew
    nk1=0

!open(31,file='test_event.dat')
!write(31,*)' fi=',fnew
!write(31,*)' tet=',tnew
!write(31,*)' dep=',zmin
!write(31,773)
!write(*,773)
773 format(' ips ','|',' sta ','|',' qua ','|','   dt   ','|','  dist  ','|','  angle')
	do i=1,nkrat
		if(ngood(i).eq.0) cycle
		write(11)ipskr(i),istkr(i),tobkr(i),trfkr(i)
        fstat=fst(istkr(i)); tstat=tst(istkr(i))
        nfle=31
        if(ipskr(i).eq.2)nfle=32
        write(nfle,*)2
        write(nfle,*)fnew,tnew
        write(nfle,*)fstat,tstat
		!write(*,*)istkr(i),ipskr(i),tobkr(i),trfkr(i)
		!write(*,'(2i4,3f8.3)')ipskr(i),istkr(i),qualkr(i),tobkr(i)-trfkr(i),diskr(i)
!write(31,772)ipskr(i),istkr(i),qualkr(i),tobkr(i)-trfkr(i),diskr(i),alkr(i)
!write(*,772)ipskr(i),istkr(i),qualkr(i),tobkr(i)-trfkr(i),diskr(i),alkr(i)
772 format(i5,'|',i5,'|',f5.1,'|',f8.3,'|',f8.2,'|',f8.2)
		dtot=dtot+abs(tobkr(i)-trfkr(i))
		ntot=ntot+1
		nk1=nk1+1
		nray=nray+1
		if(ipskr(i).eq.1) then
			nrp=nrp+1
		else
			nrs=nrs+1
		end if
	end do
!close(31)
	if(nk.ne.nk1) pause



	nztgood=nztgood+1
	!if(mod(nztgood,1).eq.0.and.index.eq.0) then
	!if(mod(nztgood,ifreq).eq.0.or.abs(fnew-27.8).lt.0.1) then
	if(mod(nztgood,ifreq).eq.0) then
!		write(*,*)' nzt=',nztgood,' nray=',nray,' np=',nrp,' ns=',nrs
!		write(*,*)' nkold=',nkold,' dispold=',dispold
		!write(*,*)' averbest=',averbest
		!write(*,341)izt,jyr,jmc,jdy,zmin,nkrat,amag,timzt
		!write(15,341)izt,jyr,jmc,jdy,zmin,nkrat,amag,timzt
341 format(i4,' y=',i4,' m=',i3,' d=',i3,' z=',f7.2,' kr=',i4,' mag=',f7.2,' t=',f12.3)
		!write(*,*)' month=',imt,' day=',idy,' hr=',ihr
		!write(*,*)' dismin 2222=',dismin
		write(*,488)izt,xold,yold,zold,aold
		write(*,489)nztgood,xmin,ymin,zmin,ank
	488	format(i6,' Old: x=',f8.2,' y=',f8.2,' z=',f8.2,' ank=',f7.2)
	489	format(i6,' New: x=',f8.2,' y=',f8.2,' z=',f8.2,' ank=',f7.2)
		dcur=dtot/ntot
		write(*,*)' nkrat=',nkrat,' nk=',nk,' ntot=',ntot
		!write(*,*)' ntot=',ntot,' dcur=',dcur
		write(*,*)
		!if(ank.lt.0.6)pause
		!stop
		!pause
	end if
	!if(ntot.gt.nraymax)goto 991

	goto 992
991 continue
close(1)
close(14)
close(21)
close(22)
close(31)
close(32)
close(33)
write(*,*)' nztgood=',nztgood
stop
end


