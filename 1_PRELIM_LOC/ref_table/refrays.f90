character*8 ar,md,line
real tref(100000),dref(100000),alref(100000),href(100000)

common/pi/pi,per
common/refmod/nrefmod,zref(600),vref(600,2)
common/ref_table/zstat,dmin,depmax,distmax,nlay,zst(20),dzst(20),zztmax

one=1.d0
pi=asin(one)*2.d0
per=pi/180.d0
rz=6371.

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the model
close(1)

write(*,*)' Computing the reference table:'
write(*,*)' ar=',ar,' md=',md		


call read_vref(ar,md)
call read_param(ar,md)

!do i=1,nrefmod
!    write(*,*)zref(i),vref(i,1),vref(i,2)
!end do

dzzt=dzst(1)
zzt=zst(1)-dzzt
izt=0
alfa1=180
alfa2=0
dalfa=-0.02
nalfa=int((alfa2-alfa1)/dalfa) + 1
write(*,*)' nalfa=',nalfa

open(11,file='../../../DATA/'//ar//'/'//md//'/data/table.dat',form='unformatted')

do ilay=1,nlay
    z1=zst(ilay)
    z2=zst(ilay+1)
    nstep=int((z2-z1)/dzst(ilay))
    dstep=(z2-z1)/nstep
    !write(*,*)' ilay=',ilay,' z1=',z1,' z2=',z2
    if(ilay.eq.nlay)nstep=nstep+1
    do izzz=1,nstep
        zzt=z1+(izzz-1)*dstep
        !zzt=15
	izt=izt+1
    !write(*,*)' ilay=',ilay,' izzz=',izzz,' zzt=',zzt
        zlow=zzt
        zup=zstat

        if(zstat.gt.zzt) then
            zlow=zstat
            zup=zzt
        end if
	
	do ips=1,2
		!write(*,*)' ips=',ips,' zlow=',zlow,' zup=',zup
		dlast=999
		nref=0
		!open(31,file='tmp.dat')
		do ial=1,nalfa
		!do al=180, 0.d0, -0.02d0
			alfa0=alfa1+(ial-1)*dalfa
            !alfa0=54.82
		!do alfa0=89.986, 0.d0, -0.002d0
        !alfa0=90

			!write(*,*)alfa0,zlow,zup,ips

			call reftrace(alfa0,zlow,zup,ips,  time,dist,hmax)
			!write(*,*)alfa0,dist,hmax,time
            !stop

			dgrad=(dist/rz)/per
		!	if (hmax.ge.hmod(nrefmod-1)) exit
                        if(hmax.lt.zlow-1.e-5) cycle
			if(dist.lt.-0.01) cycle
			dkm=dist
			if(abs(dkm-dlast).gt.dmin) then
				nref=nref+1
				tref(nref)=time
				dref(nref)=dkm
				if(nref.eq.1)dref(nref)=0
				alref(nref)=alfa0
				href(nref)=hmax
				dlast=dkm
				!write(*,*)nref,alfa0,dist,time,hmax
                !stop
				!write(31,*)nref,alfa0,dkm,time,hmax
			end if
			if (hmax.gt.depmax)exit 
			if (dist.gt.distmax)exit 
				
		end do
        !stop
       ! pause
		!close(31)
		write(11)zzt,nref
		!write(31,*)zzt,nref
		!write(*,*)' i=',izt,' z=',zzt,' ips=',ips,' nref=',nref
        !stop
		if(mod(izt,10).eq.0) write(*,*)' i=',izt,' z=',zzt,' ips=',ips,' nref=',nref
		do i=1,nref
			write(11)dref(i),tref(i),alref(i),href(i)
			!write(31,*)dref(i),tref(i),alref(i),href(i)
			!write(*,*)dref(i),tref(i),alref(i),href(i)
		end do
		!pause
	end do
    end do
end do
35 close(11)

write(*,*)' number of Z levels:',izt
stop
end 
