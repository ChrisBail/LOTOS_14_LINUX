character*8 ar,md,line
integer kodes(10,20000)
real gkode(20000)
real xvert1(20),yvert1(20),xvert2(20),yvert2(20)
character*2 ver
real xmark(200,20),ymark(200,20),smark(200,20)
integer nmark(100),ngood(500)
character*5 stacod(500),stac

common/refmod/nref,href(600),vref(2,600)
common/pi/pi,per
common/krat/nkrat,istkr(500),ipskr(500),ndrkr(500),tobkr(500),trfkr(500)
common/stations/nst,xstat(500),ystat(500),zstat(500)
common/ray/ nodes,xray(1000),yray(1000),zray(1000)

common/loc_param/res_loc1,res_loc2,w_P_S_diff

common/center/fi0,tet0
common/loc_table/krat_min,dist_max,wgs,dist_limit0,n_pwr_dist0,ncyc_av0,bad_max,res_1_km,sss_max, &
    nfreq_print,niter_loc,dx_it(10),dy_it(10),dz_it(10),res_it1(10),res_it2(10),wps_it(10)
common/ray_param/ds_ini,ds_part_min,bend_min0,bend_max0
common/general/k_re1_syn2,VPSX_key,koe,kref,key_ft1_xy2,key_true1,&
    key_flat1,key_loc


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

nxmod=0

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
close(1)

write(*,*)' Preliminary location of sources using straight lines in the starting model'
write(*,*)' ar=',ar,' md=',md

call read_param(ar,md)
call read_z_lim(ar,md)
call read_vref(ar,md)
call read_3D_mod_v(ar,md,0)
call read_topo(ar)

!xm=5; ym=5; zm=5; ips=2
!vvv=velocity (xm,ym,zm,ips)
!write(*,*)' vvv=',vvv

i=system('mkdir -p ../../../TMP_files/loc')
i=system('mkdir -p ../../../PICS/'//ar//'/'//md//'/LOC')


key_info=0
open(3,file='../../../DATA/'//ar//'/inidata/event_info.dat',status='old',err=46)
key_info=1
close(3)
46  continue


! Read the coordinates of the stations
open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
open(12,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
nst=0
33	    continue
        if(key_info.eq.0)read(1,*,end=44)fi,tet,zst
        if(key_info.eq.1)read(1,*,end=44)fi,tet,zst,stac
    call SFDEC(fi,tet,0.,X,Y,Z,fi0,tet0)
    nst=nst+1
    xstat(nst)=x
    ystat(nst)=y
    zstat(nst)=zst
        if(key_info.eq.1)stacod(nst)=stac
    write(12,*)xstat(nst),ystat(nst),zstat(nst)
    !write(*,*)xstat(nst),ystat(nst),zstat(nst)
    goto 33
44	close(12)
close(1)
write(*,*)' nst=',nst

if(k_re1_syn2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/inidata/rays.dat')
    write(*,*)' ******************************************'
    write(*,*)' REAL DATA INVERSION'
else if(k_re1_syn2.eq.2) then
	write(*,*)' ******************************************'
	write(*,*)' SYNTHETIC DATA INVERSION'
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays_syn.dat')
else
	write(*,*)' Synthetic key in MAJOR_PARAM is not defined correctly: k_re1_syn2=',k_re1_syn2
	stop
end if

open(11,file='../../../DATA/'//ar//'/'//md//'/data/rays0.dat',form='unformatted')
open(41,file='../../../DATA/'//ar//'/'//md//'/data/resid0.dat')
open(12,file='../../../DATA/'//ar//'/'//md//'/data/srce_0.dat')
open(14,file='../../../DATA/'//ar//'/'//md//'/data/srce_old.dat')
if(key_info.eq.1)then
    open(3,file='../../../DATA/'//ar//'/inidata/event_info.dat')
    open(15,file='../../../DATA/'//ar//'/'//md//'/data/info_0.dat')
    open(16,file='../../../DATA/'//ar//'/'//md//'/data/rays_full_0.dat')
end if

!write(*,*)' key_info=',key_info

nzt=0
errtot=0
nerr=0
nztgood=0
nrp=0
nrs=0
nray=0
228	continue
    read(1,*,end=229)fini,tini,zini,nkrat
    if(nkrat.eq.0) goto 228
    zzt=zini
    !write(*,*)fini,tini,zold,nkrat
    call SFDEC(fini,tini,0.,xzt,yzt,Z,fi0,tet0)
   

    xold=xzt; yold=yzt; zold=zzt

    do ikr=1,nkrat
        read(1,*)ips,ist,tobs
        !write(*,*)ikr,ips,ist,tobs
        ipskr(ikr)=ips
        istkr(ikr)=ist
        tobkr(ikr)=tobs
    end do

    xmax=xzt
    ymax=yzt
    zmax=zzt

    nzt=nzt+1
    if(koe.eq.1.and.mod(nzt,2).eq.0) goto 228
    if(koe.eq.2.and.mod(nzt,2).eq.1) goto 228

    if(nkrat.lt.krat_min) goto 228


!    if(nzt.lt.6103) goto 228
!
!
    !write(*,*)nzt,xzt,yzt,zzt,nkrat
!    write(*,*)' fini, tini=',fini,tini
!    write(*,*)' fi0, tet0=',fi0,tet0
!    call SFDEC(fini,tini,0.,xzt2,yzt2,Z,fi0,tet0)
!    write(*,*)' xzt2, yzt2=', xzt2, yzt2
    !do ikr=1,nkrat
        !write(*,*)ipskr(ikr),istkr(ikr),tobkr(ikr)
    !end do
    !stop

    do iter=1,niter_loc

        res_loc1=res_it1(iter)
        res_loc2=res_it2(iter)
        w_P_S_diff=wps_it(iter)
        dx_loc=dx_it(iter)
        dy_loc=dy_it(iter)
        dz_loc=dz_it(iter)

!write(*,*)ds_ini,ds_part_min,bend_min0,bend_max0
!write(*,*)res_loc1,res_loc2,w_P_S_diff

!        write(*,*)xmax,ymax,zmax
!        do ikr=1,nkrat
!            write(*,*)istkr(ikr),ipskr(ikr),tobkr(ikr)
!        end do
        !pause
        call goal_fun_lin(xmax,ymax,zmax, aver,goal)
        !write(*,*)' iter=',iter,' goal=',goal
        !stop
        gmax=goal

        nkode=1
        kodes(1,nkode)=0
        kodes(2,nkode)=0
        kodes(3,nkode)=0
        gkode(nkode)=gmax
        ixmax1=0
        iymax1=0
        izmax1=0

        icycle=0

        282 continue
        index=0
        do iix=1,5
            ix=ixmax1+iix-3
            !ix=14
            dx=dx_loc*ix
            do iiy=1,5
                iy=iymax1+iiy-3
                !iy=-14
                dy=dy_loc*iy
                call decsf(xmax+dx,ymax+dy,0.,fi0,tet0,fff,ttt,h)
                zlim=z_lim(fff,ttt)
                zlim_up=relief_surf(fff,ttt)
                !write(*,*)' ix=',ix,' iy=',iy
                !write(*,*)' fff=',fff,' ttt=',ttt,' zlim_up=',zlim_up
                do iiz=1,5
                    iz=izmax1+iiz-3
                    !iz=-1
                    dz=dz_loc*iz
                    if(zmax+dz.gt.zlim) cycle
                    if(zmax+dz.lt.zlim_up) cycle

                    if(nkode.ne.0) then
                        do ik=1,nkode
                            if(kodes(1,ik).eq.ix.and.kodes(2,ik).eq.iy.and.kodes(3,ik).eq.iz) goto 281
                        end do
                    end if


                    !write(*,'(3i3,3f6.1,f7.3)')ix,iy,iz,xmax+dx,ymax+dy,zmax+dz

                    call goal_fun_lin(xmax+dx,ymax+dy,zmax+dz, aver,goal)

                    !write(*,'(3i3,3f6.1,f7.3)')ix,iy,iz,xmax+dx,ymax+dy,zmax+dz,goal
                    !write(*,*)' goal=',goal

                    nkode=nkode+1
                    if(nkode.gt.2000) goto 228
                    kodes(1,nkode)=ix
                    kodes(2,nkode)=iy
                    kodes(3,nkode)=iz
                    gkode(nkode)=goal

                    if(goal.le.gmax) cycle
                    index=1
                    ixmax=ix
                    iymax=iy
                    izmax=iz
                    gmax=goal
                    !write(*,*)ixmax,iymax,izmax,' goal=',goal
281		     continue
                end do
            end do
        end do
        icycle=icycle+1

        !write(*,*)icycle,ixmax,iymax,izmax,gmax

        if(index.eq.1) then
            ixmax1=ixmax
            iymax1=iymax
            izmax1=izmax
            goto 282
        end if

        xmax=xmax+dx_loc*(ixmax1)
        ymax=ymax+dy_loc*(iymax1)
        zmax=zmax+dz_loc*(izmax1)

        !write(*,*)' after iteration:',iter
        !write(*,*)' x=',xmax,' y=',ymax,' z=',zmax,' g=',gmax


    end do

    xzt=xmax
    yzt=ymax
    zzt=zmax

    call goal_fun_lin(xold,yold,zold, aver,gold)
    call goal_fun_lin(xzt,yzt,zzt, aver,goal)

            !write(*,*)' old:',xold,yold,zold,gold
            !write(*,*)' new:',xzt,yzt,zzt,goal


 
    nbad=0
    ngood=1
    dismin=9999999
    do i=1,nkrat
        tobkr(i)=tobkr(i)-aver
        ist=istkr(i)
        xs=xstat(ist)
        ys=ystat(ist)
        zs=zstat(ist)
        dist=sqrt((xs-xzt)*(xs-xzt)+(ys-yzt)*(ys-yzt)+(zs-zzt)*(zs-zzt))
        hordist=sqrt((xs-xzt)*(xs-xzt)+(ys-yzt)*(ys-yzt))
        if(hordist.lt.dismin) dismin=hordist
        if(dist.gt.sss_max) dist=sss_max
        res_limit=dist*res_1_km
        if(ipskr(i).eq.2)res_limit=res_limit*wgs
        dt=tobkr(i)-trfkr(i)
        !write(*,*)ipskr(i),istkr(i),dt,res_limit
        if(abs(dt).lt.res_limit) cycle
        ngood(i)=0
        nbad=nbad+1
    end do
    !stop

    nk=nkrat-nbad
    abad=nbad
    akrat=nkrat
    ratio_bad=(abad/akrat)
    !write(*,*)' nk=',nk,' nkrat=',nkrat

    !write(*,*)' nbad=',nbad,' ngood=',nk
    !pause
    if(ratio_bad*100.gt.bad_max) goto 228
    if(nk.lt.krat_min) goto 228
    if(dismin.gt.dist_max) goto 228


    if(key_ft1_xy2.eq.2) then 
        write(12,*)xzt,yzt,zzt
    else
        call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
        write(12,*)fzt,tzt,zzt
    end if

    write(11)xzt,yzt,zzt,nk
    write(14,*)fini,tini,zini
    write(41,*)xzt,yzt,zzt,nk
    do i=1,nkrat
        if(ngood(i).eq.0) cycle
        write(11)ipskr(i),istkr(i),tobkr(i),trfkr(i)
        write(41,*)ipskr(i),istkr(i),tobkr(i)-trfkr(i)
        dtot=dtot+abs(tobkr(i)-trfkr(i))
        nray=nray+1
        if(ipskr(i).eq.1) nrp=nrp+1
        if(ipskr(i).eq.2) nrs=nrs+1
    end do

    nztgood=nztgood+1

    !nfreq_print=1
    if(mod(nztgood,nfreq_print).eq.0) then
        write(*,488)nzt,xold,yold,zold,gold
        write(*,489)nztgood,xzt,yzt,zzt,goal
        488	format(i6,' Old: x=',f8.2,' y=',f8.2,' z=',f8.2,' ank=',f7.2)
        489	format(i6,' New: x=',f8.2,' y=',f8.2,' z=',f8.2,' ank=',f7.2)
        dcur=dtot/nray
        write(*,*)' nkrat=',nkrat,' nk=',nk,' nray=',nray
        write(*,*)
        !stop
    end if

    goto 228
229 close(1)
close(11)
close(12)
close(14)
close(41)
if(key_info.eq.1)close(3)
if(key_info.eq.1)close(15)
if(key_info.eq.1)close(16)


write(*,*)' nztgood=',nztgood,' nray=',nray
write(*,*)' nrp=',nrp,' nrs=',nrs


!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************


1313 continue
!call visual_srces_stat(ar,md,itt)


i=system('mkdir -p ../../../PICS/'//ar//'/'//md//'/LOC')


open(1,file='../../../DATA/'//ar//'/config.txt')
read(1,*)
read(1,*) npix_h_x0,npix_h_y0
read(1,*)tick_h_x,tick_h_y
read(1,*)
read(1,*) npix_v_x0,npix_v_z0
read(1,*)tick_v_x,tick_v_z
close(1)

open(1,file='../../../DATA/'//ar//'/sethor.dat')
read(1,*)
read(1,*)
read(1,*)xmap1,xmap2,dx,ymap1,ymap2,dy
close(1)

dxmark=0; dymark=0
open(1,file='../../../DATA/'//ar//'/setver.dat')
read(1,*) nvert
do ivert=1,nvert
    read(1,*)xvert1(ivert),yvert1(ivert),xvert2(ivert),yvert2(ivert)
end do
read(1,*)dist_max
read(1,*)dxsec0
read(1,*)zmap1,zmap2
read(1,*)dsmark
read(1,*)
read(1,*,end=991,err=991)dxmark,dymark
991 close(1)


do ivert=1,nvert


    write(ver,'(i2)')ivert

    if(key_ft1_xy2.eq.1) then
        fia=xvert1(ivert)
        fib=xvert2(ivert)
        teta=yvert1(ivert)
        tetb=yvert2(ivert)
            !write(*,*)' fia=',fia,' teta=',teta
            !write(*,*)' fib=',fib,' tetb=',tetb

        call SFDEC(fia,teta,0.,xa,ya,Z,fi0,tet0)
        call SFDEC(fib,tetb,0.,xb,yb,Z,fi0,tet0)

    else

        xa=xvert1(ivert)
        ya=yvert1(ivert)
        xb=xvert2(ivert)
        yb=yvert2(ivert)

    end if




    dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist

    !write(*,*)' xa=',xa,' ya=',ya
    !write(*,*)' xb=',xb,' yb=',yb

    !write(*,*)' dist=',dist
    nxsec=dist/dxsec0+1
    dxsec=dist/(nxsec-1)

    !Draw topography on the section
    open(11,file='../../../TMP_files/loc/topo_'//ver//'.bln')
    write(11,*)nxsec
    do ix=1,nxsec
        sss=(ix-1)*dxsec
        xcur=xa+((xb-xa)/dist)*sss
        ycur=ya+((yb-ya)/dist)*sss
        if(key_ft1_xy2.eq.1) then
            call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
        else
            fff=xcur
            ttt=ycur
        end if
        depth=h_lim(fff,ttt)
        ztopo=flat_sph(key_flat1,xcur,ycur,depth)
        write(11,*)sss,-ztopo
    end do
    close(11)



    open(11,file='../../../TMP_files/loc/mark_'//ver//'.dat')
    imark=0
    do sss=0.,dist,dsmark
        x=xa+cospov*sss
        y=ya+sinpov*sss
 
        if(key_ft1_xy2.eq.1) then
            call decsf(x,y,0,fi0,tet0,fi,tet,hhh)
            x=fi; y=tet
        end if
 
        write(11,*)x,y,sss
        imark=imark+1
        xmark(ivert,imark)=x
        ymark(ivert,imark)=y
        smark(ivert,imark)=sss
    end do
    imark=imark+1
        if(key_ft1_xy2.eq.1) then
            xmark(ivert,imark)=fib
            ymark(ivert,imark)=tetb
        else
            xmark(ivert,imark)=xb
            ymark(ivert,imark)=yb
       end if
    smark(ivert,imark)=dist
    close(11)
    nmark(ivert)=imark

	
    ! Draw the position of the section on the surface (line)
    open(11,file='../../../TMP_files/loc/mark_'//ver//'.bln')
    !write(*,*)' nmark(ivert)=',nmark(ivert)
    write(11,*) imark
    do i=1,imark

        write(11,*)xmark(ivert,i),ymark(ivert,i)	!,smark(i)
        !write(*,*)xmark(ivert,i),ymark(ivert,i)	!,smark(i)
    end do
    close(11)



    if(npix_v_x0.ne.0) then
        npix_x=npix_v_x0
    else
        npix_x=int(npix_v_z0*(dist/(zmap2-zmap1)))
        npix_z=npix_v_z0
    end if

    write(*,*)' Vertical:',ver,' dist=',dist,' npix_x=',npix_x,' npix_z=',npix_z

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
    open(11,file='../../../TMP_files/loc/stat_ver'//ver//'.dat')
    201	continue
        read(1,*,end=202)xst,yst,depth
        zst = flat_sph(key_flat1,xst,yst,depth)

        xx1=(xst-xa)*cospov+(yst-ya)*sinpov
        yy1=-(xst-xa)*sinpov+(yst-ya)*cospov

        if(abs(yy1).lt.dist_max) write(11,*)xx1,-zst,yy1
        goto 201
    202 close(1)
    close(11)

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays0.dat',form='unformatted')
    open(12,file='../../../TMP_files/loc/srces_ver0'//ver//'.dat')
    if(key_true1.eq.1) open(13,file='../../../TMP_files/loc/shift_ver0'//ver//'.bln')
    203	continue
        if(key_true1.eq.1)read(1,end=204)xzt0,yzt0,zzt0
        read(1,end=204)xzt,yzt,depth,nkrat
        zzt = flat_sph(key_flat1,xzt,yzt,depth)



        if(key_true1.eq.1) xx0=(xzt0-xa)*cospov+(yzt0-ya)*sinpov
        if(key_true1.eq.1) yy0=-(xzt0-xa)*sinpov+(yzt0-ya)*cospov
        xx1=(xzt-xa)*cospov+(yzt-ya)*sinpov
        yy1=-(xzt-xa)*sinpov+(yzt-ya)*cospov

 
        do i=1,nkrat
            read(1)ips,ist,tobs,tref
        end do
        if(abs(yy1).lt.dist_max) write(12,*)xx1,-zzt

        if(key_true1.eq.1) then
            if(abs(yy1).lt.dist_max) then
                write(13,*)2
                write(13,*)xx0,-zzt0
                write(13,*)xx1,-zzt
            end if
        end if

        goto 203
    204 close(1)
    close(12)
    if(key_true1.eq.1) close(13)




end do !vertical sections


if(key_ft1_xy2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
else
    open(1,file='../../../DATA/'//ar//'/inidata/stat_xy.dat')
end if
open(11,file='../../../TMP_files/loc/stat_hor.dat')
nst=0
431	continue
    read(1,*,end=432)xst,yst,zst
    write(11,*)xst,yst
    nst=nst+1
    goto 431
432 close(1)
close(11)
write(*,*)' nst=',nst

open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays0.dat',form='unformatted')
open(12,file='../../../TMP_files/loc/srces_hor0.dat')
if(key_true1.eq.1) open(13,file='../../../TMP_files/loc/shift_hor0.bln')
433	continue
    if(key_true1.eq.1) read(1,end=434)xzt0,yzt0,zzt0
    read(1,end=434)xzt,yzt,zzt,nkrat
    do i=1,nkrat
        read(1)ips,ist,tobs,tref
    end do
    if(key_ft1_xy2.eq.1) then
      call decsf(xzt,yzt,zzt,fi0,tet0,fzt,tzt,hhh)
      call decsf(xzt0,yzt0,zzt,fi0,tet0,fzt0,tzt0,hhh)
      write(12,*)fzt,tzt
        if(key_true1.eq.1) then
            write(13,*)2
            write(13,*)fzt0,tzt0
            write(13,*)fzt,tzt
        end if
    else
        write(12,*)xzt,yzt
        if(key_true1.eq.1) then
            write(13,*)2
            write(13,*)xzt0,yzt0
            write(13,*)xzt,yzt
        end if
    end if
    goto 433
434 close(1)
close(12)
if(key_true1.eq.1) close(13)


stop
end