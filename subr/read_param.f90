subroutine read_param(ar,md)
character*8 ar,md,line
character*20 line20
character*20 scale_line,scale_vpvs,scale,scale_vp,scale_vs
integer krrr(100),kggg(100),kbbb(100)

common/flat_sf/key_flat11
common/general/key_1real_2syn,VPSX_key,koe,kref,key_ft1_xy2,key_true1,&
	key_flat1,key_loc
common/orient/nornt,ornt(10)
common/center/fi0,tet0
common/ref_table/zstat,dmin,depmax,distmax,nlay,zst(20),dzst(20),zztmax
common/loc_table/krat_min,dist_max,wgs,dist_limit0,n_pwr_dist0,ncyc_av0,bad_max,res_1_km,sss_max,ifreq0, &
	niter_loc,dx_it(10),dy_it(10),dz_it(10),res_it1(10),res_it2(10),wps_it(10)

common/ray_param/ds_ini,ds_part_min,val_bend_min,bend_max0
common/loc_param/wgs1,res_loc1,res_loc2,dist_limit1,n_pwr_dist1,ncyc_av1,w_P_S_diff
common/loc_other/stepmax,stepmin,ifreq1
common/grid/xlim1,xlim20,dxpl,ylim1,ylim20,dypl,zlim1,zlim2,dzpl,plotmin,plotmax
common/inversion/iter_lsqr,wg_vel_p,wg_vel_s,sm_hor_p,sm_hor_s,sm_ver_p,sm_ver_s,&
	rg_amp_p,rg_amp_s,wg_st_p,wg_st_s,wzt_hor,wzt_ver,wzt_time
common/mod/xx1,xx2,dxx,yy1,yy2,dyy,zz1,zz2,dzz,smaxx_3D,ismth_3D
common/visual_hor/ nlev,hlev(20),fmap1,fmap2,dfmap,tmap1,tmap2,dtmap,smaxx,ismth_h,size_hor
common/visual_ver/ nver,fia0(20),teta0(20),fib0(20),tetb0(20),dist_from_sec_event,&
	dxsec,zmin,zmax,dzsec,dsmark,dismax,ismth_v,size_x,size_z
common/scales/dv_min,dv_max,vpvs_min,vpvs_max,vp_min,vp_max,vs_min,vs_max
common/scl_names/scale_line,scale_vpvs,scale_vp,scale_vs


!******************************************************************
key_ft1_xy2=1
key_true1=0
key_flat1=1
VPSX_key=0
kref=0
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
    read(1,'(a8)',end=513)line
    if(line.eq.'GENERAL ') goto 514
end do
513 continue
write(*,*)' cannot find GENERAL INFORMATION in MAJOR_PARAM.DAT!!!'
pause
514 continue
    read(1,*)key_1real_2syn
    read(1,*)koe
	read(1,*)key_loc
	read(1,*,end=441,err=441)key_flat1      ! 1: calculations in flat model, 2: spherical velocity
441 close(1)
key_flat11=key_flat1

!******************************************************************
if(key_ft1_xy2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
    do i=1,10000
        read(1,'(a8)',end=853)line
        if(line.eq.'AREA_CEN') goto 854
    end do
    853 continue
    write(*,*)' cannot find AREA CENTER in MAJOR_PARAM.DAT!!!'
    pause
    854 read(1,*)fi0,tet0
    close(1)
else
    fi0=0
    tet0=0
end if

!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=533)line
	if(line.eq.'GRID_PAR') goto 534
end do
533 continue
write(*,*)' cannot find GRID_PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
534 continue
574 read(1,*)nornt
read(1,*)(ornt(i),i=1,nornt)
read(1,*)xlim1,xlim20,dxpl
read(1,*)ylim1,ylim20,dypl
read(1,*)zlim1,zlim2,dzpl
read(1,*)plotmin,plotmax	! maximal ray density, relative to average
close(1)
xx1=xlim1; xx2=xlim20; dxx=dxpl
yy1=ylim1; yy2=ylim20; dyy=dypl
zz1=zlim1; zz2=zlim2; dzz=dzpl
smaxx_3D=dxx*20; ismth_3D=0

!******************************************************************
zstat=0
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,*,end=593)line
	if(line.eq.'REF_PARA') goto 594
end do
593 continue
write(*,*)' cannot find REF_PARAM in MAJOR_PARAM.DAT!!!'
pause
594 continue
read(1,*)dmin		!=0.1
read(1,*)depmax		!=100.
read(1,*)distmax		!=2000.
read(1,*)nlay	
do i=1,nlay
	read(1,*)zst(i),dzst(i)
end do
zztmax=depmax		!=50.
read(1,*,end=134,err=134)zstat		! average station level
134 close(1)
zst(nlay+1)=zztmax

!******************************************************************
wgs=1.7; dist_limit0=100; n_pwr_dist0=1; ncyc_av0=30; bad_max=30; sss_max=1
res_it1=0; wps_it=2
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=943)line
	if(line.eq.'LIN_LOC_') goto 944
end do
943 continue
write(*,*)' cannot find LIN_LOC_PARAM in MAJOR_PARAM.DAT!!!'
pause
944 continue
read(1,*)krat_min
read(1,*)dist_max
read(1,*)res_1_km
read(1,*)ifreq0
read(1,*)niter_loc
do it=1,niter_loc
	read(1,*)
	read(1,*)dx_it(it),dy_it(it),dz_it(it)
	read(1,*)res_it2(it)
end do
close(1)

wgs1=wgs
!******************************************************************

dist_limit1=50; n_pwr_dist1=1; ncyc_av1=30; res_loc1=0; w_P_S_diff=2
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=583)line
	if(line.eq.'LOC_PARA') goto 584
end do
583 continue
write(*,*)' cannot find LOC_PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
584 continue
read(1,*)ds_ini
!write(*,*)ds_ini
read(1,*)ds_part_min
read(1,*)val_bend_min
read(1,*)bend_max0
read(1,*)res_loc2
read(1,*)stepmax
read(1,*)stepmin
read(1,*)ifreq1
close(1)
!ds_part_min=ds_ini*50.
!bend_max0=ds_ini*50.
!val_bend_min=ds_ini/20.
!stepmax=ds_ini*100.
!stepmin=ds_ini*2. 



!******************************************************************
wg_vel_p=1; wg_vel_s=1
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=373)line
	if(line.eq.'INVERSIO') goto 374
end do
373 continue
write(*,*)' cannot find INVERSION PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
374 continue
read(1,*)iter_lsqr
read(1,*)sm_hor_p,sm_hor_s
read(1,*)sm_ver_p,sm_ver_s
read(1,*)rg_amp_p,rg_amp_s
read(1,*)
read(1,*)wg_st_p,wg_st_s
read(1,*)wzt_hor
read(1,*)wzt_ver
read(1,*)wzt_time
close(1)
!******************************************************************

! VISUALIZATION PARAMETERS:

ismth_h=0; ismth_v=0
open(2,file='../../../DATA/'//ar//'/set_visual.dat')
read(2,*)   
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
read(2,*) smaxx
read(2,*) 
read(2,*) 
read(2,*)nver
do ii=1,nver
	read(2,*) fia0(ii),teta0(ii),fib0(ii),tetb0(ii)
end do
read(2,*) dist_from_sec_event
read(2,*) dxsec
read(2,*) zmin,zmax,dzsec
read(2,*) dsmark
read(2,*) dismax
read(2,*) 
read(2,*) 
read(2,*)size_hor
read(2,*)size_x,size_z
read(2,*)scale_line
read(2,*)dv_min,dv_max
read(2,*)scale_vpvs
read(2,*)vpvs_min,vpvs_max
read(2,*)scale_vp
read(2,*)vp_min,vp_max
read(2,*)scale_vs
read(2,*)vs_min,vs_max
close(2)

!write(*,*)' dv_max=',dv_max,' dv_min=',dv_min

return
end