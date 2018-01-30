character*4 dsaa/'DSAA'/
character*20 scale_line,scale_vpvs,scale,scale_vp,scale_vs,char_x,char_z,char_lim(4)

character*8 ar,md,line
character*2 lv, ver
character*1 ps ,rm,it,ch1
character*10 char10,ch_lim(4),ch_col(4)
character*50 grd, grdToRead, outFile, cpt
character*300 char300(10)
allocatable dvan(:,:),aprio_dv(:,:),vvv(:,:),vtmp(:,:),vabs(:,:),vabsp(:,:)
real fmark(200,20),tmark(200,20),smark(200,20)
integer nmark(100)
integer nrps(2)


common/pi/pi,per

common/general/key_1real_2syn,VPSX_key,koe,kref,key_ft1_xy2,key_true1,key_flat1
common/orient/nornt,ornt(10)
common/center/fi0,tet0
common/visual_hor/ nlev,hlev(20),fmap1,fmap2,dfmap,tmap1,tmap2,dtmap,smaxx,ismth_h
common/visual_ver/ nver,fia0(20),teta0(20),fib0(20),tetb0(20),dist_from_sec_event,&
    dxsec,zmin,zmax,dzsec,dsmark,dismax,ismth_v,size_x,size_z
common/scales/dv_min,dv_max,vpvs_min,vpvs_max,vp_min,vp_max,vs_min,vs_max

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0
rz=6371.

w_limit=0.2

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)iter		
close(1)

!md='MODEL_11'; iter=5
write(it,'(i1)')iter

aaa=-9857.567

!call real2char(aaa,4,char10)
!write(*,*)' char10=',char10


write(*,*)
write(*,*)' ***********************************************'
write(*,*)' VISUALISATION in vertical section: '
write(*,*)' ar=',ar,' md=',md,' iter=',iter

!******************************************************************
call read_param(ar,md)
call create_color_scales()
call read_vref(ar,md)
call read_topo(ar)

aaa=(dv_max-dv_min)/10.; call real2char(aaa,0,ch_col(1))
!write(*,*)' dv_max=',dv_max,' dv_min=',dv_min
!write(*,*)' aaa=',aaa,' ch_col(1)=',ch_col(1)
aaa=(vpvs_max-vpvs_min)/10.; call real2char(aaa,2,ch_col(2))
aaa=(vp_max-vp_min)/10.; call real2char(aaa,2,ch_col(3))
aaa=(vs_max-vs_min)/10.; call real2char(aaa,2,ch_col(4))



ismth=ismth_v

i=system('mkdir -p ../../../TMP_files/vert')
i=system('mkdir -p ../../../PICS/'//ar//'/'//md//'/IT'//it)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray1.dat')
read(1,*) nrps(1),nrps(2)
close(1)

ngr1=1
ngr2=nornt


add_perc=0
kod_av_bias=0
kod_apriori=0
ind_srce=1
!******************************************************************

ksec_ver=1
kps_ver=1


write(it,'(i1)')iter


rsmth=ismth+0.5


open(16,file='../../../TMP_files/info_map.txt')
write(16,*)nver
write(16,*)fmap1+fi0,fmap2+fi0
write(16,*)tmap1+tet0,tmap2+tet0

do iver=1,nver
    write(ver,'(i2)')iver

    if(key_ft1_xy2.eq.1) then
        fia=fia0(iver)
        teta=teta0(iver)
        fib=fib0(iver)
        tetb=tetb0(iver)
        call SFDEC(fia,teta,0.,xa,ya,Z,fi0,tet0)
        call SFDEC(fib,tetb,0.,xb,yb,Z,fi0,tet0)
    else
        xa=fia0(iver)
        ya=teta0(iver)
        xb=fib0(iver)
        yb=tetb0(iver)
    end if
    !write(*,*)' xa=',xa,' ya=',ya
    !write(*,*)' xb=',xb,' yb=',yb
    dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
    ch1='/'

    call real2char(0.,0,ch_lim(1))
    call real2char(dist,1,ch_lim(2))
    call real2char(-zmax,0,ch_lim(3))
    call real2char(-zmin,0,ch_lim(4))
    !write(*,*)' limits=',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4))
    !stop

    if(size_x.lt.0.0001) then
        size_xx = size_z * (dist/(zmax-zmin))
        size_zz = size_z
    else if(size_z.lt.0.0001) then
        size_zz = size_x * ((zmax-zmin)/dist)
        size_xx = size_x
    else
        size_xx = size_x
        size_zz = size_z
    end if
    write(char_x,'(f5.1)')size_xx
    if(size_xx.lt.100.) write(char_x,'(f4.1)')size_xx
    if(size_xx.lt.10.) write(char_x,'(f3.1)')size_xx
    write(char_z,'(f5.1)')size_zz
    if(size_zz.lt.100.) write(char_z,'(f4.1)')size_zz
    if(size_zz.lt.10.) write(char_z,'(f3.1)')size_zz
    write(char_lim(1),'(f5.1)')size_zz
    if(size_zz.lt.100.) write(char_z,'(f4.1)')size_zz
    if(size_zz.lt.10.) write(char_z,'(f3.1)')size_zz
    write(*,*)' dist=',dist,' size=',TRIM(char_x),ch1,TRIM(char_z)

    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist
    nxsec=dist/dxsec+1
    dxsec=dist/(nxsec-1)
    nzsec=(zmax-zmin)/dzsec+1
    dzsec=(zmax-zmin)/(nzsec-1)
    !write(*,*)' dist=',dist,' nxsec=',nxsec

    allocate (dvan(nxsec,nzsec),aprio_dv(nxsec,nzsec),vvv(nxsec,nzsec))
    allocate (vtmp(nxsec,nzsec),vabs(nxsec,nzsec),vabsp(nxsec,nzsec))
    vvv=0
    dvan=0


    open(11,file='../../../TMP_files/vert/mark_'//ver//'.dat')
    imark=0
    do sss=0.,dist,dsmark
        x=xa+cospov*sss
        y=ya+sinpov*sss
        if(key_ft1_xy2.eq.1) then
            call decsf(x,y,0.,fi0,tet0,FI,TET,h)
        else
            fi=x
            tet=y
        end if
        write(11,*)fi,tet,sss
        imark=imark+1
        fmark(iver,imark)=fi
        tmark(iver,imark)=tet
        smark(iver,imark)=sss
    end do
    imark=imark+1
    fmark(iver,imark)=fib
    tmark(iver,imark)=tetb
    smark(iver,imark)=dist
    close(11)
    nmark(iver)=imark

	
    ! Draw the position of the section on the surface (line)
    open(11,file='../../../TMP_files/vert/mark_'//ver//'.bln')
    write(11,*) imark
    do i=1,imark
	    write(11,*)fmark(iver,i),tmark(iver,i)	!,smark(i)
    end do
    close(11)

    write(16,*)'mark_',ver,'.bln'

    !Draw topography on the section
    open(11,file='../../../TMP_files/vert/topo0_'//ver//'.bln')
    open(12,file='../../../TMP_files/vert/topo_'//ver//'.bln')
    open(13,file='../../../TMP_files/vert/topo_'//ver//'.dat')
    write(11,*)nxsec
    write(12,*)nxsec
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
        write(11,*)sss,-depth
        write(12,*)sss,-ztopo
        write(13,*)sss,-ztopo
    end do
    close(11)
    close(12)
    close(13)

    if(ind_srce.ne.0) then
        ! Read the coordinates of the stations
        open(2,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
        open(12,file='../../../TMP_files/vert/stat_'//ver//'.dat')
        i=0
        nst1=0
        3   i=i+1
            read(2,*,end=4)xst,yst,depth
            zst = flat_sph(key_flat1,xst,yst,depth)

            xx1=(xst-xa)*cospov+(yst-ya)*sinpov
            yy1=-(xst-xa)*sinpov+(yst-ya)*cospov

            if(abs(yy1).lt.dist_from_sec_event) then
                nst1=nst1+1
                write(12,*)xx1,-zst
            end if
            goto 3
        4   close(2)
        close(12)



        !open(1,file='../../data/'//ar//'/inidata/events.dat')

        open(1,file='../../../DATA/'//ar//'/'//md//'/data/srces'//it//'.dat')
        open(11,file='../../../TMP_files/vert/ztr_'//ver//'.dat')
        nzt1=0
        nzt=0
    872	read(1,*,end=871)fzt,tzt,depth
            if(key_ft1_xy2.eq.1) then
                call SFDEC(fzt,tzt,0.,xzt,yzt,Z,fi0,tet0)
            else
                xzt=fzt; yzt=tzt
            end if
            zzt = flat_sph(key_flat1,xzt,yzt,depth)
            nzt=nzt+1
            xx1=(xzt-xa)*cospov+(yzt-ya)*sinpov
            yy1=-(xzt-xa)*sinpov+(yzt-ya)*cospov
            if(abs(yy1).lt.dist_from_sec_event) then
                nzt1=nzt1+1
                write(11,*)xx1,-zzt,yy1
            end if
            goto 872
        871 close(1)
        close(11)
        write(*,*)' nst1=',nst1,' nzt1=',nzt1
    end if

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
    do ips=1,2

	if(nrps(ips).eq.0) cycle

        write(ps,'(i1)')ips
        vvv=0
        dvan=0

    !write(*,*)' dist=',dist,' nxsec=',nxsec

        do igr=ngr1,ngr2
            call prepare_model_v(ar,md,ips,iter,igr)

            do ix=1,nxsec
                sss=(ix-1)*dxsec
                !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
                !sss=500
                !if(mod(ix,10).eq.0)write(*,*)' ix=',ix,' sss=',sss
                xcur=xa+((xb-xa)/dist)*sss
                ycur=ya+((yb-ya)/dist)*sss
   
                if(key_ft1_xy2.eq.1) then
                    call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
                else
                    fff=xcur; ttt=ycur
                end if



                !write(*,*)' xcur=',xcur,' ycur=',ycur
                !write(*,*)' fi=',fff,' tet=',ttt
                do iz=1,nzsec
                    zcur=zmin+(iz-1)*dzsec
                    depth = sph_flat(key_flat1,xcur,ycur,zcur)
                    vref=vrefmod(depth,ips)
                    if(igr.eq.ngr1) then
	                    aprio_dv(ix,iz)=0
                    end if
                    !zcur=15
                    !write(*,*)' avmoho=',avmoho,' z_moho=',z_moho

                    call dv_1_grid_v(fff,ttt,zcur,dismax,   dv,umn)
					
                    depth = h_lim(fff,ttt)
                    zlim_up = flat_sph(key_flat1,xcur,ycur,depth)



                    if (zcur.lt.zlim_up) umn=0

                    !write(*,*)' zcur=',zcur,' dv=',dv,' umn=',umn
                    !pause

                    dvan(ix,iz) = dvan(ix,iz) + dv*umn
                    vvv(ix,iz)=vvv(ix,iz)+umn
                end do
            end do
        end do

!        write(*,*)' nxsec=',nxsec,' nzsec=',nzsec



        !***************************************************************
        !***************************************************************
        !***************************************************************

        do iz=nzsec,1,-1
        !write(*,*)' iz=',iz
            do ix=1,nxsec
                vanom=-999
                if(vvv(ix,iz).gt.0.0001) then
	                vanom=dvan(ix,iz)/vvv(ix,iz)
                end if
                dvan(ix,iz)=vanom
            end do
        end do

!        write(*,*)' 2: nxsec=',nxsec,' nzsec=',nzsec

! smoothing:
        vtmp=dvan
        do iz=1,nzsec
	        do ix=1,nxsec
		        dvan(ix,iz)=-999
		        if(vvv(ix,iz).lt.0.01) cycle
		        vanm=0.
		        iv=0
		        do ixx=-ismth,ismth
			        if (ix+ixx.lt.1) cycle
			        if (ix+ixx.gt.nxsec) cycle
			        do izz=-ismth,ismth
				        if (iz+izz.lt.1) cycle
				        if (iz+izz.gt.nzsec) cycle
				        if(vvv(ix+ixx,iz+izz).lt.0.01) cycle
				        rr=ixx*ixx+izz*izz
				        r=sqrt(rr)
				        if(r.gt.rsmth) cycle
				        iv=iv+1
				        vanm=vanm+vtmp(ix+ixx,iz+izz)
			        end do
		        end do
		        dvan(ix,iz)=vanm/iv
	        end do
        end do

!        write(*,*)' 3: nxsec=',nxsec,' nzsec=',nzsec


        aver=0
        naver=0
        do iz=1,nzsec
            zcur=zmin+iz*dzsec
            do ix=1,nxsec
                sss=(ix-1)*dxsec
                !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
                !sss=500
                !if(mod(ix,10).eq.0)write(*,*)' ix=',ix,' sss=',sss
                xcur=xa+((xb-xa)/dist)*sss
                ycur=ya+((yb-ya)/dist)*sss
                depth = sph_flat(key_flat1,xcur,ycur,zcur)
                v0=vrefmod(depth,ips)
                vanom=-999
                vab=-999
                if(vvv(ix,iz).gt.w_limit) then
                    vanom=100*dvan(ix,iz)/v0
                    dvan(ix,iz)=vanom
                    vab=v0*(1+0.01*(vanom+aprio_dv(ix,iz)))
                    !if(ips.eq.2) then
                    !write(*,*)' vanom=',vanom,' vab=',vab
                    !pause
                    !end if
                    aver=aver+vanom
                    naver=naver+1
                end if
                vabs(ix,iz)=vab
            end do
        end do
!        write(*,*)' 4: nxsec=',nxsec,' nzsec=',nzsec
                   !if(ips.eq.2) then
                        !write(*,*)((vabs(ix,iz),ix=1,nxsec),iz=1,nzsec)
                        !stop
                    !end if

        !pause
        aver=aver/naver
        !dvan=dvan-aver

        if(ips.eq.1) then
	        vabsp=vabs
        else
	        do iz=1,nzsec
		        !zcur=zmin+iz*dzsec
		        !write(*,*)' zcur=',zcur
		        do ix=1,nxsec
			        !sss=(ix-1)*dxsec
			        !write(*,*)' sss=',sss
			        !write(*,*)' vp=',vabsp(ix,iz),' vs=',vabs(ix,iz),' vp/vs=',vabsp(ix,iz)/vabs(ix,iz)
			        if(abs(vabs(ix,iz)).gt.900.or.abs(vabsp(ix,iz)).gt.900) then
				        vtmp(ix,iz)=-999
			        else
				        vtmp(ix,iz)=vabsp(ix,iz)/vabs(ix,iz)
				        !write(*,*)' vp=',vabsp(ix,iz),' vs=',vabs(ix,iz),' vp/vs=',vabsp(ix,iz)/vabs(ix,iz)
			        end if

			        !vabs(ix,iz)=vab
		        end do
	        end do


	        open(14,file='../../../TMP_files/vert/vpvs_'//it//ver//'.xyz')
	        do iz=nzsec,1,-1
		        zcur=zmin+iz*dzsec
		        do ix=1,nxsec
			      sss=(ix-1)*dxsec
                            write(14,*)sss,-zcur,vtmp(ix,iz)
		        end do
	        end do
            close(14)

    !write(*,*)' ****************** dist=',dist,' nxsec=',nxsec


	        open(14,file='../../../TMP_files/vert/vpvs_'//it//ver//'.grd')
	        write(14,'(a4)')dsaa
	        write(14,*)nxsec,nzsec
	        write(14,*)0,dist
	        write(14,*)-zmax,-zmin
	        write(14,*)-9999,999
	        do iz=nzsec,1,-1
		        write(14,*)(vtmp(ix,iz),ix=1,nxsec)
	        end do
	        close(14)
			
grdToRead = '../../../TMP_files/vert/vpvs_'//it//ver//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/IT'//it//'/ver_vpvs'//ver//'.ps'
cpt = 'scale2.cpt'	

write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Vp/Vs ratio; Iteration: ',it,'; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'


grdToRead = '../../../TMP_files/vert/ztr_'//ver//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -Sc0.2c -G0/0/0 -O -K -B >> "',TRIM(outFile),'"'

grdToRead = '../../../TMP_files/vert/topo_'//ver//'.dat'
write(char300(3),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -A -W3,200/70/70 -O -K -B >> "',TRIM(outFile),'"'

write(char300(4),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(2)),' -By+l"Vp/Vs ratio" -O -Y-1 >> "',TRIM(outFile),'"'
write(char300(5),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'

i=system(TRIM(char300(1)))
i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))
i=system(TRIM(char300(5)))


			
        end if

        open(14,file='../../../TMP_files/vert/anom_'//it//ver//'.xyz')
        open(15,file='../../../TMP_files/vert/abs_'//it//ver//'.xyz')
        do iz=nzsec,1,-1
	        zcur=zmin+iz*dzsec
	        do ix=1,nxsec
                sss=(ix-1)*dxsec
                write(14,*)sss,-zcur,dvan(ix,iz)+add_perc+aprio_dv(ix,iz)
                write(15,*)sss,-zcur,vabs(ix,iz)
	        end do
        end do
        close(14)
        close(15)

		grd = 'TMP_files/vert/ver_'//ps//it//ver//'.grd'
        open(14,file='../../../'//grd)
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
	        write(14,*)(dvan(ix,iz)+add_perc+aprio_dv(ix,iz),ix=1,nxsec)
        end do
        close(14)
		! call plot vertical sections for vp, vs
		
grdToRead = '../../../TMP_files/vert/ver_'//ps//it//ver//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/IT'//it//'/ver_dv'//ps//ver//'.ps'
cpt = 'scale1.cpt'	
if(ips.eq.1) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"dVp; Iteration: ',it,'; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
else if(ips.eq.2) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"dVs; Iteration: ',it,'; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
end if

grdToRead = '../../../TMP_files/vert/ztr_'//ver//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -Sc0.2c -G0/0/0 -O -K -B >> "',TRIM(outFile),'"'

grdToRead = '../../../TMP_files/vert/topo_'//ver//'.dat'
write(char300(3),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -A -W3,200/70/70 -O -K -B >> "',TRIM(outFile),'"'

write(char300(4),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(1)),' -By+l"P-anomalies, %" -O -Y-1 >> "',TRIM(outFile),'"'
if(ips.eq.2) then
write(char300(4),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(1)),' -By+l"S-anomalies, %" -O -Y-1 >> "',TRIM(outFile),'"'
end if    

write(char300(5),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'

!write(*,*)TRIM(char300(1))
!write(*,*)TRIM(char300(2))
!write(*,*)TRIM(char300(3))
!write(*,*)TRIM(char300(4))
!write(*,*)TRIM(char300(5))

i=system(TRIM(char300(1)))
i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))
i=system(TRIM(char300(5)))

		
        open(14,file='../../../TMP_files/vert/abs_'//ps//it//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
	        write(14,*)(vabs(ix,iz),ix=1,nxsec)
        end do
        close(14)

grdToRead = '../../../TMP_files/vert/abs_'//ps//it//ver//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/IT'//it//'/abs_ver'//ps//ver//'.ps'
cpt = 'scale3.cpt'  
if(ips.eq.2)cpt = 'scale4.cpt'


if(ips.eq.1) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Abs Vp; Iteration: ',it,'; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
else if(ips.eq.2) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Abs Vs; Iteration: ',it,'; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
end if

grdToRead = '../../../TMP_files/vert/ztr_'//ver//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -Sc0.2c -G0/0/0 -O -K -B >> "',TRIM(outFile),'"'

grdToRead = '../../../TMP_files/vert/topo_'//ver//'.dat'
write(char300(3),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -A -W3,200/70/70 -O -K -B >> "',TRIM(outFile),'"'

write(char300(4),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(3)),' -By+l"P-abs. velocity, km/s" -O -Y-1 >> "',TRIM(outFile),'"'
if(ips.eq.2) then
write(char300(4),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(4)),' -By+l"S-abs. velocity, km/s" -O -Y-1 >> "',TRIM(outFile),'"'
end if

write(char300(5),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'

i=system(TRIM(char300(1)))
i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))
i=system(TRIM(char300(5)))
 


541	continue

    end do



5557 continue

    deallocate(dvan,aprio_dv,vvv,vtmp,vabs,vabsp)

end do
close(16)

i=system('rm scale1.cpt')
i=system('rm scale2.cpt')
i=system('rm scale3.cpt')
i=system('rm scale4.cpt')



stop
end
