! Visualization of synthetic model in VERTICAL sections
character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*20 scale_line,scale_vpvs,scale,scale_vp,scale_vs,char_x,char_z,char_lim(4)
character*1 ps,ch1
character*2 ver,lv
character*10 char10,ch_lim(4),ch_col(4)
character*50 grd, grdToRead, outFile, cpt
character*300 char300(10)

real fmark(200,20),tmark(200,20),smark(200,20)
integer nmark(100)
integer kdot_rgb(3)

common/pi/pi,per

common/general/key_1real_2syn,VPSX_key,koe,kref,key_ft1_xy2,key_true1,key_flat1
common/orient/nornt,ornt(10)
common/center/fi0,tet0
common/visual_hor/ nlev,hlev(20),fmap1,fmap2,dfmap,tmap1,tmap2,dtmap,smaxx,ismth_h
common/visual_ver/ nver,fia0(20),teta0(20),fib0(20),tetb0(20),dist_from_sec_event,&
    dxsec,zmin,zmax,dzsec,dsmark,dismax,ismth_v,size_x,size_z
common/scales/dv_min,dv_max,vpvs_min,vpvs_max,vp_min,vp_max,vs_min,vs_max

allocatable v_ini(:,:),v_abs_ini(:,:,:)

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0
ch1='/'


i=system('mkdir ../../../TMP_files/vert')


open(1,file='../../../model.dat')
read(1,'(a8)')ar	! synthetic model
read(1,'(a8)')md	! synthetic model
close(1)

write(*,*)
write(*,*)' ***********************************************'
write(*,*)' SYNTHETIC in vertical section: '
write(*,*)' ar=',ar,' md=',md

call read_param(ar,md)
call create_color_scales()
call read_vref(ar,md)
call read_topo(ar)
call read_anom(ar,md)

aaa=(dv_max-dv_min)/10.; call real2char(aaa,0,ch_col(1))
!write(*,*)' dv_max=',dv_max,' dv_min=',dv_min
!write(*,*)' aaa=',aaa,' ch_col(1)=',ch_col(1)
aaa=(vpvs_max-vpvs_min)/10.; call real2char(aaa,2,ch_col(2))
aaa=(vp_max-vp_min)/10.; call real2char(aaa,2,ch_col(3))
aaa=(vs_max-vs_min)/10.; call real2char(aaa,2,ch_col(4))


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
    npix_x = int(npix_y * dist/(zmax-zmin))
    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist
    nxsec=dist/dxsec+1
    dxsec=dist/(nxsec-1)
    nzsec=(zmax-zmin)/dzsec+1
    dzsec=(zmax-zmin)/(nzsec-1)

    call real2char(0.,0,ch_lim(1))
    call real2char(dist,1,ch_lim(2))
    call real2char(-zmax,0,ch_lim(3))
    call real2char(-zmin,0,ch_lim(4))

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
    if(key_ft1_xy2.eq.1) then
        fmark(iver,imark)=fib
        tmark(iver,imark)=tetb
    else
        fmark(iver,imark)=xb
        tmark(iver,imark)=yb
    end if
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
    open(11,file='../../../TMP_files/vert/topo_'//ver//'.bln')
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
        topo=h_lim(fff,ttt)
        write(11,*)sss,-topo
    end do
    close(11)


    allocate (v_ini(nxsec,nzsec),v_abs_ini(2,nxsec,nzsec))

    write(*,*)' section:',ver,' dist=',dist,' nxsec=',nxsec,' nzsec=',nzsec
    v_ini=0
    v_abs_ini=0
	
    do ips=1,2
        write(ps,'(i1)')ips

        do ix=1,nxsec
            sss=(ix-1)*dxsec
            !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
            !sss=40
            !if(mod(ix,10).eq.0)write(*,*)' ix=',ix,' sss=',sss
            xcur=xa+((xb-xa)/dist)*sss
            ycur=ya+((yb-ya)/dist)*sss
            if(key_ft1_xy2.eq.1) then
                call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
            else
                fff=xcur
                ttt=ycur
            end if

            relief=h_lim(fff,ttt)



            !write(*,*)' xcur=',xcur,' ycur=',ycur
            !write(*,*)' fi=',fi,' tet=',tet,' relief=',relief
            do iz=1,nzsec
                zcur=zmin+iz*dzsec
                v0=vrefmod(zcur,ips)
                dv=anomaly(xcur,ycur,zcur,ips)

                if (zcur.lt.relief) then
                    v_ini(ix,iz)=-999
                    v_abs_ini(ips,ix,iz)=-999
                else
                    v_ini(ix,iz)=dv
                    v_abs_ini(ips,ix,iz)=v0*(1+dv/100)
                end if
            end do

        end do

       open(14,file='../../../TMP_files/vert/syn_dv'//ps//ver//'.xyz')
       open(15,file='../../../TMP_files/vert/syn_abs'//ps//ver//'.xyz')
       do ix=1,nxsec
            sss=(ix-1)*dxsec
            do iz=1,nzsec
                zcur=zmin+iz*dzsec
                write(14,*)sss,-zcur,v_ini(ix,iz)
                write(15,*)sss,-zcur,v_abs_ini(ips,ix,iz)
           end do
        end do
        close(14)
        close(15)




        open(14,file='../../../TMP_files/vert/syn_dv'//ps//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
            write(14,*)(v_ini(ix,iz),ix=1,nxsec)
        end do
        close(14)


grdToRead = '../../../TMP_files/vert/syn_dv'//ps//ver//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/SYN/ver_dv'//ps//ver//'.ps'
cpt = 'scale1.cpt'  
if(ips.eq.1) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Synthetic dVp; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
else if(ips.eq.2) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Synthetic dVs; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
end if

grdToRead = '../../../TMP_files/vert/topo_'//ver//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -A -W3,200/70/70 -O -K -B >> "',TRIM(outFile),'"'

write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(1)),' -By+l"P-anomalies, %" -O -Y-1 >> "',TRIM(outFile),'"'
if(ips.eq.2) then
write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(1)),' -By+l"S-anomalies, %" -O -Y-1 >> "',TRIM(outFile),'"'
end if    

write(char300(4),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'

!write(*,*)TRIM(char300(1))
!write(*,*)TRIM(char300(2))
!write(*,*)TRIM(char300(3))
!write(*,*)TRIM(char300(4))


i=system(TRIM(char300(1)))
i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))


        open(14,file='../../../TMP_files/vert/syn_abs'//ps//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
            write(14,*)(v_abs_ini(ips,ix,iz),ix=1,nxsec)
        end do
        close(14)

grdToRead = '../../../TMP_files/vert/syn_abs'//ps//ver//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/SYN/abs_ver'//ps//ver//'.ps'
cpt = 'scale3.cpt'  
if(ips.eq.2)cpt = 'scale4.cpt'


if(ips.eq.1) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Synth Abs Vp; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
else if(ips.eq.2) then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Synth Abs Vs; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
end if

grdToRead = '../../../TMP_files/vert/topo_'//ver//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -A -W3,200/70/70 -O -K -B >> "',TRIM(outFile),'"'

write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(3)),' -By+l"P-abs. velocity, km/s" -O -Y-1 >> "',TRIM(outFile),'"'
if(ips.eq.2) then
write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(4)),' -By+l"S-abs. velocity, km/s" -O -Y-1 >> "',TRIM(outFile),'"'
end if

write(char300(4),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'

!write(*,*)TRIM(char300(1))
!write(*,*)TRIM(char300(2))
!write(*,*)TRIM(char300(3))
!write(*,*)TRIM(char300(4))

i=system(TRIM(char300(1)))
i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))


441	continue

    end do      ! ips=1,2

    v_ini=0

    !write(*,*)' only  S:',v_abs_ini(1,135,183),v_abs_ini(2,135,183),v_abs_ini(1,135,183)/v_abs_ini(2,135,183)
    !write(*,*)' P and S:',v_abs_ini(1,151,133),v_abs_ini(2,151,133),v_abs_ini(1,151,133)/v_abs_ini(2,151,133)
    !pause


    do ix=1,nxsec
        sss=(ix-1)*dxsec

        do iz=1,nzsec
            zcur=zmin+iz*dzsec
            !write(*,*)sss,zcur,v_abs_ini(1,ix,iz),v_abs_ini(2,ix,iz),v_abs_ini(1,ix,iz)/v_abs_ini(2,ix,iz)
            v_ini(ix,iz)=v_abs_ini(1,ix,iz)/v_abs_ini(2,ix,iz)
        end do
        !pause
    end do

    open(14,file='../../../TMP_files/vert/syn_vpvs'//ver//'.xyz')
    do ix=1,nxsec
        sss=(ix-1)*dxsec
        do iz=1,nzsec
            zcur=zmin+iz*dzsec
            write(14,*)sss,-zcur,v_ini(ix,iz)
        end do
    end do
    close(14)

    open(14,file='../../../TMP_files/vert/syn_vpvs'//ver//'.grd')
    write(14,'(a4)')dsaa
    write(14,*)nxsec,nzsec
    write(14,*)0,dist
    write(14,*)-zmax,-zmin
    write(14,*)-9999,999
    do iz=nzsec,1,-1
        write(14,*)(v_ini(ix,iz),ix=1,nxsec)
    end do
    close(14)

grdToRead = '../../../TMP_files/vert/syn_vpvs'//ver//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/SYN/ver_vpvs'//ver//'.ps'
cpt = 'scale2.cpt'  

write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWeSn+t"Synth Vp/Vs ratio; Section:',ver,'"',&
    ' -Bx+l"distance, km"',&
    ' -By+l"depth, km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JX',TRIM(char_x),ch1,TRIM(char_z),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'


grdToRead = '../../../TMP_files/vert/topo_'//ver//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JX -R -A -W3,200/70/70 -O -K -B >> "',TRIM(outFile),'"'

write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(2)),' -By+l"Vp/Vs ratio" -O -Y-1 >> "',TRIM(outFile),'"'
write(char300(4),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'

!write(*,*)TRIM(char300(1))
!write(*,*)TRIM(char300(2))
!write(*,*)TRIM(char300(3))
!write(*,*)TRIM(char300(4))


i=system(TRIM(char300(1)))
i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))


    442		continue

    deallocate(v_ini,v_abs_ini)

end do      ! iver=1,nver
close(16)

i=system('rm scale1.cpt')
i=system('rm scale2.cpt')
i=system('rm scale3.cpt')
i=system('rm scale4.cpt')


stop
end