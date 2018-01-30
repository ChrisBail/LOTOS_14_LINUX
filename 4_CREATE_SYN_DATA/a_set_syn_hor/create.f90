character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*1 ps,ch1
character*2 lv
character*20 scale_dv,scale_vpvs
integer line1_rgb(3),line2_rgb(3),kdot_rgb(3)

character*20 scale_line, scale_line2,char_x
character*70 grd, grdToRead, outFile, cpt
character*300 char300(10)
character*10 char10,ch_lim(4),ch_col(4),ch_hor


allocatable v_ini(:,:),v_abs(:,:,:)
common/pi/pi,per
common/keys/key_ft1_xy2_

common/general/key_1real_2syn,VPSX_key,koe,kref,key_ft1_xy2,key_true1,key_flat1
common/orient/nornt,ornt(10)
common/center/fi0,tet0
common/visual_hor/ nlev,hlev(20),fmap1,fmap2,dfmap,tmap1,tmap2,dtmap,smaxx,ismth,size_hor
common/scales/dv_min,dv_max,vpvs_min,vpvs_max,vp_min,vp_max,vs_min,vs_max
common/visual_ver/ nver,fia0(20),teta0(20),fib0(20),tetb0(20),dist_from_sec_event,&
    dxsec,zmin,zmax,dzsec,dsmark,dismax,ismth_v,size_x,size_z


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

ch1='/'

open(1,file='../../../model.dat')
read(1,'(a8)')ar	! synthetic model
read(1,'(a8)')md	! synthetic model
close(1)

i=system('mkdir ../../../TMP_files/hor')
i=system('mkdir -p ../../../PICS/'//ar//'/'//md//'/SYN')


write(*,*)' Synthetic model in horizontal sections'
write(*,*)' ar=',ar,' md=',md

call read_param(ar,md)
call create_color_scales()
call read_topo(ar)
call read_vref(ar,md)
call read_anom(ar,md)

key_ft1_xy2_=key_ft1_xy2
dfmap=dfmap/2.
dtmap=dtmap/2.
nfmap=int_best((fmap2-fmap1)/dfmap+1.)
ntmap=int_best((tmap2-tmap1)/dtmap+1.)
write(*,*)' nfmap=',nfmap,' ntmap=',ntmap

ngr1=1
ngr2=nornt

aaa=(dv_max-dv_min)/10.; call real2char(aaa,0,ch_col(1))
!write(*,*)' dv_max=',dv_max,' dv_min=',dv_min
!write(*,*)' aaa=',aaa,' ch_col(1)=',ch_col(1)
aaa=(vpvs_max-vpvs_min)/10.; call real2char(aaa,2,ch_col(2))
aaa=(vp_max-vp_min)/10.; call real2char(aaa,2,ch_col(3))
aaa=(vs_max-vs_min)/10.; call real2char(aaa,2,ch_col(4))
call real2char(size_hor,1,ch_hor)
!write(*,*)size_hor,'______',TRIM(ch_hor),'_______'

call real2char(fi0+fmap1,2,ch_lim(1))
call real2char(fi0+fmap2,2,ch_lim(2))
call real2char(tet0+tmap1,2,ch_lim(3))
call real2char(tet0+tmap2,2,ch_lim(4))




allocate(v_ini(nfmap,ntmap),v_abs(nfmap,ntmap,2))

DO ilev=1,nlev
    zzz=hlev(ilev)
    write(lv,'(i2)')ilev
    write(*,*)' ilev=',ilev,' zzz=',zzz
    v_abs=0
    do ips=1,2
        write(ps,'(i1)')ips
        vref=vrefmod(zzz,ips)

        v_ini=0

        do itet=1,ntmap
            ttt=(itet-1)*dtmap+tmap1+tet0
            !write(*,*)' itet=',itet,' ttt=',ttt
            !ttt=-7.35
            do ifi=1,nfmap
                fff=(ifi-1)*dfmap+fmap1+fi0
                relief=h_lim(fff,ttt)

                !fff=-60; ttt=13

                if (zzz.lt.relief) then
                    v_ini(ifi,itet)=-999
                    cycle
                end if
                if(key_ft1_xy2.eq.1) then
                    call SFDEC(fff,ttt,0.,xxx,yyy,Z,fi0,tet0)
                else
                    xxx=fff
                    yyy=ttt
                end if
                dv=anomaly(xxx,yyy,zzz,ips)
                !write(*,*)' fi=',fff,' tet=',ttt,' dv=',dv
                !write(*,*)' x=',xxx,' y=',yyy,' z=',zzz,' dv=',dv
                !pause

                if (zzz.lt.h_lim(fff,ttt)) dv=0
                v_ini(ifi,itet)=dv
                v_abs(ifi,itet,ips)=vref*(1+dv/100)


                !pause
                !pause
            end do
        end do
  
       open(14,file='../../../TMP_files/hor/syn'//ps//'_'//lv//'.xyz')
       do itet=1,ntmap
            ttt=(itet-1)*dtmap+tmap1+tet0
            do ifi=1,nfmap
                fff=(ifi-1)*dfmap+fmap1+fi0
                write(14,*)fff,ttt,v_ini(ifi,itet)
            end do
        end do
        close(14)


        open(14,file='../../../TMP_files/hor/syn'//ps//'_'//lv//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nfmap,ntmap
        write(14,*)fmap1+fi0,fmap2+fi0
        write(14,*)tmap1+tet0,tmap2+tet0
        write(14,*)-999,999
        do itet=1,ntmap
            write(14,*)(v_ini(ifi,itet),ifi=1,nfmap)
        end do
        close(14)


grdToRead = '../../../TMP_files/hor/syn'//ps//'_'//lv//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/SYN/syn_hor_dv'//ps//lv//'.ps'
cpt = 'scale1.cpt'  

call real2char(zzz,1,char10)
if(ips.eq.1)then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWESn+t"Synthetic dVp; Depth:',TRIM(char10),' km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JM',TRIM(ch_hor),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
else
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWESn+t"Synthetic dVs; Depth:',TRIM(char10),' km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JM',TRIM(ch_hor),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
end if

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
!i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))



    end do


! IMAGE THE VP/VS RATIO:


    open(14,file='../../../TMP_files/hor/syn_vpvs_'//lv//'.xyz')
    do itet=1,ntmap
        ttt=(itet-1)*dtmap+tmap1+tet0
        do ifi=1,nfmap
            fff=(ifi-1)*dfmap+fmap1+fi0
            write(14,*)fff,ttt,v_abs(ifi,itet,1)/v_abs(ifi,itet,2)
        end do
    end do
    close(14)



    open(14,file='../../../TMP_files/hor/syn_vpvs_'//lv//'.grd')
    write(14,'(a4)')dsaa
    write(14,*)nfmap,ntmap
    write(14,*)fmap1+fi0,fmap2+fi0
    write(14,*)tmap1+tet0,tmap2+tet0
    write(14,*)-999,999
    do itet=1,ntmap
        write(14,*)(v_abs(ifi,itet,1)/v_abs(ifi,itet,2),ifi=1,nfmap)
    end do
    close(14)

! call plot horizontal sections for vp/vs relationship
grdToRead = '../../../TMP_files/hor/syn_vpvs_'//lv//'.grd'
outFile = '../../../PICS/'//ar//'/'//md//'/SYN/syn_hor_vpvs'//lv//'.ps'
cpt = 'scale2.cpt'  

call real2char(zzz,1,char10)
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWESn+t"Synthetic Vp/Vs; Depth:',TRIM(char10),' km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JM',TRIM(ch_hor),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'

write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(2)),' -By+l"P-anomalies, %" -O -Y-1 >> "',TRIM(outFile),'"'

write(char300(4),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'


i=system(TRIM(char300(1)))
!i=system(TRIM(char300(2)))
i=system(TRIM(char300(3)))
i=system(TRIM(char300(4)))



end do

i=system('rm scale1.cpt')
i=system('rm scale2.cpt')
i=system('rm scale3.cpt')
i=system('rm scale4.cpt')


stop
end