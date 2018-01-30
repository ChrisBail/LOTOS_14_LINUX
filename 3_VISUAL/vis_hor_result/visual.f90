! HORIZONTAL !!!
! NODES !!!!!!!

character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*2 lv
character*1 ps ,rm,it,ch1
character*20 scale_line, scale_line2,char_x
character*50 grd, grdToRead, outFile, cpt
character*300 char300(10)
character*10 char10,ch_lim(4),ch_col(4),ch_hor


allocatable dvan(:,:),vvv(:,:),v1tmp(:,:),v2tmp(:,:),vabs(:,:,:)
real vaver(2,20)
integer nrps(2)
real fzzt(10000,100),tzzt(10000,100),zzzt(10000,100)
integer nzzt(100),line1_rgb(3),line2_rgb(3),kdot_rgb(3)
integer*2 izzz

common/pi/pi,per
!common/keys/key_ft1_xy2

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
rz=6371.

w_limit=0.2
igr=1
kod_av_bias=0
kod_apriori=0
ind_srce=1
ch1='/'

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)iter		
close(1)
write(it,'(i1)')iter

write(*,*)
write(*,*)' ***********************************************'
write(*,*)' VISUALISATION in horizontal sections: '
write(*,*)' ar=',ar,' md=',md,' iter=',iter


i=system('mkdir -p ../../../TMP_files/hor')
i=system('mkdir -p ../../../PICS/'//ar//'/'//md//'/IT'//it)

call read_param(ar,md)
call create_color_scales()
call read_topo(ar)
call read_vref(ar,md)

!write(*,*) nlev,hlev(1),fmap1,fmap2,dfmap,tmap1,tmap2,dtmap,smaxx,ismth,size_hor

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


rsmth=ismth+0.5
nfmap=int_best((fmap2-fmap1)/dfmap+1.)
ntmap=int_best((tmap2-tmap1)/dtmap+1.)
write(*,*)' nfmap=',nfmap,' ntmap=',ntmap
allocate(dvan(nfmap,ntmap),vvv(nfmap,ntmap),v1tmp(nfmap,ntmap),v2tmp(nfmap,ntmap))
allocate(vabs(2,nfmap,ntmap))

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray1.dat')
read(1,*) nrps(1),nrps(2)
close(1)

nzzt=0
nzt=0
   open(1,file='../../../DATA/'//ar//'/'//md//'/data/srces'//it//'.dat')
872	read(1,*,end=871)fzt,tzt,zzt
!call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
nzt=nzt+1
do ilev=1,nlev
	if(ilev.eq.1) then
		z1=-10
		z2=(hlev(1)+hlev(2))/2
	else if(ilev.eq.nlev) then
		z1=(hlev(nlev-1)+hlev(nlev))/2
		z2=hlev(nlev) + (hlev(nlev)-hlev(nlev-1))/2
	else 
		z1=(hlev(ilev-1)+hlev(ilev))/2
		z2=(hlev(ilev+1)+hlev(ilev))/2
	end if
	if((zzt-z1)*(zzt-z2).le.0) goto 995
end do
goto 872

995 continue
!write(*,*)' zzt=',zzt,' ilev=',ilev,' z1=',z1,' z2=',z2
nzzt(ilev)=nzzt(ilev)+1
fzzt(nzzt(ilev),ilev)=fzt
tzzt(nzzt(ilev),ilev)=tzt
zzzt(nzzt(ilev),ilev)=zzt

goto 872
871 close(1)
DO ilev=1,nlev
	write(lv,'(i2)')ilev
	open(11,file='../../../TMP_files/hor/ztr'//lv//'.dat')
	write(*,*)' ilev=',ilev,' nzzt=',nzzt(ilev)
	do izzt=1,nzzt(ilev)
		write(11,*)fzzt(izzt,ilev),tzzt(izzt,ilev),zzzt(izzt,ilev)
	end do
	close(11)
end do


vaver=0
DO ilev=1,nlev
	zzz=hlev(ilev)
	write(lv,'(i2)')ilev
	write(*,*)' ilev=',ilev,' zzz=',zzz

	do ips=1,2
		v0=vrefmod(zzz,ips)
		if(nrps(ips).eq.0) cycle
		write(ps,'(i1)')ips
		dvan=0
		vvv=0
		do igr=ngr1,ngr2
			call prepare_model_v(ar,md,ips,iter,igr)

			do itet=1,ntmap
				ttt=(itet-1)*dtmap+tmap1+tet0
				!write(*,*)' itet=',itet,' ttt=',ttt
				!ttt=tet0
				do ifi=1,nfmap
					fff=(ifi-1)*dfmap+fmap1+fi0
					!fff=fi0
					!fff=110.87
					dv=0
					www=0
					!write(*,*)' smaxx=',smaxx
					call dv_1_grid_v(fff,ttt,zzz,smaxx, dv,www)
					zlim_up=h_lim(fff,ttt)
					if (zzz.lt.zlim_up) www=0
					!write(*,*)' dv=',dv,' www=',www
					dvproc=100*dv/v0
					dvan(ifi,itet)=dvan(ifi,itet)+dvproc*www
					vvv(ifi,itet)=vvv(ifi,itet)+www
					!if(itet.eq.101) write(*,*)dv,www
				end do
			end do
		end do

        nonzer=0
        do ifi=1,nfmap
            do itet=1,ntmap
                vanm=-999.
                vabs(ips,ifi,itet)=-999
                if (vvv(ifi,itet).gt.w_limit) then
                    vanm=dvan(ifi,itet)/vvv(ifi,itet)
                    vabs(ips,ifi,itet)=v0*(1+0.01*vanm)
                    !write(*,*)ips,ifi,itet,vabs(ips,ifi,itet)
                    nonzer=nonzer+1
                    vaver(ips,ilev)=vaver(ips,ilev)+vabs(ips,ifi,itet)
                end if
                v1tmp(ifi,itet)=vanm
            end do
            !pause
        end do
        vaver(ips,ilev)=vaver(ips,ilev)/nonzer
        write(*,*)' ips=',ips,' vaver=',vaver(ips,ilev),' nonzer=',nonzer


		do ifi=1,nfmap
			do itet=1,ntmap
				if(vvv(ifi,itet).lt.w_limit) cycle
				vanm=0.
				iv=0
				do iff=-ismth,ismth
					if (ifi+iff.lt.1) cycle
					if (ifi+iff.gt.nfmap) cycle
					do itt=-ismth,ismth
						if (itet+itt.lt.1) cycle
						if (itet+itt.gt.ntmap) cycle
						if(vvv(ifi+iff,itet+itt).lt.w_limit) cycle
						rr=iff*iff+itt*itt
						r=sqrt(rr)
						if(r.gt.rsmth) cycle
						iv=iv+1
						vanm=vanm+v1tmp(ifi+iff,itet+itt)
					end do
				end do
				v2tmp(ifi,itet)=vanm/iv
			end do
		end do

		aver=0
		naver=0
		do ifi=1,nfmap
			do itet=1,ntmap
				vanom=-999
				if(vvv(ifi,itet).gt.w_limit) then
					vanom=v2tmp(ifi,itet)
					aver=aver+vanom
					naver=naver+1
				end if
				v2tmp(ifi,itet)=vanom
				if(kod_apriori.eq.1) then
					ttt=(itet-1)*dtmap+tmap1+tet0
					fff=(ifi-1)*dfmap+fmap1+fi0

                                    if(key_ft1_xy2.eq.1) then
                                        call SFDEC(fff,ttt,0.,xxx,yyy,Z,fi0,tet0)
                                    else
                                        xxx=fff
                                        yyy=ttt
                                    end if

					dv_aprio = vert_anom(xxx,yyy,zzz,ips)
					v2tmp(ifi,itet)=v2tmp(ifi,itet)+dv_aprio
				end if
				!if(itet.eq.101) write(*,*)vanom,vvv(ifi,itet)
			end do
		end do
		!pause
		if(naver.gt.0) then
			aver=aver/naver
		end if
		if(kod_av_bias.eq.1) v2tmp=v2tmp-aver


                open(14,file='../../../TMP_files/hor/dv'//ps//it//lv//'.xyz')
                do ifi=1,nfmap
                    fff=(ifi-1)*dfmap+fmap1+fi0
                    do itet=1,ntmap
                        ttt=(itet-1)*dtmap+tmap1+tet0
                        write(14,*)fff,ttt,v2tmp(ifi,itet)
                    end do
                end do
                close(14)
		
		grd = 'TMP_files/hor/dv'//ps//it//lv//'.grd'
		open(14,file='../../../'//grd)
		write(14,'(a4)')dsaa
		write(14,*)nfmap,ntmap
		write(14,*)fmap1+fi0,fmap2+fi0
		write(14,*)tmap1+tet0,tmap2+tet0
		write(14,*)-999,999
		do itet=1,ntmap
			write(14,*)(v2tmp(ifi,itet),ifi=1,nfmap)
		end do
		close(14)
		
		! call plot horizontal sections for vp, vs
		
		grdToRead = '../../../TMP_files/hor/dv'//ps//it//lv//'.grd'
		outFile = '../../../PICS/'//ar//'/'//md//'/IT'//it//'/hor_dv'//ps//lv//'.ps'
		cpt = 'scale1.cpt'	
		!i=system('gmt grdimage "'//TRIM(grdToRead)//'" -P -K -Ba -Jm2i -C"'//TRIM(cpt)//'" > "'//TRIM(outFile)//'"')
		!pause
		!write(char300(1),*)'gmt pscontour "',TRIM(grdToRead),'" -P -K -Ba -JM',TRIM(char_x),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'

call real2char(zzz,1,char10)
if(ips.eq.1)then
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWESn+t"dVp; Iteration: ',it,'; Depth:',TRIM(char10),' km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JM',TRIM(ch_hor),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
else
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWESn+t"dVs; Iteration: ',it,'; Depth:',TRIM(char10),' km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JM',TRIM(ch_hor),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'
end if

grdToRead = '../../../TMP_files/hor/ztr'//lv//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JM -R -Sc0.1c -G0/0/0 -O -K -B >> "',TRIM(outFile),'"'

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



	end do

	v1tmp=-999
	do ifi=1,nfmap
		do itet=1,ntmap
			if (abs(vabs(1,ifi,itet)).gt.900..or.abs(vabs(2,ifi,itet)).gt.900.) cycle
			vpvs=vabs(1,ifi,itet)/vabs(2,ifi,itet)
			v1tmp(ifi,itet)=vpvs
			!write(*,*)vabs(1,ifi,itet),vabs(2,ifi,itet),vpvs
		end do
		!pause
	end do


	open(14,file='../../../TMP_files/hor/vpvs'//it//lv//'.grd')
	write(14,'(a4)')dsaa
	write(14,*)nfmap,ntmap
	write(14,*)fmap1+fi0,fmap2+fi0
	write(14,*)tmap1+tet0,tmap2+tet0
	write(14,*)-999,999
	do itet=1,ntmap
		write(14,*)(v1tmp(ifi,itet),ifi=1,nfmap)
	end do
	close(14)

	! call plot horizontal sections for vp/vs relationship
	grdToRead = '../../../TMP_files/hor/vpvs'//it//lv//'.grd'
	outFile = '../../../PICS/'//ar//'/'//md//'/IT'//it//'/hor_vpvs'//lv//'.ps'
	cpt = 'scale2.cpt'	

call real2char(zzz,1,char10)
write(char300(1),*)'gmt grdimage "',TRIM(grdToRead),'" -K -P -Y4 -Ba',&
    ' -BWESn+t"Vp/Vs; Iteration: ',it,'; Depth:',TRIM(char10),' km"',&
    ' -R',TRIM(ch_lim(1)),ch1,TRIM(ch_lim(2)),ch1,TRIM(ch_lim(3)),ch1,TRIM(ch_lim(4)),&
    ' -JM',TRIM(ch_hor),' -C"',TRIM(cpt),'" > "',TRIM(outFile),'"'

grdToRead = '../../../TMP_files/hor/ztr'//lv//'.dat'
write(char300(2),*)'gmt psxy "',TRIM(grdToRead),'" -JM -R -Sc0.1c -G0/0/0 -O -K -B >> "',TRIM(outFile),'"'

write(char300(3),*)'gmt psscale -D5.0c/-0.8c/10c/0.5ch -C"',TRIM(cpt),&
    '" -Bx',TRIM(ch_col(2)),' -By+l"P-anomalies, %" -O -Y-1 >> "',TRIM(outFile),'"'

write(char300(4),*)'psconvert "',TRIM(outFile),'" -A -E300 -Tg'


	i=system(TRIM(char300(1)))
	i=system(TRIM(char300(2)))
	i=system(TRIM(char300(3)))
	i=system(TRIM(char300(4)))


end do

open(11,file='../../../DATA/'//ar//'/'//md//'/data/vaver'//it//'.dat')
do ilev=1,nlev
    write(11,*)hlev(ilev),vaver(1,ilev),vaver(2,ilev)
end do
close(11)

i=system('rm scale1.cpt')
i=system('rm scale2.cpt')
i=system('rm scale3.cpt')
i=system('rm scale4.cpt')

stop
end
