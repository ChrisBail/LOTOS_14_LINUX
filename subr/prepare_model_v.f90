subroutine prepare_model_v(ar,md,ips,iter,igr)

character*1 gr,ps,it
character*8 ar,md,line

common/pi/pi,per
common/center/fi0,tet0
common/orient/nornt0,ornt0(10)
common/grid_cur/nornt,ornt(4),sinal,cosal,&
		nlev,ylev(300),ntop(300),nobr,&
		xtop(20000,300), ztop(20000,300),n_pop(20000,300),&
		dv_mod(100000),vab_mod(100000)
common/limits/xlimm1,xlimm2,ylimm1,ylimm2,zlim1,zlim2
common/grid/xlim1,xlim2,dxpl,ylim1,ylim2,dypl,zlim10,zlim20,dzpl,plotmin,plotmax

nornt=nornt0
do i=1,nornt0
	ornt(i)=ornt0(i)
end do
zlim1=zlim10; zlim2=zlim20

write(ps,'(i1)')ips
write(it,'(i1)')iter
write(gr,'(i1)')igr

!******************************************************************

!write(*,*)' fi0=',fi0,' tet0=',tet0
orient=ornt(igr)
sinal=sin(orient*per)
cosal=cos(orient*per)
!write(*,*)' orient=',orient

!******************************************************************
xlimm1=xlim1+0.5*dxpl
xlimm2=xlim2-0.5*dxpl
ylimm1=ylim1+0.5*dypl
ylimm2=ylim2-0.5*dypl


open(1,file='../../../DATA/'//ar//'/'//md//'/data/levinfo'//ps//gr//'.dat')
i=0
722 i=i+1
	read(1,*,end=721)n,ylev(i)
	!write(*,*)i,n,ylev(i)
	goto 722
721 nlev=i-1
close(1)

ntop=0
open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr'//ps//gr//'.dat')

do n=1,nlev
	read(1,*,end=556)ntop(n)
	!write(*,*)' y=',ylev(n),' ntop=',ntop(n)
	do i=1,ntop(n)
		read(1,*)xtop(i,n),ztop(i,n),n_pop(i,n)
	end do
end do
556		close(1)

!open(21,file='grid2tmp.dat')
!do i=1,ntop(2)
!	x=xtop(i,2)
!	y=ytop(i,2)
!	zz=0
!	ind=0
!	call decsf(x,y,zz,ind,fi0,tet0,FI,TET,h)
!	write(21,*)fi,tet
!end do
!close(21)
!pause

open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr'//ps//gr//'.dat')
read(1,*) nobr
if(nobr.gt.100000) then
	write(*,*)' nobr=',nobr
	write(*,*)' One should increase size of dv_mod'
	pause
end if
close(1)

if(ips.eq.1) then
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/vel_p_'//it//gr//'.dat')
else
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/vel_s_'//it//gr//'.dat')
end if

!write(*,*)' nobr=',nobr
do i=1,nobr
	read(1,*)  dv_mod(i),vab_mod(i)
	!dv_mod(i)=100*dv_mod(i)/vab_mod(i)
end do


close(1)
return
end

