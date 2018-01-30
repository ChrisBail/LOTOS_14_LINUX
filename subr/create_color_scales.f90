subroutine create_color_scales()
character*20 scale_line,scale_vpvs,scale,scale_vp,scale_vs
integer krrr(100),kggg(100),kbbb(100)


common/scales/dv_min,dv_max,vpvs_min,vpvs_max,vp_min,vp_max,vs_min,vs_max
common/scl_names/scale_line,scale_vpvs,scale_vp,scale_vs


!******************************************************************
open(1,file='../../../COMMON/scales_scl/'//TRIM(scale_line))
read(1,*)a1,a2
ncol=0
276 continue
	read(1,*,end=277)irrr,iggg,ibbb
	ncol=ncol+1
	krrr(ncol)=irrr; kggg(ncol)=iggg; kbbb(ncol)=ibbb
	goto 276
277 continue
close(1)

!write(*,*)' ncol=',ncol,' dv1=',dv_min,' dv2=',dv_max
dcol=(dv_max-dv_min)/ncol
open(11,file='scale1.cpt')
write(11,*)'# blue red 10lvl color pallete table 26/08/2009'
write(11,*)'# COLOR_MODEL = RGB'
do icol=1,ncol-1
	a1=dv_min+(icol-1)*dcol
	a2=dv_min+(icol)*dcol
	!if(icol.eq.1) a1=dv_min*2.
	!if(icol.eq.ncol-1) a2=dv_max*2.
	write(11,'(f6.1,3i5,f6.1,3i5)')a1, krrr(icol),kggg(icol),kbbb(icol),&
	a2, krrr(icol+1),kggg(icol+1),kbbb(icol+1)
end do
write(11,*)'N 125 125 125'
close(11)

!******************************************************************
open(1,file='../../../COMMON/scales_scl/'//TRIM(scale_vpvs))
read(1,*)a1,a2
ncol=0
226 continue
	read(1,*,end=227)irrr,iggg,ibbb
	ncol=ncol+1
	krrr(ncol)=irrr; kggg(ncol)=iggg; kbbb(ncol)=ibbb
	goto 226
227 continue
close(1)

!write(*,*)' ncol=',ncol,' vpvs1=',vpvs_min,' vpvs2=',vpvs_max
dcol=(vpvs_max-vpvs_min)/ncol
open(11,file='scale2.cpt')
write(11,*)'# blue red 10lvl color pallete table 26/08/2009'
write(11,*)'# COLOR_MODEL = RGB'
do icol=1,ncol-1
	a1=vpvs_min+(icol-1)*dcol
	a2=vpvs_min+(icol)*dcol
	!if(icol.eq.1) a1=vpvs_min-dcol*10
	!if(icol.eq.ncol-1) a2=vpvs_max+dcol*10
	write(11,'(f5.2,3i4,f5.2,3i4)')a1, krrr(icol),kggg(icol),kbbb(icol),&
	a2, krrr(icol+1),kggg(icol+1),kbbb(icol+1)
end do
write(11,*)'N 125 125 125'
close(11)


!******************************************************************
open(1,file='../../../COMMON/scales_scl/'//TRIM(scale_vp))
read(1,*)a1,a2
ncol=0
236 continue
	read(1,*,end=237)irrr,iggg,ibbb
	ncol=ncol+1
	krrr(ncol)=irrr; kggg(ncol)=iggg; kbbb(ncol)=ibbb
	goto 236
237 continue
close(1)

!write(*,*)' ncol=',ncol,' vp1=',vp_min,' vp2=',vp_max
dcol=(vp_max-vp_min)/ncol
open(11,file='scale3.cpt')
write(11,*)'# blue red 10lvl color pallete table 26/08/2009'
write(11,*)'# COLOR_MODEL = RGB'
do icol=1,ncol-1
	a1=vp_min+(icol-1)*dcol
	a2=vp_min+(icol)*dcol
	!if(icol.eq.1) a1=vp_min-dcol*10
	!if(icol.eq.ncol-1) a2=vp_max+dcol*10.
	write(11,'(f5.2,3i4,f5.2,3i4)')a1, krrr(icol),kggg(icol),kbbb(icol),&
	a2, krrr(icol+1),kggg(icol+1),kbbb(icol+1)
end do
write(11,*)'N 125 125 125'
close(11)

!******************************************************************
open(1,file='../../../COMMON/scales_scl/'//TRIM(scale_vs))
read(1,*)a1,a2
ncol=0
246 continue
	read(1,*,end=247)irrr,iggg,ibbb
	ncol=ncol+1
	krrr(ncol)=irrr; kggg(ncol)=iggg; kbbb(ncol)=ibbb
	goto 246
247 continue
close(1)

!write(*,*)' ncol=',ncol,' vs1=',vs_min,' vs2=',vs_max
dcol=(vs_max-vs_min)/ncol
open(11,file='scale4.cpt')
write(11,*)'# blue red 10lvl color pallete table 26/08/2009'
write(11,*)'# COLOR_MODEL = RGB'
do icol=1,ncol-1
	a1=vs_min+(icol-1)*dcol
	a2=vs_min+(icol)*dcol
	!if(icol.eq.1) a1=vs_min-dcol*10.
	!if(icol.eq.ncol-1) a2=vs_max+dcol*10
	write(11,'(f6.2,3i5,f6.2,3i5)')a1, krrr(icol),kggg(icol),kbbb(icol),&
	a2, krrr(icol+1),kggg(icol+1),kbbb(icol+1)
end do
write(11,*)'N 125 125 125'
close(11)

return
end