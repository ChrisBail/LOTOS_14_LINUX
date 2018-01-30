character*8 ar,ar_all(100),md,md_all(100),line
character*1 it
integer kod_loc(100),kod_iter(100),kod_oe(100),kod_rf(100)

common/general/key_1real_2syn,VPSX_key,koe,kref,key_ft1_xy2,key_true1,&
	key_flat1,key_table1_line2
common/orient/nornt,ornt(10)


i=system('mkdir -p ../../../TMP_files/tmp')
i=system('mkdir -p ../../../TMP_files/hor')
i=system('mkdir -p ../../../TMP_files/rays')
i=system('mkdir -p ../../../TMP_files/1D_mod')
i=system('mkdir -p ../../../TMP_files/vert')

open(1, file='../../../all_areas.dat')
do i=1,4
	read(1,*)
end do
do i=1,100
	read(1,'(a8,1x,a8,1x,i1)',end=7)ar_all(i),md_all(i),kod_iter(i)
end do	
7 close(1)
n_ar=i-1
write(*,*)' n_ar=',n_ar

do iar=1,n_ar
    ar=ar_all(iar)
    md=md_all(iar)
    niter=kod_iter(iar)

    call read_param(ar,md)

     open(11,file='../../../model.dat')
    write(11,'(a8)')ar		
    write(11,'(a8)')md		
    close(11)

    i=system('mkdir ../../../DATA/'//ar//'/'//md//'/data')

    write(*,*)' DATASET:',ar,' MODEL:',md
    write(*,*)' key_ft1_xy2=',key_ft1_xy2,' kref=',kref

    !GOTO 771
    if(key_1real_2syn.eq.2) then
		write(*,*)' _________________________________________________________'
		write(*,*)' VISUALIZATION OF THE SYNTHETIC MODEL IN HORIZONTAL SECTIONS'
		i=system('../../4_CREATE_SYN_DATA/a_set_syn_hor/create.exe')
		write(*,*)' _________________________________________________________'
		write(*,*)' VISUALIZATION OF THE SYNTHETIC MODEL IN VERTICAL SECTIONS'
		i=system('../../4_CREATE_SYN_DATA/a_set_syn_ver/create.exe')
		write(*,*)' _________________________________________________________'
		write(*,*)' COMPUTE SYNTHETIC TRAVEL TIMES'
		i=system('../../4_CREATE_SYN_DATA/b_synth_times/rays.exe')
    end if

    i=system('cp ../../../DATA/'//ar//'/'//md//'/ref_start.dat ../../../DATA/'//ar//'/'//md//'/data/refmod.dat')
 
         
    if (key_table1_line2.eq.1) then
        write(*,*)' _________________________________________________________'
        write(*,*)' COMPUTING THE REFERENCE TABLE WITH THE UPDATED 1D MODEL'
        i=system('../../1_PRELIM_LOC/ref_table/refrays.exe')
         write(*,*)' _________________________________________________________'
        write(*,*)' LOCALIZATION OF SOURCES USING THE 1D REFERENCE TABLE'
       i=system('../../1_PRELIM_LOC/loc_table/locate.exe')
    else
        write(*,*)' _________________________________________________________'
        write(*,*)' PRELIMINARY LOCALIZATION OF SOURCES USING STRAIGHT RAYS'
        i=system('../../1_PRELIM_LOC/loc_straight/loc_1D.exe')
    end if

    771 continue

    ! Execute the ITERATIONS:
    do iter=1,niter	
        write(it,'(i1)')iter
        open(11,file='../../../model.dat')
        write(11,'(a8)')ar		
        write(11,'(a8)')md		
        write(11,'(i1)')iter		
        close(11)

        write(*,*)' _________________________________________________________'
        write(*,*)' LOCATE THE SOURCES USING THE 3D RAY TRACING '
        i=system('../../2_INVERS_3D/1_locate/locate.exe')

        do igr=1,nornt
            open(11,file='../../../model.dat')
            write(11,'(a8)')ar		
            write(11,'(a8)')md		
            write(11,'(i1)')iter		
            write(11,'(i1)')igr	
            close(11)

            if(iter.eq.1) then
                write(*,*)' _________________________________________________________'
                write(*,*)' COMPUTE THE RAY DENSITY '
                i=system('../../2_INVERS_3D/2_ray_density/plotray.exe')
                write(*,*)' _________________________________________________________'
                write(*,*)' DEFINE THE PARAMETERIZATION GRID '
                i=system('../../2_INVERS_3D/3_grid/grid.exe')
                i=system('../../2_INVERS_3D/4_links/add_matr.exe')
                !write(*,*)' _________________________________________________________'
                !write(*,*)' VISUALIZE THE RAY PATHS AND GRID IN HORIZONTAL AND VERTICAL SECTIONS '
                !i=system('../../3_VISUAL/_vis_ray_path/paths.exe')
            end if

            write(*,*)' _________________________________________________________'
            write(*,*)' COMPUTE THE 1ST DERIVATIVE MATRIX  '
            i=system('../../2_INVERS_3D/5_matrix/matr.exe')
            write(*,*)' _________________________________________________________'
            write(*,*)' PERFORM THE INVERSION (VP-VS SCHEME) '
            i=system('../../2_INVERS_3D/6_inversion/Invbig.exe')
        end do

        write(*,*)' _________________________________________________________'
        write(*,*)' COMPUTE THE VELOCITY FIELS IN 3D REGULAR GRID  '
        i=system('../../2_INVERS_3D/7_3D_model/mod_3D.exe')
        write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZE THE RESULT IN HORIZONTAL SECTIONS  '
        i=system('../../3_VISUAL/vis_hor_result/visual.exe')
        write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZE THE RESULT IN VERTICAL SECTIONS  '
        i=system('../../3_VISUAL/vis_ver_result/visual.exe')
    end do ! Different iterations
    write(*,*)' _________________________________________________________'
    write(*,*)' CREATING THE REPORT ABOUT THE VARIANCE REDUCTION '
    i=system('../../3_VISUAL/variance_reduction/var_red.exe')
end do	! Different areas

stop
end