# Start of the makefile
# Defining variables
objects = anom_3D_xyz_lin_v.o read_3D_mod_v.o velocity.o vrefmod.o decsf.o dispers.o int_best.o inv_svd_simple.o read_z_limit.o\
sfdec.o z_lim.o part_bending_hor.o part_bending_z.o part_ray.o straight_line.o trace_bending.o locate.o read_topo.o read_vref.o h_lim.o flat_sph.o sph_flat.o read_param.o
tmp_src = anom_3D_xyz_lin_v.f90 read_3D_mod_v.f90 velocity.f90 vrefmod.f90 decsf.f90 dispers.f90 int_best.f90 inv_svd_simple.f90 read_z_limit.f90\
sfdec.f90 z_lim.f90 part_bending_hor.f90 part_bending_z.f90 part_ray.f90 straight_line.f90 trace_bending.f90 read_topo.f90 read_vref.f90 h_lim.f90 flat_sph.f90 sph_flat.f90 read_param.f90
f90comp = gfortran
switch = -O3
# Makefile
all : locate.exe
locate.exe : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o locate.exe $(switch) $(objects)
	rm $(objects)
h_lim.f90 :	
	cp ../../subr/h_lim.f90 h_lim.f90
read_topo.f90 :	
	cp ../../subr/read_topo.f90 read_topo.f90
anom_3D_xyz_lin_v.f90 :	
	cp ../../subr/anom_3D_xyz_lin_v.f90 anom_3D_xyz_lin_v.f90
read_3D_mod_v.f90 :	
	cp ../../subr/read_3D_mod_v.f90 read_3D_mod_v.f90
velocity.f90 :	
	cp ../../subr/velocity.f90 velocity.f90
vrefmod.f90 :
	cp ../../subr/vrefmod.f90 vrefmod.f90
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90
decsf.f90 :
	cp ../../subr/decsf.f90 decsf.f90
dispers.f90 :	
	cp ../../subr/dispers.f90 dispers.f90
int_best.f90 :	
	cp ../../subr/int_best.f90 int_best.f90
inv_svd_simple.f90 :	
	cp ../../subr/inv_svd_simple.f90 inv_svd_simple.f90
read_z_limit.f90 :
	cp ../../subr/read_z_limit.f90 read_z_limit.f90
sfdec.f90 :
	cp ../../subr/sfdec.f90 sfdec.f90
z_lim.f90 :
	cp ../../subr/z_lim.f90 z_lim.f90
part_bending_hor.f90 :
	cp ../../subr/part_bending_hor.f90 part_bending_hor.f90
part_bending_z.f90 :
	cp ../../subr/part_bending_z.f90 part_bending_z.f90
part_ray.f90 :
	cp ../../subr/part_ray.f90 part_ray.f90
straight_line.f90 :
	cp ../../subr/straight_line.f90 straight_line.f90
trace_bending.f90 :
	cp ../../subr/trace_bending.f90 trace_bending.f90
flat_sph.f90 :
	cp ../../subr/flat_sph.f90 flat_sph.f90	
sph_flat.f90 :
	cp ../../subr/sph_flat.f90 sph_flat.f90	
read_param.f90 :
	cp ../../subr/read_param.f90 read_param.f90	
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean: 
	rm $(objects)
	rm $(tmp_src)
	rm locate.exe
# End of the makefile
