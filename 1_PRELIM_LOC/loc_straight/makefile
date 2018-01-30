# Start of the makefile
# Defining variables
objects = anom_3D_xyz_lin_v.o goal_fun_lin.o read_vref.o straight_line.o velocity.o vrefmod.o decsf.o flat_sph.o sph_flat.o h_lim.o int_best.o read_topo.o read_z_limit.o relief_surf.o sfdec.o z_lim.o loc_1D.o read_param.o read_3D_mod_v.o
tmp_src = anom_3D_xyz_lin_v.f90 goal_fun_lin.f90 read_vref.f90 straight_line.f90 velocity.f90 vrefmod.f90 decsf.f90 flat_sph.f90 sph_flat.f90 h_lim.f90 int_best.f90 read_topo.f90 read_z_limit.f90 relief_surf.f90 sfdec.f90 read_param.f90 z_lim.f90 read_3D_mod_v.f90 
f90comp = gfortran
switch = -O3
# Makefile
all : loc_1D.exe
loc_1D.exe : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o loc_1D.exe $(switch) $(objects)
	rm $(objects)
anom_3D_xyz_lin_v.f90 :	
	cp ../../subr/anom_3D_xyz_lin_v.f90 anom_3D_xyz_lin_v.f90
goal_fun_lin.f90 :
	cp ../../subr/goal_fun_lin.f90 goal_fun_lin.f90
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90
straight_line.f90 :
	cp ../../subr/straight_line.f90 straight_line.f90
velocity.f90 :
	cp ../../subr/velocity.f90 velocity.f90
vrefmod.f90 :	
	cp ../../subr/vrefmod.f90 vrefmod.f90
decsf.f90 :	
	cp ../../subr/decsf.f90 decsf.f90
flat_sph.f90 :	
	cp ../../subr/flat_sph.f90 flat_sph.f90
sph_flat.f90 :	
	cp ../../subr/sph_flat.f90 sph_flat.f90
h_lim.f90 :	
	cp ../../subr/h_lim.f90 h_lim.f90
int_best.f90 :
	cp ../../subr/int_best.f90 int_best.f90
read_topo.f90 :	
	cp ../../subr/read_topo.f90 read_topo.f90
read_z_limit.f90 :
	cp ../../subr/read_z_limit.f90 read_z_limit.f90
relief_surf.f90 :	
	cp ../../subr/relief_surf.f90 relief_surf.f90
sfdec.f90 :
	cp ../../subr/sfdec.f90 sfdec.f90
z_lim.f90 :	
	cp ../../subr/z_lim.f90 z_lim.f90
read_param.f90 :	
	cp ../../subr/read_param.f90 read_param.f90
read_3D_mod_v.f90 :	
	cp ../../subr/read_3D_mod_v.f90 read_3D_mod_v.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean: 
	rm $(objects)
	rm $(tmp_src)
	rm loc_1D.exe
# End of the makefile
