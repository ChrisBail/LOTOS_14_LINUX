# Start of the makefile
# Defining variables
objects = vrefmod.o read_vref.o read_3D_mod_v.o velocity.o anom_3D_xyz_lin_v.o matr.o flat_sph.o sph_flat.o int_best.o read_param.o
tmp_src = vrefmod.f90 read_vref.f90 read_3D_mod_v.f90 velocity.f90 anom_3D_xyz_lin_v.f90 flat_sph.f90 sph_flat.f90 int_best.f90 read_param.f90
f90comp = gfortran
switch = -O3
executable = matr.exe
# Makefile
all : $(executable)
$(executable) : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o $(executable) $(switch) $(objects)
	rm $(objects)
vrefmod.f90 :
	cp ../../subr/vrefmod.f90 vrefmod.f90
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90
read_3D_mod_v.f90 :	
	cp ../../subr/read_3D_mod_v.f90 read_3D_mod_v.f90
velocity.f90 :	
	cp ../../subr/velocity.f90 velocity.f90
anom_3D_xyz_lin_v.f90 :	
	cp ../../subr/anom_3D_xyz_lin_v.f90 anom_3D_xyz_lin_v.f90
flat_sph.f90 :	
	cp ../../subr/flat_sph.f90 flat_sph.f90
sph_flat.f90 :	
	cp ../../subr/sph_flat.f90 sph_flat.f90
int_best.f90 :	
	cp ../../subr/int_best.f90 int_best.f90
read_param.f90 :	
	cp ../../subr/read_param.f90 read_param.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean: 
	rm $(objects)
	rm $(tmp_src)
	rm $(executable)
# End of the makefile
