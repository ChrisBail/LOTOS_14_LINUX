# Start of the makefile
# Defining variables
objects = vrefmod.o nolet.o anom_3D_xyz_lin_v.o read_3D_mod_v.o velocity.o read_vref.o sfdec.o Invbig.o flat_sph.o sph_flat.o read_param.o
tmp_src = vrefmod.f90 nolet.f90 anom_3D_xyz_lin_v.f90 read_3D_mod_v.f90 velocity.f90 read_vref.f90 sfdec.f90 flat_sph.f90 sph_flat.f90 read_param.f90
f90comp = gfortran
switch = -O3
executable = Invbig.exe
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
nolet.f90 :
	cp ../../subr/nolet.f90 nolet.f90
sfdec.f90 :
	cp ../../subr/sfdec.f90 sfdec.f90
velocity.f90 :
	cp ../../subr/velocity.f90 velocity.f90
read_3D_mod_v.f90 :
	cp ../../subr/read_3D_mod_v.f90 read_3D_mod_v.f90
anom_3D_xyz_lin_v.f90 :
	cp ../../subr/anom_3D_xyz_lin_v.f90 anom_3D_xyz_lin_v.f90
flat_sph.f90:
	cp ../../subr/flat_sph.f90 flat_sph.f90
sph_flat.f90:
	cp ../../subr/sph_flat.f90 sph_flat.f90	
read_param.f90:
	cp ../../subr/read_param.f90 read_param.f90	
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean: 
	rm $(objects)
	rm $(tmp_src)
	rm $(executable)
# End of the makefile
