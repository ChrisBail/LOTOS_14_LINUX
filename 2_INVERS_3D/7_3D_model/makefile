# Start of the makefile
# Defining variables
objects = vrefmod.o prepare_model_v.o kram3.o dv_1_gr_xyz_v.o read_vref.o mod_3D.o read_param.o
tmp_src = vrefmod.f90 prepare_model_v.f90 kram3.f90 dv_1_gr_xyz_v.f90 read_vref.f90 read_param.f90
f90comp = gfortran
switch = -O3
executable = mod_3D.exe
# Makefile
all : $(executable)
$(executable) : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o $(executable) $(switch) $(objects)
	rm $(objects)
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90
vrefmod.f90 :
	cp ../../subr/vrefmod.f90 vrefmod.f90
prepare_model_v.f90 :
	cp ../../subr/prepare_model_v.f90 prepare_model_v.f90
kram3.f90 :
	cp ../../subr/kram3.f90 kram3.f90
dv_1_gr_xyz_v.f90 :
	cp ../../subr/dv_1_gr_xyz_v.f90 dv_1_gr_xyz_v.f90	
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
