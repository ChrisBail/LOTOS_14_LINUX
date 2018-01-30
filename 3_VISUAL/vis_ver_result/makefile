# Start of the makefile
# Defining variables
objects = dv_1_grid_v.o prepare_model_v.o decsf.o kram3.o sfdec.o vrefmod.o read_vref.o read_topo.o h_lim.o visual.o flat_sph.o sph_flat.o read_param.o real2char.o create_color_scales.o
tmp_src =  dv_1_grid_v.f90 prepare_model_v.f90 decsf.f90 kram3.f90 sfdec.f90 vrefmod.f90 read_vref.f90 read_topo.f90 h_lim.f90 flat_sph.f90 sph_flat.f90 read_param.f90 real2char.f90 create_color_scales.f90
f90comp = gfortran
switch = -O3
# Makefile
all : visual.exe
visual.exe : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o visual.exe $(switch) $(objects)
	rm $(objects)
dv_1_grid_v.f90 :
	cp ../../subr/dv_1_grid_v.f90 dv_1_grid_v.f90
prepare_model_v.f90 :	
	cp ../../subr/prepare_model_v.f90 prepare_model_v.f90
decsf.f90 :	
	cp ../../subr/decsf.f90 decsf.f90
kram3.f90 :	
	cp ../../subr/kram3.f90 kram3.f90
sfdec.f90 :	
	cp ../../subr/sfdec.f90 sfdec.f90
vrefmod.f90 :	
	cp ../../subr/vrefmod.f90 vrefmod.f90
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90
read_topo.f90 :	
	cp ../../subr/read_topo.f90 read_topo.f90
h_lim.f90 :	
	cp ../../subr/h_lim.f90 h_lim.f90
flat_sph.f90:
	cp ../../subr/flat_sph.f90 flat_sph.f90
sph_flat.f90:
	cp ../../subr/sph_flat.f90 sph_flat.f90
read_param.f90:
	cp ../../subr/read_param.f90 read_param.f90
real2char.f90:
	cp ../../subr/real2char.f90 real2char.f90
create_color_scales.f90:
	cp ../../subr/create_color_scales.f90 create_color_scales.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean: 
	rm $(objects)
	rm $(tmp_src)
	rm visual.exe
# End of the makefile 