# Start of the makefile
# Defining variables
objects = decsf.o dispers.o goal_new.o prepare_ref.o read_z_limit.o refmod_all.o sfdec.o vrefmod.o z_lim.o read_topo.o read_vref.o h_lim.o read_param.o locate.o
tmp_src = decsf.f90 dispers.f90 goal_new.f90 prepare_ref.f90 read_z_limit.f90 refmod_all.f90 sfdec.f90 vrefmod.f90 z_lim.f90 read_topo.f90 read_vref.f90 h_lim.f90 read_param.f90 
f90comp = gfortran
switch = -O3
# Makefile
all : locate.exe
locate.exe : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o locate.exe $(switch) $(objects)
	rm $(objects)
decsf.f90 :	
	cp ../../subr/decsf.f90 decsf.f90
dispers.f90 :
	cp ../../subr/dispers.f90 dispers.f90
goal_new.f90 :
	cp ../../subr/goal_new.f90 goal_new.f90
prepare_ref.f90 :	
	cp ../../subr/prepare_ref.f90 prepare_ref.f90
read_z_limit.f90 :
	cp ../../subr/read_z_limit.f90 read_z_limit.f90
refmod_all.f90 :
	cp ../../subr/refmod_all.f90 refmod_all.f90
sfdec.f90 :
	cp ../../subr/sfdec.f90 sfdec.f90
vrefmod.f90 :	
	cp ../../subr/vrefmod.f90 vrefmod.f90
z_lim.f90 :	
	cp ../../subr/z_lim.f90 z_lim.f90
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90
read_topo.f90 :	
	cp ../../subr/read_topo.f90 read_topo.f90
h_lim.f90 :	
	cp ../../subr/h_lim.f90 h_lim.f90
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
