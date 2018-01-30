# Start of the makefile
# Defining variables
objects = int_best.o decsf.o alin.o grid.o h_lim.o read_topo.o read_param.o
tmp_src = int_best.f90 decsf.f90 alin.f90 h_lim.f90 read_topo.f90 read_param.f90
f90comp = gfortran
switch = -O3
executable = grid.exe
# Makefile
all : $(executable)
$(executable) : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o $(executable) $(switch) $(objects)
	rm $(objects)
int_best.f90 :
	cp ../../subr/int_best.f90 int_best.f90
alin.f90 :
	cp ../../subr/alin.f90 alin.f90
decsf.f90 :
	cp ../../subr/decsf.f90 decsf.f90	
h_lim.f90 :
	cp ../../subr/h_lim.f90 h_lim.f90	
read_topo.f90 :
	cp ../../subr/read_topo.f90 read_topo.f90	
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
