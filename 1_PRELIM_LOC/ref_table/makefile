# Start of the makefile
# Defining variables
objects = ray_lay.o reftrace.o vrefmod.o refrays.o read_vref.o read_param.o
tmp_src = ray_lay.f90 reftrace.f90 vrefmod.f90 read_vref.f90 read_param.f90
f90comp = gfortran
switch = -O3
# Makefile
all : refrays.exe
refrays.exe : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o refrays.exe $(switch) $(objects)
	rm $(objects)
ray_lay.f90 : 
	cp ../../subr/ray_lay.f90 ray_lay.f90
reftrace.f90 :
	cp ../../subr/reftrace.f90 reftrace.f90
vrefmod.f90 :
	cp ../../subr/vrefmod.f90 vrefmod.f90  	
read_vref.f90 :
	cp ../../subr/read_vref.f90 read_vref.f90  	
read_param.f90 :
	cp ../../subr/read_param.f90 read_param.f90  	
%.o: %.f90
	$(f90comp) -Wconversion -c $(switch) $< 
# Cleaning everything
clean: 
	rm $(objects)
	rm $(tmp_src)
	rm refrays.exe
# End of the makefile 