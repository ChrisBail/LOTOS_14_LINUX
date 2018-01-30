# Start of the makefile
# Defining variables
objects = int_best.o read_param.o loop.o
tmp_src = int_best.f90 read_param.f90 
f90comp = gfortran
switch = -O3
executable = loop.exe
# Makefile
all : $(executable)
$(executable) : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o $(executable) $(switch) $(objects)
	rm $(objects)
int_best.f90 :
	cp ../../subr/int_best.f90 int_best.f90
read_param.f90 :
	cp ../../subr/read_param.f90 read_param.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean:
	rm $(executable) 
	rm $(objects)
	rm $(tmp_src)
# End of the makefile 
