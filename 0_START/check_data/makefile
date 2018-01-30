# Start of the makefile
# Defining variables
objects = int_best.o epic_dist.o check_data.o
tmp_src = int_best.f90 epic_dist.f90 
f90comp = gfortran
switch = -O3
executable = check_data.exe
# Makefile
all : $(executable)
$(executable) : $(objects)
	rm $(tmp_src)	
	$(f90comp) -o $(executable) $(switch) $(objects)
	rm $(objects)
int_best.f90 :
	cp ../../subr/int_best.f90 int_best.f90
epic_dist.f90 :
	cp ../../subr/epic_dist.f90 epic_dist.f90
%.o: %.f90
	$(f90comp) -c $(switch) $<
# Cleaning everything
clean:
	rm $(executable) 
	rm $(objects)
	rm $(tmp_src)
# End of the makefile 