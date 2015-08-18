############################
.SUFFIXES: .o .f90 .mod

F90 = gfortran
FFLAGS = -c -O3
LFLAGS = -O3

OBJS = module.o main.o aux.o scf.o

oneband: $(OBJS)
	$(F90) $(LFLAGS) -o oneband $(OBJS)

module.o: module.f90
	$(F90) $(FFLAGS) module.f90

main.o: main.f90
	$(F90) $(FFLAGS) main.f90

aux.o: aux.f90
	$(F90) $(FFLAGS) aux.f90

scf.o: scf.f90
	$(F90) $(FFLAGS) scf.f90

clean:
	rm -f oneband *.o core *.mod *~

clr:
	rm -f *.x0.* Gi.txt error.txt
