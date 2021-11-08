#FCOMP = gfortran
#OPTS  = -c -ffpe-trap=zero -fbacktrace -fimplicit-none -fcheck=all -g
FCOMP = ifort
OPTS = -c #-fpe0 -traceback -CB -g
LINKOPTS = -O3 -o
OBJS =  Modules.o main.o Util_IO.o \
		memory_mgmt.o make_grid.o IBM.o initsetup.o \
		flux_advective.o flux_diffusive.o setBC.o\
		iterations.o Poisson_Solver.o LA_lib.o Probe.o

main: $(OBJS) 
	  $(FCOMP) $(LINKOPTS) 2D_NS_IBM $(OBJS)

clean:
	rm *.o *.mod
	rm 2D_NS_IBM

.SUFFIXES: .o .F90 .f90

.f90.o:
	$(FCOMP) $(OPTS) $(DEFINES) $<

.F90.o:
	$(FCOMP) $(OPTS) $(DEFINES) $<

