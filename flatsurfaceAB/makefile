FC = gfortran -O3

OBJS = point_mod.o CelType_mod.o allocCel_mod.o basicfun_mod.o read_mod.o \
       Manipulate_Cell_mod.o dyn_mod.o Surf_Operation_mod.o force_mod.o \
       head_Cell_mod.o 0Main_FlatSurfaceAB.o

0Main_FlatSurfaceAB: $(OBJS)
	$(FC) $(OBJS) -o 0Main_FlatSurfaceAB

point_mod.o: src/point_mod.f90
	$(FC) -c src/point_mod.f90

CelType_mod.o: src/CelType_mod.f90 point_mod.o
	$(FC) -c src/CelType_mod.f90

allocCel_mod.o: src/allocCel_mod.f90 CelType_mod.o
	$(FC) -c src/allocCel_mod.f90

basicfun_mod.o: src/basicfun_mod.f90 allocCel_mod.o CelType_mod.o
	$(FC) -c src/basicfun_mod.f90

read_mod.o: src/read_mod.f90 basicfun_mod.o CelType_mod.o
	$(FC) -c src/read_mod.f90

Manipulate_Cell_mod.o: src/Manipulate_Cell_mod.f90 basicfun_mod.o CelType_mod.o allocCel_mod.o
	$(FC) -c src/Manipulate_Cell_mod.f90

dyn_mod.o: src/dyn_mod.f90 CelType_mod.o
	$(FC) -c src/dyn_mod.f90

Surf_Operation_mod.o: src/Surf_Operation_mod.f90 basicfun_mod.o CelType_mod.o allocCel_mod.o
	$(FC) -c src/Surf_Operation_mod.f90

force_mod.o: src/force_mod.f90 basicfun_mod.o CelType_mod.o point_mod.o Surf_Operation_mod.o
	$(FC) -c src/force_mod.f90

head_Cell_mod.o: src/head_Cell_mod.f90 force_mod.o allocCel_mod.o read_mod.o basicfun_mod.o \
                 Manipulate_Cell_mod.o point_mod.o dyn_mod.o Surf_Operation_mod.o
	$(FC) -c src/head_Cell_mod.f90

0Main_FlatSurfaceAB.o: 0Main_FlatSurfaceAB.f90 head_Cell_mod.o
	$(FC) -c 0Main_FlatSurfaceAB.f90

clean:
	rm -f *.mod *.o 0Main_FlatSurfaceAB
