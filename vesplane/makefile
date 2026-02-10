#########################################################
#  PROGRAM NAME (no .f90)
main = 0Main_VesiclePlane

#  SOURCE DIRECTORY FOR MODULES
SRCDIR = src

#  MODULE NAMES (no .f90)
mod1  = allocCel_mod
mod2  = basicfun_mod
mod3  = CelType_mod
mod4  = read_mod
mod5  = Manipulate_Cell_mod
mod6  = dyn_mod
mod7  = Surf_Operation_mod
mod11 = point_mod
mod12 = force_mod
mod13 = head_Cell_mod

#########################################################
#  COMPILER OPTIONS
cmplr = gfortran -O3

#########################################################
#  OBJECT FILES
objects = \
 $(mod11).o \
 $(mod3).o \
 $(mod1).o \
 $(mod2).o \
 $(mod4).o \
 $(mod5).o \
 $(mod6).o \
 $(mod7).o \
 $(mod12).o \
 $(mod13).o \
 $(main).o

#########################################################
#  LINK
$(main): $(objects)
	$(cmplr) -o $(main) $(objects)

#########################################################
#  COMPILE MODULES (from src/)
$(mod11).o : $(SRCDIR)/$(mod11).f90
	$(cmplr) -c $<

$(mod3).o  : $(SRCDIR)/$(mod3).f90
	$(cmplr) -c $<

$(mod1).o  : $(SRCDIR)/$(mod1).f90
	$(cmplr) -c $<

$(mod2).o  : $(SRCDIR)/$(mod2).f90
	$(cmplr) -c $<

$(mod4).o  : $(SRCDIR)/$(mod4).f90
	$(cmplr) -c $<

$(mod5).o  : $(SRCDIR)/$(mod5).f90
	$(cmplr) -c $<

$(mod6).o  : $(SRCDIR)/$(mod6).f90
	$(cmplr) -c $<

$(mod7).o  : $(SRCDIR)/$(mod7).f90
	$(cmplr) -c $<

$(mod12).o : $(SRCDIR)/$(mod12).f90
	$(cmplr) -c $<

$(mod13).o : $(SRCDIR)/$(mod13).f90
	$(cmplr) -c $<

#########################################################
#  COMPILE MAIN (still in current directory)
$(main).o : $(main).f90
	$(cmplr) -c $<

#########################################################
#  CLEAN
clean:
	rm -f *.mod *.o $(main)
