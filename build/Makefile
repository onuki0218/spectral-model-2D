# module load mpi-fftw
SRCDIR = ../src/
RUNDIR = ../run/
RESULTDIR = ../data/
# BUILDDIR = ../build/
# INCLUDEDIR = ../include/

OBJS = *.o
MODS = *.mod
TARGET = main.exe
# TARGET_RE = reanalysis.exe

FC = mpiifort
OPTS = -O3 -axMIC-AVX512 -parallel -qopenmp -g -assume byterecl
LIBS = -lfftw3 -lfftw3_omp
RM = rm -f

MODOBJS = mod_core.o mod_common.o mod_io.o mod_governing_equation.o misk.o
MODOBJS += mod_runge_kutta.o
MODOBJS += mod_initial_condition.o
MODOBJS += mod_log.o
MODOBJS += mod_dealias.o
# MODOBJS += mod_wave.o
MODOBJS += mod_parallel.o mod_fftw.o 
MODOBJS += mod_analysis.o

MODSRCS = $(MODOBJS:%.o=%.F90)

MAINOBJ = main.o
# MAINOBJ_RE = reanalysis.o
MAINSRC = $(MAINOBJ:%.o=%.F90)

all: $(RUNDIR) $(RESULTDIR) $(RUNDIR)$(TARGET)
re : $(RUNDIR)$(TARGET_RE)

# make target file
$(RUNDIR)$(TARGET) :  $(MAINOBJ) $(MODOBJS)
	$(FC) $(OPTS) -o $@ $(MAINOBJ) $(MODOBJS) $(LIBS)

# $(RUNDIR)$(TARGET_RE) :  $(MAINOBJ_RE) $(MODOBJS)
# 	$(FC) $(OPTS) -o $@ $(MAINOBJ_RE) $(MODOBJS) $(LIBS)

# make ofject file
%.o : $(SRCDIR)%.F90
	$(FC) $(OPTS) -o $@ -c $< $(LIBS)

# dependency of object files
main.o : mod_common.o mod_core.o

# reanalysis.o : mod_common.o mod_log.o mod_io.o mod_parallel.o  \
 mod_fftw.o mod_wave.o

mod_core.o : mod_common.o mod_io.o mod_governing_equation.o mod_runge_kutta.o  \
 mod_dealias.o mod_parallel.o mod_fftw.o mod_log.o misk.o mod_initial_condition.o

mod_io.o : mod_common.o mod_analysis.o misk.o

mod_governing_equation.o : mod_common.o mod_fftw.o misk.o

misk.o : mod_common.o

mod_runge_kutta.o : mod_common.o mod_governing_equation.o mod_fftw.o  \
 mod_dealias.o misk.o

mod_initial_condition.o : mod_common.o mod_fftw.o mod_dealias.o  \
 mod_governing_equation.o misk.o mod_analysis.o

mod_log.o : mod_common.o mod_io.o mod_analysis.o misk.o

mod_dealias.o : mod_common.o

# mod_wave.o : mod_common.o

mod_parallel.o : mod_common.o misk.o

mod_fftw.o : mod_common.o misk.o

mod_analysis.o : mod_common.o mod_fftw.o misk.o mod_governing_equation.o

# make directory
$(RUNDIR) :
	mkdir $(RUNDIR) 

$(RESULTDIR) :
	mkdir $(RESULTDIR)

# clean
clean :
	$(RM) *.mod
	$(RM) *.o
	$(RM) $(RUNDIR)$(TARGET)

cleanD :
	$(RM) $(RESULTDIR)*
