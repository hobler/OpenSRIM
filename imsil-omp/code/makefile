# ----------------------------------------------------------------------
# Makefile for the installation of a Fortran project
#
# With FC=ifort, the Intel Fortran compiler is used, otherwise gfortran.
# If necessary, the path of the Intel compiler has to be specified below.
# 
# Requires the files "dirs.mk" ("dirs_db.mk" for a debug version), 
# "objs.list" and "depend.mk". These files are included from the current
# directory, i.e. from the directory you start make from. If these files
# are not in the same directory as the makefile, go to the include file
# directory and invoke make with
#
# "cd /path/to/include_files; make -f /path/to/makefile"
#
# "dirs.mk" defines the root directory of the source code (SRCDIR), the 
# directories for the object files (OBJDIR) and the module files (MODDIR), 
# and the desired path of the executable (macro EXE).
#
# "objs.mk" defines a list of object files (macro OBJS). 
#
# "depend.mk" defines the dependencies.
#
# USAGE:
#
# "make" |
# "make all"
# "make FC=ifort" |
# "make all FC=ifort"       ... makes optimized version of the code
#
# "make DB=db" |
# "make all DB=db"
# "make DB=db FC=ifort" |
# "make all DB=db FC=ifort" ... makes debug version of the code
#
# "make clean" |
# "make clean DB=db"    ... delete object and module directories as well as 
#							the executable
#
# To all of the above invocations of make the command line variable OMP=omp
# can be added to produce a parallel version using OpenMP 
#
# In addition, the command line argument N can specify the default number of 
# threads the program will run with. With the gfortran compiler, the program
# will always be compiled for 1 thread independent of N
#
# 10-Aug-2009  G.Hobler  creation from imsil Makefile.tem
# 26-Jan-2016  G.Hobler  adpted to new code organization
# 29-Nov-2023  G.Hobler  OpenMP version (OMP command line variable)
# 20-Mar-2024  G.Hobler  Adjust compiler and linker flags
# ----------------------------------------------------------------------


# Macros ---------------------------------------------------------------

DIRS__=dirs.mk
DIRS_db_=dirs_db.mk
DIRS__omp=dirs_omp.mk
DIRS_db_omp=dirs_db_omp.mk
include ${DIRS_${DB}_${OMP}}

ifeq (${FC}, ifort)

# INTEL Fortran compiler
# This requires make to start from a shell that has invoked the ifort startup
# script so that all necessary environment variables are set
# When using eclipse, eclipse must be started from such a shell
#  override FC=ifort
override FC=ifx
LINK=${FC}

# flags for INTEL compiler
# -O3 brings no clear advantage over -O2
# processor specific optimization (-x, -ax) higher that SSE4.2 seems to be detrimental
# our processors support: pc9,10: SSE4.2, pc8,158: AVX, pc155,156:CORE-AVX512
# A previous version of IMSIL linked the MKL library. For this
#   - Add -mkl flag to FFLAGS
#   - Add -mkl -lmkl_lapack95_lp64 to LFLAGS for dynamic linking
#   - Add -mkl -static-intel -L${MKLROOT}/lib/intel64 -lmkl_lapack95_lp64 to LFLAGS for static
#     linking
# Since 2017, -assume realloc_lhs is no longer needed to conform to Fortran 2003 standard
# To put automatic and temporary arrays on heap, use -heap-arrays.
# -check all produces a lot of warning of the type "forrtl: warning (406): fort: (1): In call 
# to ABS2VECTOR, an array temporary was created for argument #1, which seem irrelevant. Avoid
# this with -check noarg_temp_created
FFLAGS__ := -O2 -ipo   # -xSSE4.2 
#FFLAGS_db_ := -O0 -g -warn all -stand f08 -traceback -check pointers -check bounds -p -fpe0
FFLAGS_db_ := -O0 -g -warn all -check all -check noarg_temp_created -stand f08 -fp-model=precise
FFLAGS_${DB}_omp := ${FFLAGS_${DB}_} -qopenmp #-static
FFLAGS_${DB}_ := ${FFLAGS_${DB}_} #-static 
# we do not specify the number of threads because it overrides the environment variable
FFLAGS := ${FFLAGS_${DB}_${OMP}} -traceback -threads -assume nobscc -module ${MODDIR} -heap-arrays 512 -fpp
ifeq (${VERSION}, AMAT)
FFLAGS:=${FFLAGS} -I${MKLROOT}/include/intel64/lp64 -I"${MKLROOT}/include"
endif

ifneq (${VERSION}, GTS)
LFLAGS_${DB}_${OMP} := -static
endif
ifeq (${VERSION}, AMAT)
LFLAGS_${DB}_${OMP} := ${LFLAGS_${DB}_${OMP}} ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -static-intel
LFLAGS_${DB}_${OMP} := ${LFLAGS_${DB}_${OMP}} -lstdc++
endif
LFLAGS_${DB}_omp := ${LFLAGS_${DB}_omp}  -threads -qopenmp -qopenmp-link=static
LFLAGS__${OMP} := ${LFLAGS__${OMP}} -ipo 
LFLAGS := ${LFLAGS_${DB}_${OMP}}

#LFLAGS__ := -ipo
#LFLAGS_db_ :=
#LFLAGS_${DB}_omp := ${LFLAGS_${DB}_} -qopenmp -qopenmp-link=static # -l:libiomp5.a
#LFLAGS_${DB}_ := ${LFLAGS_${DB}_} -static
## statically linking OpenMP is possible with -qopenmp-link=static (not recommended by Intel)
## not sure what the -Xlinker options were for (coarray?)
#LFLAGS := ${LFLAGS_${DB}_${OMP}} -threads # -Xlinker -rpath -Xlinker . # -static -L libifcoremt.a
#ifeq (${VERSION}, AMAT)
#LFLAGS:=${LFLAGS} ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
#LFLAGS:=${LFLAGS} -lstdc++
#endif

else

# GNU Fortran compiler
override FC=gfortran
LINK=${FC}

# -Ofast leads to infinite loop in rrobin_screen  
FFLAGS__ := -O3 -flto
FFLAGS_db_ := -g -O0
# use the following line in desparate situations:
FFLAGS_db_ := ${FFLAGS_db_} -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace
FFLAGS_${DB}_omp := ${FFLAGS_${DB}_} -fopenmp 
FFLAGS_${DB}_ := ${FFLAGS_${DB}_} #-static-libgfortran -static  
FFLAGS := ${FFLAGS_${DB}_${OMP}} -J ${MODDIR} -cpp

LFLAGS__ := -flto
LFLAGS_db_ :=
LFLAGS_${DB}_omp := ${LFLAGS_${DB}_} -fopenmp
LFLAGS_${DB}_ := ${LFLAGS_${DB}_} #-static-libgfortran -static  
# statically linking OpenMP still requires at runtime the shared libraries from the 
# glibc version used for linking
LFLAGS := ${LFLAGS_${DB}_${OMP}}
ifeq (${VERSION}, AMAT)
LFLAGS:=${LFLAGS} -lstdc++
endif

endif

# to add the search path for shared object files
ifeq (${VERSION}, GTS)
LFLAGS:=${LFLAGS} -lstdc++										# when using static library
#LFLAGS:=${LFLAGS} -Wl,-rpath=${HOME}/work/imsil/gts/code/lib	# when using dynamic library
endif

ifdef VERSION
FFLAGS:=${FFLAGS} -D${VERSION}
endif

include srcs.mk
include objs.mk
include libs.mk

# Suffix Rules ---------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .exe .o .f90

.f90.o:
	${FC} ${FFLAGS} -o $@ -c src/${@F:.o=.f90}

# Targets --------------------------------------------------------------

all : ${EXE}

ifeq (${DB}, db)
${EXE} : ${OBJS} ${LIBS} ${DIRS_${DB}_${OMP}} objs.mk depend.mk ${SRCDIR}/makefile
	${LINK} -o $@ ${OBJS} ${LIBS} ${LFLAGS}
else
${EXE} : ${SRCS} ${LIBS} ${DIRS_${DB}_${OMP}} srcs.mk depend.mk ${SRCDIR}/makefile
	${FC} ${FFLAGS} -o ${EXE} ${SRCS} ${LIBS} ${LFLAGS}
endif

clean :
	-rm -rf ${OBJDIR}
	-rm -rf ${MODDIR}
	-rm ${EXE}

# Dependencies of object files and compile commands --------------------

include depend.mk

# Avoid looking for a rule to make .so files ------------------------

ifeq (${VERSION}, GTS)
%.so: ;
endif

