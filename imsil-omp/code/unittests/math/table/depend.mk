${OBJDIR}/misc/const.o : ${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/const.f90

${OBJDIR}/misc/open.o : ${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/open.f90

${OBJDIR}/misc/array.o : ${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/array.f90

${OBJDIR}/math/binsearch.o : ${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/binsearch.f90

${OBJDIR}/math/solve.o : ${SRCDIR}/math/solve.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/solve.f90

${OBJDIR}/math/table.o : ${SRCDIR}/math/table.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/math/solve.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/table.f90

${OBJDIR}/unittests/math/table/test.o : ${SRCDIR}/unittests/math/table/test.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/unittests/math/table/test.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/unittests/math/table/test.f90