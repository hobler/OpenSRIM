${OBJDIR}/misc/const.o : ${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/const.f90

${OBJDIR}/math/solve.o : ${SRCDIR}/math/solve.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/solve.f90

${OBJDIR}/unittests/math/solve/test.o : ${SRCDIR}/unittests/math/solve/test.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/math/solve.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/unittests/math/solve/test.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/unittests/math/solve/test.f90