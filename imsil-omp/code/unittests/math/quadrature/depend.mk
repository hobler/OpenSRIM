${OBJDIR}/misc/const.o : ${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/const.f90

${OBJDIR}/misc/array.o : ${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/array.f90

${OBJDIR}/misc/open.o : ${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/open.f90

${OBJDIR}/math/quadrature.o : ${SRCDIR}/math/quadrature.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/quadrature.f90

${OBJDIR}/unittests/math/quadrature/test.o : ${SRCDIR}/unittests/math/quadrature/test.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/math/quadrature.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/unittests/math/quadrature/test.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/unittests/math/quadrature/test.f90