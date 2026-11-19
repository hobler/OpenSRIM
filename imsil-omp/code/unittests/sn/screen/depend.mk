${OBJDIR}/misc/const.o : ${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/const.f90

${OBJDIR}/misc/open.o : ${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/open.f90

${OBJDIR}/misc/undefine.o : ${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/undefine.f90

${OBJDIR}/misc/string.o : ${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/string.f90

${OBJDIR}/param/setup.o : ${SRCDIR}/param/setup.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/setup.f90

${OBJDIR}/math/debye.o : ${SRCDIR}/math/debye.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/debye.f90

${OBJDIR}/param/elem.o : ${SRCDIR}/param/elem.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/math/debye.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/elem.f90

${OBJDIR}/param/writetab.o : ${SRCDIR}/param/writetab.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/writetab.f90

${OBJDIR}/stack/int4.o : ${SRCDIR}/stack/int4.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/int4.f90

${OBJDIR}/stack/logical.o : ${SRCDIR}/stack/logical.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/logical.f90

${OBJDIR}/stack/real8.o : ${SRCDIR}/stack/real8.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/real8.f90

${OBJDIR}/stack/vec2.o : ${SRCDIR}/stack/vec2.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/vec2.f90

${OBJDIR}/stack/vector.o : ${SRCDIR}/stack/vector.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/vector.f90

${OBJDIR}/misc/types.o : ${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/types.f90

${OBJDIR}/stack/proj.o : ${SRCDIR}/stack/proj.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/proj.f90

${OBJDIR}/stack/stack.o : ${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/stack/int4.f90 \
	${SRCDIR}/stack/logical.f90 \
	${SRCDIR}/stack/real8.f90 \
	${SRCDIR}/stack/vec2.f90 \
	${SRCDIR}/stack/vector.f90 \
	${SRCDIR}/stack/proj.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/stack.f90

${OBJDIR}/geom/line.o : ${SRCDIR}/geom/line.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/line.f90

${OBJDIR}/math/sort.o : ${SRCDIR}/math/sort.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/sort.f90

${OBJDIR}/geom/polygon_old.o : ${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/math/sort.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/polygon_old.f90

${OBJDIR}/geom/ffpcyl.o : ${SRCDIR}/geom/ffpcyl.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/geom/line.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/ffpcyl.f90

${OBJDIR}/geom/polygon.o : ${SRCDIR}/geom/polygon.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/geom/ffpcyl.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/polygon.f90

${OBJDIR}/param/geom.o : ${SRCDIR}/param/geom.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/geom/polygon.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/geom.f90

${OBJDIR}/misc/element.o : ${SRCDIR}/misc/element.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/element.f90

${OBJDIR}/param/ion.o : ${SRCDIR}/param/ion.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/misc/element.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/ion.f90

${OBJDIR}/misc/vector.o : ${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/vector.f90

${OBJDIR}/random/haas.o : ${SRCDIR}/random/haas.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/random/haas.f90

${OBJDIR}/random/knuth.o : ${SRCDIR}/random/knuth.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/random/knuth.f90

${OBJDIR}/random/fortran.o : ${SRCDIR}/random/fortran.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/random/fortran.f90

${OBJDIR}/random/random.o : ${SRCDIR}/random/random.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/random/haas.f90 \
	${SRCDIR}/random/knuth.f90 \
	${SRCDIR}/random/fortran.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/random/random.f90

${OBJDIR}/misc/coordsys.o : ${SRCDIR}/misc/coordsys.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/random/random.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/coordsys.f90

${OBJDIR}/param/xtal.o : ${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/ion.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/xtal.f90

${OBJDIR}/param/mat.o : ${SRCDIR}/param/mat.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/coordsys.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/math/debye.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/mat.f90

${OBJDIR}/param/sn.o : ${SRCDIR}/param/sn.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/writetab.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/xtal.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/sn.f90

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

${OBJDIR}/sn/screen.o : ${SRCDIR}/sn/screen.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/sn.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/screen.f90

${OBJDIR}/unittests/sn/screen/test.o : ${SRCDIR}/unittests/sn/screen/test.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/sn/screen.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/unittests/sn/screen/test.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/unittests/sn/screen/test.f90