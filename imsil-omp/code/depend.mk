${OBJDIR}/misc/const.o : ${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/const.f90

${OBJDIR}/misc/open.o : ${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/open.f90

${OBJDIR}/misc/time.o : ${SRCDIR}/misc/time.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/time.f90

${OBJDIR}/misc/undefine.o : ${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/undefine.f90

${OBJDIR}/misc/string.o : ${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/string.f90

${OBJDIR}/math/intmath.o : ${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/intmath.f90

${OBJDIR}/param/setup.o : ${SRCDIR}/param/setup.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/math/intmath.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/setup.f90

${OBJDIR}/misc/types.o : ${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/types.f90

${OBJDIR}/misc/array.o : ${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/array.f90

${OBJDIR}/stack/int4.o : ${SRCDIR}/stack/int4.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/int4.f90

${OBJDIR}/stack/logical.o : ${SRCDIR}/stack/logical.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/logical.f90

${OBJDIR}/stack/real8.o : ${SRCDIR}/stack/real8.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/real8.f90

${OBJDIR}/stack/vec2.o : ${SRCDIR}/stack/vec2.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/vec2.f90

${OBJDIR}/stack/vector.o : ${SRCDIR}/stack/vector.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stack/vector.f90

${OBJDIR}/stack/proj.o : ${SRCDIR}/stack/proj.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
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
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/math/sort.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/polygon_old.f90

${OBJDIR}/geom/ffpcyl.o : ${SRCDIR}/geom/ffpcyl.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/geom/line.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/ffpcyl.f90

${OBJDIR}/geom/polygon.o : ${SRCDIR}/geom/polygon.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/geom/ffpcyl.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/polygon.f90

${OBJDIR}/param/geom.o : ${SRCDIR}/param/geom.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/geom/polygon.f90 \
	${SRCDIR}/param/setup.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/geom.f90

${OBJDIR}/math/debye.o : ${SRCDIR}/math/debye.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/debye.f90

${OBJDIR}/param/elem.o : ${SRCDIR}/param/elem.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/math/debye.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/elem.f90

${OBJDIR}/misc/element.o : ${SRCDIR}/misc/element.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/misc/element.f90

${OBJDIR}/param/ion.o : ${SRCDIR}/param/ion.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/misc/element.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/ion.f90

${OBJDIR}/param/xtal.o : ${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/ion.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/xtal.f90

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

${OBJDIR}/geom/coordsys.o : ${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/random/random.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/coordsys.f90

${OBJDIR}/param/mat.o : ${SRCDIR}/param/mat.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/math/debye.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/mat.f90

${OBJDIR}/param/writetab.o : ${SRCDIR}/param/writetab.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/writetab.f90

${OBJDIR}/param/dam.o : ${SRCDIR}/param/dam.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/writetab.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/dam.f90

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
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/sn.f90

${OBJDIR}/param/se.o : ${SRCDIR}/param/se.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/writetab.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/se.f90

${OBJDIR}/param/out.o : ${SRCDIR}/param/out.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/dam.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/out.f90

${OBJDIR}/geom/bndcond.o : ${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/out.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/bndcond.f90

${OBJDIR}/geom/curvicoord.o : ${SRCDIR}/geom/curvicoord.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/geom.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/curvicoord.f90

${OBJDIR}/math/binsearch.o : ${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/binsearch.f90

${OBJDIR}/target/surface1d.o : ${SRCDIR}/target/surface1d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/geom/curvicoord.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/surface1d.f90

${OBJDIR}/target/cell1d_nu.o : ${SRCDIR}/target/cell1d_nu.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/geom/curvicoord.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/target/surface1d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/cell1d_nu.f90

${OBJDIR}/math/solve.o : ${SRCDIR}/math/solve.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/solve.f90

${OBJDIR}/stat/his_adapt.o : ${SRCDIR}/stat/his_adapt.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/his_adapt.f90

${OBJDIR}/stat/his_extend.o : ${SRCDIR}/stat/his_extend.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/his_extend.f90

${OBJDIR}/stat/his1d.o : ${SRCDIR}/stat/his1d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/stat/his_adapt.f90 \
	${SRCDIR}/stat/his_extend.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/his1d.f90

${OBJDIR}/target/cell1d_reg.o : ${SRCDIR}/target/cell1d_reg.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/curvicoord.f90 \
	${SRCDIR}/math/solve.f90 \
	${SRCDIR}/stat/his1d.f90 \
	${SRCDIR}/target/cell1d_nu.f90 \
	${SRCDIR}/target/surface1d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/cell1d_reg.f90

${OBJDIR}/geom/cartgrid2d.o : ${SRCDIR}/geom/cartgrid2d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/math/sort.f90 \
	${SRCDIR}/geom/polygon_old.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/geom/cartgrid2d.f90

${OBJDIR}/stat/his2d.o : ${SRCDIR}/stat/his2d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/stat/his_adapt.f90 \
	${SRCDIR}/stat/his_extend.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/his2d.f90

${OBJDIR}/math/table.o : ${SRCDIR}/math/table.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/math/solve.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/table.f90

${OBJDIR}/target/surface2d_sdf.o : ${SRCDIR}/target/surface2d_sdf.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/cartgrid2d.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/stat/his2d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/surface2d_sdf.f90

${OBJDIR}/target/cell2d_reg.o : ${SRCDIR}/target/cell2d_reg.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/ffpcyl.f90 \
	${SRCDIR}/geom/cartgrid2d.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/geom/polygon.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/stat/his2d.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/cell2d_reg.f90

${OBJDIR}/target/cell_reg.o : ${SRCDIR}/target/cell_reg.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/target/cell1d_reg.f90 \
	${SRCDIR}/target/cell2d_reg.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/cell_reg.f90

${OBJDIR}/target/surface2d.o : ${SRCDIR}/target/surface2d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/geom/curvicoord.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/surface2d.f90

${OBJDIR}/target/grid2d_reg.o : ${SRCDIR}/target/grid2d_reg.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/geom/ffpcyl.f90 \
	${SRCDIR}/math/sort.f90 \
	${SRCDIR}/target/surface2d.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 \
	${SRCDIR}/geom/cartgrid2d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/grid2d_reg.f90

${OBJDIR}/target/grid2d_curvi.o : ${SRCDIR}/target/grid2d_curvi.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/geom/line.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/geom/ffpcyl.f90 \
	${SRCDIR}/math/sort.f90 \
	${SRCDIR}/target/surface2d.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 \
	${SRCDIR}/geom/cartgrid2d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/grid2d_curvi.f90

${OBJDIR}/target/grid2d.o : ${SRCDIR}/target/grid2d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/geom/polygon_old.f90 \
	${SRCDIR}/target/grid2d_reg.f90 \
	${SRCDIR}/target/grid2d_curvi.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/grid2d.f90

${OBJDIR}/target/cell2d_old.o : ${SRCDIR}/target/cell2d_old.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 \
	${SRCDIR}/target/grid2d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/cell2d_old.f90

${OBJDIR}/target/cell.o : ${SRCDIR}/target/cell.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/target/cell1d_nu.f90 \
	${SRCDIR}/target/cell_reg.f90 \
	${SRCDIR}/target/cell2d_old.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/cell.f90

${OBJDIR}/param/param.o : ${SRCDIR}/param/param.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/time.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/se.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/target/cell.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/param/param.f90

${OBJDIR}/sn/restffp.o : ${SRCDIR}/sn/restffp.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/restffp.f90

${OBJDIR}/xtal/ipos.o : ${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/geom/coordsys.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/ipos.f90

${OBJDIR}/xtal/oldsites.o : ${SRCDIR}/xtal/oldsites.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/xtal/ipos.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/oldsites.f90

${OBJDIR}/target/region1d.o : ${SRCDIR}/target/region1d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/geom/curvicoord.f90 \
	${SRCDIR}/target/surface1d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region1d.f90

${OBJDIR}/target/region2d_poly.o : ${SRCDIR}/target/region2d_poly.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/geom/polygon.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region2d_poly.f90

${OBJDIR}/target/region2d_sdf.o : ${SRCDIR}/target/region2d_sdf.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region2d_sdf.f90

${OBJDIR}/target/region2d.o : ${SRCDIR}/target/region2d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/geom/curvicoord.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/target/region2d_poly.f90 \
	${SRCDIR}/target/region2d_sdf.f90 \
	${SRCDIR}/target/surface2d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region2d.f90

${OBJDIR}/stat/his3d.o : ${SRCDIR}/stat/his3d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/stat/his_adapt.f90 \
	${SRCDIR}/stat/his_extend.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/his3d.f90

${OBJDIR}/target/surface3d_sdf.o : ${SRCDIR}/target/surface3d_sdf.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/stat/his3d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/surface3d_sdf.f90

${OBJDIR}/target/region3d_sdf.o : ${SRCDIR}/target/region3d_sdf.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/math/binsearch.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/target/surface3d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region3d_sdf.f90

${OBJDIR}/target/surface3d.o : ${SRCDIR}/target/surface3d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/target/surface3d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/surface3d.f90

${OBJDIR}/target/region3d.o : ${SRCDIR}/target/region3d.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/target/region3d_sdf.f90 \
	${SRCDIR}/target/surface3d.f90 \
	${SRCDIR}/target/surface3d_sdf.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region3d.f90

${OBJDIR}/target/region.o : ${SRCDIR}/target/region.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/target/region1d.f90 \
	${SRCDIR}/target/region2d.f90 \
	${SRCDIR}/target/region3d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/region.f90

${OBJDIR}/target/target.o : ${SRCDIR}/target/target.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/target/region.f90 \
	${SRCDIR}/target/cell.f90 \
	${SRCDIR}/geom/polygon_old.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/target.f90

${OBJDIR}/target/surface.o : ${SRCDIR}/target/surface.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/target/surface1d.f90 \
	${SRCDIR}/target/surface2d.f90 \
	${SRCDIR}/target/surface3d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/target/surface.f90

${OBJDIR}/casc/surface.o : ${SRCDIR}/casc/surface.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/surface.f90 \
	${SRCDIR}/target/target.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/surface.f90

${OBJDIR}/xtal/defects_diam.o : ${SRCDIR}/xtal/defects_diam.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/defects_diam.f90

${OBJDIR}/xtal/lattice_diam.o : ${SRCDIR}/xtal/lattice_diam.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/xtal.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/lattice_diam.f90

${OBJDIR}/xtal/lattice_orth.o : ${SRCDIR}/xtal/lattice_orth.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/xtal.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/lattice_orth.f90

${OBJDIR}/xtal/lattice.o : ${SRCDIR}/xtal/lattice.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/xtal/lattice_diam.f90 \
	${SRCDIR}/xtal/lattice_orth.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/lattice.f90

${OBJDIR}/eval/det.o : ${SRCDIR}/eval/det.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/stat/his1d.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/xtal/defects_diam.f90 \
	${SRCDIR}/xtal/lattice.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/det.f90

${OBJDIR}/sn/ffpcyl.o : ${SRCDIR}/sn/ffpcyl.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/ffpcyl.f90

${OBJDIR}/xtal/tvib.o : ${SRCDIR}/xtal/tvib.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/math/debye.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/tvib.f90

${OBJDIR}/ion/proj.o : ${SRCDIR}/ion/proj.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/casc/surface.f90 \
	${SRCDIR}/eval/det.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/sn/ffpcyl.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/surface.f90 \
	${SRCDIR}/target/surface.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/xtal/lattice.f90 \
	${SRCDIR}/xtal/tvib.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/ion/proj.f90

${OBJDIR}/stat/mom1d.o : ${SRCDIR}/stat/mom1d.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/mom1d.f90

${OBJDIR}/stat/mom2d.o : ${SRCDIR}/stat/mom2d.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/mom2d.f90

${OBJDIR}/stat/his_util.o : ${SRCDIR}/stat/his_util.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/stat/his1d.f90 \
	${SRCDIR}/stat/his2d.f90 \
	${SRCDIR}/stat/his3d.f90 \
	${SRCDIR}/target/surface2d.f90 \
	${SRCDIR}/target/surface3d.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/stat/his_util.f90

${OBJDIR}/eval/rrobin.o : ${SRCDIR}/eval/rrobin.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/rrobin.f90

${OBJDIR}/eval/stat.o : ${SRCDIR}/eval/stat.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/geom/polygon.f90 \
	${SRCDIR}/stat/mom1d.f90 \
	${SRCDIR}/stat/mom2d.f90 \
	${SRCDIR}/stat/his1d.f90 \
	${SRCDIR}/stat/his2d.f90 \
	${SRCDIR}/stat/his3d.f90 \
	${SRCDIR}/stat/his_util.f90 \
	${SRCDIR}/eval/rrobin.f90 \
	${SRCDIR}/ion/proj.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/misc/time.f90 \
	${SRCDIR}/eval/det.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/stat.f90

${OBJDIR}/eval/casc.o : ${SRCDIR}/eval/casc.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/eval/stat.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/xtal/lattice.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/casc.f90

${OBJDIR}/eval/dam.o : ${SRCDIR}/eval/dam.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/eval/stat.f90 \
	${SRCDIR}/target/target.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/dam.f90

${OBJDIR}/eval/flow.o : ${SRCDIR}/eval/flow.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/stack/stack.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/eval/dam.f90 \
	${SRCDIR}/eval/stat.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/target/cell.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/flow.f90

${OBJDIR}/se/lind.o : ${SRCDIR}/se/lind.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/lind.f90

${OBJDIR}/se/zbl.o : ${SRCDIR}/se/zbl.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/zbl.f90

${OBJDIR}/se/kkk.o : ${SRCDIR}/se/kkk.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/kkk.f90

${OBJDIR}/se/table.o : ${SRCDIR}/se/table.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/table.f90

${OBJDIR}/se/rand.o : ${SRCDIR}/se/rand.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/se.f90 \
	${SRCDIR}/se/lind.f90 \
	${SRCDIR}/se/zbl.f90 \
	${SRCDIR}/se/kkk.f90 \
	${SRCDIR}/se/table.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/rand.f90

${OBJDIR}/se/rbs.o : ${SRCDIR}/se/rbs.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/se/rand.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/rbs.f90

${OBJDIR}/eval/rbs.o : ${SRCDIR}/eval/rbs.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/stat/his1d.f90 \
	${SRCDIR}/se/rbs.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/rbs.f90

${OBJDIR}/eval/tra.o : ${SRCDIR}/eval/tra.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/ion/proj.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/tra.f90

${OBJDIR}/eval/eval.o : ${SRCDIR}/eval/eval.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/geom/bndcond.f90 \
	${SRCDIR}/target/region.f90 \
	${SRCDIR}/target/cell.f90 \
	${SRCDIR}/target/surface2d_sdf.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/eval/casc.f90 \
	${SRCDIR}/eval/flow.f90 \
	${SRCDIR}/eval/rbs.f90 \
	${SRCDIR}/eval/stat.f90 \
	${SRCDIR}/eval/tra.f90 \
	${SRCDIR}/eval/dam.f90 \
	${SRCDIR}/eval/det.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/eval/eval.f90

${OBJDIR}/casc/check.o : ${SRCDIR}/casc/check.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/casc/surface.f90 \
	${SRCDIR}/eval/eval.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/target/target.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/check.f90

${OBJDIR}/casc/branch.o : ${SRCDIR}/casc/branch.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/eval/eval.f90 \
	${SRCDIR}/eval/tra.f90 \
	${SRCDIR}/xtal/oldsites.f90 \
	${SRCDIR}/sn/ffpcyl.f90 \
	${SRCDIR}/sn/restffp.f90 \
	${SRCDIR}/casc/check.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/branch.f90

${OBJDIR}/sn/ndira.o : ${SRCDIR}/sn/ndira.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/random/random.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/ndira.f90

${OBJDIR}/xtal/atoms_det.o : ${SRCDIR}/xtal/atoms_det.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/eval/det.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/xtal/ipos.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/atoms_det.f90

${OBJDIR}/xtal/atoms_stat.o : ${SRCDIR}/xtal/atoms_stat.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/xtal/lattice.f90 \
	${SRCDIR}/eval/stat.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/atoms_stat.f90

${OBJDIR}/xtal/atoms.o : ${SRCDIR}/xtal/atoms.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/xtal/lattice.f90 \
	${SRCDIR}/xtal/atoms_det.f90 \
	${SRCDIR}/xtal/atoms_stat.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/xtal/atoms.f90

${OBJDIR}/casc/selncc.o : ${SRCDIR}/casc/selncc.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/surface.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/sn/ndira.f90 \
	${SRCDIR}/eval/stat.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/xtal/lattice.f90 \
	${SRCDIR}/xtal/atoms.f90 \
	${SRCDIR}/xtal/tvib.f90 \
	${SRCDIR}/xtal/oldsites.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/selncc.f90

${OBJDIR}/math/bessel.o : ${SRCDIR}/math/bessel.f90 \
	${SRCDIR}/misc/const.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/bessel.f90

${OBJDIR}/sn/ffp.o : ${SRCDIR}/sn/ffp.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/math/bessel.f90 \
	${SRCDIR}/target/target.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/ffp.f90

${OBJDIR}/casc/selnca.o : ${SRCDIR}/casc/selnca.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/sn/ffp.f90 \
	${SRCDIR}/sn/ffpcyl.f90 \
	${SRCDIR}/sn/restffp.f90 \
	${SRCDIR}/sn/ndira.f90 \
	${SRCDIR}/xtal/ipos.f90 \
	${SRCDIR}/xtal/tvib.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/surface.f90 \
	${SRCDIR}/target/target.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/selnca.f90

${OBJDIR}/casc/selncp.o : ${SRCDIR}/casc/selncp.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/geom.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/casc/selncc.f90 \
	${SRCDIR}/casc/selnca.f90 \
	${SRCDIR}/sn/restffp.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/selncp.f90

${OBJDIR}/se/firs.o : ${SRCDIR}/se/firs.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/firs.f90

${OBJDIR}/se/or.o : ${SRCDIR}/se/or.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/or.f90

${OBJDIR}/se/loc.o : ${SRCDIR}/se/loc.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 \
	${SRCDIR}/se/firs.f90 \
	${SRCDIR}/se/or.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/loc.f90

${OBJDIR}/se/straggle.o : ${SRCDIR}/se/straggle.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/straggle.f90

${OBJDIR}/se/mixed.o : ${SRCDIR}/se/mixed.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/se.f90 \
	${SRCDIR}/se/loc.f90 \
	${SRCDIR}/se/rand.f90 \
	${SRCDIR}/se/straggle.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/target/target.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/mixed.f90

${OBJDIR}/se/se.o : ${SRCDIR}/se/se.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/se.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/se/mixed.f90 \
	${SRCDIR}/target/surface.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/se/se.f90

${OBJDIR}/sn/screen.o : ${SRCDIR}/sn/screen.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/misc/undefine.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/sn.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/screen.f90

${OBJDIR}/sn/apsis.o : ${SRCDIR}/sn/apsis.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/stack/real8.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/math/table.f90 \
	${SRCDIR}/sn/screen.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/apsis.f90

${OBJDIR}/sn/cmangm.o : ${SRCDIR}/sn/cmangm.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/sn/screen.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/cmangm.f90

${OBJDIR}/sn/cmangl.o : ${SRCDIR}/sn/cmangl.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/cmangl.f90

${OBJDIR}/math/quadrature.o : ${SRCDIR}/math/quadrature.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/math/quadrature.f90

${OBJDIR}/sn/cmscatter.o : ${SRCDIR}/sn/cmscatter.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/sn/apsis.f90 \
	${SRCDIR}/sn/screen.f90 \
	${SRCDIR}/math/quadrature.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/cmscatter.f90

${OBJDIR}/sn/timeint.o : ${SRCDIR}/sn/timeint.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/misc/open.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/timeint.f90

${OBJDIR}/sn/scatter.o : ${SRCDIR}/sn/scatter.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/vector.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/elem.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/se.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/se/se.f90 \
	${SRCDIR}/sn/apsis.f90 \
	${SRCDIR}/sn/cmangm.f90 \
	${SRCDIR}/sn/cmangl.f90 \
	${SRCDIR}/sn/cmscatter.f90 \
	${SRCDIR}/sn/screen.f90 \
	${SRCDIR}/sn/timeint.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/surface.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/sn/scatter.f90

${OBJDIR}/casc/step.o : ${SRCDIR}/casc/step.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/param/xtal.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/param/dam.f90 \
	${SRCDIR}/param/sn.f90 \
	${SRCDIR}/param/out.f90 \
	${SRCDIR}/casc/selncp.f90 \
	${SRCDIR}/casc/branch.f90 \
	${SRCDIR}/sn/scatter.f90 \
	${SRCDIR}/sn/ffpcyl.f90 \
	${SRCDIR}/xtal/oldsites.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/step.f90

${OBJDIR}/casc/casc.o : ${SRCDIR}/casc/casc.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/misc/array.f90 \
	${SRCDIR}/misc/types.f90 \
	${SRCDIR}/misc/string.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/ion.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/sn/restffp.f90 \
	${SRCDIR}/xtal/oldsites.f90 \
	${SRCDIR}/ion/proj.f90 \
	${SRCDIR}/casc/branch.f90 \
	${SRCDIR}/casc/step.f90 \
	${SRCDIR}/eval/eval.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/casc/casc.f90

${OBJDIR}/main.o : ${SRCDIR}/main.f90 \
	${SRCDIR}/misc/const.f90 \
	${SRCDIR}/param/param.f90 \
	${SRCDIR}/param/setup.f90 \
	${SRCDIR}/param/mat.f90 \
	${SRCDIR}/casc/casc.f90 \
	${SRCDIR}/eval/eval.f90 \
	${SRCDIR}/geom/coordsys.f90 \
	${SRCDIR}/math/intmath.f90 \
	${SRCDIR}/random/random.f90 \
	${SRCDIR}/se/se.f90 \
	${SRCDIR}/sn/ffp.f90 \
	${SRCDIR}/sn/ffpcyl.f90 \
	${SRCDIR}/sn/restffp.f90 \
	${SRCDIR}/sn/scatter.f90 \
	${SRCDIR}/target/target.f90 \
	${SRCDIR}/target/surface.f90 \
	${SRCDIR}/xtal/lattice.f90 \
	${SRCDIR}/xtal/oldsites.f90 \
	${SRCDIR}/xtal/tvib.f90 
	${FC} ${FFLAGS} -o $@ -c ${SRCDIR}/main.f90