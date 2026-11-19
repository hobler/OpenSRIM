SRCDIR=~/git/imsil/code
EXEDIR=~/work/imsil/master-omp/code

OBJDIR=${EXEDIR}/obj_db
MODDIR=${EXEDIR}/mod_db
LIBDIR=${EXEDIR}/lib_db
EXE=${EXEDIR}/unittests/sn/cmscatter/test_db.exe

$(shell mkdir -p ${EXEDIR})
$(shell mkdir -p ${OBJDIR})
$(shell mkdir -p ${MODDIR})
$(shell mkdir -p ${EXEDIR}/unittests/sn/cmscatter)

$(shell mkdir -p ${OBJDIR}/eval)
$(shell mkdir -p ${OBJDIR}/se)
$(shell mkdir -p ${OBJDIR}/cinterface)
$(shell mkdir -p ${OBJDIR}/sn)
$(shell mkdir -p ${OBJDIR}/random)
$(shell mkdir -p ${OBJDIR}/misc)
$(shell mkdir -p ${OBJDIR}/geom)
$(shell mkdir -p ${OBJDIR}/target)
$(shell mkdir -p ${OBJDIR}/casc)
$(shell mkdir -p ${OBJDIR}/list)
$(shell mkdir -p ${OBJDIR}/math)
$(shell mkdir -p ${OBJDIR}/stack)
$(shell mkdir -p ${OBJDIR}/param)
$(shell mkdir -p ${OBJDIR}/ion)
$(shell mkdir -p ${OBJDIR}/stat)
$(shell mkdir -p ${OBJDIR}/unittests)
$(shell mkdir -p ${OBJDIR}/unittests/sn)
$(shell mkdir -p ${OBJDIR}/unittests/sn/cmscatter)
$(shell mkdir -p ${OBJDIR}/unittests/sn/apsis)
$(shell mkdir -p ${OBJDIR}/unittests/sn/screen)
$(shell mkdir -p ${OBJDIR}/unittests/math)
$(shell mkdir -p ${OBJDIR}/unittests/math/solve)
$(shell mkdir -p ${OBJDIR}/unittests/math/quadrature)
$(shell mkdir -p ${OBJDIR}/unittests/math/table)
$(shell mkdir -p ${OBJDIR}/unittests/ion)
$(shell mkdir -p ${OBJDIR}/unittests/ion/proj)
$(shell mkdir -p ${OBJDIR}/xtal)