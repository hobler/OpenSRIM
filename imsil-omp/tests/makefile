# ---------------------------------------------------------------------
# Makefile for running tests
#
# To run all tests with <n> simultaneous jobs:  make -j <n> all
# To use the parallel version of IMSIL, e.g. with 8 images: make ... t=8
#
# 04-Feb-2016  G.Hobler
# 27-Aug-2021  G.Hobler  remove -b from imsil and mimsil recipes
# ---------------------------------------------------------------------
t?=0

IMSIL=imsil
MIMSIL=mimsil

.PHONY: all

ifeq (${t}, 0)
EXE=../code/imsil.exe
else
EXE=../code/imsil_omp.exe
endif

print :
	@echo EXE=${EXE}

# Suffix rules --------------------------------------------------------

.SUFFIXES :
.SUFFIXES : .log .tem .par .out .inp
%.out : %.inp ${EXE}
	${IMSIL} -s -i -t ${t} $<

%.log : %.tem %.par ${EXE}
	${MIMSIL} -s -i -l ${<:.tem=.log} -t ${t} $<

# Targets -------------------------------------------------------------

all : | krdyn bspot arsputter asdamage asneti b2d bchannel bdamage bdiverg bse btilt btiltx hchannel highstop multilayer pchannel rbs rrobin

arsputter : arsputter/arsputter.log

asdamage : asdamage/asdamage.log

asneti : asneti/asneti.log

b2d : b2d/b2d.out b2d/b2d_pbc.out b2d/btrench.log b2d/btrench3.out b2d/btrench4.out

bchannel : bchannel/bchannel.log

bdiverg : bdiverg/bdiverg.log

bdamage : bdamage/bdamage.log

bse : bse/bse.log

bspot : bspot/bspot.log

btilt : btilt/btilt.log

btiltx : btiltx/btiltx.log

dyn2d : dyn2d/membrane.out

hchannel : hchannel/hchannel.log

highstop : highstop/highstop.log

krdyn : krdyn/krdyn.log

multilayer : multilayer/multilayer1d.log multilayer/multilayer2d.log

pchannel : pchannel/pchannel.log

rbs : rbs/rbs.log

rrobin : rrobin/rrobin.log rrobin/rrobin_screen.log rrobin/rrobin_quad.log rrobin/rrobin_coef.log rrobin/rrobin.out rrobin/rrobin_cub.out

# Dependencies --------------------------------------------------------

