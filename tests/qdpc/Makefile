L3NAME     = ssef

include ../../config/$(L3NAME).conf

L3OBJECT   = ../../dwf-$(L3NAME).o

CDEFS  = -DPREFIX=MIT_$(L3NAME)_
QMP_TOP=/home/bj/Devel/QCD/install/qmp/HEAD/single
CEXTRA = $(shell $(QMP_TOP)/bin/qmp-config --ldflags)

LIBS   = -lqdp_d3 -lqdp_int -lqdp_common \
         -lqla_d -lqla_d3 -lqla_int -lqla_random \
         -lqio \
         -llime -lm \
         $(shell $(QMP_TOP)/bin/qmp-config --libs)

sources = cg-solver.c cg-solver1.c dwf-op.c

programs = $(sources:%.c=%-$(L3NAME))

.PHONY: all clean realclean

all: $(programs)

clean:
	$(RM) $(sources:%.c=%-*.o)

realclean: clean
	$(RM) $(programs)

$(programs): % : %.o $(L3OBJECT)
	$(XCC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -o $@ $^ $(LIBS)

$(sources:%.c=%-$(L3NAME).o): %-$(L3NAME).o: %.c
	$(XCC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -c -o $@ $<

.PHONY: $(L3OBJECT)

$(L3OBJECT):
	make -C ../.. CONFIG=config/$(L3NAME).conf compile
