L3NAME     = altivecf
CONFIG     = $(L3NAME)

include ../../config/$(CONFIG).conf

L3OBJECT   = ../../dwf-$(L3NAME).o

CDEFS  = -DPREFIX=MIT_$(L3NAME)_
CEXTRA = $(shell $(QMP_TOP)/bin/qmp-config --ldflags)

LIBS   = -lqdp_d3 -lqdp_int -lqdp_common \
         -lqla_d -lqla_d3 -lqla_int -lqla_random \
         -lqio \
         -llime -lm \
         $(shell $(QMP_TOP)/bin/qmp-config --libs)

sources = cg-solver.c \
          cg-solver1.c \
          dwf-op.c \
          dwf-op1.c \
          dwf-opx.c \
          dwf-opx1.c \
          fermions.c \
          dwf-dd.c \
          dwf-dd1.c

programs = $(sources:%.c=%-$(CONFIG))

.PHONY: all clean realclean

all: $(programs)

clean:
	$(RM) $(sources:%.c=%-*.o)

realclean: clean
	$(RM) $(programs)

$(programs): % : %.o $(L3OBJECT)
	$(XCC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -o $@ $^ $(LIBS)

$(sources:%.c=%-$(CONFIG).o): %-$(CONFIG).o: %.c
	$(XCC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -c -o $@ $<

.PHONY: $(L3OBJECT)

$(L3OBJECT):
	make -C ../.. CONFIG=config/$(CONFIG).conf compile
