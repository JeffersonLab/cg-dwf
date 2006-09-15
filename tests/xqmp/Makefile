L3NAME     = altivecf
CONFIG     = $(L3NAME)
QDP_TOP    = $$HOME/QDP-C
QIO_TOP    = $(QDP_TOP)
QLA_TOP    = $(QDP_TOP)
QMP_TOP    = $(QDP_TOP)

L3OBJECT   = ../../dwf-${L3NAME}.o
CC     = gcc

CFLAGS = -Wall -g -DPREFIX=MIT_$(L3NAME)_ \
         -I$(QDP_TOP)/include -L$(QDP_TOP)/lib \
         -I$(QLA_TOP)/include -L$(QLA_TOP)/lib \
         -I$(QIO_TOP)/include -L$(QIO_TOP)/lib \
         -I../.. \
         $(shell $(QMP_TOP)/bin/qmp-config --ldflags)
LIBS   = -lqdp_d3 -lqdp_int -lqdp_common \
         -lqla_d -lqla_d3 -lqla_int -lqla_random \
         -lqio \
         -llime -lm

source = qmp-check.c qmp.c

program = qmp-check

.PHONY: all clean realclean

all: $(program)

clean:
	$(RM) $(source:%.c=%.o)

realclean: clean
	$(RM) $(program)

$(program): $(source:%.c=%.o) $(L3OBJECT)
	$(CC) $(CFLAGS) -o $(program) $(source:%.c=%.o) $(L3OBJECT) $(LIBS)

$(source:%.c=%.o): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<


.PHONY: $(L3OBJECT)

$(L3OBJECT):
	make -C ../.. CONFIG=config/$(CONFIG).conf compile
