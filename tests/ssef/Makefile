L3NAME     = ssef

include ../../config/$(L3NAME).conf

CEXTRA = $(shell $(QMP_TOP)/bin/qmp-config --ldflags)

LIBS   = -lqdp_d3 -lqdp_int -lqdp_common \
         -lqla_d -lqla_d3 -lqla_int -lqla_random \
         -lqio \
         -llime -lm \
         $(shell $(QMP_TOP)/bin/qmp-config --libs)

sources = main.c bl.c

programs = vec-test

.PHONY: all clean realclean

all: $(programs)

clean:
	$(RM) $(sources:%.c=%.o)

realclean: clean
	$(RM) $(programs)

$(programs): $(sources:%.c=%.o)
	$(XCC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -o $@ $^ $(LIBS)

main.o: main.c
	$(XCC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -c -o $@ $<

bl.o: bl.c
	$(CC) $(CFLAGS) $(CINCDIRS) $(COPT) $(CDEFS) $(CEXTRA) -c -o $@ $<
