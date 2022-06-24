CC=gcc
CFLAGS+= -m64 -Wall -O3
LDFLAGS=
LIBS= -Ilib
INC= 
ALLDEP:= $(patsubst %.h,%.o,$(wildcard lib/*.h))

.PHONY: all

all: $(ALLDEP) $(patsubst src/%.c,%.out, $(wildcard src/*.c))

%.out: src/%.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: lib/%.c makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d
	rm -rf lib/*.o lib/*.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#
