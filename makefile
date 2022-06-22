CC=gcc
CFLAGS+= -m64 -g -Wall
LDFLAGS=
LIBS=
INC= 
ALLDEP:= $(patsubst %.h,%.o,$(wildcard *.h))
ALLILP:= $(wildcard *_ILP.c)

.PHONY: test

test: $(ALLDEP) test.out partitionByLayers.out

product: CFLAGS = -O3
product: $(ALLDEP) $(patsubst %.c,%.out,$(filter-out $(patsubst %.o,%.c,$(ALLDEP)) $(ALLILP), $(wildcard *.c)))

%.out: %.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
%.o: %.c makefile
	$(CC) $(CFLAGS) -MMD -c $< -o $@

bmerVector.out: bmerVector.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) -lm

clique%.MaxID: clique%.txt
	grep -oE 'has [0-9]+' $^ | sed 's/has //' > $@; \
	awk -F' ' 'NR>4 {if(NF){sub(/[ \t]+$$/, ""); print NF, $$0}else{exit}}' $^ >> $@

%.MaxID: %.txt
	grep -oE '[0-9]+' $^ | tail -1 > $@; \
	sed -n '7 s/\r$$//p' $^ | sed 's/ /\n/g' >> $@

#gurobi make
%_ILP: %_ILP.c $(ALLDEP)
	$(CC) $(CFLAGS) -o $@ $^ -I$(INC) $(LDFLAGS) -lm

.PHONY: clean

clean:
	rm -rf *.out *.o *.dSYM *.d

.PHONY: realclean clean-backups

realclean: clean clean-backups

clean-backups:
	rm -rf *~ #*#     
