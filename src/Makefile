HTSLIB=/env/cns/bigtmp2/ONT/gitlab_nap/nanopore_assembly_pipeline/tools/samtools/htslib-1.10

C_FILE=alipile.c polished.c hash.c
OBJ = ${C_FILE:.c=.o}
MAIN_FILE=polish_consensus.c
BIN=polish_consensus

COMPILER=cc
CFLAGS += -I$(HTSLIB) -L$(HTSLIB)
CFLAGS += -O3
CFLAGS += -g


polisher: $(OBJ)
	${COMPILER} ${CFLAGS} -o $(BIN) $(MAIN_FILE) $(OBJ) -lm -lhts -lz -llzma -lbz2

clean:	
	rm -rf *.o a.out 
