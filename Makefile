CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
    CC := $(shell brew list --versions gcc | awk '{print $$2}' | cut -d'.' -f1 | awk '{print "g++-"$$1; exit}')
endif

OMPFLG=-fopenmp
HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched
OBJ_EXT=src/ExtractRNA.o
EXE_TAX=bin/PM-parallel-meta
EXE_RNA=bin/PM-extract-rna
EXE_CLT=bin/PM-plot-taxa
EXE_CLF=bin/PM-predict-func
EXE_CFN=bin/PM-predict-func-nsti
EXE_CFC=bin/PM-predict-func-contribute
EXE_TAL=bin/PM-select-taxa
EXE_FUL=bin/PM-select-func
EXE_CMP=bin/PM-comp-taxa
EXE_CMF=bin/PM-comp-func
EXE_CMC=bin/PM-comp-corr
EXE_PIP=bin/PM-pipeline
EXE_SSQ=bin/PM-split-seq
EXE_FSQ=bin/PM-format-seq
EXE_RCV=bin/PM-rare-curv
EXE_RAR=bin/PM-rand-rare
EXE_UTX=bin/PM-update-taxa

tax:$(OBJ_TAX) src/frame.cpp
	$(CC) -c -o $(OBJ_EXT) src/ExtractRNA.cpp $(HASHFLG)
	$(CC) -o $(EXE_TAX) src/frame.cpp $(OBJ_MAP) $(OBJ_EXT) $(OMPFLG) $(HASHFLG)
	$(CC) -o $(EXE_RNA) src/ExtractRNA_plus.cpp $(OBJ_EXT) $(HASHFLG)	
	$(CC) -o $(EXE_CLT) src/class_tax.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CLF) src/class_func.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CFN) src/class_func_nsti.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CFC) src/class_func_contribute.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_TAL) src/taxa_sel.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_FUL) src/func_sel.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CMP) src/comp_sam.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CMF) src/comp_sam_func.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CMC) src/comp_corr.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_PIP) src/pipeline.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_SSQ) src/split_seq.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_FSQ) src/format_seq.cpp $(HASHFLG) $(BUILDFLG)
	$(CC) -o $(EXE_RCV) src/rare_curv.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_RAR) src/rand_rare.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_UTX) src/update_taxa.cpp $(HASHFLG) $(BUILDFLG)

clean:
	rm -rf bin/PM-* src/*.o
