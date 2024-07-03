CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
    CC := $(shell brew list --versions gcc | awk '{print $$2}' | cut -d'.' -f1 | awk '{print "g++-"$$1; exit}')
endif

OMPFLG=-fopenmp
BLDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -std=c++11 -g -O3
EXE=../bin/PM-profiler

all:
	$(CC) -o $(EXE) src/main.cpp $(OMPFLG) $(BLDFLG) 

clean:
	rm -rf $(EXE)
