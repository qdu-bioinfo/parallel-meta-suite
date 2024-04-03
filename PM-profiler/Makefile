CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
	exist = $(shell if [ -e '/usr/local/bin/g++-11' ]; then echo "exist"; else echo "notexist"; fi;)
	ifeq ($(exist),exist)
		CC:=g++-11
	else
        	exist = $(shell if [ -e '/usr/local/bin/g++-10' ]; then echo "exist"; else echo "notexist"; fi;)
        	ifeq ($(exist),exist)
                	CC:=g++-10
		else
			CC:=g++-9
		endif
	endif
endif
OMPFLG=-fopenmp
BLDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse -std=c++11 -g -O3
EXE=../bin/PM-profiler

all:
	$(CC) -o $(EXE) src/main.cpp $(OMPFLG) $(BLDFLG) 

clean:
	rm -rf $(EXE)
