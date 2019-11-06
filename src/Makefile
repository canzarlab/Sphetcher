ifdef LONG
INTT = -DLONG
endif

PCC = g++
PLFLAGS = $(LFLAGS)
PCFLAGS = -O3 $(INTT)

# BNCHMRK = setCover
BNCHMRK = sphetcher

# generic toplevel common files needed by this and most timing code 
#GLOBAL_GENERIC = utils.h parseCommandLine.h IO.h parallel.h gettime.h sequence.h runTests.py 

# toplevel common files needed by the timing code
TIME_GLOBAL_REQUIRE = graph.h graphIO.h $(GLOBAL_GENERIC) 

# local common files needed by the timing code
TIME_REQUIRE = $(BNCHMRK).h 

# required files
GLOBAL_REQUIRE = utils.h parallel.h sequence.h graph.h blockRadixSort.h transpose.h
LOCAL_REQUIRE = 
OBJS = setCover.o util_functions.o
# ********************

CHECK = $(BNCHMRK)Check
TIME = $(BNCHMRK)Time
TEST_FILES = testInputs

# Make the benchmark, and make the output check code in the common directory
all : $(BNCHMRK) #testInputs  #my cd ../common; make -s $(CHECK)

$(TEST_FILES) : #cp ../common/$@ .

# Files copied (linked) from toplevel common directory
# The sort used to remove duplicates
GLOBAL = $(sort $(TIME_GLOBAL_REQUIRE) $(GLOBAL_REQUIRE))

$(GLOBAL) : #ln -s ../../common/$@ .

# Files copied (linked) from benchmark common directory
BENCH = $(TIME).C $(sort $(TIME_REQUIRE) $(BENCH_REQUIRE))

$(BENCH) : #ln -s ../common/$@ .

# Make all implementation objects
%.o : %.C $(GLOBAL_REQUIRE) $(BENCH_REQUIRE) $(LOCAL_REQUIRE)
	$(PCC) $(PCFLAGS) -c $< -o $@ -fopenmp

# Make timing code object file
$(TIME).o : $(TIME).C $(TIME_GLOBAL_REQUIRE) $(TIME_REQUIRE)
	$(PCC) $(PCFLAGS) -c $< -o $@ -fopenmp

# Make benchmark
$(BNCHMRK) : $(TIME).o $(OBJS)
	$(PCC) $(PLFLAGS) -o $@ $(TIME).o $(OBJS) -fopenmp

clean :
	rm -f $(BNCHMRK) *.o *.pyc

cleansrc :
	make -s clean
	rm -f $(GLOBAL) $(BENCH) $(TEST_FILES)

