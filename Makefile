# NOTE: needs boost, tclap, and sdsl
CXX=g++ #clang++ # g++
CPP_FLAGS=-g -m64 -std=c++0x -W -Wall -Wextra -Wpointer-arith -Wcast-qual \
					  -Wwrite-strings \
#					-Wbool-conversions -Wshift-overflow -Wliteral-conversion \
					-Werror -W -fno-strict-aliasing
BOOST_PATH=/ufrc/boucher/baharpan/VARI/3rd_party_inst/boost
DEP_PATH=3rd_party_inst
INC_PATH=-isystem $(DEP_PATH)/include -isystem $(BOOST_PATH)/include
LIB_PATH=-L$(DEP_PATH)/lib -L./ -L$(BOOST_PATH)/lib
MM_PATH=/s/chopin/k/grad/baharpan
KMC_PATH=3rd_party_src/KMC
BOOST_FLAGS= -lboost_system -lboost_filesystem

DEP_FLAGS=$(INC_PATH) $(LIB_PATH) $(BOOST_FLAGS) -isystem $(KMC_PATH)  -lsdsl -fopenmp
DEBUG_FLAGS=-pg -gstabs
NDEBUG_FLAGS= -DNDEBUG
OPT_FLAGS= -O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 
NOPT_FLAGS=-O0
NUM_COLS=64
# Using Semantic Versioning: http://semver.org/
VERSION=0.5.1
CPP_FLAGS+=-DVERSION=\"$(VERSION)\" -DNUM_COLS=$(NUM_COLS)

ifeq ($(optimise),0)
CPP_FLAGS+=$(NOPT_FLAGS)
else
CPP_FLAGS+=$(OPT_FLAGS)
endif

ifeq ($(debug),1)
CPP_FLAGS+=$(DEBUG_FLAGS)
else
CPP_FLAGS+=$(NDEBUG_FLAGS)
endif

ifeq ($(verbose),1)
CPP_FLAGS+=-DVERBOSE
endif

ifneq ($(revcomps),0)
CPP_FLAGS+=-DADD_REVCOMPS
endif

ifneq ($(dummies),0)
CPP_FLAGS+=-DALL_DUMMIES
endif

ifeq ($(varord),1)
CPP_FLAGS+=-DVAR_ORDER
endif

BUILD_REQS=debruijn_graph.hpp io.hpp io.o debug.h 
COLOR_REQS=colored_debruijn_graph.hpp io.hpp io.o debug.h 
ASSEM_REQS=debruijn_graph.hpp algorithm.hpp utility.hpp kmer.hpp uint128_t.hpp
PACK_REQS=lut.hpp debug.h io.hpp io.o sort.hpp kmer.hpp dummies.hpp bubbles.hpp reduce_color.hpp
BINARIES= reduce_color#bubbles #cosmo-pack matrix fast_matrix fast_reduce_color  cosmo-build cosmo-color cosmo-benchmark pack-color #bv2repair  transpose  uncompress#match-color # cosmo-assemble

KMC_OBJS=3rd_party_src/KMC/kmc_api/kmc_file.o 3rd_party_src/KMC/kmc_api/kmer_api.o 3rd_party_src/KMC/kmc_api/mmer.o

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(CXX) $(CPP_FLAGS) -c io.cpp  $(DEP_FLAGS) 

# TODO: Roll these all into one... "cosmo"
cosmo-pack: cosmo-pack.cpp $(PACK_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 

#bubbles: bubbles.cpp $(PACK_REQS)						
#	$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS)
#bv2repair: bv2repair.cpp $(PACK_REQS)
#	$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS)
#transpose: transpose.cpp $(PACK_REQS)
#	$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS)
#uncompress: uncompress.cpp $(PACK_REQS)
#	$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS)
#matrix: matrix.cpp $(PACK_REQS)
#		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 
#fast_matrix: fast_matrix.cpp $(PACK_REQS)
#		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 
reduce_color: reduce_color.cpp $(PACK_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 

#fast_reduce_color: fast_reduce_color.cpp $(PACK_REQS)
#	$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS)

#cosmo-build: cosmo-build.cpp $(BUILD_REQS)
#		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 

#cosmo-color: cosmo-color.cpp $(BUILD_REQS)
#		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 

#pack-color: pack-color.cpp $(BUILD_REQS)
#		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(KMC_OBJS) $(DEP_FLAGS) 

#match-color: match-color.cpp $(BUILD_REQS)
#		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(DEP_FLAGS) 

#cosmo-assemble: cosmo-assemble.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
#		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) 

#cosmo-benchmark: cosmo-benchmark.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
#		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) 

all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
