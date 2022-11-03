GUROBI_INSTALL=$(shell pwd)/build/gurobi910/linux64
BCFTOOLS_INSTALL=$(shell pwd)/build/bcftools-1.9/bcftools
SAMTOOLS_INSTALL=$(shell pwd)/build/samtools-1.12/samtools
TARGET_DIR=$(shell pwd)/build
CPPFLAGS= -g -std=c++11 -DNDEBUG -O3 

all:
	mkdir -p build
	$(CXX) $(CPPFLAGS) -D BCFTOOLSPATH=$(BCFTOOLS_INSTALL) -D SAMTOOLSPATH=$(SAMTOOLS_INSTALL) -m64 -o $(TARGET_DIR)/ilp_snp_indels -I $(GUROBI_INSTALL)/include/ -I src/ext -L  $(GUROBI_INSTALL)/lib/ src/ilp_snp_indels.cpp -lgurobi_c++ $(GUROBI_INSTALL)/lib/libgurobi91.so -lm -lz -lpthread
	$(CXX) $(CPPFLAGS) -D BCFTOOLSPATH=$(BCFTOOLS_INSTALL) -D SAMTOOLSPATH=$(SAMTOOLS_INSTALL) -m64 -o $(TARGET_DIR)/reachable_subgraph -I $(GUROBI_INSTALL)/include/ -I format-benchmark/fmt -lfmt  -L  $(GUROBI_INSTALL)/lib/ src/reachable_subgraph.cpp -lgurobi_c++ $(GUROBI_INSTALL)/lib/libgurobi91.so -lm -lz -lpthread

	@echo "check executables in build directory"


clean:
	rm -f $(TARGET_DIR)/ilp_snp_indels
