GUROBI_INSTALL=$(shell pwd)/build/gurobi910/linux64
BCFTOOLS_INSTALL=$(shell pwd)/build/bcftools-1.9/bcftools
SAMTOOLS_INSTALL=$(shell pwd)/build/samtools-1.12/samtools
IGRAPH_INSTAll=$(shell pwd)/build/igraph-0.10.2
TARGET_DIR=$(shell pwd)/build
CPPFLAGS= -g -std=c++2a -DNDEBUG -O3 

all:
	mkdir -p build
	$(CXX) $(CPPFLAGS) -D BCFTOOLSPATH=$(BCFTOOLS_INSTALL) -D SAMTOOLSPATH=$(SAMTOOLS_INSTALL) -m64 -o $(TARGET_DIR)/data_wrangler -I src/ext  src/data_wrangler.cpp -lm -lz -lpthread
	$(CXX) $(CPPFLAGS) -o $(TARGET_DIR)/main_ilp -I $(GUROBI_INSTALL)/include/ -I $(IGRAPH_INSTAl)/igraph-0.10.2/include/ -I $(IGRAPH_INSTAl)/build3/include/ -L $(GUROBI_INSTALL)/lib/ -L $(IGRAPH_INSTAl)/build3/src/ src/main.cpp -lgurobi_c++ -ligraph $(GUROBI_INSTALL)/lib/libgurobi91.so -lm -lz -lpthread -fopenmp

	@echo "check two executables in build directory: data_wrangler, and main_ilp"

clean:
	rm -f $(TARGET_DIR)/data_wrangler
	rm -f $(TARGET_DIR)/main_ilp