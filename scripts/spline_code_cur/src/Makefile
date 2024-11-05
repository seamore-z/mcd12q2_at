# Makefile for the BU Shared Computing Cluster (SCC). Includes both debug and
# release targets with compiler options suitable for development and production,
# respectively. 
#
# There are a few environment variables used here that are specific to the SCC
# module system: namely: SCC_HDF_LIB, SCC_SZIP_LIB, SCC_HDF_INCLUDE,
# SCC_SZIP_INCLUDE. The following module versions are known to be compatible:
#
# hdf/4.2.11
# szip/2.1
#
# #

# user-defined parameters
CC = gcc
CFLAGS_DEBUG = -g -fopenmp -Wall -Wextra #-Werror
CFLAGS_RELEASE = -O3 -fopenmp 
LDFLAGS = -L$(SCC_HDF_LIB) -L$(SCC_SZIP_LIB) -lmfhdf -ldf -lsz -lm 
CPPFLAGS = -I$(SCC_HDF_INCLUDE) -I$(SCC_SZIP_INCLUDE) 
DEFS =
SRC_DIR = ./
BLD_DIR = ../build
EXE_DIR = ../bin
DOC_DIR = ../doc
EXE = modis_smooth_t

# automatically-defined variables
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRC))

# clear implicit rules, protect 'clean' target
.SUFFIXES: 
.SUFFIXES: .o .c
.PHONY: help release debug c5_release c5_debug doc clean 

# help message, default target
help:
	@ echo ""
	@ echo "Usage:"
	@ echo "    make help       : Display this help message"
	@ echo "    make release    : Build optimized executable"
	@ echo "    make debug      : Build non-optimized executable"
	@ echo "    make c5_release : Build optimized executable, reads MODIS C5 data"
	@ echo "    make c5_debug   : Build non-optimized executable, reads MODIS C5 data"
	@ echo "    make clean      : Delete executable, build files, and documentation"
	@ echo ""

# build type targets
debug: CFLAGS = $(CFLAGS_DEBUG)
debug: $(EXE_DIR)/$(EXE) ;

release: CFLAGS = $(CFLAGS_RELEASE)
release: $(EXE_DIR)/$(EXE) ; 

c5_debug: CFLAGS = $(CFLAGS_DEBUG)
c5_debug: DEFS += -DC5
c5_debug: $(EXE_DIR)/$(EXE) ;

c5_release: CFLAGS = $(CFLAGS_RELEASE)
c5_release: DEFS += -DC5
c5_release: $(EXE_DIR)/$(EXE) ; 

# executable and docs
$(EXE_DIR)/$(EXE): $(OBJ) | $(EXE_DIR)  
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@
	#( cat Doxyfile ; \
	#  echo "OUTPUT_DIRECTORY=$(DOC_DIR)" ; \
	#  echo "PREDEFINED=`echo $(DEFS) | sed 's/-D//g'`" ) | doxygen -

# pattern rule for C sources, all files depend on shared header
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/modis_smooth_t.h | $(BLD_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(DEFS) -c $< -o $@

# create directories if needed
$(BLD_DIR): 
	mkdir -p $(BLD_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

# delete build files and executable
clean: 
	-rm -f $(BLD_DIR)/*
	-rm -rf $(DOC_DIR)/*
	-rm -f $(EXE_DIR)/$(EXE) 
