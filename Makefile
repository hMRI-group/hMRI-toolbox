# Makefile for hMRI-toolboxKT MEX functions
# 
# This Makefile compiles all MEX functions for the hMRI toolbox
# Usage: make (compiles all) or make <target> (compiles specific function)

# MATLAB mex command - change if using different MATLAB version
MEX = mex

# Compiler flags
CFLAGS = -O -largeArrayDims

# Source files
GETHANI_SOURCES = gethani_mex.c gethani.c aws.c
PVAWS_SOURCES = pvaws_mex.c aws.c  
PVAWSLAST_SOURCES = pvawslast_mex.c aws.c

# Targets
all: gethani pvaws pvawslast

gethani: $(GETHANI_SOURCES)
	$(MEX) $(CFLAGS) $(GETHANI_SOURCES) -output gethani

pvaws: $(PVAWS_SOURCES)
	$(MEX) $(CFLAGS) $(PVAWS_SOURCES) -output pvaws

pvawslast: $(PVAWSLAST_SOURCES)
	$(MEX) $(CFLAGS) $(PVAWSLAST_SOURCES) -output pvawslast

clean:
	rm -f *.mexa64 *.mexmaci64 *.mexmaca64 *.mexglx *.mexw64 *.mexw32

# Individual targets for each function
.PHONY: all clean gethani pvaws pvawslast