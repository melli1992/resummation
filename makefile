# Paths to libraries
LHAPDF := /home/mbeekveld/LHAPDF
LOOPTOOLS := /home/mbeekveld/LoopTools-2.15/x86_64-Linux
CUBA := /home/mbeekveld/Cuba

# Default goal of make file, this is built if you run just `make`
.DEFAULT_GOAL := all

# -- Define linker settings --
LDFLAGS := -L$(LHAPDF)/lib -L$(LOOPTOOLS)/lib64 -L$(CUBA)/lib

# -- Define compiler settings --
# Setup warnings and optimizations
CXXFLAGS := -Wall -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter -O3 -std=c++11
# Include folders for libraries
CXXFLAGS := $(CXXFLAGS) -I$(LHAPDF)/include -I$(LOOPTOOLS)/include
CXXFLAGS := $(CXXFLAGS) -I$(CUBA)/include
# Include folders for project
CXXFLAGS := $(CXXFLAGS) -I./include

# -- Object file lists for the various target binaries -- 
# Splitting this off allows easier writing of utility tasks
#  and better dealing with dependency files (see common rules)

#general objects
#OBJECTSs
RESUM_OBJS := $(RESUM_OBJS) src/parameters.o src/mellin_functions.o
#MC integrations
RESUM_OBJS := $(RESUM_OBJS) src/MC/cuba_integration.o
#PDFs
RESUM_OBJS := $(RESUM_OBJS) src/pdf/fit_coefficients.o src/pdf/deriv_pdf.o src/pdf/pf_pdf.o src/pdf/mellin_pdf.o src/pdf/cheb_pdf.o 
#resummation functions
RESUM_OBJS := $(RESUM_OBJS) src/resum/resum_functions.o src/resum/SCET_functions.o 
#kfactor files
RESUM_OBJS := $(RESUM_OBJS) src/kfactors/k_factors_higgs.o src/kfactors/k_factors_nnlo_higgs.o  
RESUM_OBJS := $(RESUM_OBJS) src/kfactors/k_factors_dy.o src/kfactors/k_factors_nnlo_dy.o
RESUM_OBJS := $(RESUM_OBJS) src/kfactors/k_factors_dihiggs.o src/kfactors/k_factors_diboson.o
RESUM_OBJS := $(RESUM_OBJS) src/kfactors/k_factors_prompt_photon.o 
#SUSY
RESUM_OBJS := $(RESUM_OBJS) src/susy/susypar.o
#ulilities
RESUM_OBJS := $(RESUM_OBJS) src/utilities/inout.o src/utilities/polygamma.o

#ttH resummation
ttH_program := $(RESUM_OBJS) 
ttH_program := $(ttH_program) src/MC/monte_carlo.o src/MC/tth_vegas.o
ttH_program := $(ttH_program) src/resum/tth_softanom.o src/resum/resum_tth.o src/kfactors/k_factors_ttH.o
ttH_program := $(ttH_program) programs/ttbarh_run.o 

#resummation for higgs and DY
DYh_program := $(RESUM_OBJS) programs/DYh.o
#resummation for diboson
diboson_program := $(RESUM_OBJS) programs/diboson.o

#trySCET
SCETtry_program := $(RESUM_OBJS) programs/try_scetfunc.o
#lumnicheck
lumni_program := $(RESUM_OBJS) programs/lumni_check.o

#PDFfits
PDF_program := $(RESUM_OBJS) programs/make_mellin_pdf.o

# -- Rules for building the main binaries --
# If not all libraries are needed for a specific binary, just remove the
# -l flag for that library here.

ttH: $(ttH_program)
	g++ -o ttH $(ttH_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -looptools -lgfortran -lboost_program_options
	    
	    
DYh: $(DYh_program)
	g++ -o DYh $(DYh_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -looptools -lgfortran -lboost_program_options

diboson: $(diboson_program)
	g++ -o diboson $(diboson_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -looptools -lgfortran -lboost_program_options
	    

SCETtry: $(SCETtry_program)
	g++ -o SCETtry $(SCETtry_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -looptools -lgfortran -lboost_program_options -lboost_math_c99
	    
	    
lumni: $(lumni_program)
	g++ -o lumni $(lumni_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -looptools -lgfortran -lboost_program_options
	    	    
PDF: $(PDF_program)
	g++ -o PDF $(PDF_program) $(CXXFLAGS) $(LDFLAGS) -lgsl -lgslcblas -lm  \
	    -lcuba -lLHAPDF -looptools -lgfortran -lboost_program_options
# -- Utility rules --

# all: builds all binaries
all: ttH DYh

# clean: remove all generated files
clean:
	rm -f resum $(RESUM_OBJS) $(RESUM_OBJS:.o=.d)
	
# Mark the utility rules as phony to instruct make that they don't produce
#  files as output. This stops trouble when you accidentally create files
#  with those names.
.PHONY: all clean

# -- Rules for producing object files from source --
# Produce object files from c++ (the -MD ensures that this also creates
#  a dependency file for the source being compiled)
# The suffixes statement tells make that this is a general rule, and that
# .cpp and .o are file extensions.
.SUFFIXES: .cpp .o
.cpp.o:
	g++ -c $< -o $@ $(CXXFLAGS)

# Include dependency files here to make changes to header files
#  cause proper recompilation of source files using them
#-include $(TRYLOOP_OBJS:.o=.d)
#-include $(RESUM_OBJS:.o=.d)
#-include $(FIXED_OBJS:.o=.d)
#-include $(MELLIN_OBJS:.o=.d)
