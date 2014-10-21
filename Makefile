# ---------------------------------------------------------------------------
#
#	Makefile for wstree
#
#	Devon Powell 
#	(with some useful bits from Jonathan Zrake)
#
#	Do not modify this file!
#	All user-set options live in Makefile.in
#
#	usage: make
#
# ---------------------------------------------------------------------------

# if there is no Makefile.in then use the template
ifneq ($(strip $(MAKEFILE_IN)),)
# use value of MAKEFILE_IN if provided on the command line
else ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = Makefile.in
else
MAKEFILE_IN = Makefile.in.template
endif
include $(MAKEFILE_IN)

# Source files
SOURCES = wstree.cpp 
COMMON = HDF_IO.hh
OBJ = $(SOURCES:.cpp=.o)
EXE = wstree

# Set up HDF5 dependencies
INC += -I$(HDF5_HOME)/include
LIB += -L$(HDF5_HOME)/lib
LDFLAGS += -lhdf5

all: $(EXE)

$(EXE): $(COMMON) $(OBJ) 
	$(CC) $(LIB) $(OBJ) -o $@ $(LDFLAGS) $(CFLAGS)

.cpp.o: $(COMMON)
	$(CC) $(INC) $(DEF) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) *~ core $(EXE)
