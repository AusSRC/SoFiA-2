
# ---------------------------------------------------------------
# SWIG Sofia Makefile
#
#        SRCS       = C source files
#        CXXSRCS    = C++ source files
#        OBJCSRCS   = Objective-C source files
#        OBJS       = Additional .o files (compiled previously)
#        INTERFACE  = SWIG interface file
#        TARGET     = Name of target module or executable
#
#----------------------------------------------------------------

SRCDIR        = src
SRCS          = $(SRCDIR)/common.c $(SRCDIR)/statistics_flt.c $(SRCDIR)/statistics_dbl.c $(SRCDIR)/Table.c $(SRCDIR)/String.c $(SRCDIR)/Stack.c $(SRCDIR)/Path.c $(SRCDIR)/Array_dbl.c $(SRCDIR)/Array_siz.c $(SRCDIR)/Map.c $(SRCDIR)/Matrix.c $(SRCDIR)/LinkerPar.c $(SRCDIR)/Parameter.c $(SRCDIR)/Source.c $(SRCDIR)/Catalog.c $(SRCDIR)/Flagger.c $(SRCDIR)/WCS.c $(SRCDIR)/Header.c $(SRCDIR)/DataCube.c
SOFIASRC      = sofia.c
INTERFACE     = sofia.i
PYINTERFACE   = $(INTERFACE:i=py)
WRAPFILE      = $(INTERFACE:.i=_wrap.c)
WRAPOBJ       = $(INTERFACE:.i=_wrap.o)
TARGET        = _sofia.so # Use this kind of target for dynamic loading

exec_prefix   = /usr/bin

CC            = gcc
CFLAGS        = -fPIC --std=c99 -Wno-incompatible-pointer-types -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O0 -g3
INCLUDES      = -I/usr/include/python3.8 -I ~/.local/lib/python3.8/site-packages/numpy/core/include
LDFLAGS       = --std=c99 --pedantic -Wall -Wextra -Wshadow -Wno-unknown-pragmas -Wno-unused-function -Wfatal-errors -O0 -g3 -shared -fopenmp
LIBS          = -lm -lwcs

# SWIG Options
#     SWIG      = location of the SWIG executable
#     SWIGOPT   = SWIG compiler options
#     SWIGCC    = Compiler used to compile the wrapper file

SWIG          = $(exec_prefix)/swig
SWIGOPT       = -python
SWIGCC        = $(CC)

# Rules for creating .o files from source.

SOFIAOBJ      = $(SOFIASRC:.c=.o)
COBJS         = $(SRCS:.c=.o)
ALLOBJS       = $(COBJS) $(SOFIAOBJ)

# Command that will be used to build the final extension.
BUILD         = $(SWIGCC)

# Build options

BUILD_LIBS    = $(LIBS) # Dynamic loading

# Compilation rules for non-SWIG components

.SUFFIXES: .c

.c.o:
	cd $(SRCDIR) && $(CC) $(CFLAGS) $(INCLUDES) -c $(<F)

# ----------------------------------------------------------------------
# Rules for building the extension
# ----------------------------------------------------------------------

all: $(TARGET)

# Convert the wrapper file into an object file

$(WRAPFILE) : $(INTERFACE)
	$(SWIG) $(SWIGOPT) $(INTERFACE)

$(WRAPOBJ) : $(WRAPFILE)
	$(SWIGCC) -c $(CCSHARED) $(CFLAGS) $(WRAPFILE) $(INCLUDES)

$(SOFIAOBJ) : $(SOFIASRC)
	$(SWIGCC) -c $(CCSHARED) $(CFLAGS) $(SOFIASRC) $(INCLUDES)
    

$(TARGET): $(WRAPOBJ) $(ALLOBJS)
	$(BUILD) $(WRAPOBJ) $(ALLOBJS) $(LDFLAGS) $(BUILD_LIBS) -o $(TARGET)

clean:
	rm -f $(COBJS) $(CXXOBJS) $(WRAPOBJ) $(WRAPFILE) $(SOFIAOBJ) $(PYINTERFACE) $(TARGET)
