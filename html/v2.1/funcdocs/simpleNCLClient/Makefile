# A simple (but heavily commented) Makefile for program "simpleNCLClient" that
#	depends on NCL. The source code for this client application should be stored
#	as "simpleNCLClient.cpp" in the same directory as this Makefile.

# This Makefile can be used after you have built and installed NCL.
# NCL_INSTALL_DIR should be the directory that you specified as the
#	--prefix argument to NCL's configure script.  By default this
#	is /usr/local in *NIX systems.
# If NCL_INSTALL_DIR is not in your shell's environment, then you'll have
#	to define it. With something like the next line.
# NCL_INSTALL_DIR = /usr/local/lib

# C preprocessor flags to include the parent of NCL's intsalled header directory.
CPP_FLAGS = -I$(NCL_INSTALL_DIR)/include

# linker flags to include the directory that hold libncl.
LD_FLAGS = -L$(NCL_INSTALL_DIR)/lib/ncl

# linker to list the ncl library as a library to include
LIBS=-lncl

# we'll use CXX is to refer to the C++ compiler
CXX=g++

# list all of the object files that are linked together to make our executable
objects = simpleNCLClient.o

# We could list headers here as DEPS so the compiler will rebuild our executable
# 	whenever the headers are changed.
DEPS = simpleNCLClient.h


# A rule for how to build simpleNCLClient.
# After the colon on the first line is a list of things that simpleNCLClient
# 	depends on.
# The second line (and any other indented lines) tell make how to build
# 	simpleNCLClient.  This is where we use our LD_FLAGS and LIBS
simpleNCLClient: $(objects)
	$(CXX) -o $@ $^ $(LD_FLAGS) $(LIBS)


# A rule for how to build any .o file.
# %.o in the target name is a makefile wildcard for any file ending in .o
# The %.cpp means that XYZ.o depends on XYZ.cpp
# 	whatever is in DEPS is also noted as a dependency
# The first step in the rule is a shell if statement that will generate an
# 	error if you do not have NCL_INSTALL_DIR defined (by uncommenting
# 	line 11 above or in your shell's env).
# The syntax $@ in the rule indicates "the target of this rule" (what the
# 	is trying to create, which is "simpleNCLClient.o" in this case)
# The $< syntax means the first dependency (the .cpp file).
# Note that this is where we use our CPP_FLAGS
%.o: %.cpp $(DEPS)
	if test -z $(NCL_INSTALL_DIR) ; then echo "NCL_INSTALL_DIR must be set. See comments in the Makefile." && exit 1; fi
	$(CXX) -c -o $@  $<	$(CPP_FLAGS)

# This target without dependencies can only be invoked by specifying it as the
# 	target ("make clean" from the command line).
# rm removes files from the filesystem
# The -f arg to rm tells rm to not complain even if the files do not exist.
# *.o is shell wildcard for any file that ends in .o
clean:
	rm -f *.o simpleNCLClient
