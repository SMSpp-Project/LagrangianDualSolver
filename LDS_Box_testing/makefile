##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for DCRBlock / MCFSolver)                              #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = LDS_Box_test

# basic directory
DIR = .

# debug switches
#SW = -g -glldb -fno-inline -std=c++17 -ferror-limit=1
# debug switches with address sanitizer and extra pedantic warning
#SW = -g3 -glldb -fno-inline -std=c++17 -ferror-limit=1 -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter
# production switches
SW = -O3 -std=c++17 -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION 
#-DNDEBUG -DCLANG_1200_0_32_27_PATCH -ferror-limit=1 

# compiler
CC = clang++ -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION

# basic CPLEX directory
libCPLEXBSCDIR = /Applications/CPLEX_Studio2211/cplex/

# lib Cplex directory
libCPLEXINCDIR = $(libCPLEXBSCDIR)lib/x86-64_osx/static_pic/

libCPLEXLIB = -L$(libCPLEXINCDIR) -lilocplex -lcplex $(libCPLEXEXTLIB)
libCPLEXINC = -I$(libCPLEXBSCDIR)include/

libOSIINC = -I/usr/local/include/coin
libOSILIB = -rpath /usr/local/lib/ -lOsiClp -lClp -lOsi -lCoinUtils -lOsiCpx

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(DIR)/$(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

SMS++SDR = ../smspp
include $(SMS++SDR)/lib/makefile-inc

LAGRANGIANDR = ./lagrangiandualsolver
BUNDLEDR = ../bundlesolver

LgDSLVSDR = $(LAGRANGIANDR)
LgDHEURSR = $(LAGRANGIANHDR)
LgDSLVODR = $(LAGRANGIANDR)/obj/
#include $(LAGRANGIANDR)/makefile-c

libNDOSDR = ../bundlesolver/NdoFiOracle
include ../bundlesolver/NdoFiOracle/lib/makefile

BNDSLVSDR = $(BUNDLEDR)
BNDSLVODR = $(BUNDLEDR)
include $(BUNDLEDR)/makefile

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(UCBckOBJ) $(SMS++OBJ)

# libraries
MLIB =  $(SMS++LIB) $(UCBckLIB) -L"/Library/gurobi1001/macos_universal2/lib" -lgurobi100

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o  
	$(CC) -o $(DIR)/$(NAME) ../milpsolver/obj/MILPSolver.o ../milpsolver/obj/CPXMILPSolver.o ../milpsolver/obj/GRBMILPSolver.o ../lagrangiandualsolver/LagrangianDualSolver.o ../bundlesolver/BundleSolver.o ../bundlesolver/ParallelBundleSolver.o ../primalproximalheur/PrimalProximalHeur.o $(libNDOOBJ) $^ $(MLIB) $(libNDOINC) $(BNDSLVINC) $(libOSIINC) $(libOSILIB) $(SW) -lpthread -lm 

#$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o  
#	$(CC) -o $(DIR)/$(NAME) ../milpsolver/obj/MILPSolver.o ../milpsolver/obj/CPXMILPSolver.o ../milpsolver/obj/GRBMILPSolver.o ../lagrangiandualsolver/LagrangianDualSolver.o ../bundlesolver/BundleSolver.o ../primalproximalheur/PrimalProximalHeur.o $(libNDOOBJ) $^ $(MLIB) $(libNDOINC) $(BNDSLVINC) $(libOSIINC) $(libOSILIB) $(SW) -lpthread -lm 


# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(SMS++INC) $(UCBckINC)

# includes
MH =    $(SMS++H) $(UCBckH)

# compile command

$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) -I../lagrangiandualsolver/ -I../bundlesolver/ -I../milpsolver/include/ -I../bundlesolver/ -I$(libNDOSDR)/OPTUtils -I../bundlesolver/NdoFiOracle/NDOSlver -I../primalproximalheur/ -I$(libNDOINC) -I$(BNDSLVINC) $(libOSIINC) $(libOSILIB) -I"/Library/gurobi1001/macos_universal2/include/" $(SW)

#$(DIR)/test.o: $(DIR)/test.cpp $(MH)
#	$(CC) -c $*.cpp -o $@ $(MINC) -I../lagrangiandualsolver/ -I../bundlesolver/ -I../milpsolver/include/ -I../bundlesolver/ -I$(libNDOSDR)/OPTUtils -I../bundlesolver/NdoFiOracle/NDOSlver -I../primalproximalheur/ -I$(libNDOINC) -I$(BNDSLVINC) $(libOSIINC) $(libOSILIB) -I"/Library/gurobi1001/macos_universal2/include/" $(SW)

############################ End of makefile #################################
