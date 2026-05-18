##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of LagrangianDualSolver                                         #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by LagrangianDualSolver, i.e., core SMS++ and the       #
#   MILPSolver module (the latter is used only by PrimalProximalHeur to      #
#   warm-start the heuristic by solving the LP relaxation of f_Block).       #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by LagrangianDualSolver.                                            #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(MILPSINC)    = the -I$( MILPSolver include directory )         #
#           $(LgDSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(LgDSLVOBJ)   = the final object(s) / library                   #
#           $(LgDSLVH)     = the .h files to include                         #
#           $(LgDSLVINC)   = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

LgDSLVOBJ = $(LgDSLVSDR)/obj/LagrangianDualSolver.o \
	$(LgDSLVSDR)/obj/PrimalProximalHeur.o

LgDSLVINC = -I$(LgDSLVSDR)/include

LgDSLVH   = $(LgDSLVSDR)/include/LagrangianDualSolver.h \
	$(LgDSLVSDR)/include/PrimalProximalHeur.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(LgDSLVOBJ) $(LgDSLVSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(LgDSLVSDR)/obj/LagrangianDualSolver.o: \
	$(LgDSLVSDR)/src/LagrangianDualSolver.cpp \
	$(LgDSLVSDR)/include/LagrangianDualSolver.h $(SMS++OBJ)
	$(CC) -c $(LgDSLVSDR)/src/LagrangianDualSolver.cpp -o $@ \
	$(LgDSLVINC) $(SMS++INC) $(SW)

$(LgDSLVSDR)/obj/PrimalProximalHeur.o: \
	$(LgDSLVSDR)/src/PrimalProximalHeur.cpp $(LgDSLVH) $(SMS++OBJ)
	$(CC) -c $(LgDSLVSDR)/src/PrimalProximalHeur.cpp -o $@ \
	$(LgDSLVINC) $(SMS++INC) $(MILPSINC) $(SW)

########################## End of makefile ###################################
