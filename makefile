##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BundleSolver                                                 #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by BundleSolver, i.e., core SMS++.                      #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by MCFSolver.                                                       #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(BNDSLVSDR)  = the directory where the source is                #
#                                                                            #
#   Output: $(BNDSLVOBJ) = the final object(s) / library                     #
#           $(BNDSLVH)   = the .h files to include                           #
#           $(BNDSLVINC) = the -I$( source directory )                       #
#                                                                            #
#                                VERSION 1.00                                #
#                               13 - 05 - 2019                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

LDSSLVOBJ = $(LDSSLVSDR)LagrangianDualSolver.o 

LDSSLVINC = -I$(LDSSLVSDR)

LDSSLVH   = $(LDSSLVSDR)LagrangianDualSolver.h 

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(LDSSLVOBJ) $(LDSSLVSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(LDSSLVSDR)LagrangianDualSolver.o: $(LDSSLVSDR)LagrangianDualSolver.cpp $(LDSSLVH) \
	$(SMS++OBJ) $(BNDSLVOBJ)
	$(CC) -c $*.cpp -o $@ $(LDSSLVINC) $(SMS++INC) $(BNDSLVINC) $(SW)

########################## End of makefile ###################################
