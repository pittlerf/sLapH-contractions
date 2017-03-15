@mainpage

This is the new version of the contraction code. 

As a starting point the main function is in @ref contract.cpp
The desired Wick diagrams and physical quantum numbers must be specified in an 
infile along with the paths to sLapH perambulators and randomvectors and if
nonzero momentum is wanted eigenvectors or vdaggerv.

In the class GlobalData the information are processed to vectors of 
self-defined structs for Quarklines, Correlators and Operators. Internally all
physical loops are replaced with auto-loops over these structs.
