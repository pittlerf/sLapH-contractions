@mainpage

This is the new version of the contraction code. 

A Wick contraction in this context means the analytic result obtained when
"integrating out" fermionic degrees of freedom from the path integral in 
lattice QCD. 

In the sLapH framework space is approximated by the span of the eigenvectors
of a Laplace operator, below referred to as eigenvector space. Introducing
random vectors and dilution in time, eigenvector and Dirac space yields a 
stochastical estimate for the propagator. 
The correlation function can be expressed as a trace of products of 
propagators, eigenvectors and field operators (@f$\Gamma@f$). 

@todo Would a simple example like a @f$\pi^+ @f$ help? 

The stochastic estimates can be brought into a form which allows saving them
to hard disk: The perambulator. This allows to reuse inverted propagators to
"build" multiple correlation functions. This code calculates correlation 
functions in exactly this manner. It featurus

- IO functions for perambulators, random vectors, eigenvectors
- An Infile to specify exactly the physical process in question
  - Different quark flavors
  - Flexible Dilution schemes
  - Every linearly independent Dirac structure
  - Different momenta
- 13 different diagrams for one- and two-meson operators pre-implemented

As a starting point the main function is in @ref contract.cpp

The issue of reusing certain parts for multiple diagrams achieved by calculating 
and saving quarklines. These hold particular combinations of perambulators and 
operators that appear in multiple diagrams.

The baseline of the program flow is
1. Read infile
2. Determine the data needed to calculate the specified diagrams: The class 
    GlobalData processes them to vectors of certain structs for Quarklines, 
    Correlators and Operators.

    Internally all physical loops are replaced with auto-loops over these 
    structs. This allows to decouple the computation from the physical context 
    and keep the code very clean where the indices are very complicated.
3. Read all necessary data. The perambulators and the random vectors are saved 
    as attributes of LapH::Perambulator and LapH::RandomVector. 
4. Build parts that are reusable. The operators are constructed 
    in LapH::OperatorsForMesons. This is the only step where the eigenvectors 
    are needed, thus they are read on the fly and discarded afterwards. The 
    quarklines are constructed in LapH::Quarklines.
    @note The operators referred to above live solely in time and eigenvector 
          space. In order to save memory if multiple Dirac structures are to be 
          calculated at once it has been factored out and is performed when 
          building a quarkline or a correlator. They are in particular not 
          identical to the quantum field operators
5. Construct the correlators. This happens in LapH::Correlators. The class has
    a private member function corresponding to and building the diagrams 
    admissible in the infiles.
