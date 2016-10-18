branchCorr
=====

Z-estimation via correlation matching for partially observed branching processes and stochastic compartmental models

Note: for hypergeometric sampling from populations with a large number of different types (or colors of balls in the urn sampling description), it is necessary to modify the dependency BiasedUrn using the following steps:
-Download the BiasedUrn package source from CRAN
-Open the file Makevars in folder src
-Modify the line PKG_CPPFLAGS= -DR_BUILD=1 -DMAXCOLORS=32 by replacing 32 with the desired maximum number of ‘colors’ or different types to sample from.
-Install the BiasedUrn package locally from this modified source file after making this change