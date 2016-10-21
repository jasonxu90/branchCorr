branchCorr
=====

Z-estimation via correlation matching for partially observed branching processes and stochastic compartmental models

## Description
This package includes functions to forward-simulate from a general class of stochastic compartmental models using the Gillespie algorithm [1], as well as capabilities to produce discretely observed datasets and partially observed datasets given a sampling distribution. The software is developed with a focus on compartmental models of hematopoiesis [2,3]. In addition to simulation, routines for parameter inference are implemented based on the method of moments or z-estimation [4]. These functions use numerical optimization to match observed correlations in the data to their corresponding model-based expressions, which depend on rate parameters. Inferring these parameters therefore amounts to a nonlinear least squares problem, seeking the parameters which minimize the residuals between empirical and analytical correlations.  

## Installation
The package can be installed directly from github using the `devtools` package, which can easily be installed using the command `install.packages("devtools")`.
See https://github.com/hadley/devtools for more details.

To install `branchCorr`, run the following:
```r
library(devtools)
install_github("jasonxu90/branchCorr")
```

## Vignette
By default, the vignette is not compiled during the installation. We provide a vignette, ``branchCorr: simulation and inference for partially observed stochastic compartmental models'', that walks through simulation and inference functions necessary to recreate all simulation studies in [2]. Examples are smaller-scale, but still require quite some computing time. Users are encouraged to run all functions in the vignette locally by reducing the number of simulations, initial restarts, etc when necessary. 

The source `.Rnw` file, as well as a `.pdf` of the vignette compiled using `knitr`, are included in the repository. To automatically compile the vignette upon package installation, instead install using the line
`install_github("jasonxu90/branchCorr", build_vignettes = TRUE)`. Note that this option will be significantly slower than installing the package alone.


## References
1.  Gillespie DT (1977) "Exact stochastic simulation of coupled chemical reactions," *Journal of Physical Chemistry,* 81(25):2340-2361.

2. To be posted 

3. Catlin, SN, Abkowitz, JL and Guttorp, P (2001) "Statistical Inference in a Two-Compartment Model for Hematopoiesis,", *Biometrics,* 57(2):546-553.

4. Pakes A and Pollard D (1989) "Simulation and the asymptotics of optimization estimators," *Econometrica: Journal of the Econometric Society,*, 57(5):1027-1057. 
