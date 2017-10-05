# bugwas
R package to test for locus and lineage associations in bacterial GWAS

This branch is under development to create a new CRAN package, containing the necessary code for GEMMA.

To install the package in R, download the source code and run
R CMD build bugwas
from the bugwas directory, which will create a file called bugwas_version.tar.gz
Ensure the GSL and LAPACK libraries are correctly installed
Then in R run
install.packages("directory-to-file/bugwas_version.tar.gz")
and
library(bugwas)

