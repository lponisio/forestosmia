#!/usr/bin/env bash

Rscript -e 'install.packages("dplyr", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("brms", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("bayesplot", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("tidybayes", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("mice", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("grid", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("gridExtra", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("scales", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("tidyr", repos="http://cran.r-project.org")'

## plotting
Rscript -e 'install.packages("ggplot2", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("lemon", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("egg", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("RColorBrewer", repos="http://cran.r-project.org")'
