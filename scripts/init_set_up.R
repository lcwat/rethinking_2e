# initial set up

# install rstan from mc-stan
install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# install cmdstanr (more modern version of rstan)
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# install dependencies including rtools before cmdstan
# https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html
install.packages(c('coda', 'mvtnorm', 'devtools', 'dagitty'))

# Then run to complete install
cmdstanr::install_cmdstan() 

# install rethinking package associated with book
devtools::install_github('rmcelreath/rethinking')


