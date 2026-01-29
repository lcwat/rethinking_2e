# Sun Nov  2 23:26:53 2025
# Rethinking
# Luke Watson

# Chapter 2

# load libraries -------------------------------------------------------------

library(tidyverse)
library(rethinking)
library(brms)

# 2.4 make the model go ------------------------------------------------------

# practice computing the posterior likelihood using the grid approximation method
# this method multiplies the prior by the likelihood of the data given a vector
# of parameter values

# this reminds me of model prediction using a grid of data values plugged into the
# model

(
  d <- tibble(
    p_grid = seq(from = 0, to = 1, length.out = 20), # define grid
    prior = 1 # prior is uniform, all values of p equally likely
  ) |> 
    mutate(
      likelihood = dbinom(6, size = 9, prob = p_grid), # compute likelihood of 
      # observing data for each value of parameter
      unstd_posterior = likelihood * prior, # compute posterior
      posterior = unstd_posterior / sum(unstd_posterior) # standardize to sum to 1
    )
)

# plot
d |> 
  ggplot(aes(x = p_grid, y = posterior)) +
  geom_line() +
  geom_point() +
  labs(
    x = 'Parameter value (p)', 
    y = 'Posterior likelihood'
  ) +
  theme(
    axis.line = element_line(linewidth = .4),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

# posterior likelihood increases as the likelihood of data given p increases since
# the prior is uninformed

# but richard encourages us to try a few other priors to see how they change the 
# calculation of the posterior

d <- d |> 
  mutate(
    step_prior = if_else(p_grid < 0.5, 0, 1), 
    exponential_prior = exp(-5*abs(p_grid - 0.5)), # mirrored exponential form
    step_unstd_posterior = likelihood * step_prior, 
    step_posterior = step_unstd_posterior / sum(step_unstd_posterior), 
    exp_unstd_posterior = likelihood * exponential_prior, 
    exponential_posterior = exp_unstd_posterior / sum(exp_unstd_posterior)
  )

# plot all three at once
d |> 
  pivot_longer(
    cols = c(posterior, step_posterior, exponential_posterior), 
    names_to = 'posterior_estimation', values_to = 'estimates'
  ) |>
  mutate(
    posterior_estimation = as.factor(posterior_estimation)
  ) |> 
  ggplot(aes(x = p_grid, y = estimates, color = posterior_estimation)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  facet_wrap(~posterior_estimation)

# how you specify the prior changes the estimate of the posterior, with a
# laplace (mirrored exponentional) leading to a more peaked posterior around
# .5, the step posterior looks similar past where the prior jumps to 1 at .5, 
# also results in a larger portion of the std posterior area devoted to those 
# values 


# quadratic approximation -------------------------------------------------

# grid approx can be very expensive and scales nonlinearly with the number of
# parameters to be estimated
# instead, can use quadratic approximation (similar to gradient descent) to 
# find the peak of the posterior (assuming gaussian) and sample around it
# to estimate the sd/shape

# use quap to estimate
globe_qa <- quap(
  alist(
    W ~ dbinom(W + L, p), # likelihood
    p ~ dunif(0, 1) # uniform prior
  ), 
  data = list(W = 6, L = 3)
)

# get summary
precis(globe_qa)


# mcmc --------------------------------------------------------------------

# the most popular alternative for posterior estimation, does not actually compute
# the posterior but simply samples from it, it returns a count of the parameter
# values which in turn correspond to the posterior probability 

n_samples <- 1000

p <- rep(NA, n_samples)

p[1] <- 0.5

W <- 24

L <- 12

# the algorithm draws parameter values from a normal distribution and then 
# computes the likelihood, choosing the newly drawn sample if the ratio of the
# new likelihood to the old is greater than 1
for(i in 2:n_samples) {
  p_new <- rnorm(1, p[i-1], 0.1)
  
  if(p_new < 0) p_new <- abs(p_new)
  
  if(p_new > 1) p_new <- 2 - p_new
  
  q0 <- dbinom(W, W+L, p[i-1])
  
  
  q1 <- dbinom(W, W+L, p_new)
  
  # if the ratio of new likelihood (q1) to old (q0) is greater than uniform
  # sample, choose that new probability else just assign the previous probability
  p[i] <- ifelse(runif(1) < q1/q0, p_new, p[i-1])
}

# plot histogram
tibble(
  x = p
) |> 
  ggplot(aes(x = x)) +
  geom_density(fill = 'black') +
  geom_vline(aes(xintercept = mean(x)), linetype = 3, color = 'lightblue2') +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal()

# compare to brmsfit
brm_1 <- brm(
  data = list(w = 24), # number successes
  family = binomial(link = 'identity'), # binomial family
  formula = w | trials(36) ~ 0 + Intercept, # one parameter model
  prior = prior(beta(1, 1), class = b, lb = 0, ub = 1), # uniform prior b/n 0 and 1
  seed = 2, 
  file = 'fits/chap_2_brm_mcmc_example.rds'
)

# view results!
print(brm_1)

# intercept here is akin to the mean/mode (in mcmc) of the posterior 
posterior_summary(brm_1) |> 
  round(digits = 2)

# many ways to summarize results with wealth of info ab posterior draws 
as_draws_df(brm_1) |> 
  mutate(n = 'n = 36') |> 
  
  ggplot() +
  geom_density(aes(x = b_Intercept), fill = 'black') +
  geom_vline(aes(xintercept = mean(b_Intercept)), linetype = 3, color = 'grey') +
  scale_x_continuous('proportion water', limits = c(0, 1)) +
  theme_minimal() +
  facet_wrap(~n)

# looks the same as the metropolis estimation with 1000 samples