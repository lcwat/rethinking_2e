# Mon Nov  3 10:44:32 2025
# Rethinking
# Luke Watson

# Chapter 3

# load libraries -------------------------------------------------------------

library(rethinking)
library(tidyverse)
library(tidybayes)
library(brms)

# 1 -----------------------------------------------------------------------

# back to grid approx for the water problem
# how many grid points would you like?
n <- 1000
n_success <- 6
n_trials  <- 9

(
  d <-
    tibble(p_grid = seq(from = 0, to = 1, length.out = n),
           # note we're still using a flat uniform prior
           prior  = 1) %>% 
    mutate(likelihood = dbinom(n_success, size = n_trials, prob = p_grid)) %>% 
    mutate(posterior = (likelihood * prior) / sum(likelihood * prior))
)

# draw samples from this posterior to get a sense of its shape and properties
n_samples <- 1e4

# take samples but weighted by their posterior probability and replace
samples <- d |> 
  slice_sample(n = n_samples, weight_by = posterior, replace = T)

samples |> 
  mutate(
    n_samples = 1:n()
  ) |> 
  ggplot(aes(x = n_samples, y = p_grid)) +
  scale_y_continuous('proportion water (p)', limits = c(0, 1)) +
  geom_point(alpha = .4, color = 'steelblue3') +
  theme_bw()

# plot the density
samples |> 
  ggplot(aes(x = p_grid)) +
  geom_density(fill = 'black') +
  scale_x_continuous('proportion water (p)', limits = c(0,1)) +
  theme_bw()

# will converge on the ideal shape of posterior with more samples
n_samples <- 1e6

samples <- d |> 
  slice_sample(n = n_samples, weight_by = posterior, replace = T)

# replot
samples |> 
  ggplot(aes(x = p_grid)) +
  geom_density(fill = 'black') +
  scale_x_continuous('proportion water (p)', limits = c(0,1)) +
  theme_bw()

# looks better!


# 2 -----------------------------------------------------------------------

# once we have samples, we can do many things with them like define how much
# of posterior probability lies above or below a certain value

# probability below and above .5
samples |> 
  count(p_grid < .5) |> 
  mutate(
    probability = n / sum(n)
  )

# 83% of the probability lies above .5, 17% below

# can also use mean to calculate the probability of this binary data! 
samples |> 
  summarize(
    prob = mean(p_grid < .5)
  )

# intervals of defined mass are like our confidence/credible intervals, they 
# are great summaries of posterior distributions as long as this distribution is
# normal. richard prefers to call them compatibility intervals as gelman and amrhein
# b/c they represent the mass of the posterior where the data and model are 
# compatible, your data and/or model will help determine how credible this is

# can visualize this below
d |> 
  ggplot(aes(x = p_grid, y = posterior)) +
  geom_line() +
  geom_area(data = d |> filter(p_grid > .4 & p_grid < .9)) +
  scale_x_continuous('probability water (p)', limits = c(0, 1)) +
  theme_bw()

# or get the quantiles from the samples
quantile(samples$p_grid)

# 75% percent of the posterior probability lies above .54
quantile(samples$p_grid, probs = c(.05, .95))

# 90% compatibility interval lies between .39 and .84

# but for asymmetrical posterior distributions, you may end up cutting off the
# region with the highest posterior density using a simple percentile interval

# instead, we can use the density with a highest density probability interval (hdpi)
# to avoid this pitfall

# compute a skewed posterior
# here we update the `dbinom()` parameters
n_success <- 3
n_trials  <- 3
n_samples <- 1e4

# update `d`
d <-
  d %>% 
  mutate(likelihood = dbinom(n_success, size = n_trials, prob = p_grid)) %>% 
  mutate(posterior  = (likelihood * prior) / sum(likelihood * prior))

# make the next part reproducible
set.seed(3)

# here's our new samples tibble
(
  samples <-
    d %>% 
    slice_sample(n = n_samples, weight_by = posterior, replace = T)
)

# percentile
(pi <- rethinking::PI(samples$p_grid, prob = .5))
# density
(hdpi <- rethinking::HPDI(samples$p_grid, prob = .5))

# view percentile
d |> 
  ggplot(aes(x = p_grid, y = posterior)) +
  geom_line() +
  geom_area(data = d |> filter(p_grid > pi[1] & p_grid < pi[2])) +
  scale_x_continuous('probability water (p)', limits = c(0, 1)) +
  theme_bw()

# view density
d |> 
  ggplot(aes(x = p_grid, y = posterior)) +
  geom_line() +
  geom_area(data = d |> filter(p_grid > hdpi[1] & p_grid < hdpi[2])) +
  scale_x_continuous('probability water (p)', limits = c(0, 1)) +
  theme_bw()

# also family of functions in tidybayes that can produce same intervals
qi(samples$p_grid, .width = .5)
mode_hdi(samples$p_grid) # this is not working for some reason




# 3 -----------------------------------------------------------------------

# posterior predictive distribution (ppd)
# can propagate forward the uncertainty in parameter estimates described 
# by our posterior distribution to then make predictions about data we could see 
# for the different values of our parameter

# essentially we will create sampling distributions, which are the familiar 
# distributions we use to make inferences in frequentist statistical models 
# (e.g., possible data given a certain parameter value)

# ppd is the average of all these distributions weighted by their posterior prob

n_draws <- 1e5

# simulate 9 draws w rbinom for given probability
simulate_binom <- function(probability) {
  set.seed(3)
  rbinom(n_draws, size = 9, prob = probability) 
}

d_small <- tibble(probability = seq(from = .1, to = .9, by = .1)) |> 
  mutate(draws = purrr::map(probability, simulate_binom)) |> 
  unnest(draws) |> # expand out column of lists
  mutate(label = str_c("p = ", probability))

head(d_small)

# plot the sampling dists
d_small |> 
  ggplot(aes(x = draws)) +
  geom_histogram(binwidth = 1, center = 0,
                 color = "grey92", linewidth = 1/10) +
  scale_x_continuous(NULL, breaks = 0:3 * 3) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(subtitle = "Sampling distributions") +
  coord_cartesian(xlim = c(0, 9)) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ label, ncol = 9) 

# create ppd by weighting sampling dists by post prob
n_samples <- 1e4

# make it reproducible
set.seed(3)

# grid approx
d <- tibble(
  prob_water = seq(from = 0, to = 1, length.out = 1000)
) |> 
  mutate(
    prior = 1, 
    likelihood = dbinom(6, size = 9, prob = prob_water), 
    posterior = likelihood * prior / sum(likelihood * prior)
  )

# sampling, weight by posterior, then simulate 9 draws from those values
# for each probability
samples <- d |>  
  slice_sample(
    n = n_samples, weight_by = posterior, replace = T
  ) |> 
  mutate(
    w = purrr::map_dbl(prob_water, rbinom, n = 1, size = 9)
  )

glimpse(samples)

# plot ppd
samples |> 
  ggplot(aes(x = w)) +
  geom_histogram(binwidth = 1, center = 0,
                 color = "grey92", linewidth = 1/10) +
  scale_x_continuous("number of water samples",
                     breaks = 0:3 * 3) +
  scale_y_continuous(NULL, breaks = NULL) +
  ggtitle("Posterior predictive distribution") +
  coord_cartesian(xlim = c(0, 9),
                  ylim = c(0, 3000)) +
  theme(panel.grid = element_blank())

# density of number waters we could expect to see in 9 tosses
# like reverse engineering the likelihood it seems, bayes theorem works
# in both directions

# richard practice --------------------------------------------------------

p_grid <- seq(from = 0, to = 1, length.out = 1000)

d <- tibble(prob_water = p_grid) |> 
  mutate(
    prior = 1, 
    likelihood = dbinom(6, size = 9, prob = prob_water), 
    raw_posterior = likelihood * prior, 
    std_posterior = raw_posterior / sum(raw_posterior)
  )

# draw samples from posterior
set.seed(3)

# find how much posterior weight lies below certain parameter value
samples <- d |> 
  slice_sample(
    n = 1e4, weight_by = std_posterior, replace = T
  )

samples |>
  summarise(
    prob_e1 = mean(prob_water < .2), 
    prob_e2 = mean(prob_water > .8), 
    prob_e3 = mean(prob_water > .2 & prob_water < .8)
  )

# find what parameter value serves as cutoff of posterior weight
quantile(samples$prob_water, .2) # 20% falls below
quantile(samples$prob_water, .8) # 20% falls above

rethinking::HPDI(samples$prob_water, prob = .66)
rethinking::PI(samples$prob_water, prob = .66)

# construct posterior from grid with new data of 8 water in 15 tosses
d <- tibble(
  prob_water = seq(from = 0, to = 1, length.out = 1000)
) |> 
  mutate(
    prior = 1, 
    likelihood = dbinom(8, size = 15, prob = prob_water), 
    posterior = likelihood * prior / sum(likelihood * prior)
  )

samples <- d |> 
  slice_sample(n = 1e4, weight_by = posterior, replace = T)

rethinking::HPDI(samples$prob_water, prob = .9)

# posterior predictive distribution to simulate new probability of observing
# 6 in 9
ppd <- samples |> 
  mutate(
    num_water = map_dbl(prob_water, rbinom, n = 1, size = 9)
  )

glimpse(ppd) 

ppd |> 
  ggplot(aes(x = num_water)) +
  geom_histogram(bins = 9, fill = 'grey30', color = 'black') +
  scale_x_continuous(breaks = seq(0, 9)) +
  theme_bw()

ppd |> 
  count(num_water == 6) |> 
  mutate(
    prop = n / sum(n)
  )

# 17% chance to see 6 water in 9 tosses given new data 8/15

# birth data
birth1 <- c(1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,0,
            
            
            0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,
            
            
            1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,
            
            
            1,0,1,1,1,0,1,1,1,1)


birth2 <- c(0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,
            
            
            1,1,1,0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
            
            
            1,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,1,
            
            
            0,0,0,1,1,1,0,0,0,0)

# 1 = boy, 0 = girl, birth_ = birth num in family (index)
(d <- tibble(
  family = seq(1, length(birth1)), 
  birth_1 = birth1, 
  birth_2 = birth2
))

# using grid approximation, compute the prob of a birth being a boy
# prob of birth being a boy is the mean of birth 1 and 2
(num_boy = sum(c(d$birth_1, d$birth_2)))

grid <- tibble(
  prob_boy = seq(0, 1, length.out = 1000)
) |> 
  mutate(
    prior = 1, # uniform again
    likelihood = dbinom(
      sum(c(d$birth_1, d$birth_2)), # observed successes (e.g., boys)
      length(d$birth_1)*2, # observed n
      prob = prob_boy
    ), 
    posterior = likelihood * prior / sum(likelihood * prior)
  )

samples <- grid |> 
  slice_sample(
    n = 1e4, weight_by = posterior, replace = T
  )

samples |> 
  ggplot(aes(x = prob_boy)) +
  geom_density(fill = 'black') +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw()
  
# find 50, 89, and 97 hpdi
for(prob in c(.50, .89, .97)) {
  print(rethinking::HPDI(samples$prob_boy, prob = prob))
}

grid |> 
  ggplot(aes(x = prob_boy, y = posterior)) +
  geom_line() +
  
  # fill in area of 97% hpdi
  geom_area(data = grid |> filter(prob_boy <= .627 & prob_boy >= .477)) +
  theme_minimal()

# posterior predictive check
ppd <- samples |> 
  mutate(
    boy = map_dbl(prob_boy, rbinom, n = 1, size = 200)
  )

glimpse(ppd)  

# does model fit data well? is 111 about the central obs of ppd?
ppd |> 
  ggplot(aes(x = boy)) +
  geom_density(fill = 'grey90', color = NA) +
  geom_vline(aes(xintercept = 111), linetype = 3, color = 'black') +
  scale_x_continuous(limits = c(0, 200)) +
  theme_bw()

# appears to fit well!

# brms (solomon) practice ----------------------------------------------------

# fit a model with binomial likelihood for data of 6 waters in 9 tosses
b_3.1 <- brm(
  data = list(w = 6), 
  family = binomial(link = 'identity'), 
  w | trials(9) ~ 0 + Intercept, 
  prior = prior(beta(1, 1), class = b, lb = 0, ub = 1), # uniform prior
  iter = 5000, warmup = 1000, 
  seed = 3, 
  file = 'fits/b_3.01'
)

# get the summary of the posterior for our parameter of interest
posterior_summary(b_3.1)['b_Intercept', ] |> 
  round(digits = 2)

# like how we've been using samples and slice sample to sample the posterior, 
# we can get some draws with fitted()
f <- fitted(
  b_3.1, 
  summary = F, # only simulated draws
  scale = 'linear' # return results in our probability metric
) |> 
  data.frame() |> 
  set_names('p')

# plot posterior draws/samples
f |> 
  ggplot(aes(x = p)) +
  geom_density(fill = 'grey40', color = NA) +
  scale_x_continuous('P(Water)', limits = 0:1) +
  scale_y_continuous(NULL, breaks = NULL) +
  theme(panel.grid = element_blank())

# can get ppd too
set.seed(3)

f <- f |> 
  mutate(
    w = rbinom(n(), size = 9, prob = p)
  )

f |> 
  ggplot(aes(x = w)) +
  geom_histogram(
    binwidth = 1, center = 0, color = 'grey40', 
    linewidth = 1/10
  ) +
  scale_x_continuous('Number water samples', breaks = 0:3*3) +
  scale_y_continuous(NULL, breaks = NULL) +
  theme(panel.grid = element_blank())
