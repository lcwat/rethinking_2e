# Fri Nov 14 14:28:06 2025
# Rethinking
# Luke Watson

# Chapter 4

# load libraries -------------------------------------------------------------

library(tidyverse)
library(rethinking)
library(showtext)
library(tidybayes)

# load data ------------------------------------------------------------------

# quant data of heights and weights of !Kung in the Kalahari desert in southern
# Africa
data("Howell1")

howell <- Howell1

glimpse(howell)

# precis summary from rethinking, includes nice little histograms
precis(howell)


# create plot theme -------------------------------------------------------

font_add_google('Bungee Shade', 'bun')
font_add_google('Noto Sans', 'ns')

showtext_auto()

chap_4_theme <- function() {
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = .75, color = 'black'), 
      axis.ticks = element_line(linewidth = .5, color = 'grey30'), 
      axis.text = element_text(family = 'ns', size = 12), 
      axis.title = element_text(color = 'dodgerblue3', family = 'bun', size = 14)
    )
}

# 4.1 ---------------------------------------------------------------------

# chapter concerns generators or underlying processes that produce distributions
# and the focus of this chapter is the Gaussian distribution. 

# this distribution is thought to be an emergent pattern of the contribution of 
# many fluctuations

# one way to get this is through the sum of independent dichotomous events like 
# arranging 100 people at the 50 yd line and having each flip a coin to determine
# whether they will step forward or backward. have them complete several tosses
# then view the resulting deviations from the starting point approximates Normal
set.seed(4)

position <- crossing(flipper = 1:100, step = 0:16) |> # two vars, create all combos
  mutate(
    deviation = map_dbl(step, ~if_else(. == 0, 0, runif(1, -1, 1))) # flip coin for each of 16 steps
  ) |> 
  group_by(flipper) |> 
  mutate(position = cumsum(deviation)) |> # sum deviations one by one
  ungroup()

# plot out additive deviations
position |> 
  ggplot(aes(x = step, y = position, group = flipper)) +
  
  geom_vline(xintercept = c(4, 8, 12, 16), linetype = 3) +
  geom_line(aes(color = flipper < 2, alpha = flipper < 2)) + # highlight 1st flipper
  
  scale_color_viridis_d(guide = 'none', option = 'magma', end = .8, begin = .2, direction = -1) + 
  scale_alpha_manual(guide = 'none', values = c(.2, 1)) +
  scale_x_continuous('flip number', breaks = 0:4 * 4)

# can see how strings of heads or tails are exceedingly rare, independent flips
# tend to cancel each other out w more flips leading to less deviation from center

# can also generate normality through multiplying small numbers (essentially like 
# adding when small enough)

# observed variables like height are product of multiplicative interactions of 
# many genes
heights <- tibble(
  height = map_dbl(1:10000, ~prod(1 + runif(12, min = 0, max = .1))) # prod of 12 'genes' contrib. small percentage
)

heights |> 
  ggplot(aes(x = height)) +
  geom_density() +
  theme_minimal()

# the smaller the contributions in product, the better the normal approx
tibble(
  big = map_dbl(1:10000, ~ prod(1 + runif(12, 0, .5))), 
  small = map_dbl(1:10000, ~ prod(1 + runif(12, 0, .01)))
) |> 
  pivot_longer(everything(), values_to = 'samples') |> 
  ggplot(aes(x = samples)) +
  geom_density(fill = 'dodgerblue2') +
  facet_wrap(~ name, scales = 'free')

# even larger contributions can be reconfigured to model additive influence through
# log transformations (recall that multiplying within log is equivalent to adding logs)
tibble(
  big = map_dbl(1:10000, ~ prod(1 + runif(12, 0, .5)))
) |> 
  mutate(log_big = log(big)) |> 
  ggplot(aes(x = log_big)) +
  geom_density(fill = 'grey20')

# very common to measure natural phenomena with logarithmic scale (richter, decibel, etc.)

# this justification of gaussian distributions is called the ontological explanation. 
# normally distributed variables add together fluctuations of contributing sources
# of variance. this is an aggregate variable that sheds all info about underlying
# contributing processes (who/what this variable listens to and how much), but like
# epicycle model, simplified models can still be extremely useful

# other explanation is the epistemological justification, where the Gaussian 
# represents a certain state of ignorance. all we know or are willing to say about
# a distribution of measures is the mean and variance, the normal provides the most
# flexibility and fewest assumptions, but is also the least informative. this is
# most consistent with the assumptions of the model or golem, which is completely
# ignorant. if you assume or want to say more about your data, need to consider 
# alternative distributions that are more constrained (e.g., like if your variable
# is categorical or discrete counts of events)

# 4.3 --------------------------------------------------------------------

# filter out children, where age confounds with height to create exponential
howell <- howell |> 
  filter(age >= 18)

# plot the distribution of heights
dens(howell$height)

howell |> 
  ggplot(aes(x = height)) +
  geom_density(fill = 'grey80', color = 'grey80') +
  chap_4_theme()

# gawking at the outcome distribution of your model isn't the best way to choose
# or justify your model. linear models don't need normality to est mean/variance


# .3.2 model and priors -----------------------------------------------------


# assume heights are independently and identically distributed, or that each value
# of h has the same probability function. untrue in most cases, which is where 
# multilevel analyses come into play

# here we want to estimate the mean and variance of a single variable like an 
# intercept only model (single t in traditional stats terms). 

#' hi ~ Normal(mu, sigma) rank rel plausibility of different height distributions
#' mu ~ Normal(178, 20) average for height likely around 178 cm
#' sigma ~ Uniform(0, 50) variance must be positive, and likely not that high

# these priors are determined from our prior knowledge of how these variables 
# tend to be distributed. like mean height being around 5'8", within +/- 1.3 ft

# plot prior
tibble(x = seq(from = 100, to = 250, by = .1)) |> 
  ggplot(aes(x = x, y = dnorm(x, mean = 178, sd = 20))) +
  geom_line() +
  labs(y = 'density') +
  chap_4_theme()

tibble(x = seq(from = -10, to = 60, by = .1)) |> 
  ggplot(aes(x = x, y = dunif(x, min = 0, max = 50))) +
  geom_line() +
  labs(y = 'density') +
  chap_4_theme()

# can simulate from both priors at once to get a prior predictive distribution
# which will help us see if our choices for priors make sense

n <- 1e4
set.seed(4)

# sample from each prior, then use those samples to draw sample for height
sim <- tibble(
  sample_mu = rnorm(n, mean = 178, sd = 20), 
  sample_sigma = runif(n, min = 0, max = 50)
) |> 
  mutate(
    height = rnorm(n, mean = sample_mu, sd = sample_sigma)
  )

# plot
sim |> 
  ggplot(aes(x = height)) +
  geom_density(fill = 'grey80') +
  scale_y_continuous(NULL, breaks = NULL) +
  chap_4_theme()

# most likely to see height around prior mean, constrained within sensible range
# (no negative heights)

# if we described a different prior, can see resulting prior predictive can be 
# non sensical
text <- tibble(
  height = 272 - 25,
  y = 0.0013,
  label = "tallest man",
  angle  = 90
)

tibble(
  sample_mu = rnorm(n, 178, 100), # absurdly large variance for heights
  sample_sigma = runif(n, 0, 50)
) |> 
  mutate(
    height = rnorm(n, sample_mu, sample_sigma)
  ) |> 
  ggplot(aes(x = height)) +
  geom_density(fill = 'dodgerblue', color = NA) +
  geom_vline(xintercept = 0, linetype = 3, color = 'black') + # mark 0
  geom_vline(xintercept = 272, color = 'red3', linetype = 3) + # mark wadlow
  geom_text(
    data = text, aes(y = y, label = label, angle = angle), 
    color = 'red4'
  ) +
  scale_y_continuous(NULL, breaks = NULL) +
  chap_4_theme()
  
# the model before seeing the data now expects a lot of people to have negative
# heights and even more to be taller than robert wadlow


# .3.3 grid approximation of posterior --------------------------------------


# use grid approximation to estimate the posterior distribution of height dist. 
# would never use this for real analysis bc coding is intractable and time sink
n <- 200

# create grid of mean and variance values to test
# crossing() is tidyverse version of expand.grid()
d_grid <- crossing(
  mu = seq(from = 140, to = 160, length.out = n), 
  sigma = seq(from = 4, to = 9, length.out = n)
)

# for each mu, there are 200 values of sigma to consider
glimpse(d_grid)

# grid approx function for likelihood
grid_function <- function(mu, sigma, d = howell) {
  # get post prob from data for this comb of mu/sigma, log then sum bc more 
  # computationally efficient than multiplying this many numbers
  dnorm(d$height, mean = mu, sd = sigma, log = T) |> 
    sum()
}

# grid approx
d_grid <- d_grid |> 
  mutate(
    log_likelihood = map2(mu, sigma, grid_function)
  ) |> 
  unnest(log_likelihood) |> 
  mutate(
    prior_mu = dnorm(mu, 178, 20, log = T), 
    prior_sigma = dunif(sigma, 0, 50, log = T)
  ) |> 
  mutate(
    product = log_likelihood + prior_mu + prior_sigma
  ) |> 
  mutate(
    probability = exp(product - max(product))
  )

# now we have posterior probabilities across all values sampled of mu and sigma
# view posterior w contour plot or raster heatmap
d_grid |> 
  ggplot(aes(x = mu, y = sigma, fill = probability)) +
  geom_raster(interpolate = T) +
  scale_fill_viridis_c(option = 'magma') +
  coord_cartesian(xlim = range(d_grid$mu), ylim = range(d_grid$sigma)) +
  chap_4_theme()

# most of the probability centers around mean of 154 w variance of 7.5

# rather than viewing posterior, usually we draw samples from it 
d_grid_samples <- d_grid |> 
  sample_n(size = 1e4, replace = T, weight = probability)

d_grid_samples |> 
  ggplot(aes(x = mu, y = sigma)) +
  geom_point(size = .9, alpha = 1/15) +
  labs(
    x = expression('mu'[samples]), 
    y = expression('sigma'[samples])
  ) +
  chap_4_theme()

# here we see the shadow of each posterior distribution cast onto the plot
# see the individual or marginal distributions
d_grid_samples |> 
  pivot_longer(mu:sigma) |> 
  ggplot(aes(x = value)) +
  geom_density(fill = 'grey40') +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(NULL) +
  chap_4_theme() +
  facet_wrap(~ name, scales = 'free', labeller = label_parsed)

# compute posterior modes and 95% highest density intervals from samples
d_grid_samples |> 
  pivot_longer(mu:sigma) |> 
  group_by(name) |> 
  mode_hdi(value) # mode hdi

# or post medians and 50% quantile interval
d_grid_samples |> 
  pivot_longer(mu:sigma) |> 
  group_by(name) |> 
  median_qi(value) # mode hdi


# .3.4 approx w brm ---------------------------------------------------------

# while great for learning and teaching, we approximate the posterior w more 
# sensible methods like brm and quap

# detach rethinking to avoid namespace clashes
detach(package:rethinking, unload = T)
library(brms)

b4.1 <- brm(
  data = howell, 
  family = gaussian, 
  height ~ 1, # only est mean and variance
  prior = c(
    prior(normal(178, 20), class = Intercept), 
    prior(uniform(0, 50), class = sigma, ub = 50)
  ), 
  iter = 2000, warmup = 1000, chains = 4, cores = 4, 
  seed = 4, 
  file = 'fits/b_4.01'
)

# or read in
b4.1 <- read_rds('fits/b_4.01.rds')

# view sampling results
plot(b4.1)

# similar marginal posteriors to grid, chains appear to have good convergence

# get the model summary
print(b4.1)

# stan-like summary
b4.1$fit


# .3.5 sampling from brm fit ------------------------------------------------

# variance covariance matrix is glue that holds together quap, like ML methods, 
# info is stored within the matrix to determine how far up that MAP hill we are

# from brm fit, we should extract samples to calculate this
post <- as_draws_df(b4.1)

head(post)

# get vcov matrix
vcov_matrix <- post |> 
  select(b_Intercept:sigma) |> 
  cov()

# variances on diag
diag(vcov_matrix)

# unlike rethinking's extract samples, brm as_draws_df() just pulls iterations
# from the sampling chains 

# can summarize these samples as marginals
post |> 
  pivot_longer(b_Intercept:sigma) |> 
  group_by(name) |> 
  summarize(
    mean = mean(value), 
    sd = sd(value), 
    `2.5%` = quantile(value, probs = .025), 
    `97.5%` = quantile(value, probs = .975)
  ) |> 
  mutate_if(is.numeric, round, digits = 2)  |> # round everything
  mutate(
    # add cute little histograms
    histospark = c(
      rethinking::histospark(post$b_Intercept), 
      rethinking::histospark(post$sigma)
    )
  )
  
# or posterior summary
posterior_summary(b4.1)

# 4.4 linear prediction ---------------------------------------------------

# most models concerned with how variables covary, now we will add a predictor
# variable to the model to see how that changes the specification of priors

# we will covary the height with weight of the kalahari foragers, not particularly
# interesting unless you have a theory about lifespan development
howell <- howell |> 
  filter(age >= 18)

# scatterplot
howell |> 
  ggplot(aes(x = weight, y = height)) +
  geom_point(shape = 1, size = 2) +
  chap_4_theme()

# regression is created when we structure the estimation equation to have the 
# mean of the Gaussian dist. estimated as a linear function of the predictor
# the golem assumes the predictor has a constant (linear) and additive relationship
# to the mean of the outcome

# golem produces posterior, which is relative plausibility of each parameter 
# combination which represent the strengths of the associations (lines)
# rank these lines by plausibility given the data P(mu, sigma, beta, alpha | Data)

#' hi ~ N(mui, sigma)  likelihood, prob of each data point
#' ui = alpha + beta(x - xbar)  beta and alpha param control how much x aff mu
#'  and is now not estimated but constructed from other estimated param and data
#' alpha ~ N(178, 20)
#' beta ~ N(0, 10)
#' sigma ~ Unif(0, 50)

# alpha and beta manipulate mu, allowing for systematic variation across cases in 
# the data (subscript i)

# centering x makes interpretation and defintion of prior easier. alpha param is 
# now the expected height for the average weight mu = alpha 
# beta is the rate of change in height expectation across weight range

# prior predictive check
set.seed(4)

n_lines <- 100

# for each simulated data point, draw from priors and assign to extreme weight
lines <- tibble(
  n = 1:n_lines, 
  a = rnorm(n_lines, 178, 20), 
  b = rnorm(n_lines, 0, 10)
) |> 
  expand_grid(weight = range(howell$weight)) |> 
  mutate(
    height = a + b * (weight - mean(howell$weight))
  )

lines

# plot lines with horizontal height ceiling/floor
lines |> 
  ggplot(aes(x = weight, y = height, group = n)) +
  geom_hline(yintercept = c(0, 272), linetype = 2, linewidth = 1/3) +
  geom_line(alpha = 1/10) +
  coord_cartesian(ylim = c(-100, 400)) +
  chap_4_theme()

# this simulation shows how our prior specification, specifically the beta is 
# creating nonsensical relationships like negative associations between h/w

# so, can try another prior that avoids negative associations like log-normal
# or logarithm of param has normal dist

tibble(b = rlnorm(1e4, 0, 1)) |> 
  ggplot(aes(x = b)) +
  geom_density(fill = 'grey60', color = NA) +
  coord_cartesian(xlim = c(0, 5)) +
  chap_4_theme()

# log normal derivation for mu/sigma we defined
mu = 0
sigma = 1

# mean
exp(mu + (sigma^2) / 2)

# sd
sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))

# find numerically from simulated draws
tibble(d = rlnorm(1e7, 0, 1)) |> 
  summarize(
    mu = mean(d), 
    sigma = sd(d)
  )

# prior predictive will show much better agreement with our prior knowledge
tibble(
  n = 1:n_lines, 
  a = rnorm(n_lines, 178, 20), 
  b = rlnorm(n_lines, 0, 1)
) |> 
  expand_grid(weight = range(howell$weight)) |> 
  mutate(
    height = a + b * (weight - mean(howell$weight))
  ) |> 
  ggplot(aes(x = weight, y = height, group = n)) +
  geom_hline(yintercept = c(0, 272), linetype = 3, linewidth = 1/3) +
  geom_line(alpha = 1/10) +
  coord_cartesian(ylim = c(-100, 400)) +
  chap_4_theme()


# .4.2 --------------------------------------------------------------------

# fit a new model with predictor to find the posterior
# create new centered variable
howell <- howell |> 
  mutate(
    weight_c = weight - mean(weight)
  )

detach(package:rethinking, unload = T)
library(brms)

b4.3 <- brm(
  data = howell, 
  height ~ 1 + weight_c, # centered weight
  prior = c(
    prior(normal(178, 20), class = Intercept), 
    prior(lognormal(0, 1), class = b), 
    prior(uniform(0, 50), class = sigma, ub = 50)
  ), 
  iter = 2000, warmup = 1000, chains = 4, cores = 4, 
  seed = 4, 
  file = 'fits/b_4.03'
)

b4.3 <- read_rds('fits/b_4.03.rds')

plot(b4.3)

# .4.3 --------------------------------------------------------------------

# interpreting the posterior, tables and simulated visuals can describe the 
# posterior, but it is hard to get the complexity of the model and data that 
# generated it. how do parameters act together to influence prediction? 

# what do parameters mean? 
# no consensus on what they mean due to different philosophies toward models, 
# probabilities, and prediction
# rethinking: posterior probabilities of param values describe the relative 
# compatibility of different states of the world with the data, according to the 
# model

# table
posterior_summary(b4.3)[1:3,] |> 
  round(2)

# always consider context of model assumptions
# if you are committed to a line, then lines with a slope around .9 are plausible

# vcov
vcov(b4.3) |> 
  round(3)

# with sigma
as_draws_df(b4.3) |> 
  select(b_Intercept:sigma) |> 
  cov() |> 
  round(3)

# pairs will plot the marginal posteriors and covariance between them (0 for most)
pairs(b4.3)

# almost always better to plot the posterior inference against actual data, 
# helps us interpret the posterior and informally check the model assumptions 
# (e.g., does linear model capture this relationship?)

# can get the line with fixef()
labels <- c(-10, 0, 10) + mean(howell$weight) |> # labels for plot in orig scale
  round(0)

howell |> 
  ggplot(aes(x = weight_c, y = height)) +
  geom_abline(
    intercept = fixef(b4.3)[1], 
    slope = fixef(b4.3)[2]
  ) +
  geom_point(shape = 1, size = 2, color = 'dodgerblue') +
  scale_x_continuous(
    'weight', 
    breaks = c(-10, 0, 10), 
    labels = labels 
  ) +
  chap_4_theme()

# but this is just a line, no uncertainty passed along
# book features example that shows how drawing lines from posterior with less and
# less data becomes more uncertain in intercept and slope

# for each parameter value, their combination produces a mu value for height
# for a particular point of weight (like 50 kg), there is a marginal distribution 
# of mu that we can get intervals for 

# posterior draws
post <- as_draws_df(b4.3)

50 - mean(howell$weight)

mu_at_50kg <- post |> 
  mutate(
    mu_at_50kg = b_Intercept + b_weight_c * 5.01, 
    .keep = 'none' # remove all other cols
  )

mu_at_50kg |> 
  ggplot(aes(x = mu_at_50kg)) +
  stat_halfeye(
    point_interval = mode_hdi, # mode point est
    .width = .95, # 95% hdi
    fill = 'dodgerblue'
  ) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(expression(mu['height | weight = 50'])) +
  chap_4_theme()

# at mu = 50, 95% of the posterior lies within the interval of 158.5 and 159.75

# use fitted to get hdi's to plot
mu <- fitted(b4.3, summary = F)

str(mu) # returns a matrix w as many rows as there was post-warm up draws and 
# as many cols as cases in the data

# use newdata to pass in custom sequence of values for predictor
weight_seq <- tibble(weight = 25:70) |> 
  mutate(weight_c = weight - mean(howell$weight))

mu <- fitted(
  b4.3, summary = F, newdata = weight_seq
) |> 
  data.frame() |> 
  set_names(25:70) |> # label cols by weights sampled
  mutate(iter = row_number())

mu <- mu |> 
  pivot_longer(
    -iter, 
    names_to = 'weight', values_to = 'height'
  ) |> 
  mutate(weight = as.numeric(weight))

# interval as point cloud, or clumps of points that are each gaussians like 
# viewing mu at 50kg from above
howell |> 
  ggplot(aes(x = weight, y = height)) +
  geom_point(data = mu |> filter(iter < 101), alpha = .05, color = 'dodgerblue') +
  coord_cartesian(xlim = c(30, 65)) +
  chap_4_theme()

# or easier, use summary and PI with ribbons
mu_summary <- fitted(b4.3, newdata = weight_seq) |> 
  data.frame() |> 
  bind_cols(weight_seq)

head(mu_summary)

howell |> 
  ggplot(aes(x = weight, y = height)) +
  geom_smooth(
    data = mu_summary, 
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), 
    stat = 'identity', 
    fill = 'grey70', color = 'black', alpha = 1, linewidth = 1/2
  ) +
  geom_point(color = 'dodgerblue', shape = 1, linewidth = 1.5, alpha = 2/3) +
  coord_cartesian(xlim = range(howell$weight)) +
  chap_4_theme()

# given that we think relationship is a line, this is the most plausible line 
# and these are the plausible bounds

# prediction intervals 
# bring in the variability from sigma param, not just mu height samples but actual 
# heights 
# instead of plotting uncertainty in line: mu = alpha + beta * weight
# plot uncertainty in predicted heights: hi ~ Normal(mu, sigma)

# so, for unique values of weight, sample from gaussian with mu and sigma for 
# that weight

# predict from brms
pred_height <- predict(b4.3, newdata = weight_seq) |> 
  data.frame() |> 
  bind_cols(weight_seq)

pred_height |> 
  slice(1:6)

# plot everything together
howell |> 
  ggplot(aes(x = weight)) +
  geom_ribbon(
    data = pred_height, 
    aes(ymin = Q2.5, ymax = Q97.5), 
    fill = 'grey83'
  ) +
  geom_smooth(
    data = mu_summary, 
    aes(y = Estimate, ymin = Q2.5, ymax = Q97.5), 
    stat = 'identity', 
    fill = 'grey70', color = 'black', alpha = 1, linewidth = 1/2
  ) +
  geom_point(
    aes(y = height), color = 'dodgerblue', shape = 1, size = 1.5, alpha = 2/3
  ) +
  coord_cartesian(xlim = range(howell$weight), ylim = range(howell$height)) +
  chap_4_theme()

# rough looking, but with more samples in the post-warm up will smooth

