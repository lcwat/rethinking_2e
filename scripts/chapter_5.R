# Tue Mar  3 10:08:48 2026
# Rethinking
# Luke Watson

# Chapter 5

# load libraries -------------------------------------------------------------

# Load
library(tidyverse)
library(tidybayes)
library(rstan)
library(patchwork)
library(posterior)
library(tigris)
library(ggdag)
library(dagitty)
library(rethinking)

# Drop grid lines
theme_set(
  theme_gray() +
    theme(panel.grid = element_blank())
)

# load data ------------------------------------------------------------------

# waffle house divorce rates
data(WaffleDivorce, package = 'rethinking')

d <- WaffleDivorce

rm(WaffleDivorce)

# standardize
d <- d |> 
  mutate(d = rethinking::standardize(Divorce),
         m = rethinking::standardize(Marriage),
         a = rethinking::standardize(MedianAgeMarriage))

glimpse(d)

# 5.1 --------------------------------------------------------------------

# spurious associations, correlation is common but doesn't imply actual causal
# links between variables 
d |>
  ggplot(aes(x = WaffleHouses/Population, y = Divorce)) +
  stat_smooth(method = "lm", formula = 'y ~ x', fullrange = TRUE, 
              alpha = 1/5, linewidth = 1/2) +
  geom_point(alpha = 1/2, size = 1.5) +
  geom_text(data = d |> 
              filter(Loc %in% c("ME", "OK", "AR", "AL", "GA", "SC", "NJ")),  
            aes(label = Loc), 
            hjust = -0.2, size = 3, vjust = -0.4) +
  scale_x_continuous("Waffle Houses per million", limits = c(0, 55)) +
  ylab("Divorce rate") +
  coord_cartesian(xlim = c(0, 50), ylim = c(5, 15))

cor(d$WaffleHouses/d$Population, d$Divorce)

# view maps
# Get the map data
d_states <- states(cb = TRUE, resolution = "20m") |>
  shift_geometry() |> # move ak and hi to below contig us
  # Add the primary data
  right_join(d |> 
               mutate(NAME = Location |> as.character()) |> 
               select(d:a, NAME),
             by = "NAME") |> 
  # Convert to the long format for faceting
  pivot_longer(cols = c("d", "m", "a"), names_to = "variable")

# chloropleths for each variable
d_states |>
  ggplot() +
  geom_sf(aes(fill = value, geometry = geometry),
          size = 0) +
  scale_fill_gradient(low = "pink", high = "red4", breaks = NULL) +
  theme_void() +
  theme(strip.text = element_text(margin = margin(0, 0, .5, 0))) +
  facet_wrap(~ variable, labeller = label_both) 

# create a simple model between age of first marriage and divorce rate

#' Di ~ N(mui, sigma)
#' mui = alpha + beta_a * Agei // three parameters
#' alpha ~ N(0, 0.2)
#' beta_a ~ N(0, 0.5)
#' sigma ~ Exp(1)

# b/c both outcome and predictor are standardized, makes specification of priors
# much easier b/c centered around 0

# beta_a represents that one sd increase in predictor relates to one sd increase 
# in outcome
sd(d$MedianAgeMarriage)

# one sd change in age at marriage (1.2 yrs) results in one sd change in divorce rate
# which seems HUGE, so we set prior to make 1/-1 unlikely for beta

# fit stan model
model_code_5.1 <- '
data {
  int<lower=1> n;
  vector[n] d;
  vector[n] a;
}
parameters {
  real b0;
  real b2; 
  real<lower=0> sigma; // set lower bound for variance
}
model {
  vector[n] mu;
  mu = b0 + b2 * a;
  d ~ normal(mu, sigma);
  b0 ~ normal(0, .2);
  b2 ~ normal(0, .5);
  sigma ~ exponential(1);
}
generated quantities {
  vector[n] log_lik;
  for(i in 1:n) log_lik[i] = normal_lpdf(d[i] | b0 + b2 * a[i], sigma);
}
'

# must make the stan data for the function
stan_data <- d |> 
  select(d, a) |> 
  compose_data()

# turns data into lists with some additional meta data about standardization
stan_data

# fit
m5.1 <- stan(
  data = stan_data, 
  model_code = model_code_5.1, 
  cores = 4, seed = 5
)

# prior fit as separate model
model_code_5.1_prior <- '
data {
  int<lower=1> n;
  vector[n] d;
  vector[n] a;
}
parameters {
  real b0;
  real b2;
  real<lower=0> sigma;
}
model {
  // The model only contains the prior
  b0 ~ normal(0, 0.2);
  b2 ~ normal(0, 0.5);
  sigma ~ exponential(1);
}
'

m5.1_prior <- stan(
  data = stan_data,
  model_code = model_code_5.1_prior,
  cores = 4, seed = 5)

# sample prior to see plausibility of prior definitions 
as_draws_df(m5.1_prior) |> 
  pivot_longer(b0:sigma) |> 
  mutate(
    name = case_when(
      name == 'b0' ~ 'b[0]', 
      name == "b2" ~ "b[2]",
      name == "sigma" ~ "sigma"
    )
  ) |> 
  ggplot(aes(x = value, y = name)) +
  stat_halfeye(point_interval = mean_qi, .width = 0.89) +
  scale_y_discrete(NULL, labels = ggplot2:::parse_safe, expand = expansion(mult = 0.1)) +
  xlab("prior") +
  coord_cartesian(xlim = c(-1, 5))

# get prior posterior of lines
set.seed(5)

as_draws_df(m5.1_prior) |> 
  slice_sample(n = 50) |> 
  
  ggplot() +
  geom_abline(aes(intercept = b0, slope = b2, group = .draw),
              alpha = 3/4, linewidth = 1/4) +
  scale_x_continuous("Median age marriage (std)", limits = c(-2, 2)) +
  scale_y_continuous("Divorce rate (std)", limits = c(-2, 2))

# plot actual posterior
p2 <- as_draws_df(m5.1) |> 
  expand_grid(a = seq(from = min(d$a), to = max(d$a), length.out = 30)) |> 
  mutate(mu = b0 + b2 * a) |> 
  
  ggplot(aes(x = a)) +
  stat_lineribbon(aes(y = mu),
                  .width = 0.89, color = "blue", fill = alpha("blue", 1/3)) +
  geom_point(data = d,
             aes(y = d),
             size = 2/3) +
  labs(x = "Median age marriage (std)",
       y = "Divorce rate (std)")

p2

# model summary
print(m5.1, pars = c("b0", "b2", "sigma"), probs = c(0.055, 0.945))

# can fit another simple model between marriage rate and divorce rate but this 
# is no good way to compare the effects of these different predictors. they could
# provide independent value or share all information or could eliminate the value
# of another. which is why we should think causally, then construct model to 
# target our effect of interest

# dags --------------------------------------------------------------------

# three indicators or observed variables, AMD, and we can specify their relationships 
# to each other with a DAG to show information flow between variables

# arrow indicates this variable listens to other

# use dagitty and ggdag to construct rough dag of:
#' A directly influences D (D listens to A)
#' M directly influences D (D listens to M)
#' A directly influences M (M listens to A)

set.seed(5)

dagify(M ~ A, D ~ A + M) |> ggdag(node_size = 8)

# pretty it up
dag_coords <- tibble(
  name = c("A", "M", "D"),
  x    = c(1, 3, 2),
  y    = c(2, 2, 1))

p1 <- dagify(
  M ~ A,
  D ~ A + M,
  coords = dag_coords) |>
  
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(alpha = 1/4, color = "firebrick", size = 10) +
  geom_dag_text(color = "firebrick") +
  geom_dag_edges(edge_color = "firebrick") +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0.1, 0.1)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0.2, 0.2)) +
  theme_void()

p1

p2 <- dagify(
  M ~ A,
  D ~ A,
  coords = dag_coords) |>
  
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(alpha = 1/4, color = "firebrick", size = 10) +
  geom_dag_text(color = "firebrick") +
  geom_dag_edges(edge_color = "firebrick") +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0.1, 0.1)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0.2, 0.2)) +
  theme_void()

p2

# view both states of world we are considering
p1 + p2

# age influences divorce rate through a direct and indirect pathways, indirecty through 
# affecting the marriage rate (higher age of marriage, fewer eligible people to 
# be married as a fact of life and death)

# if a variable only has an effect through this indirect path, it is said to be 
# mediated by that middle variable

# these models can imply that variables are independent under certain conditions
# or conditional independencies whether simply unassociated or unassociated after
# conditioning on another set of variables 

# essentially means that after learning value of one variable (Z), learning about
# your variable of interest (X) doesn't reveal any additional information about 
# outcome (Y): Y is independent of X conditional on Z

# dag arrows imply correlations, we have strong ones here so can test that assumption
d |> 
  select(d:a) |> 
  cor()

# no other assumptions of the left dag, but right dag is assumes no relationship
# between d and m (d does not listen to m)

# and b/c m listens to a, the right dag implies conditional independence of d and m
# after conditioning on a

# some code from dagitty to find your conditional independencies
dma_dag_right <- dagitty('dag{D<-A->M}')

impliedConditionalIndependencies(dma_dag_right)

# to test this implication, we need a statistical model that conditions on a to 
# see whether that renders d independent of m

#' 1. after i already know marriage rate, what additional value is there in knowing 
#'    age at marriage?
#' 2. after i already know age at marriage, what additional value is there in also 
#'    knowing marriage rate?

# and the answer lies in the parameter estimates of a model
# some practicioners refer to this as statistical control, controlling for the 
# effect of one variable while estimating effect of other, but language is a little 
# sloppy. should have a better distinction betwen statistical and experimental control


# multivariable model -----------------------------------------------------

#' di ~ N(mui, sigma) // likelihood
#' mui = alpha + bA * ai + bM * mi // add in another var with parameter
#' alpha ~ (0, .2)
#' bA ~ N(0, .5)
#' bM ~ N(0, .5)
#' sigma ~ Exp(1)

# fit model

# standata list
stan_data <- d |> 
  select(d:a) |>
  compose_data()  

names(stan_data)[4] <- 'N'

# run cmdstan model
m5.3 <- cmdstan_model('scripts/stan/model_code_5.3.stan')

m5.3 <- m5.3$sample(
  stan_data, chains = 2, parallel_chains = 2, iter_warmup = 1500, iter_sampling = 2500
)

# much faster! like the two step approach too

# summary
print(m5.3, pars = c("b0", "b1", "b2", "sigma"),  probs = c(0.055, 0.945))

print(m5.3)

# plot coef tabs
bind_rows(
  as_draws_df(m5.1) |> mutate(model = "m5.1"),
  as_draws_df(m5.3) |> mutate(model = "m5.3")
) |> 
  pivot_longer(starts_with("b")) |> 
  filter(name != "b0") |> 
  drop_na(value) |> 
  mutate(name = case_when(
    name == "b1" ~ "beta[1]",
    name == "b2" ~ "beta[2]"
  ) |> factor(levels = c("beta[2]", "beta[1]"))) |> 
  
  ggplot(aes(x = value, y = model)) +
  geom_vline(xintercept = 0, color = "white") +
  stat_pointinterval(.width = 0.89) +
  labs(x = "posterior",
       y = NULL) +
  facet_wrap(~ name, labeller = label_parsed, ncol = 1)
