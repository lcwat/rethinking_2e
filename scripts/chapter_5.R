# Tue Mar  3 10:08:48 2026
# Rethinking
# Luke Watson

# Chapter 5

# load libraries -------------------------------------------------------------

# Load
library(tidyverse)
library(tidybayes)
# library(rstan)
library(patchwork)
library(posterior)
library(tigris)
library(ggdag)
library(dagitty)
library(rethinking)
library(brms)
library(ggrepel)

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

# standardize, also adds some metadata to those vars
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
# model_code_5.1 <- '
# data {
#   int<lower=1> n;
#   vector[n] d;
#   vector[n] a;
# }
# parameters {
#   real b0;
#   real b2; 
#   real<lower=0> sigma; // set lower bound for variance
# }
# model {
#   vector[n] mu;
#   mu = b0 + b2 * a;
#   d ~ normal(mu, sigma);
#   b0 ~ normal(0, .2);
#   b2 ~ normal(0, .5);
#   sigma ~ exponential(1);
# }
# generated quantities {
#   vector[n] log_lik;
#   for(i in 1:n) log_lik[i] = normal_lpdf(d[i] | b0 + b2 * a[i], sigma);
# }
# '

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
# model_code_5.1_prior <- '
# data {
#   int<lower=1> n;
#   vector[n] d;
#   vector[n] a;
# }
# parameters {
#   real b0;
#   real b2;
#   real<lower=0> sigma;
# }
# model {
#   // The model only contains the prior
#   b0 ~ normal(0, 0.2);
#   b2 ~ normal(0, 0.5);
#   sigma ~ exponential(1);
# }
# '

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

# rename this var N to match stan code
names(stan_data)[4] <- 'N'

set.seed(5)

# run cmdstan model
m5.3 <- cmdstan_model('scripts/stan/model_code_5.3.stan') # creates executable file, 1st compilation takes sec, all following are quick!

# sample posterior
m5.3 <- m5.3$sample(
  stan_data, chains = 2, parallel_chains = 2, iter_warmup = 1500, 
  iter_sampling = 2500
)

# much faster! like the two step approach too

# summary
print(m5.3)
precis(m5.3)

# stat halfeyes for params
as_draws_df(m5.3) |> 
  select(b0:sigma) |> 
  pivot_longer(
    everything(), names_to = 'pars', values_to = 'draw'
  ) |>
  
  ggplot(aes(x = draw, y = pars)) +
  
  geom_vline(xintercept = 0, color = 'black') +
  
  stat_halfeye(fill = 'firebrick') +
  
  scale_y_discrete(
    labels = c(
      'b0' = expression(b[0]), 
      'b1' = expression(b[1]), 
      'b2' = expression(b[2]), 
      'sigma' = expression(sigma)
    )
  ) +
  
  labs(x = 'posterior', y = NULL)
  

# plot coef tab
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
  ) |> factor(levels = c("beta[2]", "beta[1]"))
  ) |> 
  
  ggplot(aes(x = value, y = model)) +
  geom_vline(xintercept = 0, color = "white") +
  stat_pointinterval(.width = 0.89) +
  labs(x = "posterior",
       y = NULL) +
  facet_wrap(~ name, labeller = label_parsed, ncol = 1)

# assuming no confounding variables, parameter b1 implies that there is no 
# direct causal path between m and d, it is spurious even though it is highly 
# correlated

# in multivariate models, have more options to interpret the posterior based on 
# combinations of parameters

# one option is to fit a model of each predictor onto one another and extract the 
# residuals, then compare that to outcome to see if any remaining info from variable 
# can relate to the outcome 

# but there are easier ways to begin to understand the implications of your parameters
# in your model either through pure visualization of the fit to the data or 
# through causal counterfactual type manipulations testing alternative data

# can ask model for predictions not observed in data, more for causal analysis to 
# probe model for answers in particular cases 

# plotting model predictions ----------------------------------------------

# run another model regressing m on a
b5.4 <- brm(
  data = d, 
  family = gaussian,
  m ~ 1 + a,
  prior = c(prior(normal(0, 0.2), class = Intercept),
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma)),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = "fits/b_5.04"
)

summary(b5.4)

# predictor residual plots plot residuals of regressed model between predictors
# to isolate effect of that regressed variable on outcome

# get residuals with brms::residuals()
r <- residuals(b5.4) |>
  # To use this in ggplot2, we need to make it a tibble or data frame
  data.frame() |> 
  bind_cols(d)

p3 <- r |> 
  ggplot(aes(x = Estimate, y = d)) +
  stat_smooth(method = "lm", fullrange = T,
              color = "firebrick4", fill = "firebrick4", 
              alpha = 1/5, linewidth = 1/2) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2, alpha = 2/3, color = "firebrick4") +
  geom_text_repel(data = r |> filter(Loc %in% c("WY", "ND", "ME", "HI", "DC")),  
                  aes(label = Loc), 
                  size = 3, seed = 5) +
  scale_x_continuous(limits = c(-2, 2)) +
  coord_cartesian(xlim = range(r$Estimate)) +
  labs(x = "Marriage rate residuals",
       y = "Divorce rate (std)") +
  theme_bw() +
  theme(panel.grid = element_blank())

p3
# plot shows how remaining information after conditioning on age results in no
# relationship to outcome of interest

# residuals are parameters, they are estimated by the model and therefore have 
# uncertainty. may have seen residualized variables in models before, this is a 
# mistake b/c that uncertainty is lost, this is more accurate description of the 
# uncertainty in the residual predictor values
r |>
  ggplot(aes(x = Estimate, y = d)) +
  stat_smooth(method = "lm", fullrange = T,
              color = "firebrick4", fill = "firebrick4", 
              alpha = 1/5, linewidth = 1/2) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  # The only change is here
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5),
                  alpha = 2/3, color = "firebrick4") +
  geom_text_repel(data = r |> filter(Loc %in% c("ID", "HI", "DC")),  
                  aes(label = Loc), 
                  size = 3, seed = 5) +
  scale_x_continuous(limits = c(-2, 3)) +
  coord_cartesian(xlim = range(r$Estimate),
                  ylim = range(d$d)) +
  labs(x = "Age at marriage residuals",
       y = "Divorce rate (std)") +
  theme_bw() +
  theme(panel.grid = element_blank())

# posterior prediction plot of predictions against 

# brm fit of m5.3, can't figure out how to do the fitted for pure cmdstan obj rn
b5.3 <- brm(
  data = d, 
  family = gaussian,
  d ~ 1 + m + a,
  prior = c(prior(normal(0, 0.2), class = Intercept),
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma)),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = "fits/b_5.03"
)

# 
f <- fitted(b5.3) |>
  data.frame() |>
  # Un-standardize the model predictions
  mutate(across(.cols = -Est.Error, .fns = \(x) x * sd(d$Divorce) + mean(d$Divorce))) |> 
  bind_cols(d) 

f |>
  ggplot(aes(x = Divorce, y = Estimate)) +
  geom_abline(linetype = 2, color = "grey50", linewidth = 0.5) +
  geom_point(size = 1.5, color = "firebrick4", alpha = 3/4) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),
                 linewidth = 1/4, color = "firebrick4") +
  geom_text(data = f |> filter(Loc %in% c("ID", "UT", "RI", "ME")),
            aes(label = Loc), 
            hjust = 1, nudge_x = - 0.25) +
  labs(x = "Observed divorce", y = "Predicted divorce") +
  theme_bw() +
  theme(panel.grid = element_blank())

# get model predictions from cmdstan model
d_pred <- d |> 
  cbind(
    m5.3_quantities |> filter(str_detect(variable, 'prediction')) |> 
      select(mean, q5, q95) |> 
      rename(pred = mean)
  )

d_pred |>
  ggplot(aes(x = Divorce, y = pred)) +
  geom_abline(linetype = 2, color = "grey50", linewidth = 0.5) +
  geom_point(size = 1.5, color = "firebrick4", alpha = 3/4) +
  geom_linerange(aes(ymin = q5, ymax = q95),
                 linewidth = 1/4, color = "firebrick4") +
  labs(x = "Observed divorce", y = "Predicted divorce") +
  theme_bw() +
  theme(panel.grid = element_blank())


# 2. masked relationships -------------------------------------------------

data("milk", package = 'rethinking')

d <- milk

rm(milk)

glimpse(d)

# depending on physiological and developmental aspects of spp, likely to see differences
# in milk production
# what extent energy content of milk, measured here by kilocalories, is related 
# to the percent of the brain mass that is neocortex

# inspect primary vars defined by mcelreath
d |> 
  select(kcal.per.g, mass, neocortex.perc) |> 
  pairs(col = 'firebrick4')

# standardize and log transform mass
d <- d |> 
  mutate(
    kcal.per.g_s     = (kcal.per.g - mean(kcal.per.g)) / sd(kcal.per.g), 
    log_mass_s       = (log(mass) - mean(log(mass))) / sd(log(mass)), 
    neocortex.perc_s = (neocortex.perc - mean(neocortex.perc, na.rm = T)) / sd(neocortex.perc, na.rm = T)
  )

# start with simple univariate regression
# expression(kcal/g[i] ~ Normal(mu[i], sigma))
# mu[i] = alpha + beta[cortex_percent]*cortex_percent[i]

b5.5_draft <- brm(
  data = d, 
  family = gaussian,
  kcal.per.g_s ~ 1 + neocortex.perc_s,
  prior = c(
    prior(normal(0, 1), class = Intercept), # slightly more permissive priors, but still informed
    prior(normal(0, 1), class = b),
    prior(exponential(1), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  sample_prior = T, # get samples from prior
  file = "fits/b_5.05_draft"
)

# got some rows that resulted in na's, model didn't return valid probability even 
# for the starting param values, some spp are missing values for neocortex

# use listwise deletion and feel bad about it, but learn better method in chap 15

# complete case data
dcc <- d |> 
  drop_na(ends_with('_s')) # remove cases with any nas in our std vars

nrow(d) - nrow(dcc)

# update to refit model with new data
b5.5_draft_cc <- update(
  b5.5_draft,
  newdata = dcc,
  seed = 5,
  file = "fits/b_5.05_draft_cc"
)

# avoids compilation just changing data

# prior predictive check!
set.seed(5)

prior_draws(b5.5_draft_cc) |> # draw from priors, req sample_prior = T in model
  slice_sample(n = 50) |> # grab random 50 preds
  rownames_to_column() |> 
  expand_grid(neocortex.perc_s = c(-2, 2)) |> # 2 values all thats needed for line
  mutate(kcal.per.g_s = Intercept + b * neocortex.perc_s) |> # calc pred
  
  ggplot(aes(x = neocortex.perc_s, y = kcal.per.g_s)) +
  geom_line(
    aes(group = rowname), alpha = 0.4, color = "firebrick"
  ) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(
    x = "neocortex percent (std)",
    y = "kilocal per g (std)",
    subtitle = "Intercept ~ dnorm(0, 1)\nb ~ dnorm(0, 1)"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) 

# prior looks bad! refit with more constrained priors
summary(b5.5_draft_cc)

detach(package:rethinking, unload = T)

# fit 
b5.5 <- brm(
  data = dcc, 
  family = gaussian,
  kcal.per.g_s ~ 1 + neocortex.perc_s,
  prior = c(prior(normal(0, 0.2), class = Intercept), # tighter priors
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma)),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  sample_prior = T,
  file = "fits/b_5.05"
)

# now check out the if the priors are more sensible
set.seed(5)
prior_draws(b5.5) |> 
  slice_sample(n = 50) |> 
  rownames_to_column() |> 
  expand_grid(neocortex.perc_s = c(-2, 2)) |> 
  mutate(kcal.per.g_s = Intercept + b * neocortex.perc_s) |> 
  
  ggplot(aes(x = neocortex.perc_s, y = kcal.per.g_s, group = rowname)) +
  geom_line(alpha = 0.4, color = "firebrick") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "neocortex percent (std)",
       y = "kilocal per g (std)",
       subtitle = "Intercept ~ dnorm(0, 0.2)\nb ~ dnorm(0, 0.5)") +
  theme_bw() +
  theme(panel.grid = element_blank())

summary(b5.5) # not much different, but errors are tighter

# create coef tab comparison for each models estimates
# Wrangle
bind_rows(as_draws_df(b5.5_draft_cc), as_draws_df(b5.5)) |> # combine draws
  select(b_Intercept:sigma) |> # select pars
  mutate(fit = rep(c("b5.5_draft", "b5.5"), each = n() / 2)) |> # create fit col
  pivot_longer(-fit) |> # expand 
  group_by(name, fit) |> 
  
  # get preds and pci
  summarise(
    mean = mean(value),
    ll = quantile(value, prob = 0.025),
    ul = quantile(value, prob = 0.975)
  ) |> 
  mutate(fit = factor(fit, levels = c("b5.5_draft", "b5.5"))) |> 
  
  # Plot
  ggplot(aes(x = mean, y = fit, xmin = ll, xmax = ul)) +
  geom_pointrange(color = "firebrick") +
  geom_hline(yintercept = 0, alpha = 1/5, color = "firebrick") +
  labs(x = "posterior", 
       y = NULL) +
  facet_wrap(~ name, ncol = 1) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank())

# plot preds over raw data
nd <- tibble(neocortex.perc_s = seq(from = -2.5, to = 2, length.out = 30))

fitted(
  b5.5, 
  newdata = nd,
  probs = c(0.025, 0.975, 0.25, 0.75) # 50% and 95% ci
) |>
  data.frame() |>
  bind_cols(nd) |> 
  
  ggplot(aes(x = neocortex.perc_s, y = Estimate)) +
  geom_ribbon(
    aes(ymin = Q2.5, ymax = Q97.5),
    alpha = 1 / 5,
    fill = "firebrick"
  ) +
  geom_smooth(
    aes(ymin = Q25, ymax = Q75),
    stat = "identity",
    alpha = 1 / 5,
    color = "firebrick4",
    fill = "firebrick4",
    linewidth = 1 / 2
  ) +
  geom_point(
    data = dcc,
    aes(x = neocortex.perc_s, y = kcal.per.g_s),
    color = "firebrick4",
    size = 2
  ) +
  coord_cartesian(
    xlim = range(dcc$neocortex.perc_s),
    ylim = range(dcc$kcal.per.g_s)
  ) +
  labs(x = "neocortex percent (std)", y = "kilocal per g (std)") +
  theme_bw() +
  theme(panel.grid = element_blank())

# now look at mass 
b5.6 <- brm(
  data = dcc, 
  family = gaussian,
  kcal.per.g_s ~ 1 + log_mass_s,
  prior = c(prior(normal(0, 0.2), class = Intercept),
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma)),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  sample_prior = T,
  file = "fits/b_5.06")

# plot pred over data
nd <- tibble(log_mass_s = seq(from = -2.5, to = 2.5, length.out = 30))

fitted(
  b5.6,
  newdata = nd,
  probs = c(0.025, 0.975, 0.25, 0.75)
) |>
  data.frame() |>
  bind_cols(nd) |>
  
  ggplot(aes(x = log_mass_s, y = Estimate)) +
  geom_ribbon(
    aes(ymin = Q2.5, ymax = Q97.5),
    alpha = 1 / 5,
    fill = "firebrick"
  ) +
  geom_smooth(
    aes(ymin = Q25, ymax = Q75),
    stat = "identity",
    alpha = 1 / 5,
    color = "firebrick4",
    fill = "firebrick4",
    linewidth = 1 / 2
  ) +
  geom_point(
    data = dcc,
    aes(y = kcal.per.g_s),
    color = "firebrick4",
    size = 2
  ) +
  coord_cartesian(
    xlim = range(dcc$log_mass_s),
    ylim = range(dcc$kcal.per.g_s)
  ) +
  labs(x = "log body mass (std)", y = "kilocal per g (std)") +
  theme_bw() +
  theme(panel.grid = element_blank())

# multivariable model
# include both additively into model
b5.7 <- brm(
  data = dcc, 
  family = gaussian,
  kcal.per.g_s ~ 1 + neocortex.perc_s + log_mass_s,
  prior = c(prior(normal(0, 0.2), class = Intercept),
            prior(normal(0, 0.5), class = b),
            prior(exponential(1), class = sigma)),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = "fits/b_5.07"
)

# compare model params from each fit to see revealed relationships of both
# neocortex and body mass with milk nutritional content
bind_cols(
  as_draws_df(b5.5) |> 
    transmute(`b5.5_beta[N]` = b_neocortex.perc_s),
  as_draws_df(b5.6) |> 
    transmute(`b5.6_beta[M]` = b_log_mass_s),
  as_draws_df(b5.7) |> 
    transmute(`b5.7_beta[N]` = b_neocortex.perc_s,
              `b5.7_beta[M]` = b_log_mass_s)
) |> 
  pivot_longer(everything()) |> 
  group_by(name) |> 
  summarise(mean = mean(value),
            ll   = quantile(value, prob = 0.025),
            ul   = quantile(value, prob = 0.975)) |> 
  separate(name, into = c("fit", "parameter"), sep = "_") |> 
  
  ggplot(aes(x = mean, y = fit, xmin = ll, xmax = ul)) +
  geom_pointrange(color = "firebrick") +
  geom_vline(xintercept = 0, alpha = 1/5, color = "firebrick") +
  ylab(NULL) +
  facet_wrap(~ parameter, ncol = 1, labeller = label_parsed) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "transparent"))

# look at the correlations again to start to unpack why these associations only appeared 
# when modelled together
library(GGally)

dcc |> 
  select(ends_with('_s')) |> 
  ggpairs()

# high correlation between predictors body mass and neocortex percentage but one
# is negatively associated with milk richness (mass) and other positively (neo)

# create some dags to see what is going on
gg_dag <- function(d) {
  d |> 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(alpha = 1/4, color = "firebrick", size = 10) +
    geom_dag_text(color = "firebrick") +
    geom_dag_edges(edge_color = "firebrick") +
    scale_x_continuous(NULL, breaks = NULL, expand = c(0.1, 0.1)) +
    scale_y_continuous(NULL, breaks = NULL, expand = c(0.2, 0.2)) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

# Left DAG
dag_coords <- tibble(
  name = c("M", "N", "K"),
  x    = c(1, 3, 2),
  y    = c(2, 2, 1))

p1 <- dagify(
  N ~ M,
  K ~ M + N,
  coords = dag_coords) |>
  gg_dag()

# Middle DAG
p2 <- dagify(
  M ~ N,
  K ~ M + N,
  coords = dag_coords) |>
  gg_dag()

# Right DAG
dag_coords <- tibble(
  name = c("M", "N", "K", "U"),
  x    = c(1, 3, 2, 2),
  y    = c(2, 2, 1, 2))
p3 <- dagify(
  M ~ U,
  N ~ U,
  K ~ M + N,
  coords = dag_coords) |>
  gg_dag() +
  geom_point(x = 2, y = 2,
             color = "firebrick4", shape = 1, size = 10, stroke = 1.25)

# combine
p1 + p2 + p3

# unfortunately for us, each of these dags imply same set of conditional independencies, 
# making it impossible to tell from the data alone which one is better inference
# this is known as markov equivalence in this lit

# helpful to understand this masking with some simulation
cases <- 100

set.seed(5)

# simulate a world where dag 1 is true
p1

d_sim <- tibble(m = rnorm(cases, mean = 0, sd = 1)) |> 
  mutate(
    n = rnorm(cases, mean = m, sd = 1), # listens to mass
    k = rnorm(cases, mean = n - m, sd = 1) # listens to n, subtracting eff of m
  )

# see what we just did
d_sim |> 
  GGally::ggpairs()

# load fit and use plug in simulated data to model
b5.7 <- read_rds('fits/b_5.07.rds')

fixef(b5.7) |> round(2)

b5.7_sim <- update(
  b5.7, 
  newdata = d_sim, # plug in new data
  formula = k ~ 1 + n + m, 
  seed = 5, 
  file = 'fits/b_5.07_sim'
)

b5.5_sim <- update(
  b5.7, 
  newdata = d_sim, # plug in new data
  formula = k ~ 1 + n, 
  seed = 5, 
  file = 'fits/b_5.05_sim'
)

b5.6_sim <- update(
  b5.7, 
  newdata = d_sim, # plug in new data
  formula = k ~ 1 + m, 
  seed = 5, 
  file = 'fits/b_5.06_sim'
)

# compare coef for each
fixef(b5.5_sim) |> round(2)
fixef(b5.6_sim) |> round(2)
fixef(b5.7_sim) |> round(2)


# 3 cats ------------------------------------------------------------------

# discrete and unordered, how does outcome change with presence or absence of 
# category? more confusing is model machinery behind how estimates are created
# with these variables

# go back to kalahari forager dataset
data('Howell1', package = 'rethinking')
d <- Howell1
rm(Howell1)

glimpse(d)

# 1 when male, 0 when female, it is an indicator variable or dummy coded by default
# turns parameter on and off depending on row value

# model of height based on sex
#' h[i] ~ N(mu[i], sigma)
#' mu[i] = alpha + beta[i] * sex[i]
#' alpha ~ N(178, 20)
#' beta[m] ~ N(0, 10) // only influences height estimate when sex[i] = 1
#' sigma ~ Uniform(0, 50)

# alpha changes meaning, now is the estimated mean height for female instead of 
# overall mean height of sample

# another complication for bayesian model is that males have more uncertainty
# than female mean estimate due to being additive of two parameters rather than one

# better approach is the index variable, contains just arbitrary integers that 
# correspond to different levels

d <- d |> 
  mutate(sex = if_else(male == 1, 2, 1))

head(d)

# update model with index var
#' h[i] ~ N(mu[i], sigma)
#' mu[i] = alpha[sex[i]] // vary alpha based on a[j]
#' alpha[j] ~ N(178, 20)  for j = 1..2 // same prior for a[1] and a[2]
#' sigma ~ Unif(0, 50)

# to do this in brms, change sex to a factor
d <- d |> 
  mutate(sex = factor(sex))

b5.8 <- brm(
  data = d, 
  family = gaussian,
  height ~ 0 + sex,
  prior = c(
    prior(normal(178, 20), class = b),
    prior(exponential(1), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = "fits/b_5.08"
)


# stan --------------------------------------------------------------------

# for index vars in stan model, a little different b/c stan doesn't allow ints
# within their vectors, instead have to use array 
# Index variable approach
# model_code_5.8 <- '
# data {
#   int<lower=1> n;
#   int<lower=1> n_sex; // var for number of levels
#   vector[n] height;
#   array[n] int sex; // array for sex ints
# }
# parameters {
#   vector[n_sex] a;
#   real<lower=0, upper=50> sigma;
# }
# model {
#   height ~ normal(a[sex], sigma); // subset mu by sex
#   
#   a ~ normal(178, 20);
#   
#   // As an alternative, this syntax works in this example, 
#   // and allows the prior to differ by level
#   // a[1] ~ normal(178, 20);
#   // a[2] ~ normal(178, 20);
#   
#   // This snytax does not work
#   // a[sex] ~ normal(178, 20);
#   
#   sigma ~ uniform(0, 50);
# }
# '
# 
# # Dummy variable approach
# model_code_5.8b <- '
# data {
#   int<lower=1> n;
#   vector[n] height;
#   vector[n] male; // just vector of continuous int values
# }
# parameters {
#   real b0;
#   real b1;
#   real<lower=0, upper=50> sigma;
# }
# model {
#   height ~ normal(b0 + b1 * male, sigma);
#   
#   b0 ~ normal(178, 20);
#   b1 ~ normal(0, 10);
#   sigma ~ uniform(0, 50);
# }
# '


# -------------------------------------------------------------------------

# view summary
print(b5.8)

# some more summary techniques harkening back to ch 2!
get_variables(b5.8)

gather_draws(b5.8, b_sex1, b_sex2, sigma) |> 
  median_qi()

# plot
b5.8 |> 
  gather_draws(b_sex1, b_sex2, sigma) |> 
  median_hdci() |> 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval()

# more importantly, we are looking for a difference between these levels, which 
# means that we should compute a contrast to pass along the uncertainty in the 
# differences between the sexes in the outcome
as_draws_df(b5.8) |> 
  mutate(diff_fm = b_sex1 - b_sex2) |> 
  pivot_longer(cols = c(b_sex1:sigma, diff_fm)) |> 
  group_by(name) |> 
  mean_qi(value, .width = 0.89)

# halfeye plot
as_draws_df(b5.8) |> 
  mutate(diff = b_sex1 - b_sex2) |> # female to male contrast, should be neg
  
  ggplot(aes(x = diff)) +
  stat_halfeye(fill = 'goldenrod') + 
  geom_vline(xintercept = 0, alpha = 2/5, color = 'orchid') +
  labs(y = 'density', x = 'contrast (female -> male)') +
  theme_bw() +
  theme(panel.grid = element_blank())

# real advantage of index approach is when there are more than 2 cats, as indicator
# approach creates those dummy variables for k-1 levels like lm() default also
# multilevel models demand them!

# lets look at clade in milk data
data('milk', package = 'rethinking')
d <- milk
rm(milk)

levels(d$clade); d |> distinct(clade) # base vs tidyverse
# 4 levels of primate clades

# for brms, can keep this factor in tact and will work just fine as index

# but still need to standardize outcome
d <- d |> 
  mutate(kcal_per_g_s = (kcal.per.g - mean(kcal.per.g)) / sd(kcal.per.g))

glimpse(d)

#' kcal[i] ~ N(mu[i], sigma)
#' mu[i] = alpha[clade[i]] // mean for each clade level
#' a[j] ~ N(0, 0.5)   for j = 1,..,4 // four levels/clades
#' sigma ~ Exp(1)

# brms model
b5.9 <- brm(
  data = d, 
  family = gaussian,
  kcal_per_g_s ~ 0 + clade, # use 0 and factor coding will be compiled into index
  prior = c(
    prior(normal(0, 0.5), class = b),
    prior(exponential(1), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = "fits/b_5.09"
)

# model in stan
# model_code_5.9 <- '
# data {
#   int<lower=1> n;
#   int<lower=1> n_clade_id;
#   vector[n] k; // kcal
#   array[n] int clade_id; // cat ints must be in array
# }
# parameters {
#   vector[n_clade_id] a; // vector of clade indexed pars
#   real<lower=0> sigma;
# }
# model {
#   k ~ normal(a[clade_id], sigma); 
#   
#   a ~ normal(0, 0.5); // sets all priors at once, much easier
#   sigma ~ exponential(1);
# }
# '

# view summary
print(b5.9)

# plot coefs with convenience functions
mcmc_plot(b5.9, variable = '^b_', regex = T)

# or ground up ggplot and tidybayes approach
b5.9 |> 
  as_draws_df() |> 
  select(starts_with('b')) |> 
  set_names(distinct(d, clade) |> arrange(clade) |> pull()) |> # set colnames
  pivot_longer(everything()) |> 
  
  # now plot
  ggplot(aes(x = value, y = reorder(name, value))) +
  geom_vline(xintercept = 0, color = 'mistyrose4', alpha = 2/5) +
  stat_pointinterval(
    point_interval = mode_hdi, .width = .89, color = 'peachpuff3', size = 1
  ) +
  labs(x = 'expected kcal (std)', y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_text(face = 'italic', hjust = 0)
  )
