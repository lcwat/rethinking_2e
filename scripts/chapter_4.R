# Fri Nov 14 14:28:06 2025
# Rethinking
# Luke Watson

# Chapter 4

# load libraries -------------------------------------------------------------

library(tidyverse)
library(rethinking)

# load data ------------------------------------------------------------------

# quant data of heights and weights of !Kung in the Kalahari desert in southern
# Africa
data("Howell1")

howell <- Howell1

glimpse(howell)

# precis summary from rethinking, includes nice little histograms
precis(howell)


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

# filter out children, where age confounds with height
howell <- howell |> 
  filter(age >= 18)
