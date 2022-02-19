Home-Field Advantage in the Summer Olympics:
================
John Williams
2/19/2022

## Introduction

As with most sports, the Olympics appear to be susceptible to home-field
advantage, the tendency for athletes competing on their own turf to
outperform their opponents. Olympic athletes from the host country end
up on the podium more often and take home more gold medals than they did
when they were competing away from home. In this vignette, we explore
the home-country advantage in the Summer Olympics by the rate of medals
per participant in the host year compared to the previous Olympics four
years earlier.

Let *i* = 1, ..., 18 denote the 18 Summer Olympic games for which we
have data, from the 1952 games in Finland to the 2020 games in Japan.
Let *Y*<sub>*i*1</sub> be the number of medals won by the host country
during Summer Olympic game *i* and *Y*<sub>*i*0</sub> be the number of
medals won by the host country in the previous Olympics. Similarly, let
*N*<sub>*i*0</sub> and *N*<sub>*i*1</sub> be the number of participates
from the country in the corresponding Olympics. For example, Italy
hosted the Olympics in 1960, thus *i* = 3. They had
*N*<sub>3, 1</sub> = 280 participants in 1960 and
*N*<sub>3, 0</sub> = 135 participants in 1956; they won
*Y*<sub>3, 1</sub> = 36 medals in 1960 and *Y*<sub>3, 0</sub> = 25
medals in 1956.

## Aggregate Analysis

First, we want to study *λ*<sub>1</sub>, the expected number of medals
per participant for a host country when competing in their host year. To
do this, we need to derive its posterior distribution. We aggregate the
data across the 18 Olympics; thus, $Y\_1 = \\sum\_{i=1}^{18}Y\_{i1}$ and
$N\_1 = \\sum\_{i=1}^{18}N\_{i1}$. Since *Y*<sub>1</sub> is a count with
mean *N*<sub>1</sub>*λ*<sub>1</sub>, a natural model for the likelihood
is the Poisson distribution:

*Y*<sub>1</sub>\|*λ*<sub>1</sub> ∼ *P**o**i**s**s**o**n*(*N*<sub>1</sub>*λ*<sub>1</sub>)

The conjugate prior for the Poisson distribution is the
*G**a**m**m**a*(*a*, *b*) distribution. We can choose an “uninformative”
Gamma distribution by selecting *a* and *b* to be sufficiently small.
For this analysis, we will choose *a* = 1 and *b* = 1. Thus our
“uninformative” prior is

*λ*<sub>1</sub> ∼ *G**a**m**m**a*(1, 1)

Since we have specified a Poisson distribution as the likelihood and and
a Gamma distribution for the prior, we know (from the derivation in
class) that the posterior distribution will also be a Gamma
distribution. Specifically, given our data is *Y*<sub>1</sub> = 1, 016
and *N*<sub>1</sub> = 7, 979, the posterior distribution is

*λ*<sub>1</sub>\|*Y*<sub>1</sub> ∼ *G**a**m**m**a*(1 + *Y*<sub>1</sub>, 1 + *N*<sub>1</sub>) = *G**a**m**m**a*(1017, 7980)

Similarly, we can find the posterior distribution of *λ*<sub>0</sub>,
the expected number of medals per participant for the host country when
competing in the previous Olympics. Since *Y*<sub>0</sub> = 682 and
*N*<sub>0</sub> = 4, 715, we have

*λ*<sub>0</sub>\|*Y*<sub>0</sub> ∼ *G**a**m**m**a*(1 + *Y*<sub>0</sub>, 1 + *N*<sub>0</sub>) = *G**a**m**m**a*(683, 4716)

We graph the resulting posterior distributions for
*λ*<sub>1</sub>\|*Y*<sub>1</sub> and *λ*<sub>0</sub>\|*Y*<sub>0</sub>
below:

<img src="C:/Users/johnw/Repos/olympics-bayesian/Home-Field-Advantage-in-the-Summer-Olympics_files/figure-gfm/aggregate-analysis-1.png" width="85%" style="display: block; margin: auto;" />

The main assumption for this Bayesian analysis is specifying the
likelihood. We chose to model *Y*<sub>1</sub>\|*λ*<sub>1</sub> and
*Y*<sub>0</sub>\|*λ*<sub>0</sub> using Poisson distributions since
*Y*<sub>1</sub> and *Y*<sub>0</sub> are both counts of the number of
medals for a particular Olympics. In these models, the parameters
*λ*<sub>1</sub> and *λ*<sub>0</sub> are the expected number of medals
per participant from the home country in the respective Olympics. The
support of *λ*<sub>1</sub> and *λ*<sub>0</sub> is all non-negative real
numbers, reinforcing our choice to model the likelihood as a Poisson
distribution.

## Hypothesis Test

Now that we have the posterior distributions for
*λ*<sub>1</sub>\|*Y*<sub>1</sub> and *λ*<sub>0</sub>\|*Y*<sub>0</sub>,
we can test the hypothesis that there is a home-country advantage in the
summer Olympics. Formally, we are testing that the expected number of
medals per participant for the host country when competing in the host
Olympics is greater than the expected number of medals per participant
for the host country when competing in the previous Olympics,
i.e. *λ*<sub>1</sub> &gt; *λ*<sub>0</sub>.

To test this hypothesis, we use Monte Carlo sampling to approximate the
probability *λ*<sub>1</sub> &gt; *λ*<sub>0</sub>. We collect 100,000
random samples of *λ*<sub>1</sub> from its posterior distribution,
*G**a**m**m**a*(1017, 7980); and we collect 100,000 random samples of
*λ*<sub>0</sub> from its posterior distribution,
*G**a**m**m**a*(683, 4716). We pair this data and calculate the
proportion of these sample pairs where
*λ*<sub>1</sub> &gt; *λ*<sub>0</sub>. The result is the approximate
probability *λ*<sub>1</sub> &gt; *λ*<sub>0</sub>. Upon completing this
process, we get *P*(*λ*<sub>1</sub> &gt; *λ*<sub>0</sub>) = 0.00491.
Thus, there is shockingly little evidence (only a 0.5% chance!) that
there is a home-country advantage in the summer Olympics.

These results are not sensitive to the prior. *Y*<sub>1</sub>,
*Y*<sub>0</sub>, *N*<sub>1</sub>, and *N*<sub>0</sub> are sufficiently
large that, for all relatively small values of *a* and *b*, a
*G**a**m**m**a*(*a*, *b*) prior will be uninformative. A
*G**a**m**m**a*(0.1, 0.1), *G**a**m**m**a*(1, 1), and
*G**a**m**m**a*(1, 10) prior all produce the same credible intervals for
*λ*<sub>1</sub> and *λ*<sub>0</sub> to at least 3 significant figures.

## Prediction

The next Olympics will be held in France in 2024. We can predict the
number of medals France will win in the 2024 Olympics, *Y*<sup>⋆</sup>.
To do this, we will need to estimate the number of participants that
France will send to the 2024 Olympics, *N*<sup>⋆</sup>. On average, host
countries send twice as many participants to the Olympics when they host
versus the previous Olympics. We can reasonably estimate that France
will send *N*<sup>⋆</sup> = 796 participants to the 2024 Olympics as
they sent 398 athletes to the previous summer games.

Now we can use Monte Carlo sampling to approximate the posterior
predictive distribution. We collect 100,000 random samples of
*λ*<sub>1</sub> from the posterior distribution,
*λ*<sub>1</sub>\|*Y*<sub>1</sub> ∼ *G**a**m**m**a*(1017, 7980). We use
those values of *λ*<sub>1</sub> to collect 1 random sample of
*Y*<sup>⋆</sup> from each of the 100,000 likelihoods of the form
*Y*<sup>⋆</sup>\|*λ*<sub>1</sub> ∼ *P**o**i**s**s**o**n*(*N*<sup>⋆</sup>*λ*<sub>1</sub>) = *P**o**i**s**s**o**n*(796*λ*<sub>1</sub>).
These 100,000 random samples of *Y*<sup>⋆</sup> approximate the
posterior predictive distribution, which we use to calculate the point
estimate and 95% credible interval for *Y*<sup>⋆</sup>:

<img src="C:/Users/johnw/Repos/olympics-bayesian/Home-Field-Advantage-in-the-Summer-Olympics_files/figure-gfm/prediction-1.png" width="85%" style="display: block; margin: auto;" />

If France sends 796 participants to the 2024 Olympics, we predict they
will earn 101 medals. Of course, this is just a point estimate. We can
add more certainty to our prediction by providing a range for the number
of medals France will win: there is a 95% probability that France will
win between 81 and 123 medals.

## Country-Specific Analysis

In this section we conduct a Bayesian analysis separately by country
(combining the data across the two Olympics for Australia, Japan, and
the United States). The likelihood selected for *Y*<sub>*i*1</sub> and
*Y*<sub>*i*0</sub>–the expected number of medals per participant during
the host Olympics and previous Olympics, respectively, for each
country–has the form
*Y*<sub>*i**j*</sub>\|*λ*<sub>*j*</sub> ∼ *P**o**i**s**s**o**n*(*N*<sub>*i**j*</sub>*λ*<sub>*j*</sub>)
where *i* = 1, ..., 15 and *j* = 0, 1. We choose the prior to be
*λ*<sub>*j*</sub> ∼ *G**a**m**m**a*(0.1, 0.1). Having described the
likelihood and the prior, we use the data to calculate the posterior
distributions for *λ*<sub>1</sub> and *λ*<sub>0</sub> for each county,
summarized in Table 1.

| Country       | Posterior *λ*<sub>*i*1</sub>\|*Y*<sub>*i*1</sub> | Posterior *λ*<sub>*i*0</sub>\|*Y*<sub>*i*0</sub> | *P*(*λ*<sub>*i*1</sub>/*λ*<sub>*i*0</sub> &gt; 1) |
|---------------|--------------------------------------------------|--------------------------------------------------|---------------------------------------------------|
| Finland       |  ∼ *G**a**m**m**a*(22.1, 258.1)                  |  ∼ *G**a**m**m**a*(24.1, 129.1)                  | 0.00453                                           |
| Australia     |  ∼ *G**a**m**m**a*(92.1, 911.1)                  |  ∼ *G**a**m**m**a*(52.1, 498.1)                  | 0.4294                                            |
| Italy         |  ∼ *G**a**m**m**a*(36.1, 280.1)                  |  ∼ *G**a**m**m**a*(25.1, 135.1)                  | 0.08583                                           |
| Japan         |  ∼ *G**a**m**m**a*(80.1, 949.1)                  |  ∼ *G**a**m**m**a*(59.1, 557.1)                  | 0.09417                                           |
| Mexico        |  ∼ *G**a**m**m**a*(9.1, 275.1)                   |  ∼ *G**a**m**m**a*(1.1, 94.1)                    | 0.91923                                           |
| West Germany  |  ∼ *G**a**m**m**a*(40.1, 423.1)                  |  ∼ *G**a**m**m**a*(26.1, 275.1)                  | 0.50774                                           |
| Canada        |  ∼ *G**a**m**m**a*(11.1, 385.1)                  |  ∼ *G**a**m**m**a*(5.1, 208.1)                   | 0.64136                                           |
| Soviet Union  |  ∼ *G**a**m**m**a*(195.1, 489.1)                 |  ∼ *G**a**m**m**a*(125.1, 410.1)                 | 0.99105                                           |
| United States |  ∼ *G**a**m**m**a*(275.1, 1169.1)                |  ∼ *G**a**m**m**a*(202.1, 941.1)                 | 0.84129                                           |
| South Korea   |  ∼ *G**a**m**m**a*(33.1, 401.1)                  |  ∼ *G**a**m**m**a*(19.1, 175.1)                  | 0.17696                                           |
| Spain         |  ∼ *G**a**m**m**a*(22.1, 422.1)                  |  ∼ *G**a**m**m**a*(4.1, 229.1)                   | 0.98997                                           |
| Greece        |  ∼ *G**a**m**m**a*(16.1, 426.1)                  |  ∼ *G**a**m**m**a*(13.1, 140.1)                  | 0.00954                                           |
| China         |  ∼ *G**a**m**m**a*(100.1, 599.1)                 |  ∼ *G**a**m**m**a*(63.1, 384.1)                  | 0.54966                                           |
| Great Britain |  ∼ *G**a**m**m**a*(65.1, 530.1)                  |  ∼ *G**a**m**m**a*(47.1, 304.1)                  | 0.11868                                           |
| Brazil        |  ∼ *G**a**m**m**a*(19.1, 462.1)                  |  ∼ *G**a**m**m**a*(17.1, 236.1)                  | 0.05015                                           |

Country-Specific Analysis

Is there evidence that the home-country advantage differs by country? We
can answer this question by studying the distributions of the ratios
*λ*<sub>*i*1</sub>/*λ*<sub>*i*0</sub>, which can be approximated using
Monte Carlo sampling. We will look at the probability that
*λ*<sub>*i*1</sub>/*λ*<sub>*i*0</sub> &gt; 1 for each country. If these
probabilities vary wildly, we have evidence that the home-country
advantage differs by country. Looking at the last column in Table 1, we
see these probabilities are spread out across the support. Thus, the
analysis shows home-country advantage does differ by country.
Specifically, we find evidence that Mexico, Soviet Union, and Spain
likely had home-country advantage
(*P*(*λ*<sub>*i*1</sub>/*λ*<sub>*i*0</sub> &gt; 1) close to 1); Finland,
Italy, Japan, Greece, and Brazil likely did not
(*P*(*λ*<sub>*i*1</sub>/*λ*<sub>*i*0</sub> &gt; 1) close to 0).

## Conclusions

Predominantly, host countries win more medals in the host Olympics than
four years earlier. Many attribute this to home-country advantage. We
did not discover any substantial evidence for this in our analysis. In
fact, we found evidence to the contrary: the rate of medals earned per
participant in the host Olympics is most likely lower than the rate of
medals earned four years earlier, not higher. There have only been a few
countries where we found evidence that a home-country advantage may have
contributed to the increased medal count during the host Olympics.

Due to limitations in the data, we were unable to consider some
questions that could determine if home-country advantage exists. Do a
higher proportion of participants from the host country win gold medals
in the host Olympics than in the previous Olympics? How does the rate of
medals earned per participant in the host country during the host
Olympics compare to the same rate during the next Olympics? Also, we do
not know how medals/participants in team events were calculated in the
data. Were medals awarded to teams counted as one or multiple? Were
participants in a team event counted as one or multiple. Answers to
these questions could revise our calculations of the rate of medals
earned per participant and change the results of the analysis.

## Appendix: R Code

``` r
#####  SETUP  #####

# Set global options for code chunks
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                      out.width = "85%", fig.align = "center")

# Load necessary packages
library(tidyverse)

# Download data
medals <- read_csv("C://ST540/Medals.csv")

#####  AGGREGATE ANALYSIS  #####

# For reproducibility
set.seed(316187263)

# For Monte Carlo sampling
S = 100000

# Aggregate data across all years
Y0 <- sum(medals$`MEDALS WON DURING PREVIOUS OLYMPICS`)
Y1 <- sum(medals$`MEDALS WON DURING HOST YEAR`)
N0 <- sum(medals$`PARTICIPATING ATHLETES DURING PREVIOUS OLYMPICS`)
N1 <- sum(medals$`PARTICIPATING ATHLETES DURING HOST YEAR`)

# Set values for the parameters in the prior
a = 1; b = 1

# Collect S random samples of lambda1 from its posterior distribution
lambda1 <- rgamma(S, Y1+a, N1+b)

# Collect S random samples of lambda0 from its posterior distribution
lambda0 <- rgamma(S, Y0+a, N0+b)

# Set grid for plotting
x <- seq(0.05, 0.25, 0.0001)

# Calculate the posterior densities for lambda1 and lambda0 at each grid value
posterior.lambda1 <- dgamma(x, Y1+a, N1+b)
posterior.lambda0 <- dgamma(x, Y0+a, N0+b)

# Create plot of the posterior distributions for lambda1 and lambda0
plot(x, posterior.lambda1, type = "l", col = "blue", lwd = 2,
     main = expression(paste("Posterior distributions for ", lambda[1],
                             "|", Y[1], " and ", lambda[0], "|", Y[0])),
     xlab = expression(lambda), ylab = "density")
lines(x, posterior.lambda0, col = "dark orange", lwd = 2)
legend("topleft", legend = c(expression(paste(lambda[1], "|", Y[1])),
                              expression(paste(lambda[0], "|", Y[0]))),
       lwd = 2, col = c("blue", "dark orange"))
text(0.18, 63, pos = 4, cex = 1.05,
     labels = expr(paste(mu[lambda[1]], " = ", !!round(mean(lambda1), 3))))
text(0.18, 55, pos = 4, cex = 1.05,
     labels = expr(paste(sigma[lambda[1]], " = ", !!round(sd(lambda1), 3))))
text(0.18, 42, pos = 4, cex = 1.05,
     labels = expr(paste(mu[lambda[0]], " = ", !!round(mean(lambda0), 3))))
text(0.18, 34, pos = 4, cex = 1.05,
     labels = expr(paste(sigma[lambda[0]], " = ", !!round(sd(lambda0), 3))))

#####  PREDICTION  #####

# For reproducibility
set.seed(316187263)

# Enter the data
Nstar <- 796

# Collect one random sample of Ystar from each of the S=100,000 likelihoods
Ystar <- rpois(S, Nstar*lambda1)

# Initialize the vector of probability mass values
pmf <- test <- rep(0, 100)

# Determine the probability mass of Ystar at each integer value from 51 to 150
for (i in 1:100){
  test <- (Ystar == i+50)
  pmf[i] <- mean(test)
}

# Set the grid for plotting
grid <- seq(51, 150, by = 1)

# Create plot of the posterior predictive distribution for Ystar
plot(grid, pmf, type="l", col = "blue", lwd = 2, xlab = "Y*", ylab = "density",
     main = "Posterior Predictive Distribution for Y*")
abline(v = round(mean(Ystar), 0), col = "darkorange", lwd = 2)
abline(v = quantile(Ystar, probs = 0.025), col = "red", lwd = 2)
abline(v = quantile(Ystar, probs = 0.975), col = "red", lwd = 2)
text(101, 0.005, labels = "mean = 101",
     col = "darkorange3", pos = 4, cex = 0.8)
text(81, 0.013, labels = "5th quantile = 81",
     col = "red3", pos = 2, cex = 0.8)
text(123, 0.016, labels = "95th quantile = 123",
     col = "red3", pos = 4, cex = 0.8)

#####  COUNTRY-SPECIFIC ANALYSIS  #####

# For reproducibility
set.seed(316187263)

# Collect S random samples of lambda1 for each country from its posterior
# distribution; do the same for lambda0 for each country; then calculate the
# probability that lamda1/lambda0 is greater than 1 for each country

lambda1.Finland <- rgamma(S, 22.1,258.1)
lambda0.Finland <- rgamma(S, 24.1,129.1)
# mean(lambda1.Finland/lambda0.Finland > 1)

lambda1.Australia <- rgamma(S, 92.1,911.1)
lambda0.Australia <- rgamma(S, 52.1,498.1)
#mean(lambda1.Australia/lambda0.Australia > 1)

lambda1.Italy <- rgamma(S, 36.1,280.1)
lambda0.Italy <- rgamma(S, 25.1,135.1)
#mean(lambda1.Italy/lambda0.Italy > 1)

lambda1.Japan <- rgamma(S, 80.1,949.1)
lambda0.Japan <- rgamma(S, 59.1,557.1)
#mean(lambda1.Japan/lambda0.Japan > 1)

lambda1.Mexico <- rgamma(S, 9.1,275.1)
lambda0.Mexico <- rgamma(S, 1.1,94.1)
#mean(lambda1.Mexico/lambda0.Mexico > 1)

lambda1.WGermany <- rgamma(S, 40.1,423.1)
lambda0.WGermany <- rgamma(S, 26.1,275.1)
#mean(lambda1.WGermany/lambda0.WGermany > 1)

lambda1.Canada <- rgamma(S, 11.1,385.1)
lambda0.Canada <- rgamma(S, 5.1,208.1)
#mean(lambda1.Canada/lambda0.Canada > 1)

lambda1.SovietUnion <- rgamma(S, 195.1,489.1)
lambda0.SovietUnion <- rgamma(S, 125.1,410.1)
#mean(lambda1.SovietUnion/lambda0.SovietUnion > 1)

lambda1.US <- rgamma(S, 275.1,1169.1)
lambda0.US <- rgamma(S, 202.1,941.1)
#mean(lambda1.US/lambda0.US > 1)

lambda1.SKorea <- rgamma(S, 33.1,401.1)
lambda0.SKorea <- rgamma(S, 19.1,175.1)
#mean(lambda1.SKorea/lambda0.SKorea > 1)

lambda1.Spain <- rgamma(S, 22.1,422.1)
lambda0.Spain <- rgamma(S, 4.1,229.1)
#mean(lambda1.Spain/lambda0.Spain > 1)

lambda1.Greece <- rgamma(S, 16.1,426.1)
lambda0.Greece <- rgamma(S, 13.1,140.1)
#mean(lambda1.Greece/lambda0.Greece > 1)

lambda1.China <- rgamma(S, 100.1,599.1)
lambda0.China <- rgamma(S, 63.1,384.1)
#mean(lambda1.China/lambda0.China > 1)

lambda1.GreatBritain <- rgamma(S, 65.1,530.1)
lambda0.GreatBritain <- rgamma(S, 47.1,304.1)
#mean(lambda1.GreatBritain/lambda0.GreatBritain > 1)

lambda1.Brazil <- rgamma(S, 19.1,462.1)
lambda0.Brazil <- rgamma(S, 17.1,236.1)
#mean(lambda1.Brazil/lambda0.Brazil > 1)
```
