Comparing Staggered DiD
================
Florian M. Hollenbach
2021-12-23

# Comparing different staggered Difference-in-Differences Estimators

    ## Loading required package: fixest

    ## ── Attaching packages ────────────────────────────────── tidyverse 1.3.1.9000 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.1.1     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## ## See bit.ly/panelview4r for more info.
    ## ## Report bugs -> yiqingxu@stanford.edu.

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

First, we create a simulated data set, with staggered treatments,
heterogeneous and dynamic treatment effects. The code is based on the
simulations in Baker, Larcker, and Wang (2021), except we have a
never-treated group and decreased the number of units.[^1]

``` r
# Data 6 - Multiple Treatment Periods and Dynamic Treatment Effects --------------
make_data6 <- function(...) {
  
  # Fixed Effects ------------------------------------------------
  # unit fixed effects
  unit <- tibble(
    unit = 1:200, 
    unit_fe = rnorm(200, 0, 0.5),
    # generate state
    state = sample(1:50, 200, replace = TRUE),
    # generate treatment groups
    group = case_when(
      state %in% 1:10 ~ 1989,
      state %in% 11:20 ~ 1998,
      state %in% 21:35 ~ NA_real_, ## never treated
      state %in% 36:50 ~ 2005
    ),
    # avg yearly treatment effects by group
    hat_gamma = case_when(
      is.na(group) ~ 0, ## never treated
      group == 1989 ~ .5,
      group == 1998 ~ .3,
      group == 2005 ~ .1
    )) %>%
    # generate unit specific yearly treatment effects 
    rowwise() %>% 
    mutate(gamma = if_else(is.na(group) == TRUE, 0, rnorm(1, hat_gamma, .2))) %>% 
    ungroup()
  
  # year fixed effects 
  year <- tibble(
    year = 1980:2015,
    year_fe = rnorm(36, 0, 0.5))
  
  # full interaction of unit X year 
  crossing(unit, year) %>% 
    # make error term and get treatment indicators and treatment effects
    mutate(error = rnorm(nrow(.), 0, 0.5),
           treat = ifelse(year >= group & is.na(group)==F, 1, 0), # 0 for ## never treated
           tau = ifelse(treat == 1 & is.na(group)==F, gamma, 0)) %>% # 0 for ## never treated
    # calculate the dep variable
    group_by(unit) %>% 
    mutate(cumtau = cumsum(tau)) %>% 
    mutate(dep_var = unit_fe + year_fe + cumtau + error)
}

# make data
#and treatment group variable for CSA
data <- make_data6() %>% 
  as_tibble() %>%
  mutate(group_CSA = if_else(is.na(group), 0, group), # CSA wants never treated cohort variable to be 0
         group = if_else(is.na(group), 10000, group), # never treated cohort variable 10000 for fixest
         time_to_treatment = ifelse(group != 10000, year - group, -1000)) # set time to treatment to -1000 for fixest
```

We can plot the data and treatment status using Licheng Liu and Yiqing
Xu’s awesome `panelView` package (Liu and Xu 2021).

``` r
panelView(dep_var ~ treat, data = data, index = c("unit","year"), xlab = "Year", ylab = "Unit", axis.lab.gap = 5)
```

![](ComparingDiD_files/figure-gfm/plot%20treatment-1.png)<!-- --> Next,
we create the stacked data set, once again following the code by [Andrew
Baker](https://github.com/andrewchbaker/DiD_Codes).

``` r
### for stacking
groups <- data %>% 
  filter(group != 10000) %>% 
  pull(group) %>% 
  unique()

### create stacked data
getdata <- function(i) {
  
  #keep what we need
  data %>% 
    # keep treated units and all units not treated within -5 to 5
    filter(group == i | group > i + 5) %>% 
    # keep just year -5 to 5
    filter(year >= i - 5 & year <= i + 5) %>%
    # create an indicator for the dataset
    mutate(df = i) %>% 
    mutate(time_to_treatment = year - group) %>% 
    # make dummies
    mutate(time_to_treatment = if_else(group == i, time_to_treatment, 0))
}
stacked_data <- map_df(groups, getdata) %>% 
  mutate(bracket_df = paste(state,df))
```

Now we can move on to estimating the different models. First, the
standard two-way fixed effects model with dynamic event time estimates.
We estimate the model using the `fixest` package (Bergé 2018) and
extract the dynamic event time estimates.

``` r
twfe <- data %>% 
  do(broom::tidy(feols(dep_var ~ + i(time_to_treatment, ref = c(-1, -1000)) | unit + year, 
                       data = ., cluster = ~state), conf.int = TRUE)) %>% 
  mutate(t =  as.double(str_replace_all(term, c("time_to_treatment::" = "", ":treated" = "")))) %>% 
  filter(t >= -10 & t <=10) %>% 
  select(t, estimate, conf.low, conf.high) %>% 
  # add in data for year -1
  bind_rows(tibble(t = -1, estimate = 0, 
                   conf.low = 0, conf.high = 0
  )) %>% 
  mutate(method = "TWFE")
```

Next, the same model but on the stacked data. Following Baker, Larcker,
and Wang (2021), we cluster standard errors at the unit×dataset
interaction.

``` r
stacked <- stacked_data %>% 
  # fit the model
  do(broom::tidy(feols(dep_var ~ i(time_to_treatment, ref = c(-1, -1000)) | unit^df + year^df, data = ., cluster = "bracket_df"),
                 conf.int = TRUE)) %>% 
  filter(!(term %in% c("x1"))) %>% 
  mutate(t =  as.double(str_replace(term, "time_to_treatment::", ""))) %>% 
  filter(t >= -10 & t <=10) %>% 
  select(t, estimate, conf.low, conf.high) %>% 
  # add in data for year -1
  bind_rows(tibble(t = -1, estimate = 0, 
                   conf.low = 0, conf.high = 0
  )) %>% 
  mutate(method = "Stacked")
```

We continue using the `fixest` package and its `sunab` function to
estimate the dynamic effects using the Sun & Abraham method (Sun and
Abraham 2021).

``` r
SA <- data %>% 
  do(broom::tidy(feols(dep_var ~ sunab(group, year) | unit + year, data = .,
                 cluster = ~ state))) %>% 
  mutate(t =  as.double(str_replace(term, "year::", "")),
         conf.low = estimate - (qnorm(0.975)*std.error),
         conf.high = estimate + (qnorm(0.975)*std.error)) %>% 
  filter(t >= -10 & t <=10) %>% 
  select(t, estimate, conf.low, conf.high) %>% 
  mutate(method = "Sun & Abraham")
```

The next model to estimate is the doubly-robust estimator developed by
Callaway and Sant’Anna (2021b) and available in the `did` package
(Callaway and Sant’Anna 2021a). We use the `not-yet-treated` as the
control group, standard errors are clustered at the treatment level
(state).

``` r
csa.est<- att_gt(yname= 'dep_var',
             tname= 'year',
             idname = 'unit',
             gname = 'group_CSA',
             clustervars = 'state',
             est_method = 'dr',
             control_group = 'not-yet-treated',
             data = data) 

CSA <- aggte(csa.est, type = "dynamic", na.rm = TRUE) %>% 
  tidy() %>% 
  rename(t = event.time) %>% 
  filter(t >= -10 & t <=10) %>% 
  select(t, estimate, conf.low, conf.high) %>% 
  mutate(method = "CSA")
```

Now we use the `didimputation` package written by Kyle Butts
\[butts.2021.didimputation\] based on the paper by Borusyak, Jaravel,
and Spiess (2021).

``` r
did_imp <- did_imputation(data = data, yname = "dep_var", gname = "group_CSA",
                          tname = "year", idname = "unit", 
                          horizon=TRUE, pretrends = -10:-1) 
coef_imp <- did_imp %>% 
  select(t = term, estimate, std.error) %>%
  mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    t = as.numeric(t)
  ) %>%
  mutate(method = "DID Imputation") %>% 
  select(c(t, estimate, conf.low, conf.high, method)) %>% 
  filter(t >= -10 & t <= 10)
```

Next, we add the augmented synthetic control estimates for staggered
adoption (Ben-Michael, Feller, and Rothstein 2021) using the `augsynth`
package provided by Ben-Michael (2021).

``` r
asyn_res <- multisynth(dep_var ~ treat,
                   unit, 
                   year, 
                   data)

asyn <- summary(asyn_res)$att %>% 
  filter(Time > -10 & Time < 10 & (Level == 'Average')) %>%
  rename(t = Time, estimate = Estimate, conf.low = lower_bound, conf.high = upper_bound) %>% 
  mutate(method = "Aug. Synth") %>% 
  select(c(t, estimate, conf.low, conf.high, method))
```

Lastly, we use the `fect` package by Liu et al. (2021) and estimate a
counterfactul estimator chosen via cross-validation (Liu, Wang, and Xu
2021).

``` r
fect.res <- data %>% 
  fect(dep_var ~ treat, data = ., 
       index = c("unit","year"), 
       method = "both",
       CV = TRUE, 
       se = TRUE, 
       nboots = 500, 
       parallel = TRUE, 
       cv.treat = FALSE)
```

    ## Parallel computing ...
    ## Cross-validating ... 
    ## Criterion: Mean Squared Prediction Error
    ## Interactive fixed effects model...
    ## 
    ##  r = 0; sigma2 = 0.24222; IC = -0.99712; PC = 0.23017; MSPE = 0.24709; GMSPE = 0.06618; Moment = 0.07233; MSPTATT = 0.00220; MSE = 0.23251*
    ##  r = 1; sigma2 = 0.23302; IC = -0.61861; PC = 0.27625; MSPE = 0.27944; GMSPE = 0.07793; Moment = 0.08554; MSPTATT = 0.00207; MSE = 0.20686
    ##  r = 2; sigma2 = 0.22663; IC = -0.23278; PC = 0.32208; MSPE = 0.33906; GMSPE = 0.08744; Moment = 0.07140; MSPTATT = 0.00211; MSE = 0.18861
    ##  r = 3; sigma2 = 0.21913; IC = 0.14362; PC = 0.36315; MSPE = 0.37461; GMSPE = 0.09406; Moment = 0.06142; MSPTATT = 0.00115; MSE = 0.16827
    ##  r = 4; sigma2 = 0.21324; IC = 0.52281; PC = 0.40381; MSPE = 0.47464; GMSPE = 0.11494; Moment = 0.08073; MSPTATT = 0.00103; MSE = 0.15085
    ##  r = 5; sigma2 = 0.20694; IC = 0.89573; PC = 0.44091; MSPE = 0.63652; GMSPE = 0.13871; Moment = 0.09578; MSPTATT = 0.00096; MSE = 0.13438
    ## 
    ##  r* = 0
    ## 
    ## Matrix completion method...
    ## 
    ##  lambda.norm = 1.00000; MSPE = 0.24709; GMSPE = 0.06618; Moment = 0.07233; MSPTATT = 0.00220; MSE = 0.23251*
    ##  lambda.norm = 0.42170; MSPE = 0.26313; GMSPE = 0.07681; Moment = 0.07163; MSPTATT = 0.00066; MSE = 0.07581
    ##  lambda.norm = 0.17783; MSPE = 0.26848; GMSPE = 0.07548; Moment = 0.07174; MSPTATT = 0.00012; MSE = 0.01397
    ##  lambda.norm = 0.07499; MSPE = 0.26538; GMSPE = 0.07532; Moment = 0.07183; MSPTATT = 0.00002; MSE = 0.00249
    ## 
    ##  lambda.norm* = 1
    ## 
    ## 
    ## 
    ##  Recommended method through cross-validation: ife
    ## 
    ## Bootstrapping for uncertainties ... 500 runs
    ## Cannot use full pre-treatment periods. The first period is removed.
    ## Call:
    ## fect.formula(formula = dep_var ~ treat, data = ., index = c("unit", 
    ##     "year"), CV = TRUE, cv.treat = FALSE, method = "both", se = TRUE, 
    ##     nboots = 500, parallel = TRUE)
    ## 
    ## ATT:
    ##                             ATT   S.E. CI.lower CI.upper p.value
    ## Tr obs equally weighted   3.916 0.3202    3.288    4.543       0
    ## Tr units equally weighted 2.972 0.2689    2.445    3.499       0

``` r
fect <- fect.res$est.att %>% 
  as_tibble() %>% 
  mutate(t = as.double(rownames(fect.res$est.att))) %>% 
  filter(t > -10 & t < 10) %>% 
  mutate(method = "FECT") %>% 
  rename(estimate = ATT, conf.low = CI.lower, conf.high = CI.upper) %>% 
  select(c(t, estimate, conf.low, conf.high, method))
```

``` r
coefs <- bind_rows(twfe, stacked, CSA, SA, coef_imp, asyn, fect) 

plot <- coefs %>% 
  ggplot(aes(x = t, y = estimate, color = method)) + 
  geom_point(aes(x = t, y = estimate), position = position_dodge2(width = 0.75), size = 1) +
  geom_linerange(aes(x = t, ymin = conf.low, ymax = conf.high), position = position_dodge2(width = 0.75), size = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = .25, alpha = 0.75) + 
  geom_vline(xintercept = -0.5, linetype = "dashed", size = .25) +
  scale_color_manual(name="Estimation Method", values= met.brewer("Hokusai1", 7, "discrete")) +
  theme_clean() + theme(legend.position= 'bottom') +
  labs(title = 'Event Time Estimates', y="ATT", x = "Relative Time") + 
  guides(col = guide_legend(nrow = 3)) 
plot
```

![](ComparingDiD_files/figure-gfm/plot-1.png)<!-- -->

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-baker.2021.how" class="csl-entry">

Baker, Andrew, David F. Larcker, and Charles C. Y Wang. 2021. “How Much
Should We Trust Staggered Difference-in-Differences Estimates?”
<http://dx.doi.org/10.2139/ssrn.3794018>.

</div>

<div id="ref-eli.2021.augsynth" class="csl-entry">

Ben-Michael, Eli. 2021. *Augsynth: The Augmented Synthetic Control
Method*.

</div>

<div id="ref-ben-michael.2021.augmented" class="csl-entry">

Ben-Michael, Eli, Avi Feller, and Jesse Rothstein. 2021. “The Augmented
Synthetic Control Method.” *Journal of the American Statistical
Association* 116 (536): 1789–1803.
<https://doi.org/10.1080/01621459.2021.1929245>.

</div>

<div id="ref-berge.2018.efficient" class="csl-entry">

Bergé, Laurent. 2018. “Efficient Estimation of Maximum Likelihood Models
with Multiple Fixed-Effects: The R Package FENmlm.” *CREA Discussion
Papers*, no. 13.

</div>

<div id="ref-borusyak.2021.revisiting" class="csl-entry">

Borusyak, Kirill, Xavier Jaravel, and Jann Spiess. 2021. “Revisiting
Event Study Designs: Robust and Efficient Estimation.”
<https://www.dropbox.com/s/y92mmyndlbkufo1/Draft_RobustAndEfficient.pdf?raw=1>.

</div>

<div id="ref-callaway.2021.did" class="csl-entry">

Callaway, Brantly, and Pedro H. C. Sant’Anna. 2021a. “Did: Difference in
Differences.” <https://bcallaway11.github.io/did/>.

</div>

<div id="ref-callaway.2021.difference" class="csl-entry">

———. 2021b. “Difference-in-Differences with Multiple Time Periods.”
*Journal of Econometrics*.
<https://doi.org/10.1016/j.jeconom.2020.12.001>.

</div>

<div id="ref-liu.2021.practical" class="csl-entry">

Liu, Licheng, Ye Wang, and Yiqing Xu. 2021. “A Practical Guide to
Counterfactual Estimators for Causal Inference with Time-Series
Cross-Sectional Data.” <http://dx.doi.org/10.2139/ssrn.3555463>.

</div>

<div id="ref-liu.2021.fect" class="csl-entry">

Liu, Licheng, Ye Wang, Yiqing Xu, and Ziyi Liu. 2021. *Fect: Fixed
Effects Counterfactuals*.
<https://yiqingxu.org/packages/fect/fect.html>.

</div>

<div id="ref-liu.2021.panelview" class="csl-entry">

Liu, Licheng, and Yiqing Xu. 2021. *panelView: Visualizing Panel Data*.
<https://yiqingxu.org/packages/panelView/panelView.html>.

</div>

<div id="ref-sun.2021.estimating" class="csl-entry">

Sun, Liyang, and Sarah Abraham. 2021. “Estimating Dynamic Treatment
Effects in Event Studies with Heterogeneous Treatment Effects.” *Journal
of Econometrics* 225: 175–99.
https://doi.org/<https://doi.org/10.1016/j.jeconom.2020.09.006>.

</div>

</div>

[^1]: Thanks to Andrew Baker for sharing his code
    [here](https://github.com/andrewchbaker/DiD_Codes).
