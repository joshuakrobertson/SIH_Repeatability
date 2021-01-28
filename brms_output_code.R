# Compute Canada Analyses

# Acute eye temperature

library("easypackages")
library("devtools")
install_github("joshuakrobertson/brmsMethods")
libraries(
  "bayesboot", "BayesFactor", "bayesplot", "bigleaf",
  "brms", "brmsMethods", "emmeans", "gridExtra",
  "itsadug", "latex2exp", "loo", "mgcv",
  "rstan", "rstanarm", "Thermimage", "tidybayes", "tidyverse",
  "viridis"
)
rstan_options(auto_write=TRUE)
options(mc.cores=4)

my.theme <- theme(
  panel.grid.minor = element_blank(),
  axis.title = element_text(size = 14, family = "Noto Sans"),
  axis.text = element_text(size = 12, colour = "black", family = "Noto Sans"),
  axis.title.y = element_text(vjust = 0.5), panel.grid.major =
    element_line(colour = "grey75"), legend.title = element_text(
    size = 12,
    colour = "black", family = "Noto Sans"
  ), legend.text = element_text(
    size = 12,
    colour = "black", family = "Noto Sans"
  )
)

Third_Bound = read.csv("/home/joshk/git_repositories/BCCH_IndVar/Data/BCCH_Dom_Scramble.csv")
Third_Bound$Treatment = factor(Third_Bound$Treatment, ordered = TRUE)
Third_Bound$Treat_Day = paste(Third_Bound$Treatment, Third_Bound$Date.of.Photo, sep = "_")
Third_Bound$Treat_Day = factor(Third_Bound$Treat_Day)
Third_Bound$Maximum.Eye.Temp[c(which(Third_Bound$Bird.ID == "YAO1"))] = NA
Third_Bound$Sex = factor(Third_Bound$Sex)

prior_acute <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "Sex"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  set_prior("gamma(4, 2)", class = "b", coef = "sAmb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(60, 2)", class = "Intercept")
)

Acute_Full_DS <- brm(bf(
  Maximum.Eye.Temp ~ Treatment + Sex +
    s(Amb.Temp, k = 4, bs = "cr") +
    s(Amb.Temp, k = 4, bs = "cr", by = Treatment, m = 1) +
    s(Timeline, k = 5, bs = "cc") +
    s(Timeline, k = 5, bs = "cc", by = Treatment, m = 1) +
    t2(Timeline, Amb.Temp,
      by = Treatment, k = c(5, 4),
      bs = c("cc", "fs"), m = 1
    ) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re"), full = TRUE) +
    s(Timeline, Bird_Scramble, bs = "re") +
    s(Timeline, by = Treatment, Bird_Scramble, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird_Scramble),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian", iter = 10000,
warmup = 1000, thin = 10, control = list(
  adapt_delta = 0.98,
  max_treedepth = 14
), prior = prior_acute,
knots = list(Timeline = c(-1200, 0, 1200, 2400, 3600)),
file = "/home/joshk/git_repositories/SIH_Repeatability/brms_outputs/Acute_HT_DomS.Rds")




# Acute heat transfer

dim = 0.011
Area = ((1.1/2)*(1.0/2)*pi)*0.0001
Area = Area*2 # For two eyes

qrad_Clean = c()
for (i in 1:nrow(Third_Bound)) {
  qrad_Clean[i] = Area * (5.67 * 10^-8) * 0.95 * 0.95 *
    ((Third_Bound$Maximum.Eye.Temp[i] + 273.15)^4 -
(Third_Bound$Maximum.Eye.Temp[i] + 273.15)^4)
}

Kt_clean = c()
v_clean = c()
alpha_clean = c()
Gr_clean = c()
Nu_clean = c()
Hc_clean = c()
qconv_clean = c()

for (i in 1:nrow(Third_Bound)) {
  Kt_clean[i] = airtconductivity(Third_Bound$Amb.Temp[i])
  v_clean[i] = kinematic.viscosity(Third_Bound$Amb.Temp[i], 101.325)
  alpha_clean[i] = 1 / (Third_Bound$Amb.Temp[i] + 273.15)
}

for (i in 1:nrow(Third_Bound)) {
  Gr_clean[i] = ((alpha_clean[i]) * (9.81) * (dim)^3 *
                   (Third_Bound$Maximum.Eye.Temp[i] - Third_Bound$Amb.Temp[i])) / (v_clean[i])^2
}

for (i in 1:nrow(Third_Bound)) {
  Nu_clean[i] = sign(Gr_clean[i]) * 0.50 * (abs(Gr_clean[i]))^0.25
}

for (i in 1:nrow(Third_Bound)) {
  Hc_clean[i] = Nu_clean[i] * (Kt_clean[i] / dim)
}

for (i in 1:nrow(Third_Bound)) {
  qconv_clean[i] = Area * (Hc_clean[i]) *
(Third_Bound$Maximum.Eye.Temp[i] - Third_Bound$Amb.Temp[i])
}

qtot_clean = c()

for (i in 1:nrow(Third_Bound)) {
  qtot_clean[i] = qconv_clean[i] + qrad_Clean[i]
}

Third_Bound$qtot = qtot_clean

# Converting W to mW

Third_Bound$mW = Third_Bound$qtot*1000

prior_acute <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  set_prior("normal(-5, 5)", class = "b", coef = "sAmb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(50, 2)", class = "Intercept")
)

Acute_Full_HT_DS <- brm(bf(
  mW ~ Treatment + Sex +
    s(Amb.Temp, k = 4, bs = "cr") +
    s(Amb.Temp, k = 4, bs = "cr", by = Treatment, m = 1) +
    s(Timeline, k = 5, bs = "cc") +
    s(Timeline, k = 5, bs = "cc", by = Treatment, m = 1) +
    t2(Timeline, Amb.Temp,
      by = Treatment, k = c(5, 4),
      bs = c("cc", "fs"), m = 1
    ) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re"), full = TRUE) +
    s(Timeline, Bird_Scramble, bs = "re") +
    s(Timeline, by = Treatment, Bird_Scramble, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird_Scramble),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian", iter = 10000,
warmup = 1000, thin = 10, control = list(
  adapt_delta = 0.98,
  max_treedepth = 14
), prior = prior_acute,
knots = list(Timeline = c(-1200, 0, 1200, 2400, 3600)),
file = "/home/joshk/git_repositories/SIH_Repeatability/brms_outputs/Acute_HT_DomS.Rds")

# Checking posteriors 

pp_check(Acute_Full_HT_DS, nsamples = 100) # Good, but sbutle straying on the low end. Checking across treatments with violin

yrep_Acute <- posterior_predict(Acute_Full_HT_DS, nsamples = 100)

Acute_y <- Third_Bound %>%
  dplyr::select(
    mW, Treatment, Sex,
    Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo, Timeline
  ) %>%
  na.omit(.) %>%
  pull(mW)

Acute_Group <- Third_Bound %>%
  dplyr::select(
    mW, Treatment, Sex,
    Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo, Timeline
  ) %>%
  na.omit(.) %>%
  pull(Treatment)

ppc_violin_grouped(Acute_y, yrep_Acute, group = Acute_Group, 
                   probs = c(.05, .95), alpha = 0.05, y_draw = "points")

# Appears okay. Checking fit by treatment

ppc_scatter_avg_grouped(Acute_y, yrep_Acute, group = Acute_Group)

# Good. 

stan_ac(Acute_Full_HT_DS$fit) # Autocorrelation low
mcmc_neff(neff_ratio(Acute_Full_HT_DS), size = 2)
mcmc_rhat(rhat(Acute_Full_HT_DS))

# Ratio of effective sample size to sample size is slightly low, but passable (> 0.6). Rhat values okay.

Acute_Full_HT_DS$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  group_by(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

# No clear trend, but some possible extreme values.

Acute_Full_HT_DS$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

# Looks great.

Acute_Full_HT_DS$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = z_residual)) +
  geom_density()

# Good. Checking residuals by predictors.

{
  p1 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual, fill = Treatment)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Timeline, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p4 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual, fill = Sex)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p5 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Hour, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p6 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual, fill = Pen)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p7 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Direction, y = z_residual, fill = Direction)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p8 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird_Scramble, y = z_residual, fill = Bird_Scramble)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2, pch = 21,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)

# Taking a closer look at bird identities

print(p8) # SEs of most cross zero, and others draw close. Acceptable. 

# Visualising trends

conditional_smooths(Acute_Full_HT_DS) # Treatment by ambient temperature interaction weak, but reasonable perameterisation. 

# Running original model and comparing repeatabilities.

Acute_Full <- brm(bf(
  Maximum.Eye.Temp ~ Treatment + Sex +
    s(Amb.Temp, k = 4, bs = "cr") +
    s(Amb.Temp, k = 4, bs = "cr", by = Treatment, m = 1) +
    s(Timeline, k = 5, bs = "cc") +
    s(Timeline, k = 5, bs = "cc", by = Treatment, m = 1) +
    t2(Timeline, Amb.Temp,
      by = Treatment, k = c(5, 4),
      bs = c("cc", "fs"), m = 1
    ) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re"), full = TRUE) +
    s(Timeline, Bird.ID, bs = "re") +
    s(Timeline, by = Treatment, Bird.ID, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird.ID),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian", iter = 10000,
warmup = 1000, thin = 10, control = list(
  adapt_delta = 0.98,
  max_treedepth = 14
), prior = prior_acute,
knots = list(Timeline = c(-1200, 0, 1200, 2400, 3600)),
file =
  "/home/joshk/Desktop/CC_Results/Final/Acute_OP_Comp.Rds"
)

Rep_Treats_Acute <- rbind(
  data.frame(Acute_Full %>%
    spread_draws(
      sds_sTimelineBird.ID_1,
      sds_sTimelineBird.IDTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>%
    mutate(
      "Mod_Type" = "True Model",
      "C_Sigma" = exp(b_sigma_Intercept),
      "T_Sigma" = exp(b_sigma_Treatment.L + b_sigma_Intercept),
    ) %>%
    rename(
      "T_Bird" = sds_sTimelineBird.IDTreatmentStress_1,
      "C_Bird" = sds_sTimelineBird.ID_1
    ) %>%
    mutate(
      "C_Repeatability" = C_Bird^2 / (C_Bird^2 + C_Sigma^2),
      "T_Repeatability" = T_Bird^2 / (T_Bird^2 + T_Sigma^2)
    )),
  data.frame(Acute_Full_HT_DS %>%
    spread_draws(
      sds_sTimelineBird_Scramble_1,
      sds_sTimelineBird_ScrambleTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>%
    mutate(
      "Mod_Type" = "Null Model               ",
      "C_Sigma" = exp(b_sigma_Intercept),
      "T_Sigma" = exp(b_sigma_Treatment.L + b_sigma_Intercept),
    ) %>%
    rename(
      "T_Bird" = sds_sTimelineBird_ScrambleTreatmentStress_1,
      "C_Bird" = sds_sTimelineBird_Scramble_1
    ) %>%
    mutate(
      "C_Repeatability" = C_Bird^2 / (C_Bird^2 + C_Sigma^2),
      "T_Repeatability" = T_Bird^2 / (T_Bird^2 + T_Sigma^2)
    ))
)

S_Rep_Acute <- Rep_Treats_Acute %>% ggplot(aes(x = T_Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  #scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
  scale_fill_manual(values = c("black", "firebrick4"), name = NULL) + 
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Slopes") +
  my.theme +
  guides(fill = guide_legend(nrow = 1))

C_Rep_Acute <- Rep_Treats_Acute %>% ggplot(aes(x = C_Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  scale_fill_manual(values = c("black", "firebrick4"), name = NULL) + 
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Control Slopes") +
  my.theme +
  guides(fill = guide_legend(nrow = 1)) +
  theme(
    legend.text = element_text(size = 10, family = "Noto Sans"),
    legend.spacing.x = unit(0.4, "cm")
  )

grid.arrange(S_Rep_Acute, C_Rep_Acute, nrow = 2)

# Subtle increases among true models, but little. Saving plot then runnning formal comparison

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_Repeat_New.jpeg",
  S_Rep_Acute,
  height = 6.0, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_Repeat_New.pdf",
       S_Rep_Acute,
       height = 6.0, width = 6.0, dpi = 800
)

Rep_Treats_Acute %>%
  group_by(Mod_Type) %>%
  summarise(
    "Mean" = mean(T_Repeatability),
    "UCL" = quantile(T_Repeatability, 0.975, type = 8),
    "LCL" = quantile(T_Repeatability, 0.025, type = 8)
  )

Wide_Acute_Rep <- Rep_Treats_Acute %>%
  mutate("Sample" = paste(.chain, .iteration, sep = "_")) %>%
  dplyr::select(Sample, Mod_Type, T_Repeatability) %>%
  spread(Mod_Type, T_Repeatability) %>% 
  as.data.frame(.)
colnames(Wide_Acute_Rep)[c(2,3)] = c("Null", "True")

DF_Acute_Rep <- build_hdf(vars = list(Wide_Acute_Rep$True, Wide_Acute_Rep$Null), priors = list(rbeta(nrow(Wide_Acute_Rep), 1, 4), rbeta(nrow(Wide_Acute_Rep), 1, 4)), names = c("True", "Null"))
df_acute_hyp = hypothesis_df("True > Null", DF_Acute_Rep, class = "b", alpha = 0.05)
plot(df_acute_hyp)
print(df_acute_hyp)

# Still significant, but not by much.

##########################################################################################
## Chronic response ######################################################################
##########################################################################################

Third_Bound$s.Amb.Temp <- (Third_Bound$Amb.Temp - mean(Third_Bound$Amb.Temp, na.rm = T)) /
  (2 * sd(Third_Bound$Amb.Temp, na.rm = T))

# Setting prior

prior_chronic <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  # set_prior("gamma(2, 0.5)", class = "b", coef = "ss.Amb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(60, 2)", class = "Intercept")
)

CC_Full_DS <- brm(bf(
  Maximum.Eye.Temp ~ Treatment + Sex +
    s(s.Amb.Temp, k = 4, bs = "cr") +
    s(s.Amb.Temp, by = Treatment, k = 4, bs = "cr", m = 1) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re")) +
    s(s.Amb.Temp, Bird_Scramble, bs = "re") +
    s(s.Amb.Temp, by = Treatment, Bird_Scramble, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird_Scramble),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian",
iter = 10000, warmup = 1000, thin = 10,
control = list(adapt_delta = 0.98, max_treedepth = 14),
prior = prior_chronic,
file = "/home/joshk/git_repositories/SIH_Repeatability/brms_outputs/Chronic_Regular_DomScramble.Rds")

pp_check(CC_Full_DS, nsamples = 100) # Surprising skewing above 35 degrees. Running prior predictive check
View(prior_summary(CC_Full_DS))

# Note that flat priors are replaced with wide cauchy priors and a flat beta

prior_chronic_check <- c(
  set_prior("cauchy(0,2.5)", class = "b", dpar = "sigma"),
  set_prior("cauchy(0,2.5)", class = "b", coef = "ss.Amb.Temp_1"), 
  set_prior("cauchy(0,2.5)", class = "b", coef = "ss.Amb.Temp:TreatmentStress_1"),
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  set_prior("gamma(60, 2)", class = "Intercept"),
  set_prior("beta(1, 1)", class = "ar")
)

CC_Prior <- brm(bf(
  Maximum.Eye.Temp ~ Treatment + Sex +
    s(s.Amb.Temp, k = 4, bs = "cr") +
    s(s.Amb.Temp, by = Treatment, k = 4, bs = "cr", m = 1) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re")) +
    s(s.Amb.Temp, Bird_Scramble, bs = "re") +
    s(s.Amb.Temp, by = Treatment, Bird_Scramble, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird_Scramble),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 3, chains = 3,
seed = 100, refresh = 0, family = "gaussian",
iter = 10000, warmup = 1000, thin = 10,
control = list(adapt_delta = 0.95, max_treedepth = 13),
sample_prior = "only", prior = prior_chronic_check)

Prior_Preds <- posterior_predict(CC_Prior, nsamples = 100)
Prior_Preds <- cbind(data.frame("Rep" = c(1:100)), as.data.frame(Prior_Preds)) %>% 
  pivot_longer(!Rep, names_to = "Sample_Number", values_to = "Estimate") %>%
  mutate(Rep = factor(Rep))

ggplot(Prior_Preds, aes(x = Estimate, fill = Rep)) + 
  geom_density(adjust = 1, colour = "black", alpha = 0.5) + 
  theme_classic() + 
  scale_fill_viridis_d() + 
  xlab("Eye Temperature") + ylab("Count") + 
  theme(legend.position = "none") + 
  xlim(c(-20,60))

# Adding treatment into the mix

Prior_Dist <- Prior_Preds %>% mutate("Treatment" = 
  rep(Third_Bound %>% 
    drop_na(Maximum.Eye.Temp, Treatment, Sex,
      s.Amb.Temp, Hour, Direction, Bird_Scramble,
      Pen, Date.of.Photo) %>% 
      pull(Treatment), 100)
  ) %>% ggplot(aes(x = Estimate, fill = Rep)) + 
  geom_density(adjust = 1, colour = "black", alpha = 0.5) + 
  theme_classic() + 
  scale_fill_viridis_d() + 
  xlab("Eye Temperature") + ylab("Count") + 
  theme(legend.position = "none") + 
  xlim(c(-20,60)) + facet_wrap(~Treatment)

ggsave("/home/joshk/Desktop/Prior_predictive_check_dist.jpeg",
  Prior_Dist,
  height = 4.5, width = 9.0, dpi = 800
)

# Appears reasonable, although perhaps slightly overfit. Collapsing and plotting predicted values against real values

Prior_Preds %>% group_by(Sample_Number) %>% 
  summarise("Prior_Eye" = median(Estimate, na.rm = T)) %>% 
  ggplot(aes(x = Prior_Eye)) + 
  geom_density(adjust = 1, colour = "black", alpha = 0.5) + 
  theme_classic() + 
  scale_fill_viridis_d() + 
  xlab("Eye Temperature") + ylab("Count") + 
  theme(legend.position = "none")

Check_Dat <- Third_Bound %>% 
  drop_na(Maximum.Eye.Temp, Treatment, Sex,
    s.Amb.Temp, Hour, Direction, Bird_Scramble,
    Pen, Date.of.Photo) %>% 
    mutate("Prior_Eye" = 
      Prior_Preds %>% group_by(Sample_Number) %>% 
      summarise("Prior_Eye" = median(Estimate, na.rm = T)) %>% 
      pull(Prior_Eye)
    )  

Prior_check = ggplot(Check_Dat, aes(x = Prior_Eye, y = Maximum.Eye.Temp, 
  fill = Treatment, linetype = Treatment)) + 
  geom_point(size = 3, pch = 21, colour = "black", alpha = 0.5) + 
  theme_classic() + xlab("Prior Estimate") + ylab("True Eye Temp") + 
  scale_fill_viridis_d(begin = 0.3, end = 0.6)

ggsave("/home/joshk/Desktop/Prior_predictive_check.jpeg",
  Prior_check,
  height = 6.0, width = 7.0, dpi = 800
)

Prior_check_split = ggplot(Check_Dat, aes(x = Prior_Eye, y = Maximum.Eye.Temp, 
  fill = Treatment, linetype = Treatment)) + 
  geom_point(size = 3, pch = 21, colour = "black", alpha = 0.5) + 
  geom_smooth(method = "loess", colour = "black", size = 1, se = FALSE) + 
  theme_classic() + xlab("Prior Estimate") + ylab("True Eye Temp") + 
  scale_fill_viridis_d(begin = 0.3, end = 0.6) + facet_wrap(~Treatment)

ggsave("/home/joshk/Desktop/Prior_predictive_check_split.jpeg",
  Prior_check_split,
  height = 4.5, width = 9.0, dpi = 800
)

Birds <- Check_Dat %>% group_by(Treatment, Bird.ID) %>% 
  summarise("True" = var(Maximum.Eye.Temp, na.rm = T),
    "Prior" = var(Prior_Eye, na.rm = T)) %>% 
    pivot_longer(!(Treatment | Bird.ID), names_to = "Type", values_to = "Variance") %>%
    ggplot(aes(x = Type, y = Variance, fill = Treatment)) + 
    stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", size = 1, colour = "black", width = 0.2, position = position_dodge(width = 0.5)) +
    stat_summary(geom = "point", fun = "mean", size = 6, pch = 21, colour = "black", position = position_dodge(width = 0.5)) + 
    scale_fill_viridis_d(begin = 0.2, end = 0.6) + 
    theme_classic() + 
    xlab("Sample Type") + ylab("Mean Variance Around Eye Temp (Across Birds)")

Days <- Check_Dat %>% group_by(Treatment, Date.of.Photo) %>% 
  summarise("True" = var(Maximum.Eye.Temp, na.rm = T),
    "Prior" = var(Prior_Eye, na.rm = T)) %>% 
    pivot_longer(!(Treatment | Date.of.Photo), 
      names_to = "Type", values_to = "Variance") %>%
    ggplot(aes(x = Type, y = Variance, fill = Treatment)) + 
    stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", size = 1, colour = "black", width = 0.2, position = position_dodge(width = 0.5)) +
    stat_summary(geom = "point", fun = "mean", size = 6, pch = 21, colour = "black", position = position_dodge(width = 0.5)) + 
    scale_fill_viridis_d(begin = 0.2, end = 0.6) + 
    theme_classic() + 
    xlab("Sample Type") + ylab("Mean Variance Around Eye Temp (Across Days)") + 
    ylim(c(1,3))

Out <- arrangeGrob(grobs = list(Birds, Days), nrow = 1)

ggsave("/home/joshk/Desktop/Var_check_split.jpeg",
  Out, height = 4.5, width = 10.0, dpi = 800
)

# Appears okay, although quite tight. Trying across all iterations, then checking distributions

Prior_Preds %>% mutate("Treatment" = 
  rep(Third_Bound %>% 
    drop_na(Maximum.Eye.Temp, Treatment, Sex,
      s.Amb.Temp, Hour, Direction, Bird_Scramble,
      Pen, Date.of.Photo) %>% 
      pull(Treatment), 100),
  "Bird.ID" = rep(Third_Bound %>% 
    drop_na(Maximum.Eye.Temp, Treatment, Sex,
      s.Amb.Temp, Hour, Direction, Bird_Scramble,
      Pen, Date.of.Photo) %>% 
      pull(Bird.ID), 100),
  "True" = rep(Third_Bound %>% 
    drop_na(Maximum.Eye.Temp, Treatment, Sex,
      s.Amb.Temp, Hour, Direction, Bird_Scramble,
      Pen, Date.of.Photo) %>% 
      pull(Maximum.Eye.Temp), 100),  
  ) %>% rename("Prior" = Estimate) %>%
    pivot_longer(!(Rep | Sample_Number | Treatment | Bird.ID),
    names_to = "Type", values_to = "Eye.Temp") %>%
  group_by(Rep, Type, Treatment, Bird.ID) %>% 
  summarise("Variance" = var(Eye.Temp, na.rm = T)) %>% 
  pivot_wider(names_from = Type, values_from = Variance) %>% 
  group_by(Rep, Treatment) %>% 
  mutate(Rep = factor(Rep)) %>% 
  summarise("Variance_Ratio" = median(True/Prior, na.rm = T)) %>%
  mutate(Variance_Ratio = as.numeric(Variance_Ratio)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Variance_Ratio, fill = Treatment)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_classic() + theme(legend.position = "none") + 
  ylab("Density")

hist(Test$Variance_Ratio)




yrep_Acute <- posterior_predict(Acute_Full_HT_DS, nsamples = 100)

Acute_y <- Third_Bound %>%
  dplyr::select(
    mW, Treatment, Sex,
    Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo, Timeline
  ) %>%
  na.omit(.) %>%
  pull(mW)

Acute_Group <- Third_Bound %>%
  dplyr::select(
    mW, Treatment, Sex,
    Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo, Timeline
  ) %>%
  na.omit(.) %>%
  pull(Treatment)

ppc_violin_grouped(Acute_y, yrep_Acute, group = Acute_Group, 
                   probs = c(.05, .95), alpha = 0.05, y_draw = "points")

# Appears okay. Checking fit by treatment

ppc_scatter_avg_grouped(Acute_y, yrep_Acute, group = Acute_Group)

# Good. 

stan_ac(Acute_Full_HT_DS$fit) # Autocorrelation low
mcmc_neff(neff_ratio(Acute_Full_HT_DS), size = 2)
mcmc_rhat(rhat(Acute_Full_HT_DS))

# Ratio of effective sample size to sample size is slightly low, but passable (> 0.6). Rhat values okay.

Acute_Full_HT_DS$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  group_by(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

# No clear trend, but some possible extreme values.

Acute_Full_HT_DS$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

# Looks great.

Acute_Full_HT_DS$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT_DS) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = z_residual)) +
  geom_density()

# Good. Checking residuals by predictors.

{
  p1 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual, fill = Treatment)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Timeline, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p4 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual, fill = Sex)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p5 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Hour, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p6 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual, fill = Pen)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p7 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Direction, y = z_residual, fill = Direction)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p8 <- Acute_Full_HT_DS$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline, 
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT_DS) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird_Scramble, y = z_residual, fill = Bird_Scramble)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2, pch = 21,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)

# Taking a closer look at bird identities

print(p8) # SEs of most cross zero, and others draw close. Acceptable. 

# Visualising trends

conditional_smooths(Acute_Full_HT_DS)















Rep_Samples <- rbind(
  data.frame(CC_Full %>%
    spread_draws(
      sds_ss.Amb.TempBird.IDTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>%
    mutate(
      "Mod_Type" = "True Model",
      "Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L)
    ) %>%
    rename("Bird_ID" = sds_ss.Amb.TempBird.IDTreatmentStress_1) %>%
    mutate("Repeatability" = Bird_ID^2 / (Bird_ID^2 + Sigma^2))),
  data.frame(CC_Full_DS %>%
    spread_draws(
      sds_ss.Amb.TempBird_ScrambleTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>%
    mutate(
      "Mod_Type" = "Null Model               ",
      "Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L)
    ) %>%
    rename("Bird_ID" = sds_ss.Amb.TempBird_ScrambleTreatmentStress_1) %>%
    mutate("Repeatability" = Bird_ID^2 / (Bird_ID^2 + Sigma^2)))
)

Rep_Samples %>% ggplot(aes(x = Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  #scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
  scale_fill_manual(name = NULL, values = c("black", "firebrick4")) + 
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Slopes") +
  my.theme

Rep_Treats <- rbind(
  data.frame(CC_Full %>%
    spread_draws(
      sds_ss.Amb.TempBird.ID_1,
      sds_ss.Amb.TempBird.IDTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>%
    mutate(
      "Mod_Type" = "True\nModel  ",
      "C_Sigma" = exp(b_sigma_Intercept),
      "T_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L),
    ) %>%
    rename(
      "T_Bird" = sds_ss.Amb.TempBird.IDTreatmentStress_1,
      "C_Bird" = sds_ss.Amb.TempBird.ID_1
    ) %>%
    mutate(
      "C_Repeatability" = C_Bird^2 / (C_Bird^2 + C_Sigma^2),
      "T_Repeatability" = T_Bird^2 / (T_Bird^2 + T_Sigma^2)
    )),
  data.frame(CC_Full_DS %>%
    spread_draws(
      sds_ss.Amb.TempBird_Scramble_1,
      sds_ss.Amb.TempBird_ScrambleTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>%
    mutate(
      "Mod_Type" = "Null\nModel  ",
      "C_Sigma" = exp(b_sigma_Intercept),
      "T_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L),
    ) %>%
    rename(
      "T_Bird" = sds_ss.Amb.TempBird_ScrambleTreatmentStress_1,
      "C_Bird" = sds_ss.Amb.TempBird_Scramble_1
    ) %>%
    mutate(
      "C_Repeatability" = C_Bird^2 / (C_Bird^2 + C_Sigma^2),
      "T_Repeatability" = T_Bird^2 / (T_Bird^2 + T_Sigma^2)
    ))
)

Wide_Chronic_Rep <- Rep_Treats %>%
  mutate("Sample" = paste(.chain, .iteration, sep = "_")) %>%
  dplyr::select(Sample, Mod_Type, T_Repeatability) %>%
  spread(Mod_Type, T_Repeatability) %>%
  rename("True" = "True\nModel  ", "Null" = "Null\nModel  ") %>%
  as.data.frame(.)

DF <- build_hdf(vars = list(Wide_Chronic_Rep$True, Wide_Chronic_Rep$Null),
  priors = list(rbeta(nrow(Wide_Chronic_Rep), 1, 4), rbeta(nrow(Wide_Chronic_Rep), 1, 4)),
  names = c("True", "Null"))

df_hyp = hypothesis_df("True - Null > 0", DF, class = "b", alpha = 0.05)
df_hyp_cons = hypothesis_df("True - Null > 0.1", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)
plot(df_hyp_cons)
print(df_hyp_cons)





# Chronic heat-loss

prior_chronic <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  # set_prior("gamma(2, 0.5)", class = "b", coef = "ss.Amb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(50, 2)", class = "Intercept")
)

CC_Full <- brm(bf(
  mW ~ Treatment + Sex +
    s(s.Amb.Temp, k = 4, bs = "cr") +
    s(s.Amb.Temp, by = Treatment, k = 4, bs = "cr", m = 1) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re")) +
    s(s.Amb.Temp, Bird_Scramble, bs = "re") +
    s(s.Amb.Temp, by = Treatment, Bird_Scramble, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird_Scramble),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian",
iter = 10000, warmup = 1000, thin = 10,
control = list(adapt_delta = 0.98, max_treedepth = 14),
prior = prior_chronic,
file = "/home/joshk/git_repositories/SIH_Repeatability/brms_outputs/Chronic_HT_DomScramble.Rds")



