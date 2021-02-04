# Loading in packages

library("easypackages")
library("devtools")
install_github("joshuakrobertson/brmsMethods")
libraries(
  "bayesboot", "BayesFactor", "bayesplot",
  "brms", "brmsMethods", "emmeans", "gridExtra",
  "itsadug", "latex2exp", "loo", "mgcv",
  "rstan", "rstanarm", "tidybayes", "tidyverse",
  "viridis"
)

# Adding plot theme. Note that this theme relies on use of Noto sans font
# pulled from the "showtext" package. This packages is loaded in below
# when required.

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

# Reading in, levelling, and randomising data for null models. Note that
# randomisation is restricted to within treatments.

Third_Bound <- read.csv("/home/joshk/git_repositories/BCCH_IndVar/Data/BCCH_Thermal_Data_Control_Matched_All.csv")
Third_Bound <- Third_Bound %>%
  mutate(
    Stress_Applied = factor(Third_Bound$Stress_Applied),
    Treatment = factor(Third_Bound$Treatment, ordered = T),
    Date.of.Photo = factor(Third_Bound$Date.of.Photo),
    Direction = ifelse(Pen == "NW" | Pen == "SW", "West", "East")
  ) %>%
  mutate(Direction = factor(Direction))

Third_Bound <- Third_Bound %>%
  mutate(
    "RowID" = 1:nrow(Third_Bound),
    "Bird_Present" = ifelse(is.na(Bird.ID), "NO", "YES"),
    "Group" = ifelse(Pen == "NW" | Pen == "SE", "G1", "G2")
  ) %>%
  mutate(Bird_Present = factor(Bird_Present), Group = factor(Group))

Third_Bound_Split <- split(Third_Bound, f = Third_Bound$Group)

for (i in 1:length(Third_Bound_Split)) {
  Third_Bound_Split[[i]] <- split(Third_Bound_Split[[i]], f = Third_Bound_Split[[i]]$Bird_Present)
}

# Setting seed for randomisation

set.seed(1001)

# Randomly sampling bird identities and assigning them to data points

for (j in 1:2) {
  Bird_Scramble <- c()
  for (i in 1:nrow(Third_Bound_Split[[j]][[2]])) {
    Bird_Scramble[i] <- paste(sample(unique(Third_Bound_Split[[j]][[2]]$Bird.ID), 1))
  }
  Third_Bound_Split[[j]][[2]]$Bird_Scramble <- factor(Bird_Scramble)
}

# Re-joining data and arranging

for (i in 1:2) {
  Third_Bound_Split[[i]] <- bind_rows(Third_Bound_Split[[i]])
}

Third_Bound_Final <- bind_rows(Third_Bound_Split) %>%
  arrange(RowID) %>%
  dplyr::select(-c(Bird_Present, Group, RowID)) %>%
  mutate(Bird_Scramble = factor(Bird_Scramble))

write.csv(Third_Bound_Final, "/home/joshk/git_repositories/BCCH_IndVar/Data/BCCH_Final_Broad.csv", row.names = FALSE)

# Re-loading and adjusting data. Note that ambient temperature is centred and scaled by
# 2 standard deviations for chronic analyses.

Third_Bound <- read.csv("/home/joshk/git_repositories/BCCH_IndVar/Data/BCCH_Final_Broad.csv")
Third_Bound$Treatment <- factor(Third_Bound$Treatment, ordered = TRUE)
Third_Bound$Treat_Day <- paste(Third_Bound$Treatment, Third_Bound$Date.of.Photo, sep = "_")
Third_Bound$Treat_Day <- factor(Third_Bound$Treat_Day)
Third_Bound$s.Amb.Temp <- (Third_Bound$Amb.Temp - mean(Third_Bound$Amb.Temp, na.rm = T)) /
  (2 * sd(Third_Bound$Amb.Temp, na.rm = T))

# Assessing number of observations per individual

Third_Bound %>%
  drop_na(Maximum.Eye.Temp) %>%
  group_by(Bird.ID) %>%
  summarise(Count = n()) %>%
  as.data.frame(.)

# Removing YAO1 due to very low number of observations.

Third_Bound$Maximum.Eye.Temp[c(which(Third_Bound$Bird.ID == "YAO1"))] <- NA

# Assigning chronic prior

prior_chronic <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  # set_prior("gamma(2, 0.5)", class = "b", coef = "ss.Amb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(60, 2)", class = "Intercept")
)

# Running chronic models

CC_Full <- brm(bf(
  Maximum.Eye.Temp ~ Treatment + Sex +
    s(s.Amb.Temp, k = 4, bs = "cr") +
    s(s.Amb.Temp, by = Treatment, k = 4, bs = "cr", m = 1) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re")) +
    s(s.Amb.Temp, Bird.ID, bs = "re") +
    s(s.Amb.Temp, by = Treatment, Bird.ID, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird.ID),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian",
iter = 10000, warmup = 1000, thin = 10,
control = list(adapt_delta = 0.98, max_treedepth = 14),
prior = prior_chronic,
file =
  "/home/joshk/Desktop/CC_Results/Final/Chronic_Corrected.Rds"
)

# Calculating descriptive statistics

Descriptive_Chronic = expand.grid(
  "s.Amb.Temp" = c(min(CC_Full$data$s.Amb.Temp), max(CC_Full$data$s.Amb.Temp)),
  "Hour" = c(8, 12, 16),
  "Sex" = na.omit(unique(CC_Full$data$Sex)),
  "Direction" = na.omit(unique(CC_Full$data$Direction)),
  "Treatment" = na.omit(unique(CC_Full$data$Treatment)),
  "Bird.ID" = na.omit(unique(CC_Full$data$Bird.ID)),
  "Treat_Day" = na.omit(unique(CC_Full$data$Treat_Day))[1]
)

Global_Pred_Chronic = posterior_predict(CC_Full, newdata = Descriptive_Chronic, re_formula = NA,
  nsamples = 1000, summary = FALSE)

Descriptive_Chronic %>% mutate("Post" = colMeans(Global_Pred_Chronic)) %>% 
  group_split(Treatment) %>% 
  reduce(left_join, by = c("s.Amb.Temp", "Hour", "Sex", "Direction", "Bird.ID", "Treat_Day")) %>% 
  rename("Control_Temp" = Post.x, "Stress_Temp" = Post.y) %>% 
  dplyr::select(-c(Treatment.x, Treatment.y)) %>% 
  mutate("Diff" = Control_Temp - Stress_Temp) %>% 
  group_by(s.Amb.Temp) %>% 
  summarise("Mean_Diff" = mean(Diff), "LCL" = quantile(Diff, 0.025, type = 8),
    "UCL" = quantile(Diff, 0.975, type = 8),
    "SD_Diff" = sd(Diff))

# Pulling effective sample sizes from summary

Effs_Chronic = rstan::summary(CC_Full$fit)$summary %>%
  as.data.frame(.) %>%
  dplyr::select(n_eff) %>%
  filter(rownames(.) %in% c(
    "b_Intercept", "b_Treatment.L", "b_SexMale",
    "sds_ss.Amb.Temp_1",
    "sds_ss.Amb.TempTreatmentStress_1",
    "sds_t2HourDirection_1",
    "sds_t2HourDirection_2",
    "sds_ss.Amb.TempBird.ID_1",
    "sds_ss.Amb.TempBird.IDTreatmentStress_1",
    "sd_Bird.ID__Intercept",
    "sd_Date.of.Photo__Intercept",
    "sd_Pen__Intercept",
    "b_sigma_Intercept",
    "b_sigma_Treatment.L"
  )) %>%
  mutate("Term" = gsub("_[[:digit:]]", "", rownames(.))) %>%
  group_by(Term) %>%
  summarise("Neff" = sum(n_eff))

Effs_Chronic

# Evaluating estimates and CIs for sigma parameters on original scale

CC_Full %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("C_Sigma" = exp(b_sigma_Intercept),
    "S_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L)) %>% 
    summarise("C_Est" = mean(C_Sigma), 
    "C_LCL" = quantile(C_Sigma, 0.025, type = 8),
    "C_UCL" = quantile(C_Sigma, 0.975, type = 8),
    "S_Est" = mean(S_Sigma), 
    "S_LCL" = quantile(S_Sigma, 0.025, type = 8),
    "S_UCL" = quantile(S_Sigma, 0.975, type = 8))

# Assessing AC and effective sample size:sample size ratio

stan_ac(CC_Full$fit)
mcmc_neff(neff_ratio(CC_Full), size = 2)

# Reasonable. All Neff/N > 0.75. Rhats?

mcmc_rhat(rhat(CC_Full))

# Great. Checking posterior and residual distributions.

pp_check(CC_Full, nsamples = 100, stat = "mean")
pp_check(CC_Full, nsamples = 100, stat = "median")
yrep <- posterior_predict(CC_Full, nsamples = 100)

y <- Third_Bound %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    s.Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo
  ) %>%
  na.omit(.) %>%
  pull(Maximum.Eye.Temp)

Group <- Third_Bound %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    s.Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo
  ) %>%
  na.omit(.) %>%
  pull(Treatment)

ppc_violin_grouped(y, yrep, group = Group, probs = c(.05, .95), alpha = 0.05, y_draw = "points")

# Some high ends missing in posterior samples, but otherwise a good fit.

ppc_scatter_avg_grouped(y, yrep, group = Group)

# Good, and slightly tighter fit for the stress group.

CC_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full) %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(CC_Full) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

CC_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

# Looks good.

CC_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(CC_Full) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Some left tailing, but appears reasonable. Seeing if
# trending manifests acrosss a given predictor.

{
  p1 <- CC_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- CC_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = s.Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- CC_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- CC_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Hour, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p5 <- CC_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- CC_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird.ID, y = z_residual, fill = Bird.ID)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Distributions look quite reasonable. Testing across response and fitted values.

F_by_R <- CC_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(
    x = Maximum.Eye.Temp, y = z_residual,
    fill = Treatment, linetype = Treatment
  )) +
  geom_point(size = 2, pch = 21, colour = "black") +
  scale_fill_viridis_d(begin = 0.3, end = 0.6)

F_by_P <- CC_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full) %>%
  summarise(
    Pred = .prediction,
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(
    x = Pred, y = z_residual,
    fill = Treatment, linetype = Treatment
  )) +
  geom_point(size = 2, pch = 21, colour = "black") +
  scale_fill_viridis_d(begin = 0.3, end = 0.6)

grid.arrange(F_by_R, F_by_P, nrow = 1)

# Checking fit

CC_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full, type = "response") %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = Maximum.Eye.Temp)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Tight fit.
# Calculating leave-one-out ("loo") residuals

loo_R2(CC_Full)

# Loo-R2 = 0.89. Testing differences in sigma between treatment types
# First double-checking direction of potential difference.

summary(CC_Full)
CC_Full %>% spread_draws(b_sigma_Intercept, b_sigma_Treatment.L) %>%
  rename("Int_S" = b_sigma_Intercept, "T_S" = b_sigma_Treatment.L) %>%
  summarise("Mean_Int_S" = mean(exp(Int_S)), "LCL_Int" = quantile(exp(Int_S), 0.025, type = 8),
  "UCL_Int" = quantile(exp(Int_S), 0.975, type = 8), "Mean_T_S" = mean(exp(Int_S + T_S)),
  "LCL_TS" = quantile(exp(Int_S + T_S), 0.025, type = 8), 
  "UCL_TS" = quantile(exp(Int_S + T_S), 0.975, type = 8))

CC_Full %>% spread_draws(b_sigma_Intercept, b_sigma_Treatment.L) %>%
  rename("Int_S" = b_sigma_Intercept, "T_S" = b_sigma_Treatment.L) %>% 
  mutate("Sigma_Difference" = Int_S - T_S) %>% 
  ggplot(aes(x = Sigma_Difference)) + 
  geom_density(colour = "black", fill = "cornflowerblue", alpha = 0.5, adjust = 3) + 
  xlim(c(-0.1,0.5)) + geom_vline(xintercept = 0, size = 1, colour = "firebrick") + 
  theme_bw()

# Good. Running hypothesis test:

Chronic_Hyp = hypothesis(CC_Full, "exp(sigma_Intercept) > exp(sigma_Intercept + sigma_Treatment.L)", alpha = 0.05)
plot(Chronic_Hyp)
print(Chronic_Hyp, digits = 3)

# Strong evidence for a reduction in variance among treatment birds; K = 48.32. Using custom priors.

sigma_chronic <- CC_Full %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("Control_Sigma" = exp(b_sigma_Intercept),
    "Stress_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L))

DF_Sig_Chronic <- build_hdf(vars = list(sigma_chronic$Control_Sigma, sigma_chronic$Stress_Sigma),
  priors = list(rnorm(nrow(sigma_chronic), 0, 0.25), rnorm(nrow(sigma_chronic), 0, 0.25)),
  names = c("Control", "Stress"))

df_sigma_chronic_hyp = hypothesis_df("Stress = Control", DF_Sig_Chronic, class = "b", alpha = 0.05)

plot(df_sigma_chronic_hyp)
print(df_sigma_chronic_hyp)

df_sigma_chronic_hyp = hypothesis_df("Stress < Control", DF_Sig_Chronic, class = "b", alpha = 0.05)

plot(df_sigma_chronic_hyp)
print(df_sigma_chronic_hyp)

# Saving supplemental plot 

sFig_2 = plot(df_sigma_chronic_hyp, plot = F, theme = theme_get())[[1]]

sFig_2_Final = sFig_2 + theme_bw() + ylab("Density") + xlab(TeX('$\\sigma^2_{Control} - \\sigma^2_{Stress}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/SFig2.jpeg",
  sFig_2_Final,
  height = 6.0, width = 7.0, dpi = 800
)

# Comparing repeatabilities between treatment groups using non-linear hypothesis test.

Treat_hyp_test_Chronic = c("(sds_ss.Amb.TempBird.IDTreatmentStress_1)^2/((sds_ss.Amb.TempBird.IDTreatmentStress_1)^2 + (exp(b_sigma_Treatment.L))^2) - (sds_ss.Amb.TempBird.ID_1)^2/((sds_ss.Amb.TempBird.ID_1)^2 + (exp(b_sigma_Intercept))^2) > 0.1")
THT_Chronic_Out = hypothesis(CC_Full, Treat_hyp_test_Chronic, alpha = 0.05, class = NULL)
plot(THT_Chronic_Out)
print(THT_Chronic_Out, digits = 3)

# And using custom priors

repeat_chronic_custom <- CC_Full %>%
    spread_draws(
      sds_ss.Amb.TempBird.ID_1,
      sds_ss.Amb.TempBird.IDTreatmentStress_1,
      b_sigma_Intercept,
      b_sigma_Treatment.L
     ) %>% 
    mutate("Control_Repeat" = sds_ss.Amb.TempBird.ID_1^2/(exp(b_sigma_Intercept)^2 + sds_ss.Amb.TempBird.ID_1^2),
    "Stress_Repeat" = sds_ss.Amb.TempBird.IDTreatmentStress_1^2/(exp(b_sigma_Intercept + b_sigma_Treatment.L)^2 + sds_ss.Amb.TempBird.IDTreatmentStress_1^2))

repeat_chronic_custom %>% summarise("R_Stress" = mean(Stress_Repeat), "S_LCL" = quantile(Stress_Repeat, 0.025, type = 8),
"S_UCL" = quantile(Stress_Repeat, 0.975, type = 8), "R_Control" = mean(Control_Repeat), "C_LCL" = quantile(Control_Repeat, 0.025, type = 8), "C_UCL" = quantile(Control_Repeat, 0.975, type = 8))

DF_chronic_rep <- build_hdf(vars = list(repeat_chronic_custom$Control_Repeat, repeat_chronic_custom$Stress_Repeat),
  priors = list(rbeta(nrow(repeat_chronic_custom), 1, 4), rbeta(nrow(repeat_chronic_custom), 1, 4)),
  names = c("Control", "Stress"))

df_chronic_rep_hyp = hypothesis_df("Stress = Control", DF_chronic_rep, class = "b", alpha = 0.05)
plot(df_chronic_rep_hyp)
print(df_chronic_rep_hyp)

df_chronic_rep_hyp = hypothesis_df("Stress > Control", DF_chronic_rep, class = "b", alpha = 0.05)
plot(df_chronic_rep_hyp)
print(df_chronic_rep_hyp)

# Again, saving supplemental plot

sFig_4 = plot(df_chronic_rep_hyp, plot = F, theme = theme_get())[[1]]

sFig_4_Final = sFig_4 + theme_bw() + ylab("Density") + xlab(TeX('$R_{Stress} - R_{Control}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/SFig4.jpeg",
  sFig_4_Final,
  height = 6.0, width = 7.0, dpi = 800
)

# Strong support for increased repeatability within stress exposed treatment (K = 13.63)

# Running null model.

CC_Null <- brm(bf(
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
file = "/home/joshk/Desktop/CC_Results/Final/Chronic_Null_Corrected.Rds"
)

# Again, checking residuals and fit. Note that residual patterns are unlikely to be
# clean owing to randomisation of bird identities.

stan_ac(CC_Null$fit)
mcmc_neff(neff_ratio(CC_Null), size = 2)
mcmc_rhat(rhat(CC_Null))
conditional_smooths(CC_Null)

# Reasonable. Minimum Neff/N falling slightly below 0.75, but Rhats suggest strong mixing.

CC_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Null) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

CC_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(CC_Null) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Reasonably normal.
# Across predictors:

{
  p1 <- CC_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- CC_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = s.Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- CC_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- CC_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Hour, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p5 <- CC_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- CC_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird_Scramble, y = z_residual, fill = Bird_Scramble)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Variance in quantile residuals per bird is considerably larger. Fit?

CC_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Null, type = "response") %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = Maximum.Eye.Temp)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

loo_R2(CC_Null)

# Loo-R2 is expectedly lower, but by litle. Presumably repeatability of slopes is low.
# Comparing repeatability estimates from each model.

get_variables(CC_Full)
get_variables(CC_Null)

# Note repeatability is calculated for the stress-exposure group
# alone here.

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
  data.frame(CC_Null %>%
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

library("showtext")

Rep_Plot <- Rep_Samples %>% ggplot(aes(x = Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  #scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
  scale_fill_manual(name = NULL, values = c("black", "firebrick4")) + 
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Slopes") +
  my.theme

Rep_Plot
showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Slopes_OP.jpeg", Rep_Plot,
  height = 6.0, width = 6.0, dpi = 800,
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Slopes_OP.pdf", Rep_Plot,
       height = 6.0, width = 6.0, dpi = 800,
)

## Plotting within each treatment type

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
  data.frame(CC_Null %>%
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

S_Rep <- Rep_Treats %>% ggplot(aes(x = T_Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Stress Slopes") +
  my.theme +
  guides(fill = guide_legend(nrow = 1)) +
  theme(
    legend.text = element_text(size = 10, family = "Noto Sans"),
    legend.spacing.x = unit(0.4, "cm")
  )

C_Rep <- Rep_Treats %>% ggplot(aes(x = C_Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
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

grid.arrange(S_Rep, C_Rep, nrow = 2)

# Stress slopes more repeatable, suggesting possibility for selection on these responses?
# Testing differences in repeatability witih treatment groups

BTT <- ttestBF(
  formula = T_Repeatability ~ Mod_Type, data = Rep_Treats,
  iterations = 10000, rscale = "ultrawide", posterior = FALSE
)

Post_Est <- posterior(BTT, iterations = 1000)
plot(Post_Est[, "mu"])
plot(acf(Post_Est[, "mu"]))

# Chains relatively stable and autocorrelation negligable. Assessing Bayes factor

summary(BTT)

# BF > 100. Strong evidence for alternative over null. Quantifying
# then using Savage-Dickey approach to calculating Bayes Factor.

Rep_Treats %>%
  group_by(Mod_Type) %>%
  summarise(
    "Mean" = mean(T_Repeatability),
    "UCL" = quantile(T_Repeatability, 0.975, type = 8),
    "LCL" = quantile(T_Repeatability, 0.025, type = 8)
  )

# No overlap of credible intervals.

Wide_Chronic_Rep <- Rep_Treats %>%
  mutate("Sample" = paste(.chain, .iteration, sep = "_")) %>%
  dplyr::select(Sample, Mod_Type, T_Repeatability) %>%
  spread(Mod_Type, T_Repeatability) %>%
  rename("True" = "True\nModel  ", "Null" = "Null\nModel  ") %>%
  as.data.frame(.)

hypothesis(Wide_Chronic_Rep, "True = Null", alpha = 0.05)
Keep = hypothesis(Wide_Chronic_Rep, "True - Null > 0.25", alpha = 0.05)
plot(hypothesis(Wide_Chronic_Rep, "True - Null > 0.25", alpha = 0.05))

# And with custom hypothesis function to load in priors 

DF <- build_hdf(vars = list(Wide_Chronic_Rep$True, Wide_Chronic_Rep$Null),
  priors = list(rbeta(nrow(Wide_Chronic_Rep), 1, 4), rbeta(nrow(Wide_Chronic_Rep), 1, 4)),
  names = c("True", "Null"))

df_hyp = hypothesis_df("True - Null > 0", DF, class = "b", alpha = 0.05)
df_hyp_cons = hypothesis_df("True - Null > 0.1", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)
plot(df_hyp_cons)
print(df_hyp_cons)

# Again, Bayes factor is extremely high (in this case, infinite)!
# And lastly, using a model in brms where variance is unequal between groups and prior is default.

Robust_Comp_data <- Rep_Treats %>%
  dplyr::select(Mod_Type, T_Repeatability) %>%
  mutate(Model = gsub("\\n.*", "", Mod_Type)) %>%
  dplyr::select(-Mod_Type)

# Checking priors first, but with temporary prior on sigma by treatment

robust_comp <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student,
  data = Robust_Comp_data, sample_prior = "only",
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.05)", class = "nu"),
    set_prior("exponential(2)", class = "b", dpar = "sigma", lb = 0.001)
  )
)

y_est <- posterior_predict(robust_comp)
prior_pc = data.frame("Type" = c(rep("Prior", length(colMeans(y_est))),
  rep("True", length(colMeans(y_est)))), 
  "Y" = c(colMeans(y_est), Robust_Comp_data$T_Repeatability))

ggplot(prior_pc, aes(x = Y, fill = Type)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.4, end = 0.6)

## Priors are very wide, but given that on sigma will drop, unlikely to be an issue.

robust_comp <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student, sample_prior = TRUE,
  data = Robust_Comp_data,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ),
  file = "/home/joshk/Desktop/CC_Results/Final/Robust_Comp_Chronic.Rds"
)

pp_check(robust_comp)
mcmc_neff(neff_ratio(robust_comp), size = 2)
mcmc_rhat(rhat(robust_comp))

# Clean coverage, Neff/N and Rhats.

plot(robust_comp)
summary(robust_comp)
plot(hypothesis(robust_comp, "ModelTrue > ModelNull"))
print(hypothesis(robust_comp, "ModelTrue > ModelNull"), digits = 3)
print(hypothesis(robust_comp, "ModelTrue - ModelNull > 0.25"), digits = 3)

### Pulling out qualitative differences in means and SDs

tidy_MCMC_Chronic <- tidyMCMC(robust_comp, conf.int = TRUE, conf.level = 0.95, 
  estimate.method = "median", conf.method = "HPDinterval") %>% 
  mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

Post_Dif_Chronic <- posterior_samples(robust_comp) %>% 
  mutate_at(vars(contains("sigma")), funs(exp)) %>% 
  mutate(nu = log10(nu)) %>% 
  mutate(diff_means = b_ModelTrue - b_ModelNull,
         diff_sigma = b_sigma_ModelTrue - b_sigma_ModelNull) %>% 
  mutate(cohen_d = diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)/2),
         cles = dnorm(diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)), 0, 1))

Out_Chronic <- tidyMCMC(Post_Dif_Chronic, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")
Out_Chronic

# Lastly, as permutation test 

set.seed(200) 

simu <- 1000
res <- numeric(simu) 

for (i in 1:simu) {
    perm <- sample(nrow(Robust_Comp_data))
    bdat <- transform(Robust_Comp_data, Model = Model[perm])
    res[i] <- mean(bdat[bdat$Model=="True", "T_Repeatability"])-
        mean(bdat[bdat$Model=="Null", "T_Repeatability"])
}

obs <- mean(Robust_Comp_data[Robust_Comp_data$Model == "True", "T_Repeatability"])-
    mean(Robust_Comp_data[Robust_Comp_data$Model == "Null", "T_Repeatability"])

ggplot(as.data.frame(res), aes(x = res)) + 
  geom_density(adjust = 3, fill = "mediumorchid", alpha = 0.5) + 
  geom_vline(xintercept = obs, colour = "black")

mean(abs(res)>=abs(obs))

# All results pointing in the same direction. 

# Extracting slopes from urban and rural birds

get_variables(CC_Full)

Slopes <- CC_Full %>%
  spread_draws(s_ss.Amb.TempBird.IDTreatmentStress_1[ID]) %>%
  rename("Coef" = s_ss.Amb.TempBird.IDTreatmentStress_1) %>%
  as.data.frame(.)

Birds <- data.frame(
  "Bird_ID" = rownames(ranef(CC_Full)$Bird.ID),
  "ID" = seq(1, 19, 1)
)
Collapsed <- left_join(Slopes, Birds, by = c("ID")) %>%
  dplyr::select(-ID) %>%
  group_by(Bird_ID) %>%
  summarise("Slope" = mean(Coef)) %>%
  as.data.frame(.)

# Binding in locale and treatment order

U_R <- c()
Locale <- c()

for (i in 1:nrow(Collapsed)) {
  if (is.na(Collapsed$Bird_ID[i])) {
    Locale[i] <- NA
    U_R[i] <- NA
  } else if (Collapsed$Bird_ID[i] == "BdABd") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdAO1") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "ABlBl") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "AYY") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "YAO1") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "OOA") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdABl1") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "YAO2") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "ABlO") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "ABdBd") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "YAR") {
    Locale[i] <- "Cambridge"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "RAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "ABlR") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "ARO") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdABl2") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BlAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdAO2") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "AOR") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "YAY") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "AYBd") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else {
    Locale[i] <- NA
    U_R[i] <- NA
  }
}
Collapsed$U_R <- factor(U_R)
Collapsed$Locale <- factor(Locale)
Orders <- Third_Bound %>%
  dplyr::select("Bird_ID" = Bird.ID, Pen) %>%
  na.omit(.) %>%
  distinct(.) %>%
  mutate("Order" = ifelse(Pen == "NE" | Pen == "SW", "SR", "RS")) %>%
  dplyr::select(-Pen)
Collapsed <- left_join(Collapsed, Orders, by = c("Bird_ID")) %>%
  mutate(Order = factor(Order), U_R = factor(U_R))

# Summarising and Running "ANOVA"

Collapsed %>%
  group_by(U_R) %>%
  summarise(
    Count = n(),
    Mean = mean(Slope),
    LCL = quantile(Slope, 0.025, type = 8),
    UCL = quantile(Slope, 0.975, type = 8)
  )

# Considerable overlap.

aovBF <- anovaBF(Slope ~ Order * U_R,
  data = Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(aovBF)

# Minimal evidence for order effect. Dropping and re-running.

aovBF <- anovaBF(Slope ~ U_R,
  data = Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(BayesFactor::posterior(aovBF, iterations = 10000)[, "mu"])
plot(acf(BayesFactor::posterior(aovBF, iterations = 10000)[, "mu"]))

# Chain appears relatively stable, although with some early peaks. No autocorrelation evident. Summarising

summary(aovBF)
summary(BayesFactor::posterior(aovBF, iterations = 10000)[, "mu"])

# No evidence for differences between urban and rural individuals. Sitting quite close to zero! Again, approaching using Savage-Dickey method. First checking priors.

Slope_Contrast <- brm(Slope ~ 0 + U_R + (1 | Locale),
  prior = c(
    set_prior("normal(0, 5)", class = "b"),
    set_prior("gamma(1, 0.5)", class = "sd"),
    set_prior("exponential(0.05)", class = "sigma")
  ), data = Collapsed, cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  sample_prior = "only")

y_est <- posterior_predict(Slope_Contrast)
prior_pc = data.frame("Type" = c(rep("Prior", length(colMeans(y_est))),
  rep("True", length(colMeans(y_est)))), 
  "Y" = c(colMeans(y_est), Collapsed$Slope))

ggplot(prior_pc, aes(x = Y, fill = Type)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.4, end = 0.6)

# Priors seem reasonably appropriate, however, true response is clearly noisy.

Slope_Contrast <- brm(Slope ~ 0 + U_R + (1 | Locale),
  prior = c(
    set_prior("normal(0, 5)", class = "b"),
    set_prior("gamma(1, 0.5)", class = "sd"),
    set_prior("exponential(0.05)", class = "sigma")
  ), data = Collapsed, cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  sample_prior = TRUE,
  file = "/home/joshk/Desktop/CC_Results/Final/Chronic_UR_Contrast_New.Rds"
)

Chronic_Hypoth <- hypothesis(Slope_Contrast, c(
  "U_RRural = U_RUrban"
))
plot(Chronic_Hypoth)
1 / Chronic_Hypoth$hypothesis$Evid.Ratio

# Results corroborate reasonably. Plotting.

UR_Plot <- Collapsed %>%
  rename("Ecotype" = U_R) %>%
  group_by(Ecotype) %>%
  summarise(
    "Coefficient" = mean(Slope),
    "lcl" = quantile(Slope, 0.025, type = 8),
    "ucl" = quantile(Slope, 0.975, type = 8)
  ) %>%
  ggplot(aes(x = Ecotype, y = Coefficient, fill = Ecotype)) +
  geom_errorbar(aes(x = Ecotype, ymin = lcl, ymax = ucl),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 6, pch = 21, colour = "black") +
  #scale_fill_viridis_d(begin = 0.2, end = 0.5) +
  scale_fill_manual(values = c("wheat4","cornsilk")) + 
  theme_bw() +
  my.theme +
  theme(legend.position = "none") +
  ylab("Reaction Norm Slope")

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Chronic_OP_New.jpeg",
  UR_Plot,
  height = 5.67, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Chronic_OP_New.pdf",
       UR_Plot,
       height = 5.67, width = 6.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Plotting individual variability

Simp_Dat_Chronic <- expand.grid(
  "s.Amb.Temp" = seq(min(CC_Full$data$s.Amb.Temp, na.rm = TRUE),
    max(CC_Full$data$s.Amb.Temp, na.rm = TRUE),
    length.out = 25
  ),
  "Hour" = c(12),
  "Date.of.Photo" = na.omit(unique(CC_Full$data$Date.of.Photo))[30],
  "Treatment" = na.omit(unique(CC_Full$data$Treatment)),
  "Pen" = na.omit(unique(CC_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(CC_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(CC_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(CC_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(CC_Full$data$Treat_Day))[c(14, 70)]
)

Third_Bound$s.Amb.Temp <- (Third_Bound$Amb.Temp - mean(Third_Bound$Amb.Temp, na.rm = T)) /
  (2 * sd(Third_Bound$Amb.Temp, na.rm = T))

Chronic_Plot_Dat <- Simp_Dat_Chronic %>%
  add_predicted_draws(CC_Full, n = 1000, scale = "response") %>%
  mutate(
    "Group_ID" = paste(Bird.ID, Treatment, sep = "_"),
    "Amb.Temp" = (s.Amb.Temp * (2 * sd(Third_Bound$Amb.Temp, na.rm = T))) +
      mean(Third_Bound$Amb.Temp, na.rm = T)
  ) %>%
  group_by(Group_ID, Bird.ID, Treatment, Amb.Temp) %>%
  summarise("Pred" = mean(.prediction))

Chronic_Curves <- Chronic_Plot_Dat %>%
  ggplot(aes(x = Amb.Temp, y = Pred, colour = Treatment)) +
  geom_line(aes(group = Group_ID)) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_OP.jpeg",
  Chronic_Curves,
  height = 6.0, width = 7.0, dpi = 800
)

# With labels for later interpretation

Chronic_Curves <- Chronic_Plot_Dat %>%
  ggplot(aes(x = Amb.Temp, y = Pred, linetype = Treatment, colour = Bird.ID)) +
  geom_line(aes(group = Group_ID)) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme

#### Acute models

prior_acute <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "Sex"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  set_prior("gamma(4, 2)", class = "b", coef = "sAmb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(60, 2)", class = "Intercept")
)

# And running full model.

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

# Checking posteriors

pp_check(Acute_Full, nsamples = 100)

# Nice ovelay but some straying around 36 degrees and above. Minimal, but 
# assessing with violoin plot.

yrep_Acute <- posterior_predict(Acute_Full, nsamples = 100)

Acute_y <- Third_Bound %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo, Timeline
  ) %>%
  na.omit(.) %>%
  pull(Maximum.Eye.Temp)

Acute_Group <- Third_Bound %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex,
    Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo, Timeline
  ) %>%
  na.omit(.) %>%
  pull(Treatment)

ppc_violin_grouped(Acute_y, yrep_Acute, group = Acute_Group, 
                   probs = c(.05, .95), alpha = 0.05, y_draw = "points")

# Yes, subtle underestimation at high end of stress group. Seeing if this 
# manifests in a scatter plot.

ppc_scatter_avg_grouped(Acute_y, yrep_Acute, group = Acute_Group)

# No, it does not seem so; underestimation at high end appears negligable. Checking effective sample sizes

effectiveSize(Acute_Full)
plot(Acute_Full)
summary(Acute_Full)
conditional_smooths(Acute_Full)
get_variables(Acute_Full)

Mean_Tensor_Coefs <- Acute_Full %>%
  spread_draws(
    sds_t2TimelineAmb.TempTreatmentStress_1,
    sds_t2TimelineAmb.TempTreatmentStress_2,
    sds_t2TimelineAmb.TempTreatmentStress_3,
    sds_t2HourDirection_1,
    sds_t2HourDirection_2,
    sds_t2HourDirection_3,
    b_sigma_Intercept,
    b_sigma_Treatment.L
  ) %>%
  rename(
    "TA_Tens_1" = sds_t2TimelineAmb.TempTreatmentStress_1,
    "TA_Tens_2" = sds_t2TimelineAmb.TempTreatmentStress_2,
    "TA_Tens_3" = sds_t2TimelineAmb.TempTreatmentStress_3,
    "Hour_Tens_1" = sds_t2HourDirection_1,
    "Hour_Tens_2" = sds_t2HourDirection_2,
    "Hour_Tens_3" = sds_t2HourDirection_3
  ) %>%
  summarise(
    "Mean_TAE" = mean(c(TA_Tens_1, TA_Tens_2, TA_Tens_3)),
    "TA_LCL" = quantile(c(TA_Tens_1, TA_Tens_2, TA_Tens_3),
      0.025,
      type = 8
    ),
    "TA_UCL" = quantile(c(TA_Tens_1, TA_Tens_2, TA_Tens_3),
      0.975,
      type = 8
    ),
    "Mean_HourE" = mean(c(Hour_Tens_1, Hour_Tens_2, Hour_Tens_3)),
    "Hour_LCL" = quantile(c(Hour_Tens_1, Hour_Tens_2, Hour_Tens_3),
      0.025,
      type = 8
    ),
    "Hour_UCL" = quantile(c(Hour_Tens_1, Hour_Tens_2, Hour_Tens_3),
      0.975,
      type = 8
    )
  )

Mean_Tensor_Coefs

Effs = rstan::summary(Acute_Full$fit)$summary %>%
  as.data.frame(.) %>%
  dplyr::select(n_eff) %>%
  filter(rownames(.) %in% c(
    "b_Intercept", "b_Treatment.L", "b_SexMale",
    "sds_sAmb.Temp_1", "sds_sTimeline_1",
    "sds_sAmb.TempTreatmentStress_1",
    "sds_sTimelineTreatmentStress_1",
    "sds_t2TimelineAmb.TempTreatmentStress_1",
    "sds_t2TimelineAmb.TempTreatmentStress_2",
    "sds_t2TimelineAmb.TempTreatmentStress_3",
    "sds_t2HourDirection_1",
    "sds_t2HourDirection_2",
    "sds_t2HourDirection_3",
    "sds_sTimelineBird.ID_1",
    "sds_sTimelineBird.IDTreatmentStress_1",
    "sd_Bird.ID__Intercept",
    "sd_Date.of.Photo__Intercept",
    "sd_Pen__Intercept",
    "b_sigma_Intercept",
    "b_sigma_Treatment.L"
  )) %>%
  mutate("Term" = gsub("_[[:digit:]]", "", rownames(.))) %>%
  group_by(Term) %>%
  summarise("Neff" = sum(n_eff))

Effs

# Validating model.
stan_ac(Acute_Full$fit)
mcmc_neff(neff_ratio(Acute_Full), size = 2)

# Neff/N ratios > 0.75. No clear autocorrelation. Checking Rhats.
mcmc_rhat(rhat(Acute_Full))

# Great - all near 1. Chains are well mixed. Checking residual distributions.

Acute_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full) %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(Acute_Full) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

Acute_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

# Reasonbly well fit and no clear pattern across fitted values.

Acute_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(Acute_Full) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Reasonable, but indicative of some clear outliers. Producing Cleveland dotplot then
# plotting residuals by predictor to assess where outliers are biologically meaningful or
# concerning.

P_Dat <- Acute_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  )

ggplot(P_Dat, aes(x = 1:nrow(P_Dat), y = z_residual)) +
  geom_point(colour = "black") +
  ylim(c(-4, 4))

# Note infinite values on each end. Some trailing on low end, perhaps due to
# low Ta. Assessing.

{
  p1 <- Acute_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- Acute_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- Acute_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- Acute_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Timeline, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p5 <- Acute_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- Acute_Full$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird.ID, y = z_residual, fill = Bird.ID)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Fair distributions. No apparent trends or heteroskedasticity across predictors, although peculiar ordering
# among individuals. Assessing model fit

Acute_Full$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full, scale = "response") %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = Maximum.Eye.Temp)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Again, good fit. Calculating loo-R2

loo_R2(Acute_Full)

# Quickly assessing descriptive statistics.

Descriptive = expand.grid(
  "Amb.Temp" = c(min(Acute_Full$data$Amb.Temp), max(Acute_Full$data$Amb.Temp)),
  "Timeline" = c(-1200, 1200),
  "Hour" = c(8, 12, 16),
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[1]
)

Global_Pred = posterior_predict(Acute_Full, newdata = Descriptive, re_formula = NA,
  nsamples = 1000, summary = FALSE)

Descriptive %>% mutate("Post" = colMeans(Global_Pred)) %>% 
  filter(Treatment == "Stress", 
  Amb.Temp %in% c(min(Acute_Full$data$Amb.Temp), max(Acute_Full$data$Amb.Temp))) %>% 
  mutate("Measure_Time" = ifelse(Timeline == -1200, "Start", "End")) %>% 
  dplyr::select(-Timeline) %>% 
  group_split(Measure_Time) %>% 
  reduce(left_join, by = c("Amb.Temp", "Hour", "Sex", "Direction", "Treatment", "Bird.ID", "Treat_Day")) %>% 
  rename("End_Temp" = Post.x, "Start_Temp" = Post.y) %>% 
  dplyr::select(-c(Measure_Time.x, Measure_Time.y)) %>% 
  mutate("Diff" = End_Temp - Start_Temp) %>% 
  group_by(Amb.Temp) %>% 
  summarise("Mean_Diff" = mean(Diff), "LCL" = quantile(Diff, 0.025, type = 8),
    "UCL" = quantile(Diff, 0.975, type = 8),
    "SD_Diff" = sd(Diff))

# Evaluating means and CIs for sigma parameters on original scale.

Acute_Full %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("C_Sigma" = exp(b_sigma_Intercept),
    "S_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L)) %>% 
    summarise("C_Est" = mean(C_Sigma), 
    "C_LCL" = quantile(C_Sigma, 0.025, type = 8),
    "C_UCL" = quantile(C_Sigma, 0.975, type = 8),
    "S_Est" = mean(S_Sigma), 
    "S_LCL" = quantile(S_Sigma, 0.025, type = 8),
    "S_UCL" = quantile(S_Sigma, 0.975, type = 8))

# Poorer fit than chronic model; quite surprising! Fit still reasonable (0.85).
# Comparing residual variance between treatment groups.

hyp <- c(
  "exp(sigma_Intercept + sigma_Treatment.L) < exp(sigma_Intercept)",
  "exp(sigma_Intercept + sigma_Treatment.L) - exp(sigma_Intercept) = 0"
)
plot(hypothesis(Acute_Full, hyp))
hypothesis(Acute_Full, hyp)

# Using custom sigma priors

sigma_priors <- Acute_Full %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("Control_Sigma" = exp(b_sigma_Intercept),
    "Stress_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L))

DF <- build_hdf(vars = list(sigma_priors$Control_Sigma, sigma_priors$Stress_Sigma),
  priors = list(rnorm(nrow(sigma_priors), 0, 0.25), rnorm(nrow(sigma_priors), 0, 0.25)),
  names = c("Control", "Stress"))

df_hyp = hypothesis_df("Stress = Control", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)

df_hyp = hypothesis_df("Stress < Control", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)

# Saving supplemental plot 

sFig_1 = plot(df_hyp, plot = F, theme = theme_get())[[1]]

sFig_1_Final = sFig_1 + theme_bw() + ylab("Density") + xlab(TeX('$\\sigma^2_{Control} - \\sigma^2_{Stress}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/SFig1.jpeg",
  sFig_1_Final,
  height = 6.0, width = 7.0, dpi = 800
)

# Stress-exposure treatment variance lower than controls; strong support (K = 72.47).
# Similar to chronic model, comparing repeatabilities between treatment groups using non-linear hypothesis test.

get_variables(Acute_Full)
Treat_hyp_test_Acute = c("(sds_sTimelineBird.IDTreatmentStress_1)^2/((sds_sTimelineBird.IDTreatmentStress_1)^2 + (exp(b_sigma_Treatment.L))^2) - (sds_sTimelineBird.ID_1)^2/((sds_sTimelineBird.ID_1)^2 + (exp(b_sigma_Intercept))^2) > 0.1")
THT_Acute_Out = hypothesis(Acute_Full, Treat_hyp_test_Acute, alpha = 0.05, class = NULL)
plot(THT_Acute_Out)
print(THT_Acute_Out, digits = 3)

# And again, using custom function to load in priors.

repeat_acute_custom <- Acute_Full %>%
    spread_draws(
      sds_sTimelineBird.ID_1,
      sds_sTimelineBird.IDTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
     ) %>% 
    mutate("Control_Repeat" = sds_sTimelineBird.ID_1^2/(exp(b_sigma_Intercept)^2 + sds_sTimelineBird.ID_1^2),
    "Stress_Repeat" = sds_sTimelineBird.IDTreatmentStress_1^2/(exp(b_sigma_Intercept + b_sigma_Treatment.L)^2 + sds_sTimelineBird.IDTreatmentStress_1^2))

repeat_acute_custom %>% summarise("R_Stress" = mean(Stress_Repeat), "S_LCL" = quantile(Stress_Repeat, 0.025, type = 8),
"S_UCL" = quantile(Stress_Repeat, 0.975, type = 8), "R_Control" = mean(Control_Repeat), "C_LCL" = quantile(Control_Repeat, 0.025, type = 8), "C_UCL" = quantile(Control_Repeat, 0.975, type = 8))

DF <- build_hdf(vars = list(repeat_acute_custom$Control_Repeat, repeat_acute_custom$Stress_Repeat),
  priors = list(rbeta(nrow(repeat_acute_custom), 1, 4), rbeta(nrow(repeat_acute_custom), 1, 4)),
  names = c("Control", "Stress"))

df_hyp = hypothesis_df("Stress = Control", DF, class = "b", alpha = 0.05)
plot(df_hyp)
print(df_hyp)

df_hyp = hypothesis_df("Stress > Control", DF, class = "b", alpha = 0.05)
plot(df_hyp)
print(df_hyp)

# Producing supplemental figure

sFig_3 = plot(df_hyp, plot = F, theme = theme_get())[[1]]

sFig_3_Final = sFig_3 + theme_bw() + ylab("Density") + xlab(TeX('$R_{Stress} - R_{Control}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/SFig3.jpeg",
  sFig_3_Final,
  height = 6.0, width = 7.0, dpi = 800
)

# Negligable support for a difference between groups (K = 2.06)
# Running null model

Acute_Null <- brm(bf(
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
file =
  "/home/joshk/Desktop/CC_Results/Final/Acute_OP_Comp_Null.Rds"
)

stan_ac(Acute_Null$fit)
mcmc_neff(neff_ratio(Acute_Null), size = 2)
mcmc_rhat(rhat(Acute_Null))

# Good, Neff/N > 0.67 for all and Rhats tight to 1.

conditional_smooths(Acute_Null)

# Trends remain similar to full model. Plotting residuals and exploring fit.

Acute_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

Acute_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null) %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(Acute_Null) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

Acute_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(Acute_Null) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Fairly uniform but wide. Again, producing Cleveland dotplot then plotting residuals by predictors.

P_Dat <- Acute_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null) %>%
  summarise(
    p_residual = mean(.prediction < Maximum.Eye.Temp),
    z_residual = qnorm(p_residual)
  )

ggplot(P_Dat, aes(x = 1:nrow(P_Dat), y = z_residual)) +
  geom_point(colour = "black") +
  ylim(c(-4, 4))

# Again, some trailing. Checking across predictors.

{
  p1 <- Acute_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- Acute_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- Acute_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- Acute_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = Timeline, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p5 <- Acute_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- Acute_Null$data %>%
    dplyr::select(
      Maximum.Eye.Temp, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null) %>%
    summarise(
      p_residual = mean(.prediction < Maximum.Eye.Temp),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird_Scramble, y = z_residual, fill = Bird_Scramble)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Homoskedastic across predictors. Model fit?

Acute_Null$data %>%
  dplyr::select(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null, scale = "response") %>%
  group_by(
    Maximum.Eye.Temp, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = Maximum.Eye.Temp)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Tight fit despite scrambling bird identities. Checking loo-R2.

loo_R2(Acute_Null)

# Highly similar to full model. Here, R2 = 0.85.
# Testing repeatability exclusively within stress-exposed groups.

get_variables(Acute_Full)
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
  data.frame(Acute_Null %>%
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
  scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
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

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_Repeat_New.jpeg",
  S_Rep_Acute,
  height = 6.0, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_Repeat_New.pdf",
       S_Rep_Acute,
       height = 6.0, width = 6.0, dpi = 800
)

# Very wide and uncertain with stress slopes only slightly exceeding those of control slopes.
# Comparing repeatability estimates by Bayesian t-test and Savage-Dickey non-linear hypothesis test.

BTT_Acute <- ttestBF(
  formula = T_Repeatability ~ Mod_Type, data = Rep_Treats_Acute,
  iterations = 10000, rscale = "ultrawide", posterior = FALSE
)

Post_Acute <- posterior(BTT_Acute, iterations = 1000)
plot(Post_Acute[, "mu"])
plot(acf(Post_Acute[, "mu"]))

# Chains a but noisy, but range is small.

summary(BTT_Acute)

# BF strikingly high. Testing using Savage-Dickey non-linear hypothesis test after
# qualitative comparison.

Rep_Treats_Acute %>%
  group_by(Mod_Type) %>%
  summarise(
    "Mean" = mean(T_Repeatability),
    "UCL" = quantile(T_Repeatability, 0.975, type = 8),
    "LCL" = quantile(T_Repeatability, 0.025, type = 8)
  )

# Some overlap, but not striking.

Wide_Acute_Rep <- Rep_Treats_Acute %>%
  mutate("Sample" = paste(.chain, .iteration, sep = "_")) %>%
  dplyr::select(Sample, Mod_Type, T_Repeatability) %>%
  spread(Mod_Type, T_Repeatability) %>%
  rename("True" = "True\nModel  ", "Null" = "Null\nModel  ") %>%
  as.data.frame(.)

h_test_Acute <- hypothesis(Wide_Acute_Rep, c("True > Null", "True - Null > 0.25"))
plot(h_test_Acute)
print(h_test_Acute, digits = 3)

# Using custom function to apply priors 

DF_Acute_Rep <- build_hdf(vars = list(Wide_Acute_Rep$True, Wide_Acute_Rep$Null), priors = list(rbeta(nrow(Wide_Acute_Rep), 1, 4), rbeta(nrow(Wide_Acute_Rep), 1, 4)), names = c("True", "Null"))

df_acute_hyp = hypothesis_df("True > Null", DF_Acute_Rep, class = "b", alpha = 0.05)

plot(df_acute_hyp)
print(df_acute_hyp)

# Much more logical results. Low support (K = 0.074 at 0.25 level). Note that variances between 
# groups are clearly not equal, however. Modeling with unequal variance and loose but mildly informative priors.

Robust_Comp_data_Acute <- Rep_Treats_Acute %>%
  dplyr::select(Mod_Type, T_Repeatability) %>%
  mutate(Model = gsub("\\n.*", "", Mod_Type)) %>%
  dplyr::select(-Mod_Type)

# Again, checking priors

robust_comp_acute <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student,
  data = Robust_Comp_data_Acute,
  sample_prior = "only",
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu"),
    set_prior("exponential(2)", class = "b", dpar = "sigma", lb = 0.001)
  )#,
  #file = "/home/joshk/Desktop/CC_Results/Final/Robust_Comp_Acute.Rds"
)

y_est <- posterior_predict(robust_comp_acute)
prior_pc = data.frame("Type" = c(rep("Prior", length(colMeans(y_est))),
  rep("True", length(colMeans(y_est)))), 
  "Y" = c(colMeans(y_est), Robust_Comp_data_Acute$T_Repeatability))

ggplot(prior_pc, aes(x = Y, fill = Type)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.4, end = 0.6)

# Priors are very wide, but removal of sigma prior is likely to help. Running model.

robust_comp_acute <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student,
  data = Robust_Comp_data_Acute,
  sample_prior = TRUE,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ),
  file = "/home/joshk/Desktop/CC_Results/Final/Robust_Comp_Acute.Rds"
)

pp_check(robust_comp_acute)
summary(robust_comp_acute)
mcmc_neff(neff_ratio(robust_comp_acute), size = 2)
mcmc_rhat(rhat(robust_comp_acute))

# Neff/N > 0.75 and Rhats near 1. Note that posterior is predicting a greater divergence between
# groups than is true, however, proximitey of peaks may render this seperation negligable.

# Comparing coefficients with Savage-Dickey

plot(hypothesis(robust_comp_acute, c("ModelTrue > ModelNull", "ModelTrue - ModelNull > 0.25")))
print(hypothesis(robust_comp_acute, c("ModelTrue > ModelNull", "ModelTrue - ModelNull > 0.25")))

# No difference according to conservative measure.
### Pulling out qualitative differences in means and SDs

tidy_MCMC_Acute <- tidyMCMC(robust_comp_acute, conf.int = TRUE, conf.level = 0.95, 
  estimate.method = "median", conf.method = "HPDinterval") %>% 
  mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

Post_Dif_Acute <- posterior_samples(robust_comp_acute) %>% 
  mutate_at(vars(contains("sigma")), funs(exp)) %>% 
  mutate(nu = log10(nu)) %>% 
  mutate(diff_means = b_ModelTrue - b_ModelNull,
         diff_sigma = b_sigma_ModelTrue - b_sigma_ModelNull) %>% 
  mutate(cohen_d = diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)/2),
         cles = dnorm(diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)), 0, 1))

Out_Acute <- tidyMCMC(Post_Dif_Acute, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")
Out_Acute

# Lastly, as permutation test 

set.seed(200) 

simu <- 1000
res <- numeric(simu) 

for (i in 1:simu) {
    perm <- sample(nrow(Robust_Comp_data_Acute))
    bdat <- transform(Robust_Comp_data_Acute, Model = Model[perm])
    res[i] <- mean(bdat[bdat$Model=="True", "T_Repeatability"])-
        mean(bdat[bdat$Model=="Null", "T_Repeatability"])
}

obs <- mean(Robust_Comp_data_Acute[Robust_Comp_data_Acute$Model == "True", "T_Repeatability"])-
    mean(Robust_Comp_data_Acute[Robust_Comp_data_Acute$Model == "Null", "T_Repeatability"])

ggplot(as.data.frame(res), aes(x = res)) + 
  geom_density(adjust = 3, fill = "mediumorchid", alpha = 0.5) + 
  geom_vline(xintercept = obs, colour = "black")

mean(abs(res)>=abs(obs))

# Again, clear differences between posterior means, however, sample size is likely a problem
# (that is, being inflated).

# Pulling out slopes and comparing between urban and rural birds.

Acute_Slopes <- Acute_Full %>%
  spread_draws(s_sTimelineBird.IDTreatmentStress_1[ID]) %>%
  rename("Coef" = s_sTimelineBird.IDTreatmentStress_1) %>%
  as.data.frame(.)

Acute_Birds <- data.frame(
  "Bird_ID" = rownames(ranef(Acute_Full)$Bird.ID),
  "ID" = seq(1, 19, 1)
)

Acute_Collapsed <- left_join(Acute_Slopes, Acute_Birds, by = c("ID")) %>%
  dplyr::select(-ID) %>%
  group_by(Bird_ID) %>%
  summarise("Slope" = mean(Coef)) %>%
  as.data.frame(.)

# Binding in locale and treatment order

U_R <- c()
Locale <- c()

for (i in 1:nrow(Acute_Collapsed)) {
  if (is.na(Acute_Collapsed$Bird_ID[i])) {
    Locale[i] <- NA
    U_R[i] <- NA
  } else if (Acute_Collapsed$Bird_ID[i] == "BdABd") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdAO1") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABlBl") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "AYY") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAO1") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "OOA") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdABl1") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAO2") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABlO") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABdBd") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAR") {
    Locale[i] <- "Cambridge"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "RAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABlR") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "ARO") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdABl2") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BlAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdAO2") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "AOR") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAY") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "AYBd") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else {
    Locale[i] <- NA
    U_R[i] <- NA
  }
}
Acute_Collapsed$U_R <- factor(U_R)
Acute_Collapsed$Locale <- factor(Locale)

A_Orders <- Third_Bound %>%
  dplyr::select("Bird_ID" = Bird.ID, Pen) %>%
  na.omit(.) %>%
  distinct(.) %>%
  mutate("Order" = ifelse(Pen == "NE" | Pen == "SW", "SR", "RS")) %>%
  dplyr::select(-Pen)
Acute_Collapsed <- left_join(Acute_Collapsed, A_Orders, by = c("Bird_ID")) %>%
  mutate(Order = factor(Order), U_R = factor(U_R))

# Running formal comparison by "ANOVA" and non-linear hypothesis test (similar to above for
# acute responses).

Acute_Collapsed %>%
  group_by(U_R) %>%
  summarise(
    Mean = mean(Slope),
    LCL = quantile(Slope, 0.025, type = 8),
    UCL = quantile(Slope, 0.975, type = 8)
  )

# Note high degree of overlap.

Acute_aovBF <- anovaBF(Slope ~ Order * U_R,
  data = Acute_Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(Acute_aovBF)

# Again, minimal evidence for order effect. Dropping and re-running.

Acute_aovBF <- anovaBF(Slope ~ U_R,
  data = Acute_Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(BayesFactor::posterior(Acute_aovBF, iterations = 10000)[, "mu"])
plot(acf(BayesFactor::posterior(Acute_aovBF, iterations = 10000)[, "mu"]))

# Chain is stable and peak looks clear. Note that trace is hovering surprisingly high. Summarising.

summary(Acute_aovBF)
summary(BayesFactor::posterior(Acute_aovBF, iterations = 10000)[, "mu"])

# Bayes factor is quite low (0.249). Compiling model for non-linear hypothesis test.

Slope_Contrast_Acute <- brm(Slope ~ 0 + U_R + (1 | Locale),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("gamma(1,1)", class = "sd")
  ), data = Acute_Collapsed, cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  sample_prior = TRUE
)

Acute_Hypoth <- hypothesis(Slope_Contrast_Acute, c(
  "U_RRural = U_RUrban"
))
plot(Acute_Hypoth)
1 / Acute_Hypoth$hypothesis$Evid.Ratio

# Results are weaker by non-linear hypothesis test and appear quite poor by plot.
# Plotting global trends.

UR_Plot_Acute <- Acute_Collapsed %>%
  rename("Ecotype" = U_R) %>%
  group_by(Ecotype) %>%
  summarise(
    "Coefficient" = mean(Slope),
    "lcl" = quantile(Slope, 0.025, type = 8),
    "ucl" = quantile(Slope, 0.975, type = 8)
  ) %>%
  ggplot(aes(x = Ecotype, y = Coefficient, fill = Ecotype)) +
  geom_errorbar(aes(x = Ecotype, ymin = lcl, ymax = ucl),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 6, pch = 21, colour = "black") +
  #scale_fill_viridis_d(begin = 0.2, end = 0.5) +
  scale_fill_manual(values = c("wheat4","cornsilk")) + 
  theme_bw() +
  my.theme +
  theme(legend.position = "none") +
  ylab("Reaction Norm Slope")

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Acute_New.jpeg",
  UR_Plot_Acute,
  height = 5.67, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Acute_New.pdf",
       UR_Plot_Acute,
       height = 5.67, width = 6.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Plotting global slopes.

Simp_Dat_Acute <- expand.grid(
  "Amb.Temp" = c(5, 20, 35),
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = c(12),
  "Date.of.Photo" = na.omit(unique(Acute_Full$data$Date.of.Photo))[30],
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  "Pen" = na.omit(unique(Acute_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[70]
)

facet_names <- c(`5` = "5°C", `20` = "20°C", `35` = "35°C")

Simp_Curves_Filled <- Simp_Dat_Acute %>%
  add_predicted_draws(Acute_Full, n = 1000, scale = "response") %>%
  mutate("Group_ID" = paste(Bird.ID, Treatment, sep = "_")) %>%
  group_by(Group_ID, Bird.ID, Treatment, Amb.Temp, Timeline) %>%
  summarise("Pred" = mean(.prediction))

SC_F_Plot <- ggplot(Simp_Curves_Filled, aes(x = Timeline, y = Pred, colour = Treatment)) +
  geom_line(aes(group = Group_ID), alpha = 0.7) +
  annotate("rect",
    xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
    fill = "grey30", alpha = 0.8
  ) +
  facet_wrap(~Amb.Temp, labeller = as_labeller(facet_names), nrow = 3) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Time Post Stressor (s)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  theme(legend.position = "bottom", strip.text.x = element_text(size = 12, family = "Noto Sans")) +
  my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Stacked_Acute_Comp.jpeg", SC_F_Plot,
  height = 13, width = 5.5, dpi = 800
)

# Surface plot

Surface_Plot <- expand.grid(
  "Maximum.Eye.Temp" = NA,
  "Date.of.Photo" = unique(Acute_Full$data$Date.of.Photo)[30],
  "Pen" = unique(Acute_Full$data$Pen)[c(1, 2)],
  "Amb.Temp" = seq(2.5, 37.5, by = 1),
  "Treatment" = unique(Acute_Full$data$Treatment)[2],
  "Bird.ID" = unique(Acute_Full$data$Bird.ID)[c(3, 15)],
  "Hour" = 12,
  "Timeline" = seq(-1200, 3600, 180),
  "Sex" = unique(Acute_Full$data$Sex),
  "Direction" = unique(Acute_Full$data$Direction),
  "Treat_Day" = unique(Acute_Full$data$Treat_Day)[40]
) %>%
  as.data.frame(.) %>%
  add_predicted_draws(Acute_Full, n = 500, scale = "response") %>%
  group_by(Amb.Temp, Timeline) %>%
  summarise("yvar" = mean(.prediction)) %>%
  as.data.frame(.)

# Pulling out legend

Leg_Lower <- Surface_Plot %>% ggplot(aes(x = Amb.Temp, y = Timeline, z = yvar)) +
  geom_raster(aes(fill = yvar)) +
  scale_fill_gradient(
    low = viridis::viridis(n = 10)[4],
    high = "ivory", breaks = c(23.5, 30.5, 37.5), 
    name = "Eye Region\nTemperature\n(°C)"
  ) + 
  my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Final/Lower_Leg_New.jpeg", 
Leg_Lower, height = 6, width = 8.5, dpi = 800)

extractLegend <- function(xplot) {
  grobs <- ggplot_gtable(ggplot_build(xplot))
  g_title <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[g_title]]
}

ScaleLeg <- extractLegend(Leg_Lower)

min(Surface_Plot$yvar)
max(Surface_Plot$yvar)

OP_min <- 23.0
OP_max <- 38.0

OP_mybreaks <- seq(OP_min, OP_max, length.out = 15)
OP_colours <- function(x) {
  colours <- colorRampPalette(c(viridis::viridis(n = 10)[4], "ivory"))(14)
  colours[1:x]
}

OP_breaklabel <- function(x) {
  labels <- paste0(OP_mybreaks[1:14], "-", OP_mybreaks[2:15])
  labels[1:x]
}

Surface_No_Legend <- Surface_Plot %>%
  ggplot(aes(x = Amb.Temp, y = Timeline, z = yvar)) +
  geom_contour_filled(breaks = OP_mybreaks, show.legend = TRUE) +
  xlab("Ambient Temperature (°C)") +
  ylab("Time Post Stress Exposure (s)") +
  my.theme +
  scale_fill_manual(palette = OP_colours, values = OP_breaklabel(14), name = "Eye Region\nTemperature (°C)", drop = TRUE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank()
  )

replace_Grob <- ggplot_gtable(ggplot_build(Surface_No_Legend))
rep_new <- which(sapply(replace_Grob$grobs, function(x) x$name) == "guide-box")
replace_Grob$grobs[[rep_new]] <- ScaleLeg

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_New.jpeg", replace_Grob,
       height = 7, width = 8.5, dpi = 800
)

# And with timeline on the x axis

Surface_No_Legend_Rotated <- Surface_Plot %>%
  ggplot(aes(x = Timeline, y = Amb.Temp, z = yvar)) +
  geom_contour_filled(breaks = OP_mybreaks, show.legend = TRUE) +
  xlab("Time Post Stress Exposure (s)") +
  ylab("Ambient Temperature (°C)") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(trans = ~., name = "", breaks = c(14,22,30), labels = c(" -","TNZ"," -"))) + 
  coord_cartesian(ylim = c(2.5, 37.5), xlim = c(-1200, 3600), clip="off") +
  annotate("segment", x=3600, y=23, xend=3600, yend=29,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=3600, y=21, xend=3600, yend=15,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) + 
  scale_fill_manual(palette = OP_colours, values = OP_breaklabel(14), name = "Eye Region\nTemperature (°C)", drop = TRUE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_text(hjust = -0.6, size = 16),
    axis.title.y.left = element_text(vjust = 5),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + my.theme

cow_legend <- cowplot::get_legend(Leg_Lower)
cow_grid <- cowplot::plot_grid(plotlist = list(Surface_No_Legend_Rotated), ncol = 1)
to_save = cowplot::plot_grid(cow_grid, cow_legend, ncol = 2, rel_widths = c(12,1.5), align = "hv")

to_save

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed.jpeg", to_save,
       height = 7, width = 9, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed.pdf", to_save,
       height = 7, width = 9, dpi = 800
)

# With grey rectangle over stress exposure period.

Surface_No_Legend_Rotated <- Surface_Plot %>%
  ggplot(aes(x = Timeline, y = Amb.Temp, z = yvar)) +
  geom_contour_filled(breaks = OP_mybreaks, show.legend = TRUE) +
  xlab("Time Post Stress Exposure (s)") +
  ylab("Ambient Temperature (°C)") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(trans = ~., name = "", breaks = c(14,22,30), labels = c(" -","TNZ"," -"))) + 
  coord_cartesian(ylim = c(2.5, 37.5), xlim = c(-1200, 3600), clip="off") +
  annotate("segment", x=3600, y=23, xend=3600, yend=29,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=3600, y=21, xend=3600, yend=15,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) + 
  annotate("rect",
           xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
           fill = "grey30", alpha = 0.8
  ) + 
  scale_fill_manual(palette = OP_colours, values = OP_breaklabel(14), name = "Eye Region\nTemperature (°C)", drop = TRUE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_text(hjust = -0.6, size = 16),
    axis.title.y.left = element_text(vjust = 5),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + my.theme

cow_legend <- cowplot::get_legend(Leg_Lower)
cow_grid <- cowplot::plot_grid(plotlist = list(Surface_No_Legend_Rotated), ncol = 1)
to_save = cowplot::plot_grid(cow_grid, cow_legend, ncol = 2, rel_widths = c(12,1.5), align = "hv")

to_save

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed_Rectangle.jpeg", to_save,
       height = 7, width = 9, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed_Rectangle.pdf", to_save,
       height = 7, width = 9, dpi = 800
)

# Finally, re-running formal comparisons between repeatabilities of treatment groups. Note that 
# the below comparisons are directed towards differences in means and not distributions. Such comparisons
# are likely to be quite liberal.

Rep_Chronic_HT <- Rep_Treats %>%
  mutate(Model = gsub("\n.*", "", Rep_Treats$Mod_Type)) %>%
  rename(
    "Control" = C_Repeatability,
    "Stress" = T_Repeatability
  ) %>%
  filter(Model == "True") %>%
  dplyr::select(.chain, .iteration, Control, Stress) %>%
  gather(Treatment, Repeat, Control:Stress, factor_key = TRUE)

Rep_Acute_HT <- Rep_Treats_Acute %>%
  mutate(Model = gsub("\n.*", "", Rep_Treats$Mod_Type)) %>%
  rename(
    "Control" = C_Repeatability,
    "Stress" = T_Repeatability
  ) %>%
  filter(Model == "True") %>%
  dplyr::select(.chain, .iteration, Control, Stress) %>%
  gather(Treatment, Repeat, Control:Stress, factor_key = TRUE)

{
Treat_Comp_Robust_Acute <- brm(
  bf(Repeat ~ 0 + Treatment, sigma ~ 0 + Treatment),
  family = student,
  data = Rep_Acute_HT,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101, 
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ), 
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  file = "/home/joshk/Desktop/CC_Results/Final/Treat_Comp_Robust_Acute.Rds"
) 
  
Treat_Comp_Robust_Chronic <- brm(
  bf(Repeat ~ 0 + Treatment, sigma ~ 0 + Treatment),
  family = student,
  data = Rep_Chronic_HT,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ),
  file = "/home/joshk/Desktop/CC_Results/Final/Treat_Compare_Robust_Chronic.Rds"
)
}
# Checking outcomes

plot(Treat_Comp_Robust_Acute)
plot(Treat_Comp_Robust_Chronic)
pp_check(Treat_Comp_Robust_Acute, nsamples = 500)
pp_check(Treat_Comp_Robust_Chronic, nsamples = 500)

# Strikingly high nu parameters, but posterior samples fit well. Again, these estimates should be liberal.

stan_ac(Treat_Comp_Robust_Acute$fit)
mcmc_neff(neff_ratio(Treat_Comp_Robust_Acute), size = 2)
mcmc_rhat(rhat(Treat_Comp_Robust_Acute))

stan_ac(Treat_Comp_Robust_Chronic$fit)
mcmc_neff(neff_ratio(Treat_Comp_Robust_Chronic), size = 2)
mcmc_rhat(rhat(Treat_Comp_Robust_Chronic))

# Fair chain convergence. Plotting fit against raw.

Pred_A <- Rep_Acute_HT %>%
  add_predicted_draws(Treat_Comp_Robust_Acute, class = "response") %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(.prediction),
    "UCL" = quantile(.prediction, 0.975, type = 8),
    "LCL" = quantile(.prediction, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

Raw_A <- Rep_Acute_HT %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(Repeat),
    "UCL" = quantile(Repeat, 0.975, type = 8),
    "LCL" = quantile(Repeat, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Raw")

grid.arrange(Pred_A, Raw_A, nrow = 2)

# Acceptably comparable. Checking whether heteroskedasticity was captured.

Rep_Acute_HT %>%
  add_residual_draws(Treat_Comp_Robust_Acute, class = "response") %>%
  ggplot(aes(x = Treatment, y = .residual, fill = Treatment)) +
  geom_boxplot(colour = "black") +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

# Seems loosely corrected and acceptable for now. Chronic model?

Pred_C <- Rep_Chronic_HT %>%
  add_predicted_draws(Treat_Comp_Robust_Chronic, class = "response") %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(.prediction),
    "UCL" = quantile(.prediction, 0.975, type = 8),
    "LCL" = quantile(.prediction, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

Raw_C <- Rep_Chronic_HT %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(Repeat),
    "UCL" = quantile(Repeat, 0.975, type = 8),
    "LCL" = quantile(Repeat, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Raw")

grid.arrange(Pred_C, Raw_C, nrow = 2)

# Great. Again, checking heteroskedasticity.

Rep_Chronic_HT %>%
  add_residual_draws(Treat_Comp_Robust_Chronic, class = "response") %>%
  ggplot(aes(x = Treatment, y = .residual, fill = Treatment)) +
  geom_boxplot(colour = "black") +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

# No concern.
# Looking at model summaries

summary(Treat_Comp_Robust_Acute)
summary(Treat_Comp_Robust_Chronic)

hypothesis(Treat_Comp_Robust_Acute, "TreatmentStress - TreatmentControl > 0.1", alpha = 0.05)
hypothesis(Treat_Comp_Robust_Chronic, "TreatmentStress - TreatmentControl > 0.1", alpha = 0.05)

# KEEP 

#robust_comp_acute$model = gsub("real Intercept_sigma;","real<lower=0.001> Intercept_sigma;",robust_comp_acute$model)
#new_standata = standata(robust_comp_acute)
#robust_comp_acute$fit = stan(model_code = robust_comp_acute$model, data = new_standata,
#                             chains = 4, iter = 10000, warmup = 1000, thin = 10,
#                             seed = 101, control = list(adapt_delta = 0.98, max_treedepth = 14))

## Test raw data point plot

Simp_Dat_Acute <- expand.grid(
  "Amb.Temp" = c(3,5,7,9,11,13),
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = 12,
  #"Date.of.Photo" = na.omit(unique(Acute_Full$data$Date.of.Photo))[1],
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  #"Pen" = na.omit(unique(Acute_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[1]
)

Pred = posterior_predict(Acute_Full, newdata = Simp_Dat_Acute, nsamples = 1000, re_formula = ~ (1|Bird.ID))
Simp_Dat_Acute = Simp_Dat_Acute %>% mutate("Pred" = colMeans(Pred), "Temp_Group" = "Low") 

Simp_Dat_Acute %>%
  mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_")) %>%
  group_by(Bird.ID, Timeline, Treatment, Treat_Bird) %>% 
  summarise("Maximum.Eye.Temp" = mean(Pred)) %>% 
  ggplot(aes(x = Timeline, y = Maximum.Eye.Temp, colour = Treatment, fill = Treatment)) + 
  geom_line(size = 1, aes(group = Treat_Bird)) + 
  theme_bw() + 
  scale_fill_viridis_d(begin = 0.3, end = 0.6) + 
  scale_colour_viridis_d(begin = 0.3, end = 0.6)

Plot_Dat = Third_Bound %>% drop_na(Maximum.Eye.Temp) %>% 
mutate("Treat_Bird" = paste(Treatment, Bird_Scramble, sep = "_")) %>% 
filter(Amb.Temp < 14)

Base_Plot + stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 60, mapping = aes(x = Timeline, y = Maximum.Eye.Temp,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2)

## Good! Expanding to include data points for other temperature regions.   

Mid_Temp <- expand.grid(
  "Amb.Temp" = c(14,16,18,20,22,24,26,28,30),
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = 12,
  #"Date.of.Photo" = na.omit(unique(Acute_Full$data$Date.of.Photo))[1],
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  #"Pen" = na.omit(unique(Acute_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[1]
)
Mid_Pred = posterior_predict(Acute_Full, newdata = Mid_Temp, nsamples = 1000, re_formula =  ~ (1|Bird.ID))
Mid_Temp = Mid_Temp %>% mutate("Pred" = colMeans(Mid_Pred), "Temp_Group" = "Mid")

High_Temp <- expand.grid(
  "Amb.Temp" = c(31,33,35,37,38.5),
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = 12,
  #"Date.of.Photo" = na.omit(unique(Acute_Full$data$Date.of.Photo))[1],
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  #"Pen" = na.omit(unique(Acute_Full$data$Pen))[c(1,2)],
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[1]
)

#High_Temp = left_join(High_Temp, as.data.frame(Acute_Full$data %>% 
#  dplyr::select(Bird.ID, Pen) %>% 
#  distinct() %>% 
#  drop_na(Bird.ID)), by = c("Bird.ID"))

High_Pred = posterior_predict(Acute_Full, newdata = High_Temp, nsamples = 500, re_formula = ~ (1|Bird.ID))
High_Temp = High_Temp %>% mutate("Pred" = colMeans(High_Pred), "Temp_Group" = "High")

Grouped = rbind(Simp_Dat_Acute, Mid_Temp, High_Temp) %>% 
  mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_")) %>%
  group_by(Temp_Group, Bird.ID, Timeline, Treatment, Treat_Bird) %>% 
  summarise("Maximum.Eye.Temp" = mean(Pred))

write.csv(Grouped, "/home/joshk/git_repositories/BCCH_IndVar/AcutePlotDat.csv", row.names = FALSE)
Grouped = read.csv("/home/joshk/git_repositories/BCCH_IndVar/AcutePlotDat.csv")
# Global plot

Plot_Dat = Third_Bound %>% drop_na(Maximum.Eye.Temp) %>% 
mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_"),
"Temp_Group" = ifelse(Amb.Temp < 14, "Low", ifelse(Amb.Temp > 30, "High", "Mid")))
Plot_Dat$Temp_Group = factor(Plot_Dat$Temp_Group, levels = c("Low", "Mid", "High"))

Temp_Groups = c(`Low` = "< TNZ", `Mid` = "TNZ", `High` = "> TNZ")
Grouped = Grouped %>% mutate(Temp_Group = factor(Temp_Group, levels = c("Low", "Mid", "High")))

ggplot(Grouped, aes(x = Timeline, y = Maximum.Eye.Temp, colour = Treatment, fill = Treatment)) + 
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 60, mapping = aes(x = Timeline, y = Maximum.Eye.Temp,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.8) + 
  geom_line(size = 1, aes(group = Treat_Bird)) + 
  theme_bw() + 
  scale_fill_viridis_d(begin = 0.3, end = 0.6) + 
  scale_colour_viridis_d(begin = 0.3, end = 0.6) + 
  facet_wrap(~Temp_Group, labeller = as_labeller(Temp_Groups))

Stacked_Acute_Dots = ggplot(Grouped, aes(x = Timeline, y = Maximum.Eye.Temp, colour = Treatment, fill = Treatment)) + 
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 180, mapping = aes(x = Timeline, y = Maximum.Eye.Temp,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(size = 1, aes(group = Treat_Bird), alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  facet_wrap(~Temp_Group, ncol = 1, labeller = as_labeller(Temp_Groups)) + 
  ylab("Maximum Eye Temperature (°C)") + 
  xlab("Time Post Stess Exposure (s)") + 
    annotate("rect",
    xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
    fill = "grey30", alpha = 0.8
  ) + my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Stacked_Acute_Dots.jpeg", Stacked_Acute_Dots,
  height = 13, width = 5.5, dpi = 800
)

# Wide format

Stacked_Acute_Dots = ggplot(Grouped, aes(x = Timeline, y = Maximum.Eye.Temp, colour = Treatment, fill = Treatment)) + 
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
                   binwidth = 180, mapping = aes(x = Timeline, y = Maximum.Eye.Temp,
                                                 group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(size = 1, aes(group = Treat_Bird), alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  facet_wrap(~Temp_Group, ncol = 3, labeller = as_labeller(Temp_Groups)) + 
  ylab("Maximum Eye Temperature (°C)") + 
  xlab("Time Post Stess Exposure (s)") + 
  annotate("rect",
           xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
           fill = "grey30", alpha = 0.8
  ) + my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Wide_Acute_Dots.jpeg", Stacked_Acute_Dots,
       height = 4.25, width = 10, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Wide_Acute_Dots.pdf", Stacked_Acute_Dots,
       height = 4.25, width = 10, dpi = 800
)

# Plotting chronic responses with dots.

Descriptive_Chronic = expand.grid(
  "s.Amb.Temp" = c(min(CC_Full$data$s.Amb.Temp), max(CC_Full$data$s.Amb.Temp)),
  "Hour" = c(8, 12, 16),
  "Sex" = na.omit(unique(CC_Full$data$Sex)),
  "Pen" = na.omit(unique(CC_Full$data$Pen)),
  "Date.of.Photo" = na.omit(unique(CC_Full$data$Date.of.Photo))[1],
  "Direction" = na.omit(unique(CC_Full$data$Direction)),
  "Treatment" = na.omit(unique(CC_Full$data$Treatment)),
  "Bird.ID" = na.omit(unique(CC_Full$data$Bird.ID)),
  "Treat_Day" = na.omit(unique(CC_Full$data$Treat_Day))[1]
)

Descriptive_Chronic_Pred = posterior_predict(CC_Full, newdata = Descriptive_Chronic, re_formula = NA,
  nsamples = 1000, summary = FALSE)
Descriptive_Chronic = Descriptive_Chronic %>% mutate("Pred" = colMeans(Descriptive_Chronic_Pred))

Chronic_Grouped = Descriptive_Chronic %>% 
  mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_"),
  "Amb.Temp" = (s.Amb.Temp * (2 * sd(Third_Bound$Amb.Temp, na.rm = T))) +
      mean(Third_Bound$Amb.Temp, na.rm = T)) %>% 
  group_by(Treatment, Bird.ID, Amb.Temp, Treat_Bird) %>% 
  summarise("Maximum.Eye.Temp" = mean(Pred))

Chronic_Curves_Dots <- 
  ggplot(Chronic_Grouped, aes(x = Amb.Temp, y = Maximum.Eye.Temp, colour = Treatment)) +
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 1, mapping = aes(x = Amb.Temp, y = Maximum.Eye.Temp,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_Dots.jpeg",
  Chronic_Curves_Dots,
  height = 6.0, width = 7.0, dpi = 800
)

Simp_Dat_Chronic <- expand.grid(
  "s.Amb.Temp" = seq(min(CC_Full$data$s.Amb.Temp, na.rm = TRUE),
    max(CC_Full$data$s.Amb.Temp, na.rm = TRUE),
    length.out = 25
  ),
  "Hour" = c(12),
  "Date.of.Photo" = na.omit(unique(CC_Full$data$Date.of.Photo))[30],
  "Treatment" = na.omit(unique(CC_Full$data$Treatment)),
  "Pen" = na.omit(unique(CC_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(CC_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(CC_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(CC_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(CC_Full$data$Treat_Day))[c(14, 70)]
)

Descriptive_Chronic_Pred = posterior_predict(CC_Full, newdata = Simp_Dat_Chronic,
  nsamples = 1000, summary = FALSE)

Simp_Dat_Chronic = Simp_Dat_Chronic %>% 
  mutate("Pred" = colMeans(Descriptive_Chronic_Pred),
  "Treat_Bird" = paste(Bird.ID, Treatment, sep = "_")
  ) %>%
  group_by(Treat_Bird, Bird.ID, Treatment, s.Amb.Temp) %>%
  summarise("Maximum.Eye.Temp" = mean(Pred))

Chronic_Curves_Dots <- 
  ggplot(Simp_Dat_Chronic, aes(x = s.Amb.Temp, y = Maximum.Eye.Temp, colour = Treatment)) +
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 0.01, mapping = aes(x = s.Amb.Temp, y = Maximum.Eye.Temp,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  #annotate(geom = "rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf,
  #         fill = "grey70", alpha = 0.5, colour = "black") + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme  

Chronic_Curves_Dots

# Averaged dots

Chronic_Curves_Bare <- 
  ggplot(Simp_Dat_Chronic, aes(x = s.Amb.Temp, y = Maximum.Eye.Temp, colour = Treatment)) +
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
                   binwidth = 0.416667, mapping = aes(x = s.Amb.Temp, y = Maximum.Eye.Temp,
                                                  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  annotate(geom = "rect", xmin = -0.70139, xmax = 0.40972, ymin = -Inf, ymax = Inf,
           fill = "grey30", alpha = 0.5, colour = "black") + 
  annotate(geom = "text", x = -0.14583, y = 24, label = "TNZ", size = 8, family = "Noto Sans") + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-1.5, 1, length.out = 7), 
                     labels = seq(2.5, 38.5, by = 6)) + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_TNZ_Lab.jpeg",
       Chronic_Curves_Bare,
       height = 6.0, width = 7.0, dpi = 800
)

# No dots 

Chronic_Curves_Nodots <- 
  ggplot(Simp_Dat_Chronic, aes(x = s.Amb.Temp, y = Maximum.Eye.Temp, colour = Treatment)) +
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  annotate(geom = "rect", xmin = -0.70139, xmax = 0.40972, ymin = -Inf, ymax = Inf,
           fill = "grey30", alpha = 0.5, colour = "black") + 
  annotate(geom = "text", x = -0.14583, y = 27, label = "TNZ", size = 8, family = "Noto Sans") + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-1.5, 1, length.out = 7), 
                     labels = seq(2.5, 38.5, by = 6)) + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_TNZ_Lab_NoDots.jpeg",
       Chronic_Curves_Nodots,
       height = 6.0, width = 7.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_TNZ_Lab_NoDots.pdf",
       Chronic_Curves_Nodots,
       height = 6.0, width = 7.0, dpi = 800
)

# Drawing distributions for methods figure

require(latex2exp)

Sample = rnorm(1e+4, mean = 0, sd = 1)

Rep_Curve = ggplot(as.data.frame(Sample), aes(x = Sample)) +
  geom_density(adjust = 3, size = 3, colour = "firebrick4", alpha = 0.5) + 
  annotate("text", x = 0.03, y = 0.1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram1.jpg",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
)
 
ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram1.pdf",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

# With colour 

Rep_Curve = ggplot(as.data.frame(Sample), aes(x = Sample)) +
  geom_density(adjust = 3, size = 1, colour = "black", fill = "firebrick4", alpha = 0.5) + 
  annotate("text", x = 0.03, y = 0.1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour1.jpg",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour1.pdf",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

# Adding null curve

Sample = rnorm(1e+4, mean = 1, sd = 0.6)
Beta_Sample = rbeta(1e+4, 1, 6)

Combined_Curve = ggplot(as.data.frame(Sample), aes(x = Sample + 0.5)) +
  geom_density(adjust = 3, size = 3, colour = "firebrick4", alpha = 0.5) + 
  geom_density(data = as.data.frame(Beta_Sample), aes(x = Beta_Sample),
               colour = "black", alpha = 0.5,
               size = 3, adjust = 3) + 
  xlim(c(0,4)) + 
  annotate("text", x = 1.3, y = 2, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{N}} \\cdot (\\mu_{1_{N}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("text", x = 3, y = 1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("segment", x = 1, y = 1.8, xend = 0.5, yend = 1, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) + 
  annotate("segment", x = 2.4, y = 0.8, xend = 2.2, yend = 0.5, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) +
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram2.jpg",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram2.pdf",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

# Again, with colour 

Combined_Curve = ggplot(as.data.frame(Sample), aes(x = Sample + 0.5)) +
  geom_density(adjust = 3, size = 1, colour = "black", fill = "firebrick4", alpha = 0.5) + 
  geom_density(data = as.data.frame(Beta_Sample), aes(x = Beta_Sample),
               fill = "black", colour = "black", alpha = 0.5,
               size = 1, adjust = 3) + 
  xlim(c(0,4)) + 
  annotate("text", x = 1.3, y = 2, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{N}} \\cdot (\\mu_{1_{N}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("text", x = 3, y = 1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("segment", x = 1, y = 1.8, xend = 0.5, yend = 1, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) + 
  annotate("segment", x = 2.4, y = 0.8, xend = 2.2, yend = 0.5, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) +
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour2.jpg",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour2.pdf",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

##########################################################################
########### Assessing patterns of heat transfer across birds #############
##########################################################################

# Calculating dry heat transfer across both eyes 

{
require(Thermimage)
require(bigleaf)

dim = 0.011
Area = ((1.1/2)*(1.0/2)*pi)*0.0001
Area = Area*2 # For two eyes 

qrad_Clean = c()
for (i in 1:nrow(Third_Bound)) {
  qrad_Clean[i] = Area * (5.67 * 10^-8) * 0.95 * 0.95 *
    ((Third_Bound$Maximum.Eye.Temp[i] + 273.15)^4 - (Third_Bound$Maximum.Eye.Temp[i] + 273.15)^4)
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
  qconv_clean[i] = Area * (Hc_clean[i]) * (Third_Bound$Maximum.Eye.Temp[i] - Third_Bound$Amb.Temp[i])
}

qtot_clean = c()

for (i in 1:nrow(Third_Bound)) {
  qtot_clean[i] = qconv_clean[i] + qrad_Clean[i]
}

Third_Bound$qtot = qtot_clean
}

# Converting W to mW and extending to two eye.

Third_Bound$mW = Third_Bound$qtot*1000

# Assigning chronic prior

prior_chronic <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  # set_prior("gamma(2, 0.5)", class = "b", coef = "ss.Amb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(50, 2)", class = "Intercept")
)

# Running chronic models

CC_Full_HT <- brm(bf(
  mW ~ Treatment + Sex +
    s(s.Amb.Temp, k = 4, bs = "cr") +
    s(s.Amb.Temp, by = Treatment, k = 4, bs = "cr", m = 1) +
    t2(Hour, Direction, k = 4, bs = c("tp", "re")) +
    s(s.Amb.Temp, Bird.ID, bs = "re") +
    s(s.Amb.Temp, by = Treatment, Bird.ID, bs = "re", m = 1) +
    (1 | Pen) + (1 | Date.of.Photo) + (1 | Bird.ID),
  sigma ~ Treatment
) +
  cor_ar(~ 1 | Treat_Day),
data = Third_Bound, cores = 4, chains = 4,
seed = 100, refresh = 0, family = "gaussian",
iter = 10000, warmup = 1000, thin = 10,
control = list(adapt_delta = 0.98, max_treedepth = 14),
prior = prior_chronic,
file =
  "/home/joshk/Desktop/Chronic_HT.Rds"
)

# Calculating descriptive statistics

Descriptive_Chronic = expand.grid(
  "s.Amb.Temp" = c(min(CC_Full_HT$data$s.Amb.Temp), max(CC_Full_HT$data$s.Amb.Temp)),
  "Hour" = c(8, 12, 16),
  "Sex" = na.omit(unique(CC_Full_HT$data$Sex)),
  "Direction" = na.omit(unique(CC_Full_HT$data$Direction)),
  "Treatment" = na.omit(unique(CC_Full_HT$data$Treatment)),
  "Bird.ID" = na.omit(unique(CC_Full_HT$data$Bird.ID)),
  "Treat_Day" = na.omit(unique(CC_Full_HT$data$Treat_Day))[1]
)

Global_Pred_Chronic = posterior_predict(CC_Full_HT, newdata = Descriptive_Chronic, re_formula = NA, nsamples = 1000, summary = FALSE)

Descriptive_Chronic %>% mutate("Post" = colMeans(Global_Pred_Chronic)) %>% 
  group_split(Treatment) %>% 
  reduce(left_join, by = c("s.Amb.Temp", "Hour", "Sex", "Direction", "Bird.ID", "Treat_Day")) %>% 
  rename("Control_Temp" = Post.x, "Stress_Temp" = Post.y) %>% 
  dplyr::select(-c(Treatment.x, Treatment.y)) %>% 
  mutate("Diff" = Control_Temp - Stress_Temp) %>% 
  group_by(s.Amb.Temp) %>% 
  summarise("Mean_Diff" = mean(Diff), "LCL" = quantile(Diff, 0.025, type = 8),
    "UCL" = quantile(Diff, 0.975, type = 8),
    "SD_Diff" = sd(Diff))

# Pulling effective sample sizes from summary

Effs_Chronic = rstan::summary(CC_Full_HT$fit)$summary %>%
  as.data.frame(.) %>%
  dplyr::select(n_eff) %>%
  filter(rownames(.) %in% c(
    "b_Intercept", "b_Treatment.L", "b_SexMale",
    "sds_ss.Amb.Temp_1",
    "sds_ss.Amb.TempTreatmentStress_1",
    "sds_t2HourDirection_1",
    "sds_t2HourDirection_2",
    "sds_ss.Amb.TempBird.ID_1",
    "sds_ss.Amb.TempBird.IDTreatmentStress_1",
    "sd_Bird.ID__Intercept",
    "sd_Date.of.Photo__Intercept",
    "sd_Pen__Intercept",
    "b_sigma_Intercept",
    "b_sigma_Treatment.L"
  )) %>%
  mutate("Term" = gsub("_[[:digit:]]", "", rownames(.))) %>%
  group_by(Term) %>%
  summarise("Neff" = sum(n_eff))

Effs_Chronic

# Averaging effects of tensor

CC_Full_HT %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("C_Sigma" = exp(b_sigma_Intercept),
    "S_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L)) %>% 
    summarise("C_Est" = mean(C_Sigma), 
    "C_LCL" = quantile(C_Sigma, 0.025, type = 8),
    "C_UCL" = quantile(C_Sigma, 0.975, type = 8),
    "S_Est" = mean(S_Sigma), 
    "S_LCL" = quantile(S_Sigma, 0.025, type = 8),
    "S_UCL" = quantile(S_Sigma, 0.975, type = 8))

# Evaluating estimates and CIs for sigma parameters on original scale

CC_Full_HT %>%
    spread_draws(
      sds_t2HourDirection_1, sds_t2HourDirection_2
    ) %>% 
    summarise(Mean_Tens = mean(c(sds_t2HourDirection_1, sds_t2HourDirection_2)),
    Lower_CI = quantile(c(sds_t2HourDirection_1, sds_t2HourDirection_2), 0.025, type  = 8),
    Upper_CI = quantile(c(sds_t2HourDirection_1, sds_t2HourDirection_2), 0.975, type  = 8))

# Assessing AC and effective sample size:sample size ratio

stan_ac(CC_Full_HT$fit)
mcmc_neff(neff_ratio(CC_Full_HT), size = 2)

# Reasonable. All Neff/N > 0.75. Rhats?

mcmc_rhat(rhat(CC_Full_HT))

# Great. Checking posterior and residual distributions.

pp_check(CC_Full_HT, nsamples = 100, stat = "mean")
pp_check(CC_Full_HT, nsamples = 100, stat = "median")
yrep <- posterior_predict(CC_Full_HT, nsamples = 100)

y <- Third_Bound %>%
  dplyr::select(
    mW, Treatment, Sex,
    s.Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo
  ) %>%
  na.omit(.) %>%
  pull(mW)

Group <- Third_Bound %>%
  dplyr::select(
    mW, Treatment, Sex,
    s.Amb.Temp, Hour, Direction, Bird.ID,
    Pen, Date.of.Photo
  ) %>%
  na.omit(.) %>%
  pull(Treatment)

ppc_violin_grouped(y, yrep, group = Group, probs = c(.05, .95), alpha = 0.05, y_draw = "points")

# Very reasonable fits.

ppc_scatter_avg_grouped(y, yrep, group = Group)

# Good, and slightly tighter fit for the stress group.

CC_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full_HT) %>%
  group_by(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(CC_Full_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

CC_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

# Looks good.

CC_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(CC_Full_HT) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Some subtle sway at low end, but appears reasonable. Seeing if
# trending manifests acrosss a given predictor.

{
  p1 <- CC_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- CC_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = s.Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- CC_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- CC_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full_HT) %>%
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

  p5 <- CC_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- CC_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird.ID, y = z_residual, fill = Bird.ID)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Distributions look quite fair. Testing across response and fitted values.

F_by_R <- CC_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(
    x = mW, y = z_residual,
    fill = Treatment, linetype = Treatment
  )) +
  geom_point(size = 2, pch = 21, colour = "black") +
  scale_fill_viridis_d(begin = 0.3, end = 0.6)

F_by_P <- CC_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full_HT) %>%
  summarise(
    Pred = .prediction,
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(
    x = Pred, y = z_residual,
    fill = Treatment, linetype = Treatment
  )) +
  geom_point(size = 2, pch = 21, colour = "black") +
  scale_fill_viridis_d(begin = 0.3, end = 0.6)

grid.arrange(F_by_R, F_by_P, nrow = 1)

# Checking fit

CC_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Full_HT, type = "response") %>%
  group_by(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = mW)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Tight fit.
# Calculating leave-one-out ("loo") residuals

loo_R2(CC_Full_HT)

# Loo-R2 = 0.93. Very tight - ambient temperature expectedly dictating the majority of this effect. Testing differences in sigma between treatment types. First double-checking direction of potential difference.

summary(CC_Full_HT)
CC_Full_HT %>% spread_draws(b_sigma_Intercept, b_sigma_Treatment.L) %>%
  rename("Int_S" = b_sigma_Intercept, "T_S" = b_sigma_Treatment.L) %>%
  summarise("Mean_Int_S" = mean(exp(Int_S)), "LCL_Int" = quantile(exp(Int_S), 0.025, type = 8),
  "UCL_Int" = quantile(exp(Int_S), 0.975, type = 8), "Mean_T_S" = mean(exp(Int_S + T_S)),
  "LCL_TS" = quantile(exp(Int_S + T_S), 0.025, type = 8), 
  "UCL_TS" = quantile(exp(Int_S + T_S), 0.975, type = 8))

CC_Full_HT %>% spread_draws(b_sigma_Intercept, b_sigma_Treatment.L) %>%
  rename("Int_S" = b_sigma_Intercept, "T_S" = b_sigma_Treatment.L) %>% 
  mutate("Sigma_Difference" = Int_S - T_S) %>% 
  ggplot(aes(x = Sigma_Difference)) + 
  geom_density(colour = "black", fill = "cornflowerblue", alpha = 0.5, adjust = 3) + 
  xlim(c(0,0.9)) + geom_vline(xintercept = 0, size = 1, colour = "firebrick") + 
  theme_bw()

# Good. Running hypothesis test:

Chronic_Hyp = hypothesis(CC_Full_HT, "exp(sigma_Intercept) > exp(sigma_Intercept + sigma_Treatment.L)", alpha = 0.05)
plot(Chronic_Hyp)
print(Chronic_Hyp, digits = 3)

# Moderate evidence for a reduction in variance among treatment birds; K = 9.619. Using custom priors.

sigma_chronic <- CC_Full_HT %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("Control_Sigma" = exp(b_sigma_Intercept),
    "Stress_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L))

DF_Sig_Chronic <- build_hdf(vars = list(sigma_chronic$Control_Sigma, sigma_chronic$Stress_Sigma),
  priors = list(rnorm(nrow(sigma_chronic), 0, 0.25), rnorm(nrow(sigma_chronic), 0, 0.25)),
  names = c("Control", "Stress"))

df_sigma_chronic_hyp = hypothesis_df("Stress = Control", DF_Sig_Chronic, class = "b", alpha = 0.05)

plot(df_sigma_chronic_hyp)
print(df_sigma_chronic_hyp)

df_sigma_chronic_hyp = hypothesis_df("Stress < Control", DF_Sig_Chronic, class = "b", alpha = 0.05)

plot(df_sigma_chronic_hyp)
print(df_sigma_chronic_hyp)

# Again, evidence ratio remains at 9.62.

# Saving supplemental plot 

sFig_2 = plot(df_sigma_chronic_hyp, plot = F, theme = theme_get())[[1]]

sFig_2_Final = sFig_2 + theme_bw() + ylab("Density") + xlab(TeX('$\\sigma^2_{Control} - \\sigma^2_{Stress}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/SFig2_HT.jpeg",
  sFig_2_Final,
  height = 6.0, width = 7.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Comparing repeatabilities between treatment groups using non-linear hypothesis test.

Treat_hyp_test_Chronic = c("(sds_ss.Amb.TempBird.IDTreatmentStress_1)^2/((sds_ss.Amb.TempBird.IDTreatmentStress_1)^2 + (exp(b_sigma_Treatment.L))^2) - (sds_ss.Amb.TempBird.ID_1)^2/((sds_ss.Amb.TempBird.ID_1)^2 + (exp(b_sigma_Intercept))^2) > 0.1")
THT_Chronic_Out = hypothesis(CC_Full, Treat_hyp_test_Chronic, alpha = 0.05, class = NULL)
plot(THT_Chronic_Out)
print(THT_Chronic_Out, digits = 3)

# And using custom priors

repeat_chronic_custom <- CC_Full_HT %>%
    spread_draws(
      sds_ss.Amb.TempBird.ID_1,
      sds_ss.Amb.TempBird.IDTreatmentStress_1,
      b_sigma_Intercept,
      b_sigma_Treatment.L
     ) %>% 
    mutate("Control_Repeat" = sds_ss.Amb.TempBird.ID_1^2/(exp(b_sigma_Intercept)^2 + sds_ss.Amb.TempBird.ID_1^2),
    "Stress_Repeat" = sds_ss.Amb.TempBird.IDTreatmentStress_1^2/(exp(b_sigma_Intercept + b_sigma_Treatment.L)^2 + sds_ss.Amb.TempBird.IDTreatmentStress_1^2))

repeat_chronic_custom %>% summarise("R_Stress" = mean(Stress_Repeat), "S_LCL" = quantile(Stress_Repeat, 0.025, type = 8),
"S_UCL" = quantile(Stress_Repeat, 0.975, type = 8), "R_Control" = mean(Control_Repeat), "C_LCL" = quantile(Control_Repeat, 0.025, type = 8), "C_UCL" = quantile(Control_Repeat, 0.975, type = 8))

DF_chronic_rep <- build_hdf(vars = list(repeat_chronic_custom$Control_Repeat, repeat_chronic_custom$Stress_Repeat),
  priors = list(rbeta(nrow(repeat_chronic_custom), 1, 4), rbeta(nrow(repeat_chronic_custom), 1, 4)),
  names = c("Control", "Stress"))

df_chronic_rep_hyp = hypothesis_df("Stress = Control", DF_chronic_rep, class = "b", alpha = 0.05)
plot(df_chronic_rep_hyp)
print(df_chronic_rep_hyp)

df_chronic_rep_hyp = hypothesis_df("Stress > Control", DF_chronic_rep, class = "b", alpha = 0.05)
plot(df_chronic_rep_hyp)
print(df_chronic_rep_hyp)

# Again, moderate evidence for an increase in repeatability among stress exposure treatments. Saving supplemental plot

sFig_4 = plot(df_chronic_rep_hyp, plot = F, theme = theme_get())[[1]]

sFig_4_Final = sFig_4 + theme_bw() + ylab("Density") + xlab(TeX('$R_{Stress} - R_{Control}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/SFig4_HT.jpeg",
  sFig_4_Final,
  height = 6.0, width = 7.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Running null model.

CC_Null_HT <- brm(bf(
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
file = "/home/joshk/Desktop/Chronic_HT_Null.Rds"
)

# Again, checking residuals and fit. Note that residual patterns are unlikely to be
# clean owing to randomisation of bird identities.

stan_ac(CC_Null_HT$fit)
mcmc_neff(neff_ratio(CC_Null_HT), size = 2)
mcmc_rhat(rhat(CC_Null_HT))
conditional_smooths(CC_Null_HT)

# Reasonable. Minimum Neff/N falling slightly below 0.75, but Rhats suggest strong mixing.

CC_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Null_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

CC_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(CC_Null_HT) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Reasonably normal - quite surprisingly. Across predictors:

{
  p1 <- CC_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- CC_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(
      x = s.Amb.Temp, y = z_residual,
      fill = Treatment, linetype = Treatment
    )) +
    geom_point(size = 2, pch = 21, colour = "black") +
    geom_smooth(method = "loess", colour = "black") +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p3 <- CC_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- CC_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null_HT) %>%
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

  p5 <- CC_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- CC_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex,
      Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(CC_Null_HT) %>%
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
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Largely homoskedastic. Fit?

CC_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(CC_Null_HT, type = "response") %>%
  group_by(
    mW, Treatment, Sex,
    Hour, s.Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = mW)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

loo_R2(CC_Null_HT)

# Loo-R2 is expectedly lower, but by litle. Presumably repeatability of slopes is low.
# Comparing repeatability estimates from each model.

get_variables(CC_Null_HT)
get_variables(CC_Null_HT)

# Note repeatability is calculated for the stress-exposure group
# alone here.

Rep_Samples <- rbind(
  data.frame(CC_Full_HT %>%
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
  data.frame(CC_Null_HT %>%
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

# Quickly summarising

Rep_Samples %>% group_by(Mod_Type) %>% 
  summarise(Rep = mean(Repeatability))

Rep_Plot <- Rep_Samples %>% ggplot(aes(x = Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  #scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
  scale_fill_manual(name = NULL, values = c("black", "firebrick4")) + 
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Slopes") +
  my.theme

Rep_Plot
showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Slopes_HT.jpeg", Rep_Plot,
  height = 6.0, width = 6.0, dpi = 800,
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Slopes_HT.pdf", Rep_Plot,
       height = 6.0, width = 6.0, dpi = 800,
)
showtext_auto(enable = FALSE)

## Plotting within each treatment type

Rep_Treats <- rbind(
  data.frame(CC_Full_HT %>%
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
  data.frame(CC_Null_HT %>%
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

S_Rep <- Rep_Treats %>% ggplot(aes(x = T_Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
  scale_y_continuous(trans = "sqrt", name = "Posterior Density") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Repeatability of Stress Slopes") +
  my.theme +
  guides(fill = guide_legend(nrow = 1)) +
  theme(
    legend.text = element_text(size = 10, family = "Noto Sans"),
    legend.spacing.x = unit(0.4, "cm")
  )

C_Rep <- Rep_Treats %>% ggplot(aes(x = C_Repeatability, fill = Mod_Type)) +
  geom_density(alpha = 0.5, adjust = 5) +
  xlim(0, 1) +
  scale_fill_viridis_d(begin = 0.3, end = 0.7, name = NULL) +
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

grid.arrange(S_Rep, C_Rep, nrow = 2)

# Stress slopes more repeatable, suggesting possibility for selection on these responses?
# Testing differences in repeatability witih treatment groups

BTT <- ttestBF(
  formula = T_Repeatability ~ Mod_Type, data = Rep_Treats,
  iterations = 10000, rscale = "ultrawide", posterior = FALSE
)

Post_Est <- posterior(BTT, iterations = 1000)
plot(Post_Est[, "mu"])
plot(acf(Post_Est[, "mu"]))

# Chains relatively stable and autocorrelation negligable. Assessing Bayes factor

summary(BTT)

# BF > 100. Strong evidence for alternative over null. Quantifying
# then using Savage-Dickey approach to calculating Bayes Factor.

Rep_Treats %>%
  group_by(Mod_Type) %>%
  summarise(
    "Mean" = mean(T_Repeatability),
    "UCL" = quantile(T_Repeatability, 0.975, type = 8),
    "LCL" = quantile(T_Repeatability, 0.025, type = 8)
  )

# No overlap of credible intervals.

Wide_Chronic_Rep <- Rep_Treats %>%
  mutate("Sample" = paste(.chain, .iteration, sep = "_")) %>%
  dplyr::select(Sample, Mod_Type, T_Repeatability) %>%
  spread(Mod_Type, T_Repeatability) %>%
  rename("True" = "True\nModel  ", "Null" = "Null\nModel  ") %>%
  as.data.frame(.)

hypothesis(Wide_Chronic_Rep, "True = Null", alpha = 0.05)
Keep = hypothesis(Wide_Chronic_Rep, "True - Null > 0.25", alpha = 0.05)
plot(hypothesis(Wide_Chronic_Rep, "True - Null > 0.25", alpha = 0.05))

# And with custom hypothesis function to load in priors 

DF <- build_hdf(vars = list(Wide_Chronic_Rep$True, Wide_Chronic_Rep$Null),
  priors = list(rbeta(nrow(Wide_Chronic_Rep), 1, 4), rbeta(nrow(Wide_Chronic_Rep), 1, 4)),
  names = c("True", "Null"))

df_hyp = hypothesis_df("True - Null > 0", DF, class = "b", alpha = 0.05)
df_hyp_cons = hypothesis_df("True - Null > 0.1", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)
plot(df_hyp_cons)
print(df_hyp_cons)

# Again, Bayes factor is extremely high (in this case, infinite)!
# And lastly, using a model in brms where variance is unequal between groups and prior is default.

Robust_Comp_data <- Rep_Treats %>%
  dplyr::select(Mod_Type, T_Repeatability) %>%
  mutate(Model = gsub("\\n.*", "", Mod_Type)) %>%
  dplyr::select(-Mod_Type)

# Checking priors first, but with temporary prior on sigma by treatment

robust_comp <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student,
  data = Robust_Comp_data, sample_prior = "only",
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.05)", class = "nu"),
    set_prior("exponential(2)", class = "b", dpar = "sigma", lb = 0.001)
  )
)

y_est <- posterior_predict(robust_comp)
prior_pc = data.frame("Type" = c(rep("Prior", length(colMeans(y_est))),
  rep("True", length(colMeans(y_est)))), 
  "Y" = c(colMeans(y_est), Robust_Comp_data$T_Repeatability))

ggplot(prior_pc, aes(x = Y, fill = Type)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.4, end = 0.6)

## Priors are very wide, but given that on sigma will drop, unlikely to be an issue.

robust_comp <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student, sample_prior = TRUE,
  data = Robust_Comp_data,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ),
  file = "/home/joshk/Desktop/CC_Results/Final/Robust_Comp_Chronic_HT.Rds"
)

pp_check(robust_comp)
mcmc_neff(neff_ratio(robust_comp), size = 2)
mcmc_rhat(rhat(robust_comp))

# Clean coverage, Neff/N and Rhats.

plot(robust_comp)
summary(robust_comp)
plot(hypothesis(robust_comp, "ModelTrue > ModelNull"))
print(hypothesis(robust_comp, "ModelTrue > ModelNull"), digits = 3)
print(hypothesis(robust_comp, "ModelTrue - ModelNull > 0.25"), digits = 3)

### Pulling out qualitative differences in means and SDs

tidy_MCMC_Chronic <- tidyMCMC(robust_comp, conf.int = TRUE, conf.level = 0.95, 
  estimate.method = "median", conf.method = "HPDinterval") %>% 
  mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

Post_Dif_Chronic <- posterior_samples(robust_comp) %>% 
  mutate_at(vars(contains("sigma")), funs(exp)) %>% 
  mutate(nu = log10(nu)) %>% 
  mutate(diff_means = b_ModelTrue - b_ModelNull,
         diff_sigma = b_sigma_ModelTrue - b_sigma_ModelNull) %>% 
  mutate(cohen_d = diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)/2),
         cles = dnorm(diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)), 0, 1))

Out_Chronic <- tidyMCMC(Post_Dif_Chronic, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")
Out_Chronic

# Lastly, as permutation test 

set.seed(200) 

simu <- 1000
res <- numeric(simu) 

for (i in 1:simu) {
    perm <- sample(nrow(Robust_Comp_data))
    bdat <- transform(Robust_Comp_data, Model = Model[perm])
    res[i] <- mean(bdat[bdat$Model=="True", "T_Repeatability"])-
        mean(bdat[bdat$Model=="Null", "T_Repeatability"])
}

obs <- mean(Robust_Comp_data[Robust_Comp_data$Model == "True", "T_Repeatability"])-
    mean(Robust_Comp_data[Robust_Comp_data$Model == "Null", "T_Repeatability"])

ggplot(as.data.frame(res), aes(x = res)) + 
  geom_density(adjust = 3, fill = "mediumorchid", alpha = 0.5) + 
  geom_vline(xintercept = obs, colour = "black")

mean(abs(res)>=abs(obs))

# All results pointing in the same direction. 

# Extracting slopes from urban and rural birds

get_variables(CC_Full_HT)

Slopes <- CC_Full_HT %>%
  spread_draws(s_ss.Amb.TempBird.IDTreatmentStress_1[ID]) %>%
  rename("Coef" = s_ss.Amb.TempBird.IDTreatmentStress_1) %>%
  as.data.frame(.)

Birds <- data.frame(
  "Bird_ID" = rownames(ranef(CC_Full_HT)$Bird.ID),
  "ID" = seq(1, 19, 1)
)
Collapsed <- left_join(Slopes, Birds, by = c("ID")) %>%
  dplyr::select(-ID) %>%
  group_by(Bird_ID) %>%
  summarise("Slope" = mean(Coef)) %>%
  as.data.frame(.)

# Binding in locale and treatment order

U_R <- c()
Locale <- c()

for (i in 1:nrow(Collapsed)) {
  if (is.na(Collapsed$Bird_ID[i])) {
    Locale[i] <- NA
    U_R[i] <- NA
  } else if (Collapsed$Bird_ID[i] == "BdABd") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdAO1") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "ABlBl") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "AYY") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "YAO1") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "OOA") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdABl1") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "YAO2") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "ABlO") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "ABdBd") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "YAR") {
    Locale[i] <- "Cambridge"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "RAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "ABlR") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "ARO") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdABl2") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BlAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Collapsed$Bird_ID[i] == "BdAO2") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "AOR") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "YAY") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Collapsed$Bird_ID[i] == "AYBd") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else {
    Locale[i] <- NA
    U_R[i] <- NA
  }
}
Collapsed$U_R <- factor(U_R)
Collapsed$Locale <- factor(Locale)
Orders <- Third_Bound %>%
  dplyr::select("Bird_ID" = Bird.ID, Pen) %>%
  na.omit(.) %>%
  distinct(.) %>%
  mutate("Order" = ifelse(Pen == "NE" | Pen == "SW", "SR", "RS")) %>%
  dplyr::select(-Pen)
Collapsed <- left_join(Collapsed, Orders, by = c("Bird_ID")) %>%
  mutate(Order = factor(Order), U_R = factor(U_R))

# Summarising and Running "ANOVA"

Collapsed %>%
  group_by(U_R) %>%
  summarise(
    Mean = mean(Slope),
    LCL = quantile(Slope, 0.025, type = 8),
    UCL = quantile(Slope, 0.975, type = 8)
  )

# Considerable overlap.

aovBF <- anovaBF(Slope ~ Order * U_R,
  data = Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(aovBF)
summary(aovBF)

# Order may matter here, but only in the presence of an interaction (with large error). Rationale for this appears rather unclear. Removing order and test.

aovBF <- anovaBF(Slope ~ U_R,
  data = Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(BayesFactor::posterior(aovBF, iterations = 10000)[, "mu"])
plot(acf(BayesFactor::posterior(aovBF, iterations = 10000)[, "mu"]))

# Chain appears stable. No autocorrelation evident. Summarising

summary(aovBF)
summary(BayesFactor::posterior(aovBF, iterations = 10000)[, "mu"])

# No evidence for differences between urban and rural individuals. Sitting quite close to zero!
# Again, approaching using Savage-Dickey method. First checking priors.

Slope_Contrast <- brm(Slope ~ 0 + U_R + (1 | Locale),
  prior = c(
    set_prior("normal(0, 5)", class = "b"),
    set_prior("gamma(1, 0.5)", class = "sd"),
    set_prior("exponential(0.05)", class = "sigma")
  ), data = Collapsed, cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  sample_prior = "only")

y_est <- posterior_predict(Slope_Contrast)
prior_pc = data.frame("Type" = c(rep("Prior", length(colMeans(y_est))),
  rep("True", length(colMeans(y_est)))), 
  "Y" = c(colMeans(y_est), Collapsed$Slope))

ggplot(prior_pc, aes(x = Y, fill = Type)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.4, end = 0.6)

# Priors seem reasonably appropriate, however, true response is clearly noisy.

Slope_Contrast <- brm(Slope ~ 0 + U_R + (1 | Locale),
  prior = c(
    set_prior("normal(0, 5)", class = "b"),
    set_prior("gamma(1, 0.5)", class = "sd"),
    set_prior("exponential(0.05)", class = "sigma")
  ), data = Collapsed, cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  sample_prior = TRUE,
  file = "/home/joshk/Desktop/CC_Results/Final/Chronic_UR_Contrast.Rds"
)

Chronic_Hypoth <- hypothesis(Slope_Contrast, c(
  "U_RRural = U_RUrban"
))
plot(Chronic_Hypoth)
1 / Chronic_Hypoth$hypothesis$Evid.Ratio

# Again, results corroborate. Plotting.

UR_Plot <- Collapsed %>%
  rename("Ecotype" = U_R) %>%
  group_by(Ecotype) %>%
  summarise(
    "Coefficient" = mean(Slope),
    "lcl" = quantile(Slope, 0.025, type = 8),
    "ucl" = quantile(Slope, 0.975, type = 8)
  ) %>%
  ggplot(aes(x = Ecotype, y = Coefficient, fill = Ecotype)) +
  geom_errorbar(aes(x = Ecotype, ymin = lcl, ymax = ucl),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 6, pch = 21, colour = "black") +
  #scale_fill_viridis_d(begin = 0.2, end = 0.5) +
  scale_fill_manual(values = c("wheat4","cornsilk")) + 
  theme_bw() +
  my.theme +
  theme(legend.position = "none") +
  ylab("Reaction Norm Slope")

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Chronic_HT.jpeg",
  UR_Plot,
  height = 5.67, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Chronic_HT.pdf",
       UR_Plot,
       height = 5.67, width = 6.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Plotting individual variability

Simp_Dat_Chronic <- expand.grid(
  "s.Amb.Temp" = seq(min(CC_Full_HT$data$s.Amb.Temp, na.rm = TRUE),
    max(CC_Full_HT$data$s.Amb.Temp, na.rm = TRUE),
    length.out = 25
  ),
  "Hour" = c(12),
  "Date.of.Photo" = na.omit(unique(CC_Full_HT$data$Date.of.Photo))[30],
  "Treatment" = na.omit(unique(CC_Full_HT$data$Treatment)),
  "Pen" = na.omit(unique(CC_Full_HT$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(CC_Full_HT$data$Sex)),
  "Bird.ID" = na.omit(unique(CC_Full_HT$data$Bird.ID)),
  "Direction" = na.omit(unique(CC_Full_HT$data$Direction)),
  "Treat_Day" = na.omit(unique(CC_Full_HT$data$Treat_Day))[c(14, 70)]
)

Third_Bound$s.Amb.Temp <- (Third_Bound$Amb.Temp - mean(Third_Bound$Amb.Temp, na.rm = T)) /
  (2 * sd(Third_Bound$Amb.Temp, na.rm = T))

Chronic_Plot_Dat <- Simp_Dat_Chronic %>%
  add_predicted_draws(CC_Full_HT, n = 1000, scale = "response") %>%
  mutate(
    "Group_ID" = paste(Bird.ID, Treatment, sep = "_"),
    "Amb.Temp" = (s.Amb.Temp * (2 * sd(Third_Bound$Amb.Temp, na.rm = T))) +
      mean(Third_Bound$Amb.Temp, na.rm = T)
  ) %>%
  group_by(Group_ID, Bird.ID, Treatment, Amb.Temp) %>%
  summarise("Pred" = mean(.prediction))

Chronic_Curves <- Chronic_Plot_Dat %>%
  ggplot(aes(x = Amb.Temp, y = Pred, colour = Treatment)) +
  geom_line(aes(group = Group_ID)) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Dry Heat Loss (mW)") +
  theme_bw() +
  my.theme

Chronic_Curves

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_OP_HT.jpeg",
  Chronic_Curves,
  height = 6.0, width = 7.0, dpi = 800
)
showtext_auto(enable = FALSE)

#### Moving on to acute models

prior_acute <- c(
  set_prior("normal(0, 2.5)", class = "b", coef = "SexMale"),
  set_prior("normal(-1, 2.5)", class = "b", coef = "Treatment.L"),
  set_prior("normal(-5, 5)", class = "b", coef = "sAmb.Temp_1"),
  set_prior("gamma(2, 0.5)", class = "sds"),
  set_prior("gamma(50, 2)", class = "Intercept")
)

# And running full model.

Acute_Full_HT <- brm(bf(
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
  "/home/joshk/Desktop/Acute_HT.Rds"
)

# Checking posteriors

pp_check(Acute_Full_HT, nsamples = 100)

# Nice ovelay but some straying around 9 mW. Minimal, but assessing with violoin plot.

yrep_Acute <- posterior_predict(Acute_Full_HT, nsamples = 100)

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

# Distributions reasonably fitting. Checking scatter by treatments

ppc_scatter_avg_grouped(Acute_y, yrep_Acute, group = Acute_Group) + 
geom_smooth(method = "lm", colour = "black")
summary(lm(Acute_y ~ Acute_Full$data$mW))

# Good; 1:1 relationship between estimated and true heat transfer. Checking effective sample sizes

effectiveSize(Acute_Full_HT)
plot(Acute_Full_HT)
summary(Acute_Full_HT)
conditional_smooths(Acute_Full_HT)
get_variables(Acute_Full_HT)

Mean_Tensor_Coefs <- Acute_Full_HT %>%
  spread_draws(
    sds_t2TimelineAmb.TempTreatmentStress_1,
    sds_t2TimelineAmb.TempTreatmentStress_2,
    sds_t2TimelineAmb.TempTreatmentStress_3,
    sds_t2HourDirection_1,
    sds_t2HourDirection_2,
    sds_t2HourDirection_3,
    b_sigma_Intercept,
    b_sigma_Treatment.L
  ) %>%
  rename(
    "TA_Tens_1" = sds_t2TimelineAmb.TempTreatmentStress_1,
    "TA_Tens_2" = sds_t2TimelineAmb.TempTreatmentStress_2,
    "TA_Tens_3" = sds_t2TimelineAmb.TempTreatmentStress_3,
    "Hour_Tens_1" = sds_t2HourDirection_1,
    "Hour_Tens_2" = sds_t2HourDirection_2,
    "Hour_Tens_3" = sds_t2HourDirection_3
  ) %>%
  summarise(
    "Mean_TAE" = mean(c(TA_Tens_1, TA_Tens_2, TA_Tens_3)),
    "TA_LCL" = quantile(c(TA_Tens_1, TA_Tens_2, TA_Tens_3),
      0.025,
      type = 8
    ),
    "TA_UCL" = quantile(c(TA_Tens_1, TA_Tens_2, TA_Tens_3),
      0.975,
      type = 8
    ),
    "Mean_HourE" = mean(c(Hour_Tens_1, Hour_Tens_2, Hour_Tens_3)),
    "Hour_LCL" = quantile(c(Hour_Tens_1, Hour_Tens_2, Hour_Tens_3),
      0.025,
      type = 8
    ),
    "Hour_UCL" = quantile(c(Hour_Tens_1, Hour_Tens_2, Hour_Tens_3),
      0.975,
      type = 8
    )
  )

Mean_Tensor_Coefs

Effs = rstan::summary(Acute_Full_HT$fit)$summary %>%
  as.data.frame(.) %>%
  dplyr::select(n_eff) %>%
  filter(rownames(.) %in% c(
    "b_Intercept", "b_Treatment.L", "b_SexMale",
    "sds_sAmb.Temp_1", "sds_sTimeline_1",
    "sds_sAmb.TempTreatmentStress_1",
    "sds_sTimelineTreatmentStress_1",
    "sds_t2TimelineAmb.TempTreatmentStress_1",
    "sds_t2TimelineAmb.TempTreatmentStress_2",
    "sds_t2TimelineAmb.TempTreatmentStress_3",
    "sds_t2HourDirection_1",
    "sds_t2HourDirection_2",
    "sds_t2HourDirection_3",
    "sds_sTimelineBird.ID_1",
    "sds_sTimelineBird.IDTreatmentStress_1",
    "sd_Bird.ID__Intercept",
    "sd_Date.of.Photo__Intercept",
    "sd_Pen__Intercept",
    "b_sigma_Intercept",
    "b_sigma_Treatment.L"
  )) %>%
  mutate("Term" = gsub("_[[:digit:]]", "", rownames(.))) %>%
  group_by(Term) %>%
  summarise("Neff" = sum(n_eff))

Effs

# Validating model.
stan_ac(Acute_Full_HT$fit)
mcmc_neff(neff_ratio(Acute_Full_HT), size = 2)

# Neff/N ratios > 0.75. No clear autocorrelation. Checking Rhats.
mcmc_rhat(rhat(Acute_Full_HT))

# Great - all near 1. Chains are well mixed. Checking residual distributions.

Acute_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT) %>%
  group_by(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(Acute_Full_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

Acute_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

# Reasonbly well fit and no clear pattern across fitted values. Checking residuals with density plot

Acute_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(Acute_Full_HT) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Good, but possible outliers contributing to wide left tail. Producing Cleveland dotplot then
# plotting residuals by predictor to assess where outliers are biologically meaningful or
# concerning.

P_Dat <- Acute_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  )

ggplot(P_Dat, aes(x = 1:nrow(P_Dat), y = z_residual)) +
  geom_point(colour = "black") +
  ylim(c(-4, 4))

# Note infinite values on each end. Some trailing on low and high ends, perhaps due to
# low Ta. Assessing patterns across predictors.

{
  p1 <- Acute_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- Acute_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT) %>%
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

  p3 <- Acute_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- Acute_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT) %>%
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

  p5 <- Acute_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- Acute_Full_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird.ID,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Full_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Bird.ID, y = z_residual, fill = Bird.ID)) +
    stat_summary(
      geom = "errorbar", fun.data = "mean_se", size = 1,
      colour = "black", width = 0.2
    ) +
    stat_summary(
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# No apparent trends or heteroskedasticity across predictors, although peculiar ordering
# among individuals. Assessing model fit

Acute_Full_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Full_HT, scale = "response") %>%
  group_by(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird.ID,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = mW)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Again, good fit. Calculating loo-R2

loo_R2(Acute_Full_HT)

# Again, quite high. 93.5%. Quickly assessing descriptive statistics.

Descriptive = expand.grid(
  "Amb.Temp" = c(min(Acute_Full_HT$data$Amb.Temp), max(Acute_Full_HT$data$Amb.Temp)),
  "Timeline" = c(-1200, 1200),
  "Hour" = c(8, 12, 16),
  "Sex" = na.omit(unique(Acute_Full_HT$data$Sex)),
  "Direction" = na.omit(unique(Acute_Full_HT$data$Direction)),
  "Treatment" = na.omit(unique(Acute_Full_HT$data$Treatment)),
  "Bird.ID" = na.omit(unique(Acute_Full_HT$data$Bird.ID)),
  "Treat_Day" = na.omit(unique(Acute_Full_HT$data$Treat_Day))[1]
)

Global_Pred = posterior_predict(Acute_Full_HT, newdata = Descriptive, re_formula = NA,
  nsamples = 1000, summary = FALSE)

Descriptive %>% mutate("Post" = colMeans(Global_Pred)) %>% 
  filter(Treatment == "Stress", 
  Amb.Temp %in% c(min(Acute_Full_HT$data$Amb.Temp), max(Acute_Full_HT$data$Amb.Temp))) %>% 
  mutate("Measure_Time" = ifelse(Timeline == -1200, "Start", "End")) %>% 
  dplyr::select(-Timeline) %>% 
  group_split(Measure_Time) %>% 
  reduce(left_join, by = c("Amb.Temp", "Hour", "Sex", "Direction", "Treatment", "Bird.ID", "Treat_Day")) %>% 
  rename("End_Temp" = Post.x, "Start_Temp" = Post.y) %>% 
  dplyr::select(-c(Measure_Time.x, Measure_Time.y)) %>% 
  mutate("Diff" = End_Temp - Start_Temp) %>% 
  group_by(Amb.Temp) %>% 
  summarise("Mean_Diff" = mean(Diff), "LCL" = quantile(Diff, 0.025, type = 8),
    "UCL" = quantile(Diff, 0.975, type = 8),
    "SD_Diff" = sd(Diff))

# Evaluating means and CIs for sigma parameters on original scale.

Acute_Full_HT %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("C_Sigma" = exp(b_sigma_Intercept),
    "S_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L)) %>% 
    summarise("C_Est" = mean(C_Sigma), 
    "C_LCL" = quantile(C_Sigma, 0.025, type = 8),
    "C_UCL" = quantile(C_Sigma, 0.975, type = 8),
    "S_Est" = mean(S_Sigma), 
    "S_LCL" = quantile(S_Sigma, 0.025, type = 8),
    "S_UCL" = quantile(S_Sigma, 0.975, type = 8))

# Comparing residual variance between treatment groups.

hyp <- c(
  "exp(sigma_Intercept + sigma_Treatment.L) < exp(sigma_Intercept)",
  "exp(sigma_Intercept + sigma_Treatment.L) - exp(sigma_Intercept) = 0"
)
plot(hypothesis(Acute_Full_HT, hyp))
hypothesis(Acute_Full_HT, hyp)

# Using custom sigma priors

sigma_priors <- Acute_Full_HT %>%
    spread_draws(
      b_sigma_Intercept, b_sigma_Treatment.L
    ) %>% 
    mutate("Control_Sigma" = exp(b_sigma_Intercept),
    "Stress_Sigma" = exp(b_sigma_Intercept + b_sigma_Treatment.L))

DF <- build_hdf(vars = list(sigma_priors$Control_Sigma, sigma_priors$Stress_Sigma),
  priors = list(rnorm(nrow(sigma_priors), 0, 0.25), rnorm(nrow(sigma_priors), 0, 0.25)),
  names = c("Control", "Stress"))

df_hyp = hypothesis_df("Stress = Control", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)

df_hyp = hypothesis_df("Stress < Control", DF, class = "b", alpha = 0.05)

plot(df_hyp)
print(df_hyp)

# Moderate support for reduced variance in stressed group. Similar to chronic results.

# Saving supplemental plot 

sFig_1 = plot(df_hyp, plot = F, theme = theme_get())[[1]]

sFig_1_Final = sFig_1 + theme_bw() + ylab("Density") + xlab(TeX('$\\sigma^2_{Control} - \\sigma^2_{Stress}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_SFig1_HT.jpeg",
  sFig_1_Final,
  height = 6.0, width = 7.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Similar to chronic model, comparing repeatabilities between treatment groups using non-linear hypothesis test.

get_variables(Acute_Full_HT)
Treat_hyp_test_Acute = c("(sds_sTimelineBird.IDTreatmentStress_1)^2/((sds_sTimelineBird.IDTreatmentStress_1)^2 + (exp(b_sigma_Treatment.L))^2) - (sds_sTimelineBird.ID_1)^2/((sds_sTimelineBird.ID_1)^2 + (exp(b_sigma_Intercept))^2) > 0.1")
THT_Acute_Out = hypothesis(Acute_Full_HT, Treat_hyp_test_Acute, alpha = 0.05, class = NULL)
plot(THT_Acute_Out)
print(THT_Acute_Out, digits = 3)

# And again, using custom function to load in priors.

repeat_acute_custom <- Acute_Full_HT %>%
    spread_draws(
      sds_sTimelineBird.ID_1,
      sds_sTimelineBird.IDTreatmentStress_1,
      b_sigma_Intercept, b_sigma_Treatment.L
     ) %>% 
    mutate("Control_Repeat" = sds_sTimelineBird.ID_1^2/(exp(b_sigma_Intercept)^2 + sds_sTimelineBird.ID_1^2),
    "Stress_Repeat" = sds_sTimelineBird.IDTreatmentStress_1^2/(exp(b_sigma_Intercept + b_sigma_Treatment.L)^2 + sds_sTimelineBird.IDTreatmentStress_1^2))

repeat_acute_custom %>% summarise("R_Stress" = mean(Stress_Repeat), "S_LCL" = quantile(Stress_Repeat, 0.025, type = 8),
"S_UCL" = quantile(Stress_Repeat, 0.975, type = 8), "R_Control" = mean(Control_Repeat), "C_LCL" = quantile(Control_Repeat, 0.025, type = 8), "C_UCL" = quantile(Control_Repeat, 0.975, type = 8))

DF <- build_hdf(vars = list(repeat_acute_custom$Control_Repeat, repeat_acute_custom$Stress_Repeat),
  priors = list(rbeta(nrow(repeat_acute_custom), 1, 4), rbeta(nrow(repeat_acute_custom), 1, 4)),
  names = c("Control", "Stress"))

df_hyp = hypothesis_df("Stress = Control", DF, class = "b", alpha = 0.05)
plot(df_hyp)
print(df_hyp)

df_hyp = hypothesis_df("Stress > Control", DF, class = "b", alpha = 0.05)
plot(df_hyp)
print(df_hyp)

# This time, mild support for increased repeatability in stressed groups.
# Producing supplemental figure

sFig_3 = plot(df_hyp, plot = F, theme = theme_get())[[1]]

sFig_3_Final = sFig_3 + theme_bw() + ylab("Density") + xlab(TeX('$R_{Stress} - R_{Control}$')) + 
my.theme + theme(strip.text.x = element_blank())

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_SFig3_HT.jpeg",
  sFig_3_Final,
  height = 6.0, width = 7.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Running null model

Acute_Null_HT <- brm(bf(
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
file =
  "/home/joshk/Desktop/Acute_HT_Null.Rds"
)

stan_ac(Acute_Null_HT$fit)
mcmc_neff(neff_ratio(Acute_Null_HT), size = 2)
mcmc_rhat(rhat(Acute_Null_HT))

# Good, Neff/N > 0.75 for all and Rhats tight to 1.

conditional_smooths(Acute_Null_HT)

# Trends remain similar to full model. Plotting residuals and exploring fit.

Acute_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(sample = z_residual)) +
  geom_qq() +
  geom_abline()

Acute_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null_HT) %>%
  group_by(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise("Fit" = mean(.prediction)) %>%
  add_predicted_draws(Acute_Null_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = Fit, y = z_residual)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "cornflowerblue") + 
  theme_bw()

Acute_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_residual_draws(Acute_Null_HT) %>%
  ggplot(aes(x = .residual)) +
  geom_density()

# Fairly normal but wide (similar to true model). Again, producing Cleveland dotplot then plotting residuals by predictors.

P_Dat <- Acute_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null_HT) %>%
  summarise(
    p_residual = mean(.prediction < mW),
    z_residual = qnorm(p_residual)
  )

ggplot(P_Dat, aes(x = 1:nrow(P_Dat), y = z_residual)) +
  geom_point(colour = "black") +
  ylim(c(-4, 4))

# Again, some trailing. Checking across predictors.

{
  p1 <- Acute_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Treatment, y = z_residual)) +
    geom_boxplot() +
    scale_fill_viridis_d(begin = 0.3, end = 0.6)

  p2 <- Acute_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null_HT) %>%
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

  p3 <- Acute_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Sex, y = z_residual)) +
    geom_boxplot()

  p4 <- Acute_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null_HT) %>%
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

  p5 <- Acute_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null_HT) %>%
    summarise(
      p_residual = mean(.prediction < mW),
      z_residual = qnorm(p_residual)
    ) %>%
    ggplot(aes(x = Pen, y = z_residual)) +
    geom_boxplot()

  p6 <- Acute_Null_HT$data %>%
    dplyr::select(
      mW, Treatment, Sex, Timeline,
      Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
      Date.of.Photo, Treat_Day
    ) %>%
    na.omit(.) %>%
    add_predicted_draws(Acute_Null_HT) %>%
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
      geom = "point", fun.y = "mean", size = 2,
      colour = "black"
    ) +
    scale_fill_viridis_d()
}

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)

# Homoskedastic across predictors. Model fit?

Acute_Null_HT$data %>%
  dplyr::select(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  na.omit(.) %>%
  add_predicted_draws(Acute_Null_HT, scale = "response") %>%
  group_by(
    mW, Treatment, Sex, Timeline,
    Hour, Amb.Temp, Direction, Pen, Bird_Scramble,
    Date.of.Photo, Treat_Day
  ) %>%
  summarise(Pred = mean(.prediction)) %>%
  ggplot(aes(x = Pred, y = mW)) +
  geom_point(size = 2, pch = 21, colour = "black", fill = "slateblue", alpha = 0.5) +
  geom_smooth(method = "lm", size = 1, colour = "black", linetype = "dashed") +
  theme_bw()

# Tight fit despite scrambling bird identities. Checking loo-R2.

loo_R2(Acute_Null_HT)

# Highly similar to full model. Here, R2 = 0.93.
# Testing repeatability exclusively within stress-exposed groups.

get_variables(Acute_Null_HT)
Rep_Treats_Acute <- rbind(
  data.frame(Acute_Full_HT %>%
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
  data.frame(Acute_Null_HT %>%
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
  guides(fill = guide_legend(nrow = 1))

grid.arrange(S_Rep_Acute, C_Rep_Acute, nrow = 2)

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_Repeat_HT.jpeg",
  S_Rep_Acute,
  height = 6.0, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Acute_Repeat_HT.pdf",
       S_Rep_Acute,
       height = 6.0, width = 6.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Very wide and uncertain with stress slopes only slightly exceeding those of control slopes.
# Comparing repeatability estimates by Bayesian t-test and Savage-Dickey non-linear hypothesis test.

BTT_Acute <- ttestBF(
  formula = T_Repeatability ~ Mod_Type, data = Rep_Treats_Acute,
  iterations = 10000, rscale = "ultrawide", posterior = FALSE
)

Post_Acute <- posterior(BTT_Acute, iterations = 1000)
plot(Post_Acute[, "mu"])
plot(acf(Post_Acute[, "mu"]))

# Chains a but noisy, but range is small.

summary(BTT_Acute)

# BF strikingly high. Testing using Savage-Dickey non-linear hypothesis test after
# qualitative comparison.

Rep_Treats_Acute %>%
  group_by(Mod_Type) %>%
  summarise(
    "Mean" = mean(T_Repeatability),
    "UCL" = quantile(T_Repeatability, 0.975, type = 8),
    "LCL" = quantile(T_Repeatability, 0.025, type = 8)
  )

# Some overlap, but not striking.

Wide_Acute_Rep <- Rep_Treats_Acute %>%
  mutate("Sample" = paste(.chain, .iteration, sep = "_")) %>%
  dplyr::select(Sample, Mod_Type, T_Repeatability) %>%
  spread(Mod_Type, T_Repeatability) %>% 
  rename("True" = "True Model", "Null" = "Null Model               ") %>%
  as.data.frame(.)

h_test_Acute <- hypothesis(Wide_Acute_Rep, c("True > Null", "True - Null > 0.25"))
plot(h_test_Acute)
print(h_test_Acute, digits = 3)

# Using custom function to apply priors 

DF_Acute_Rep <- build_hdf(vars = list(Wide_Acute_Rep$True, Wide_Acute_Rep$Null), priors = list(rbeta(nrow(Wide_Acute_Rep), 1, 4), rbeta(nrow(Wide_Acute_Rep), 1, 4)), names = c("True", "Null"))

df_acute_hyp = hypothesis_df("True - Null > 0", DF_Acute_Rep, class = "b", alpha = 0.05)

plot(df_acute_hyp)
print(df_acute_hyp)

# Much more logical results. Low support (K = 0.98 at 0.25 level). Note that variances between 
# groups are clearly not equal, however. Modeling with unequal variance and loose but mildly informative priors.

Robust_Comp_data_Acute <- Rep_Treats_Acute %>%
  dplyr::select(Mod_Type, T_Repeatability) %>%
  mutate(Model = gsub("\\n.*", "", Mod_Type)) %>%
  dplyr::select(-Mod_Type)

# Again, checking priors

robust_comp_acute <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student,
  data = Robust_Comp_data_Acute,
  sample_prior = "only",
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu"),
    set_prior("exponential(2)", class = "b", dpar = "sigma", lb = 0.001)
  )#,
  #file = "/home/joshk/Desktop/CC_Results/Final/Robust_Comp_Acute.Rds"
)

y_est <- posterior_predict(robust_comp_acute)
prior_pc = data.frame("Type" = c(rep("Prior", length(colMeans(y_est))),
  rep("True", length(colMeans(y_est)))), 
  "Y" = c(colMeans(y_est), Robust_Comp_data_Acute$T_Repeatability))

ggplot(prior_pc, aes(x = Y, fill = Type)) + 
  geom_density(colour = "black", alpha = 0.5) + 
  scale_fill_viridis_d(begin = 0.4, end = 0.6)

# Priors are very wide, but removal of sigma prior is likely to help. Running model.

robust_comp_acute <- brm(
  bf(T_Repeatability ~ 0 + Model, sigma ~ 0 + Model),
  family = student,
  data = Robust_Comp_data_Acute,
  sample_prior = TRUE,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ),
  file = "/home/joshk/Desktop/CC_Results/Final/Robust_Comp_Acute_HT.Rds"
)

pp_check(robust_comp_acute)
summary(robust_comp_acute)
mcmc_neff(neff_ratio(robust_comp_acute), size = 2)
mcmc_rhat(rhat(robust_comp_acute))

# Neff/N > 0.75 and Rhats near 1. Note that posterior is predicting a greater divergence between
# groups than is true, however, proximitey of peaks may render this seperation negligable.

# Comparing coefficients with Savage-Dickey

plot(hypothesis(robust_comp_acute, c("ModelTrue > ModelNull", "ModelTrue - ModelNull > 0.25")))
print(hypothesis(robust_comp_acute, c("ModelTrue > ModelNull", "ModelTrue - ModelNull > 0.25")))

# No difference according to conservative measure.
### Pulling out qualitative differences in means and SDs

tidy_MCMC_Acute <- tidyMCMC(robust_comp_acute, conf.int = TRUE, conf.level = 0.95, 
  estimate.method = "median", conf.method = "HPDinterval") %>% 
  mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

Post_Dif_Acute <- posterior_samples(robust_comp_acute) %>% 
  mutate_at(vars(contains("sigma")), funs(exp)) %>% 
  mutate(nu = log10(nu)) %>% 
  mutate(diff_means = b_ModelTrue - b_ModelNull,
         diff_sigma = b_sigma_ModelTrue - b_sigma_ModelNull) %>% 
  mutate(cohen_d = diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)/2),
         cles = dnorm(diff_means / sqrt((b_sigma_ModelNull + b_sigma_ModelTrue)), 0, 1))

Out_Acute <- tidyMCMC(Post_Dif_Acute, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval")
Out_Acute

# Lastly, as permutation test 

set.seed(200) 

simu <- 1000
res <- numeric(simu) 

for (i in 1:simu) {
    perm <- sample(nrow(Robust_Comp_data_Acute))
    bdat <- transform(Robust_Comp_data_Acute, Model = Model[perm])
    res[i] <- mean(bdat[bdat$Model=="True", "T_Repeatability"])-
        mean(bdat[bdat$Model=="Null", "T_Repeatability"])
}

obs <- mean(Robust_Comp_data_Acute[Robust_Comp_data_Acute$Model == "True", "T_Repeatability"])-
    mean(Robust_Comp_data_Acute[Robust_Comp_data_Acute$Model == "Null", "T_Repeatability"])

ggplot(as.data.frame(res), aes(x = res)) + 
  geom_density(adjust = 3, fill = "mediumorchid", alpha = 0.5) + 
  geom_vline(xintercept = obs, colour = "black")

mean(abs(res)>=abs(obs))

# Again, clear differences between posterior means, however, sample size is likely a problem
# (that is, being inflated).

# Pulling out slopes and comparing between urban and rural birds.

Acute_Slopes <- Acute_Full_HT %>%
  spread_draws(s_sTimelineBird.IDTreatmentStress_1[ID]) %>%
  rename("Coef" = s_sTimelineBird.IDTreatmentStress_1) %>%
  as.data.frame(.)

Acute_Birds <- data.frame(
  "Bird_ID" = rownames(ranef(Acute_Full_HT)$Bird.ID),
  "ID" = seq(1, 19, 1)
)

Acute_Collapsed <- left_join(Acute_Slopes, Acute_Birds, by = c("ID")) %>%
  dplyr::select(-ID) %>%
  group_by(Bird_ID) %>%
  summarise("Slope" = mean(Coef)) %>%
  as.data.frame(.)

# Binding in locale and treatment order

U_R <- c()
Locale <- c()

for (i in 1:nrow(Acute_Collapsed)) {
  if (is.na(Acute_Collapsed$Bird_ID[i])) {
    Locale[i] <- NA
    U_R[i] <- NA
  } else if (Acute_Collapsed$Bird_ID[i] == "BdABd") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdAO1") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABlBl") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "AYY") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAO1") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "OOA") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdABl1") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAO2") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABlO") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABdBd") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAR") {
    Locale[i] <- "Cambridge"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "RAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "ABlR") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "ARO") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdABl2") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BlAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Acute_Collapsed$Bird_ID[i] == "BdAO2") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "AOR") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "YAY") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Acute_Collapsed$Bird_ID[i] == "AYBd") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else {
    Locale[i] <- NA
    U_R[i] <- NA
  }
}
Acute_Collapsed$U_R <- factor(U_R)
Acute_Collapsed$Locale <- factor(Locale)

A_Orders <- Third_Bound %>%
  dplyr::select("Bird_ID" = Bird.ID, Pen) %>%
  na.omit(.) %>%
  distinct(.) %>%
  mutate("Order" = ifelse(Pen == "NE" | Pen == "SW", "SR", "RS")) %>%
  dplyr::select(-Pen)
Acute_Collapsed <- left_join(Acute_Collapsed, A_Orders, by = c("Bird_ID")) %>%
  mutate(Order = factor(Order), U_R = factor(U_R))

# Running formal comparison by "ANOVA" and non-linear hypothesis test (similar to above for
# acute responses).

Acute_Collapsed %>%
  group_by(U_R) %>%
  summarise(
    Mean = mean(Slope),
    LCL = quantile(Slope, 0.025, type = 8),
    UCL = quantile(Slope, 0.975, type = 8)
  )

# Note high degree of overlap.

Acute_aovBF <- anovaBF(Slope ~ Order * U_R,
  data = Acute_Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(Acute_aovBF)

# Again, minimal evidence for order effect. Dropping and re-running.

Acute_aovBF <- anovaBF(Slope ~ U_R,
  data = Acute_Collapsed,
  whichRandom = "Locale", iterations = 10000,
  rscaleFixed = "ultrawide"
)

plot(BayesFactor::posterior(Acute_aovBF, iterations = 10000)[, "mu"])
plot(acf(BayesFactor::posterior(Acute_aovBF, iterations = 10000)[, "mu"]))

# Chain is stable and peak looks clear. Note that trace is hovering surprisingly high. Summarising.

summary(Acute_aovBF)
summary(BayesFactor::posterior(Acute_aovBF, iterations = 10000)[, "mu"])

# Bayes factor is quite low (0.249). Compiling model for non-linear hypothesis test.

Slope_Contrast_Acute <- brm(Slope ~ 0 + U_R + (1 | Locale),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("gamma(1,1)", class = "sd")
  ), data = Acute_Collapsed, cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14),
  sample_prior = TRUE
)

Acute_Hypoth <- hypothesis(Slope_Contrast_Acute, c(
  "U_RRural = U_RUrban"
))
plot(Acute_Hypoth)
1 / Acute_Hypoth$hypothesis$Evid.Ratio

# Results are weaker by non-linear hypothesis test and appear quite poor by plot.
# Plotting global trends.

UR_Plot_Acute <- Acute_Collapsed %>%
  rename("Ecotype" = U_R) %>%
  group_by(Ecotype) %>%
  summarise(
    "Coefficient" = mean(Slope),
    "lcl" = quantile(Slope, 0.025, type = 8),
    "ucl" = quantile(Slope, 0.975, type = 8)
  ) %>%
  ggplot(aes(x = Ecotype, y = Coefficient, fill = Ecotype)) +
  geom_errorbar(aes(x = Ecotype, ymin = lcl, ymax = ucl),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 6, pch = 21, colour = "black") +
  #scale_fill_viridis_d(begin = 0.2, end = 0.5) +
  scale_fill_manual(values = c("wheat4","cornsilk")) + 
  theme_bw() +
  my.theme +
  theme(legend.position = "none") +
  ylab("Reaction Norm Slope")

UR_Plot_Acute

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Acute_HT.jpeg",
  UR_Plot_Acute,
  height = 5.67, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/OP/UR_Slopes_Acute_HT.pdf",
       UR_Plot_Acute,
       height = 5.67, width = 6.0, dpi = 800
)
showtext_auto(enable = FALSE)

# Plotting global slopes.

Simp_Dat_Acute <- expand.grid(
  "Amb.Temp" = c(5, 20, 35),
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = c(12),
  "Date.of.Photo" = na.omit(unique(Acute_Full$data$Date.of.Photo))[30],
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  "Pen" = na.omit(unique(Acute_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[70]
)

facet_names <- c(`5` = "5°C", `20` = "20°C", `35` = "35°C")

Simp_Curves_Filled <- Simp_Dat_Acute %>%
  add_predicted_draws(Acute_Full, n = 1000, scale = "response") %>%
  mutate("Group_ID" = paste(Bird.ID, Treatment, sep = "_")) %>%
  group_by(Group_ID, Bird.ID, Treatment, Amb.Temp, Timeline) %>%
  summarise("Pred" = mean(.prediction))

SC_F_Plot <- ggplot(Simp_Curves_Filled, aes(x = Timeline, y = Pred, colour = Treatment)) +
  geom_line(aes(group = Group_ID), alpha = 0.7) +
  annotate("rect",
    xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
    fill = "grey30", alpha = 0.8
  ) +
  facet_wrap(~Amb.Temp, labeller = as_labeller(facet_names), nrow = 3) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Time Post Stressor (s)") +
  ylab("Heat Transfer (mW)") +
  theme_bw() +
  theme(legend.position = "bottom", strip.text.x = element_text(size = 12, family = "Noto Sans")) +
  my.theme

SC_F_Plot
showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Stacked_Acute_Comp_HT.jpeg", SC_F_Plot,
  height = 13, width = 5.5, dpi = 800
)

showtext_auto(enable = FALSE)

# Surface plot

Surface_Plot <- expand.grid(
  "Maximum.Eye.Temp" = NA,
  "Date.of.Photo" = unique(Acute_Full$data$Date.of.Photo)[30],
  "Pen" = unique(Acute_Full$data$Pen)[c(1, 2)],
  "Amb.Temp" = seq(2.5, 37.5, by = 1),
  "Treatment" = unique(Acute_Full$data$Treatment)[2],
  "Bird.ID" = unique(Acute_Full$data$Bird.ID)[c(3, 15)],
  "Hour" = 12,
  "Timeline" = seq(-1200, 3600, 180),
  "Sex" = unique(Acute_Full$data$Sex),
  "Direction" = unique(Acute_Full$data$Direction),
  "Treat_Day" = unique(Acute_Full$data$Treat_Day)[40]
) %>%
  as.data.frame(.) %>%
  add_predicted_draws(Acute_Full, n = 500, scale = "response") %>%
  group_by(Amb.Temp, Timeline) %>%
  summarise("yvar" = mean(.prediction)) %>%
  as.data.frame(.)

# Pulling out legend

range(Surface_Plot$yvar)
(45 + 3)/2 -3

Leg_Lower <- Surface_Plot %>% ggplot(aes(x = Amb.Temp, y = Timeline, z = yvar)) +
  geom_raster(aes(fill = yvar)) +
  scale_fill_gradient(
    low = viridis::viridis(n = 10)[4],
    high = "ivory", breaks = c(-3.0, 21, 45.0), 
    name = "Heat Transfer\n(mW)"
  ) + 
  my.theme

Leg_Lower

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Lower_Leg_HT.jpeg", 
Leg_Lower, height = 6, width = 8.5, dpi = 800)

extractLegend <- function(xplot) {
  grobs <- ggplot_gtable(ggplot_build(xplot))
  g_title <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[g_title]]
}

ScaleLeg <- extractLegend(Leg_Lower)

min(Surface_Plot$yvar)
max(Surface_Plot$yvar)

OP_min <- -5
OP_max <- 50

OP_mybreaks <- seq(OP_min, OP_max, length.out = 15)
OP_colours <- function(x) {
  colours <- colorRampPalette(c(viridis::viridis(n = 10)[4], "ivory"))(14)
  colours[1:x]
}

OP_breaklabel <- function(x) {
  labels <- paste0(OP_mybreaks[1:14], "-", OP_mybreaks[2:15])
  labels[1:x]
}

Surface_No_Legend <- Surface_Plot %>%
  ggplot(aes(x = Amb.Temp, y = Timeline, z = yvar)) +
  geom_contour_filled(breaks = OP_mybreaks, show.legend = TRUE) +
  xlab("Ambient Temperature (°C)") +
  ylab("Time Post Stress Exposure (s)") +
  my.theme +
  scale_fill_manual(palette = OP_colours, values = OP_breaklabel(14), name = "Heat Transfer\n(mW)", drop = TRUE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank()
  )

replace_Grob <- ggplot_gtable(ggplot_build(Surface_No_Legend))
rep_new <- which(sapply(replace_Grob$grobs, function(x) x$name) == "guide-box")
replace_Grob$grobs[[rep_new]] <- ScaleLeg

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_HT.jpeg", replace_Grob,
       height = 7, width = 8.5, dpi = 800
)
showtext_auto(enable = FALSE)

# And with timeline on the x axis

Surface_No_Legend_Rotated <- Surface_Plot %>%
  ggplot(aes(x = Timeline, y = Amb.Temp, z = yvar)) +
  geom_contour_filled(breaks = OP_mybreaks, show.legend = TRUE) +
  xlab("Time Post Stress Exposure (s)") +
  ylab("Ambient Temperature (°C)") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(trans = ~., name = "", breaks = c(14,22,30), labels = c(" -","TNZ"," -"))) + 
  coord_cartesian(ylim = c(2.5, 37.5), xlim = c(-1200, 3600), clip="off") +
  annotate("segment", x=3600, y=23, xend=3600, yend=29,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=3600, y=21, xend=3600, yend=15,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) + 
  scale_fill_manual(palette = OP_colours, values = OP_breaklabel(14), name = "Heat Transfer\n(mW)", drop = TRUE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_text(hjust = -0.6, size = 16),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + my.theme

cow_legend <- cowplot::get_legend(Leg_Lower)
cow_grid <- cowplot::plot_grid(plotlist = list(Surface_No_Legend_Rotated), ncol = 1)
to_save = cowplot::plot_grid(cow_grid, cow_legend, ncol = 2, rel_widths = c(12,1.5), align = "hv")

to_save

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed_HT.jpeg", to_save,
       height = 7, width = 9, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed_HT.pdf", to_save,
       height = 7, width = 9, dpi = 800
)

# With grey rectangle over stress exposure period.

Surface_No_Legend_Rotated <- Surface_Plot %>%
  ggplot(aes(x = Timeline, y = Amb.Temp, z = yvar)) +
  geom_contour_filled(breaks = OP_mybreaks, show.legend = TRUE) +
  xlab("Time Post Stress Exposure (s)") +
  ylab("Ambient Temperature (°C)") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(trans = ~., name = "", breaks = c(14,22,30), labels = c(" -","TNZ"," -"))) + 
  coord_cartesian(ylim = c(2.5, 37.5), xlim = c(-1200, 3600), clip="off") +
  annotate("segment", x=3600, y=23, xend=3600, yend=29,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=3600, y=21, xend=3600, yend=15,
           col="black", arrow=arrow(length=unit(0.3, "cm"))) + 
  annotate("rect",
           xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
           fill = "grey30", alpha = 0.8
  ) + 
  scale_fill_manual(palette = OP_colours, values = OP_breaklabel(14), name = "Heat Transfer\n(mW)", drop = TRUE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_text(hjust = -0.6, size = 16),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + my.theme

cow_legend <- cowplot::get_legend(Leg_Lower)
cow_grid <- cowplot::plot_grid(plotlist = list(Surface_No_Legend_Rotated), ncol = 1)
to_save = cowplot::plot_grid(cow_grid, cow_legend, ncol = 2, rel_widths = c(12,1.5), align = "hv")

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed_Rectangle_HT.jpeg", to_save,
       height = 7, width = 9, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Contour_Reversed_Rectangle_HT.pdf", to_save,
       height = 7, width = 9, dpi = 800
)
showtext_auto(enable = FALSE)

# Finally, re-running formal comparisons between repeatabilities of treatment groups. Note that 
# the below comparisons are directed towards differences in means and not distributions. Such comparisons
# are likely to be quite liberal.

Rep_Chronic_HT <- Rep_Treats %>%
  mutate(Model = gsub("\n.*", "", Rep_Treats$Mod_Type)) %>%
  rename(
    "Control" = C_Repeatability,
    "Stress" = T_Repeatability
  ) %>%
  filter(Model == "True") %>%
  dplyr::select(.chain, .iteration, Control, Stress) %>%
  gather(Treatment, Repeat, Control:Stress, factor_key = TRUE)

Rep_Acute_HT <- Rep_Treats_Acute %>%
  mutate(Model = gsub("\n.*", "", Rep_Treats$Mod_Type)) %>%
  rename(
    "Control" = C_Repeatability,
    "Stress" = T_Repeatability
  ) %>%
  filter(Model == "True") %>%
  dplyr::select(.chain, .iteration, Control, Stress) %>%
  gather(Treatment, Repeat, Control:Stress, factor_key = TRUE)

{
Treat_Comp_Robust_Acute <- brm(
  bf(Repeat ~ 0 + Treatment, sigma ~ 0 + Treatment),
  family = student,
  data = Rep_Acute_HT,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101, 
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ), 
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  file = "/home/joshk/Desktop/CC_Results/Final/Treat_Comp_Robust_Acute.Rds"
) 
  
Treat_Comp_Robust_Chronic <- brm(
  bf(Repeat ~ 0 + Treatment, sigma ~ 0 + Treatment),
  family = student,
  data = Rep_Chronic_HT,
  chains = 4, cores = 4, iter = 10000, warmup = 1000,
  thin = 10, seed = 101,
  prior = c(
    set_prior("beta(1, 3)", class = "b", lb = 0.001, ub = 1),
    set_prior("exponential(0.1)", class = "nu")
  ),
  file = "/home/joshk/Desktop/CC_Results/Final/Treat_Compare_Robust_Chronic.Rds"
)
}
# Checking outcomes

plot(Treat_Comp_Robust_Acute)
plot(Treat_Comp_Robust_Chronic)
pp_check(Treat_Comp_Robust_Acute, nsamples = 500)
pp_check(Treat_Comp_Robust_Chronic, nsamples = 500)

# Strikingly high nu parameters, but posterior samples fit well. Again, these estimates should be liberal.

stan_ac(Treat_Comp_Robust_Acute$fit)
mcmc_neff(neff_ratio(Treat_Comp_Robust_Acute), size = 2)
mcmc_rhat(rhat(Treat_Comp_Robust_Acute))

stan_ac(Treat_Comp_Robust_Chronic$fit)
mcmc_neff(neff_ratio(Treat_Comp_Robust_Chronic), size = 2)
mcmc_rhat(rhat(Treat_Comp_Robust_Chronic))

# Fair chain convergence. Plotting fit against raw.

Pred_A <- Rep_Acute_HT %>%
  add_predicted_draws(Treat_Comp_Robust_Acute, class = "response") %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(.prediction),
    "UCL" = quantile(.prediction, 0.975, type = 8),
    "LCL" = quantile(.prediction, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

Raw_A <- Rep_Acute_HT %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(Repeat),
    "UCL" = quantile(Repeat, 0.975, type = 8),
    "LCL" = quantile(Repeat, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Raw")

grid.arrange(Pred_A, Raw_A, nrow = 2)

# Acceptably comparable. Checking whether heteroskedasticity was captured.

Rep_Acute_HT %>%
  add_residual_draws(Treat_Comp_Robust_Acute, class = "response") %>%
  ggplot(aes(x = Treatment, y = .residual, fill = Treatment)) +
  geom_boxplot(colour = "black") +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

# Seems loosely corrected and acceptable for now. Chronic model?

Pred_C <- Rep_Chronic_HT %>%
  add_predicted_draws(Treat_Comp_Robust_Chronic, class = "response") %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(.prediction),
    "UCL" = quantile(.prediction, 0.975, type = 8),
    "LCL" = quantile(.prediction, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

Raw_C <- Rep_Chronic_HT %>%
  group_by(Treatment) %>%
  summarise(
    "Pred" = mean(Repeat),
    "UCL" = quantile(Repeat, 0.975, type = 8),
    "LCL" = quantile(Repeat, 0.025, type = 8)
  ) %>%
  ggplot(aes(x = Treatment, y = Pred, fill = Treatment)) +
  geom_errorbar(aes(x = Treatment, ymin = LCL, ymax = UCL),
    colour = "black", size = 1, width = 0.3
  ) +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Raw")

grid.arrange(Pred_C, Raw_C, nrow = 2)

# Great. Again, checking heteroskedasticity.

Rep_Chronic_HT %>%
  add_residual_draws(Treat_Comp_Robust_Chronic, class = "response") %>%
  ggplot(aes(x = Treatment, y = .residual, fill = Treatment)) +
  geom_boxplot(colour = "black") +
  geom_point(size = 4, pch = 21, colour = "black") +
  scale_fill_manual(values = c("slateblue", "navajowhite")) +
  theme_bw() +
  annotate("text", x = 1, y = 0.8, label = "Predicted")

# No concern.
# Looking at model summaries

summary(Treat_Comp_Robust_Acute)
summary(Treat_Comp_Robust_Chronic)

hypothesis(Treat_Comp_Robust_Acute, "TreatmentStress - TreatmentControl > 0.1", alpha = 0.05)
hypothesis(Treat_Comp_Robust_Chronic, "TreatmentStress - TreatmentControl > 0.1", alpha = 0.05)

# KEEP 

#robust_comp_acute$model = gsub("real Intercept_sigma;","real<lower=0.001> Intercept_sigma;",robust_comp_acute$model)
#new_standata = standata(robust_comp_acute)
#robust_comp_acute$fit = stan(model_code = robust_comp_acute$model, data = new_standata,
#                             chains = 4, iter = 10000, warmup = 1000, thin = 10,
#                             seed = 101, control = list(adapt_delta = 0.98, max_treedepth = 14))

## Test raw data point plot
# Note that below, mean temperature observed within each temperature grouping is used for prediction purposes.

Mean_Cold = mean(subset(Acute_Full$data, Amb.Temp <= 13.5)$Amb.Temp)
Modal_Date = subset(Acute_Full$data, Amb.Temp <= 13.5) %>% 
    group_by(Date.of.Photo) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(Count)) %>% 
    pull(Date.of.Photo) %>% head(1)

Simp_Dat_Acute <- expand.grid(
  "Amb.Temp" = Mean_Cold,
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = c(8,12,16),
  "Date.of.Photo" = Modal_Date,
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  "Pen" = na.omit(unique(Acute_Full$data$Pen)),
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[30]
)

Pred = posterior_predict(Acute_Full, newdata = Simp_Dat_Acute, nsamples = 1000, type = "response")
Simp_Dat_Acute = Simp_Dat_Acute %>% mutate("Pred" = colMeans(Pred), "Temp_Group" = "Low") 

Base_Plot = Simp_Dat_Acute %>%
  mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_")) %>%
  group_by(Bird.ID, Timeline, Treatment, Treat_Bird) %>% 
  summarise("mW" = mean(Pred)) %>% 
  ggplot(aes(x = Timeline, y = mW, colour = Treatment, fill = Treatment)) + 
  geom_line(size = 1, aes(group = Treat_Bird)) + 
  theme_bw() + 
  scale_fill_viridis_d(begin = 0.3, end = 0.6) + 
  scale_colour_viridis_d(begin = 0.3, end = 0.6)

PDat = Third_Bound %>% 
    mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_")) %>%
    filter(Amb.Temp <= 13.5)

Base_Plot + stat_summary_bin(data = PDat, geom = "point", fun = "mean",
  binwidth = 120, mapping = aes(x = Timeline, y = mW,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2)

## Good! Expanding to include data points for other temperature regions.   

Mean_Mid = mean(subset(Acute_Full$data, Amb.Temp >= 14 & Amb.Temp <= 29.5)$Amb.Temp)
Modal_Date_2 = subset(Acute_Full$data, Amb.Temp >= 14 & Amb.Temp <= 29.5) %>% 
    group_by(Date.of.Photo) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(Count)) %>% 
    pull(Date.of.Photo) %>% head(1)

Mid_Temp <- expand.grid(
  "Amb.Temp" = Mean_Mid,
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = c(8,12,16),
  "Date.of.Photo" = Modal_Date_2,
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  "Pen" = na.omit(unique(Acute_Full$data$Pen)),
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[30]
)
Mid_Pred = posterior_predict(Acute_Full, newdata = Mid_Temp, nsamples = 1000, type = "response")
Mid_Temp = Mid_Temp %>% mutate("Pred" = colMeans(Mid_Pred), "Temp_Group" = "Mid")

Mean_Hot = mean(subset(Acute_Full$data, Amb.Temp >= 30.5)$Amb.Temp)
Modal_Date_3 = subset(Acute_Full$data, Amb.Temp >= 30.5) %>% 
    group_by(Date.of.Photo) %>% 
    summarise(Count = n()) %>% 
    arrange(desc(Count)) %>% 
    pull(Date.of.Photo) %>% head(1)

High_Temp <- expand.grid(
  "Amb.Temp" = Mean_Hot,
  "Timeline" = seq(-1200, 3600, by = 180),
  "Hour" = c(8,12,16),
  "Date.of.Photo" = Modal_Date_3,
  "Treatment" = na.omit(unique(Acute_Full$data$Treatment)),
  "Pen" = na.omit(unique(Acute_Full$data$Pen)),
  "Sex" = na.omit(unique(Acute_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(Acute_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(Acute_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(Acute_Full$data$Treat_Day))[30]
)

High_Pred = posterior_predict(Acute_Full, newdata = High_Temp, nsamples = 1000, type = "response")
High_Temp = High_Temp %>% mutate("Pred" = colMeans(High_Pred), "Temp_Group" = "High")

Grouped = rbind(Simp_Dat_Acute, Mid_Temp, High_Temp) %>% 
  mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_")) %>%
  group_by(Temp_Group, Bird.ID, Timeline, Treatment, Treat_Bird) %>% 
  summarise("mW" = mean(Pred))

write.csv(Grouped, "/home/joshk/git_repositories/BCCH_IndVar/AcutePlotDat_HT.csv", row.names = FALSE)
Grouped = read.csv("/home/joshk/git_repositories/BCCH_IndVar/AcutePlotDat_HT.csv")

# Global plot

Plot_Dat = Third_Bound %>% drop_na(mW) %>% 
mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_"),
"Temp_Group" = ifelse(Amb.Temp <= 13.5, "Low", ifelse(Amb.Temp >= 30.5, "High", "Mid")))
Plot_Dat$Temp_Group = factor(Plot_Dat$Temp_Group, levels = c("Low", "Mid", "High"))

Temp_Groups = c(`Low` = "< TNZ", `Mid` = "TNZ", `High` = "> TNZ")
Grouped = Grouped %>% mutate(Temp_Group = factor(Temp_Group, levels = c("Low", "Mid", "High")))

ggplot(Grouped, aes(x = Timeline, y = mW, colour = Treatment, fill = Treatment)) + 
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 60, mapping = aes(x = Timeline, y = mW,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.8) + 
  geom_line(size = 1, aes(group = Treat_Bird)) + 
  theme_bw() + 
  scale_fill_viridis_d(begin = 0.3, end = 0.6) + 
  scale_colour_viridis_d(begin = 0.3, end = 0.6) + 
  facet_wrap(~Temp_Group, labeller = as_labeller(Temp_Groups))

# Slightly shifted fit, but acceptable.

Stacked_Acute_Dots = ggplot(Grouped, aes(x = Timeline, y = mW, colour = Treatment, fill = Treatment)) + 
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 180, mapping = aes(x = Timeline, y = mW,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(size = 1, aes(group = Treat_Bird), alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  facet_wrap(~Temp_Group, ncol = 1, labeller = as_labeller(Temp_Groups)) + 
  ylab("Heat Transfer Rate (mW)") + 
  xlab("Time Post Stess Exposure (s)") + 
    annotate("rect",
    xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
    fill = "grey30", alpha = 0.8
  ) + my.theme

Stacked_Acute_Dots
showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Stacked_Acute_Dots_HT.jpeg", Stacked_Acute_Dots,
  height = 13, width = 5.5, dpi = 800
)

# Wide format

Stacked_Acute_Dots = ggplot(Grouped, aes(x = Timeline, y = mW, colour = Treatment, fill = Treatment)) + 
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
                   binwidth = 180, mapping = aes(x = Timeline, y = mW,
                                                 group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(size = 1, aes(group = Treat_Bird), alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
  facet_wrap(~Temp_Group, ncol = 3, labeller = as_labeller(Temp_Groups)) + 
  ylab("Heat Transfer Rate (mW)") + 
  xlab("Time Post Stess Exposure (s)") + 
  annotate("rect",
           xmin = 0, xmax = 1200, ymin = -Inf, ymax = Inf, colour = "black",
           fill = "grey30", alpha = 0.8
  ) + my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Wide_Acute_Dots_HT.jpeg", Stacked_Acute_Dots,
       height = 4.25, width = 10, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Wide_Acute_Dots_HT.pdf", Stacked_Acute_Dots,
       height = 4.25, width = 10, dpi = 800
)

# Plotting chronic responses with dots.

Descriptive_Chronic = expand.grid(
  "s.Amb.Temp" = c(min(CC_Full$data$s.Amb.Temp), max(CC_Full$data$s.Amb.Temp)),
  "Hour" = c(8, 12, 16),
  "Sex" = na.omit(unique(CC_Full$data$Sex)),
  "Pen" = na.omit(unique(CC_Full$data$Pen)),
  "Date.of.Photo" = na.omit(unique(CC_Full$data$Date.of.Photo))[1],
  "Direction" = na.omit(unique(CC_Full$data$Direction)),
  "Treatment" = na.omit(unique(CC_Full$data$Treatment)),
  "Bird.ID" = na.omit(unique(CC_Full$data$Bird.ID)),
  "Treat_Day" = na.omit(unique(CC_Full$data$Treat_Day))[1]
)

Descriptive_Chronic_Pred = posterior_predict(CC_Full, newdata = Descriptive_Chronic, re_formula = NA,
  nsamples = 1000, summary = FALSE)
Descriptive_Chronic = Descriptive_Chronic %>% mutate("Pred" = colMeans(Descriptive_Chronic_Pred))

Chronic_Grouped = Descriptive_Chronic %>% 
  mutate("Treat_Bird" = paste(Treatment, Bird.ID, sep = "_"),
  "Amb.Temp" = (s.Amb.Temp * (2 * sd(Third_Bound$Amb.Temp, na.rm = T))) +
      mean(Third_Bound$Amb.Temp, na.rm = T)) %>% 
  group_by(Treatment, Bird.ID, Amb.Temp, Treat_Bird) %>% 
  summarise(mW = mean(Pred))

Chronic_Curves_Dots <- 
  ggplot(Chronic_Grouped, aes(x = Amb.Temp, y = mW, colour = Treatment)) +
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 1, mapping = aes(x = Amb.Temp, y = mW,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Heat Transfer Rate (mW)") +
  theme_bw() +
  my.theme

Chronic_Curves_Dots

# Widening ambient temperature spectrum

Simp_Dat_Chronic <- expand.grid(
  "s.Amb.Temp" = seq(min(CC_Full$data$s.Amb.Temp, na.rm = TRUE),
    max(CC_Full$data$s.Amb.Temp, na.rm = TRUE),
    length.out = 25
  ),
  "Hour" = c(12),
  "Date.of.Photo" = na.omit(unique(CC_Full$data$Date.of.Photo))[30],
  "Treatment" = na.omit(unique(CC_Full$data$Treatment)),
  "Pen" = na.omit(unique(CC_Full$data$Pen))[c(1, 2)],
  "Sex" = na.omit(unique(CC_Full$data$Sex)),
  "Bird.ID" = na.omit(unique(CC_Full$data$Bird.ID)),
  "Direction" = na.omit(unique(CC_Full$data$Direction)),
  "Treat_Day" = na.omit(unique(CC_Full$data$Treat_Day))[c(14, 70)]
)

Descriptive_Chronic_Pred = posterior_predict(CC_Full, newdata = Simp_Dat_Chronic,
  nsamples = 1000, summary = FALSE)

Simp_Dat_Chronic = Simp_Dat_Chronic %>% 
  mutate("Pred" = colMeans(Descriptive_Chronic_Pred),
  "Treat_Bird" = paste(Bird.ID, Treatment, sep = "_")
  ) %>%
  group_by(Treat_Bird, Bird.ID, Treatment, s.Amb.Temp) %>%
  summarise("mW" = mean(Pred))

Chronic_Curves_Dots <- 
  ggplot(Simp_Dat_Chronic, aes(x = s.Amb.Temp, y = mW, colour = Treatment)) +
  stat_summary_bin(data = Plot_Dat, geom = "point", fun = "mean",
  binwidth = 0.01, mapping = aes(x = s.Amb.Temp, y = mW,
  group = Treat_Bird, fill = Treatment), pch = 21, size = 2, alpha = 0.3) + 
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  annotate(geom = "rect", xmin = -0.70139, xmax = 0.40972, ymin = -Inf, ymax = Inf,
           fill = "grey30", alpha = 0.5, colour = "black") + 
  annotate(geom = "text", x = -0.14583, y = -5, label = "TNZ", size = 8, family = "Noto Sans") + 
  #annotate(geom = "rect", xmin = 14, xmax = 30, ymin = -Inf, ymax = Inf,
  #         fill = "grey70", alpha = 0.5, colour = "black") + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Heat Transfer Rate (mW)") +
  theme_bw() +
  my.theme  

Chronic_Curves_Dots

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_TNZ_HT.jpeg",
       Chronic_Curves_Dots,
       height = 6.0, width = 7.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_TNZ_HT.pdf",
       Chronic_Curves_Dots,
       height = 6.0, width = 7.0, dpi = 800
)
showtext_auto(enable = FALSE)

# No dots 

Chronic_Curves_Nodots <- 
  ggplot(Simp_Dat_Chronic, aes(x = s.Amb.Temp, y = mW, colour = Treatment)) +
  geom_line(aes(group = Treat_Bird), alpha = 0.7) +
  annotate(geom = "rect", xmin = -0.70139, xmax = 0.40972, ymin = -Inf, ymax = Inf,
           fill = "grey30", alpha = 0.5, colour = "black") + 
  annotate(geom = "text", x = -0.14583, y = -5, label = "TNZ", size = 8, family = "Noto Sans") + 
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-1.5, 1, length.out = 7), 
                     labels = seq(2.5, 38.5, by = 6)) + 
  scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Heat Transfer Rate (mW)") +
  theme_bw() +
  my.theme

showtext_opts(dpi = 800)
showtext_auto()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_HT_NoDots.jpeg",
       Chronic_Curves_Nodots,
       height = 6.0, width = 7.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Chronic_Curves_HT_NoDots.pdf",
       Chronic_Curves_Nodots,
       height = 6.0, width = 7.0, dpi = 800
)

showtext_auto(enable = FALSE)

# Drawing distributions for methods figure

require(latex2exp)

Sample = rnorm(1e+4, mean = 0, sd = 1)

Rep_Curve = ggplot(as.data.frame(Sample), aes(x = Sample)) +
  geom_density(adjust = 3, size = 3, colour = "firebrick4", alpha = 0.5) + 
  annotate("text", x = 0.03, y = 0.1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram1.jpg",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
)
 
ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram1.pdf",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

# With colour 

Rep_Curve = ggplot(as.data.frame(Sample), aes(x = Sample)) +
  geom_density(adjust = 3, size = 1, colour = "black", fill = "firebrick4", alpha = 0.5) + 
  annotate("text", x = 0.03, y = 0.1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour1.jpg",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
)

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour1.pdf",
       Rep_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

# Adding null curve

Sample = rnorm(1e+4, mean = 1, sd = 0.6)
Beta_Sample = rbeta(1e+4, 1, 6)

Combined_Curve = ggplot(as.data.frame(Sample), aes(x = Sample + 0.5)) +
  geom_density(adjust = 3, size = 3, colour = "firebrick4", alpha = 0.5) + 
  geom_density(data = as.data.frame(Beta_Sample), aes(x = Beta_Sample),
               colour = "black", alpha = 0.5,
               size = 3, adjust = 3) + 
  xlim(c(0,4)) + 
  annotate("text", x = 1.3, y = 2, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{N}} \\cdot (\\mu_{1_{N}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("text", x = 3, y = 1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("segment", x = 1, y = 1.8, xend = 0.5, yend = 1, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) + 
  annotate("segment", x = 2.4, y = 0.8, xend = 2.2, yend = 0.5, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) +
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram2.jpg",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram2.pdf",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

# Again, with colour 

Combined_Curve = ggplot(as.data.frame(Sample), aes(x = Sample + 0.5)) +
  geom_density(adjust = 3, size = 1, colour = "black", fill = "firebrick4", alpha = 0.5) + 
  geom_density(data = as.data.frame(Beta_Sample), aes(x = Beta_Sample),
               fill = "black", colour = "black", alpha = 0.5,
               size = 1, adjust = 3) + 
  xlim(c(0,4)) + 
  annotate("text", x = 1.3, y = 2, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{N}} \\cdot (\\mu_{1_{N}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("text", x = 3, y = 1, size = 8, family = "Noto Sans",
           label = TeX('$\\mu_{1_{T}} \\cdot (\\mu_{1_{T}}+\\sigma_{_{T}}^{2})^{-1}$')) + 
  annotate("segment", x = 1, y = 1.8, xend = 0.5, yend = 1, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) + 
  annotate("segment", x = 2.4, y = 0.8, xend = 2.2, yend = 0.5, size = 1,
           arrow = arrow(length = unit(0.1, "inches"))) +
  theme_void()

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour2.jpg",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

ggsave("/home/joshk/git_repositories/BCCH_IndVar/Figures/Curve_Diagram_Colour2.pdf",
       Combined_Curve,
       height = 6.0, width = 6.0, dpi = 800
) 

##########################################################################
### Testing evidence of divergent traits between urban and rural birds ###
##########################################################################

get_variables(CC_Full)

Comp_Slopes <- CC_Full %>%
  spread_draws(s_ss.Amb.TempBird.IDTreatmentStress_1[ID],
  s_ss.Amb.TempBird.ID_1[ID]) %>%
  rename("Stress_B" = s_ss.Amb.TempBird.IDTreatmentStress_1,
  "Control_B" = s_ss.Amb.TempBird.ID_1) %>%
  as.data.frame(.)

Birds <- data.frame(
  "Bird_ID" = rownames(ranef(CC_Full)$Bird.ID),
  "ID" = seq(1, 19, 1)
)

Comp_Slopes <- left_join(Comp_Slopes, Birds, by = c("ID")) %>%
  dplyr::select(-ID) %>%
  group_by(Bird_ID) %>%
  summarise("Stress_B" = mean(Stress_B),
  "Control_B" = mean(Control_B)) %>%
  as.data.frame(.)

U_R <- c()
Locale <- c()

for (i in 1:nrow(Comp_Slopes)) {
  if (is.na(Comp_Slopes$Bird_ID[i])) {
    Locale[i] <- NA
    U_R[i] <- NA
  } else if (Comp_Slopes$Bird_ID[i] == "BdABd") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "BdAO1") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "ABlBl") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "AYY") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "YAO1") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "OOA") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "BdABl1") {
    Locale[i] <- "Ruthven.Park"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "YAO2") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "ABlO") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "ABdBd") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "YAR") {
    Locale[i] <- "Cambridge"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "RAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "ABlR") {
    Locale[i] <- "Brantford"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "ARO") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "BdABl2") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "BlAR") {
    Locale[i] <- "Corwhin"
    U_R[i] <- "Rural"
  } else if (Comp_Slopes$Bird_ID[i] == "BdAO2") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "AOR") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "YAY") {
    Locale[i] <- "Guelph"
    U_R[i] <- "Urban"
  } else if (Comp_Slopes$Bird_ID[i] == "AYBd") {
    Locale[i] <- "Erin"
    U_R[i] <- "Rural"
  } else {
    Locale[i] <- NA
    U_R[i] <- NA
  }
}
Comp_Slopes$U_R <- factor(U_R)
Comp_Slopes$Locale <- factor(Locale)

Orders <- Third_Bound %>%
  dplyr::select("Bird_ID" = Bird.ID, Pen) %>%
  na.omit(.) %>%
  distinct(.) %>%
  mutate("Order" = ifelse(Pen == "NE" | Pen == "SW", "SR", "RS")) %>%
  dplyr::select(-Pen)

Comp_Slopes <- left_join(Comp_Slopes, Orders, by = c("Bird_ID")) %>%
  mutate(Order = factor(Order), U_R = factor(U_R))

Model = brm(Stress_B ~ Control_B*U_R + (1|Locale), data = Comp_Slopes,
  cores = 4, chains = 4,
  seed = 100, refresh = 0, family = "gaussian",
  iter = 10000, warmup = 1000, thin = 10,
  control = list(adapt_delta = 0.98, max_treedepth = 14))

plot(Model)
pp_check(Model, nsamples = 500)

stan_ac(Model$fit)
mcmc_neff(neff_ratio(Model), size = 2)
mcmc_rhat(rhat(Model))

require(modelr)

Comp_Slopes %>%
  group_by(U_R) %>%
  data_grid(Control_B = seq_range(Control_B, n = 101)) %>%
  slice(rep(row_number(),6)) %>% 
  mutate("Locale" = c(rep("Guelph", 101), rep("Ruthven.Park", 101),
      rep("Brantford", 101), rep("Corwhin", 101),
      rep("Erin", 101), rep("Cambridge", 101))
  ) %>% 
  mutate(Locale = factor(Locale)) %>%
  add_fitted_draws(Model, n = 100) %>%
  ggplot(aes(x = Control_B, y = Stress_B, color = ordered(U_R))) +
  geom_line(aes(y = .value, group = paste(U_R, Locale, .draw)), alpha = .1) +
  geom_point(data = Comp_Slopes) +
  scale_color_brewer(palette = "Dark2", name = "Ecotype") + 
  theme_bw() + my.theme + theme(panel.grid.major = element_blank()) + 
  xlab("Control Slope") + 
  ylab("Stress Slope")

# Interesting difference in slopes, although messy. Adjusting axes to be absolute and not relative.

CC_Full %>% 
  spread_draws(bs_ss.Amb.Temp_1,
  `bs_ss.Amb.Temp:TreatmentStress_1`,
  b_Intercept) %>% 
  rename("Control_B" = bs_ss.Amb.Temp_1,
  "Stress_B" = `bs_ss.Amb.Temp:TreatmentStress_1`,
  "Intercept" = b_Intercept) %>%
  summarise(Control_B = mean(Control_B),
  Stress_B = mean(Stress_B),
  Intercept = mean(Intercept))

Control_Data = data.frame("Temp" = rep(seq_range(c(-1.2,1.2), 100)),
  "Treatment" = rep("Control", 100)) %>% 
  mutate(Eye_Temp = Temp*2.58 + 32.9)

Stress_Data = data.frame("Temp" = rep(seq_range(c(-1.2,1.2), 100)),
  "Treatment" = rep("Stress", 100)) %>% 
  mutate(Eye_Temp = Temp*2.58*1.14 + 32.9)

All = rbind(Control_Data, Stress_Data)

ggplot(All, aes(x = Temp, y = Eye_Temp, linetype = Treatment)) + 
  geom_line()

# Good. Generic linear responses prooduce. Next step is adding individual-level adjustments.

Temp = rep(seq_range(c(-1.2,1.2), 100))
Eye_Temps_Control = vector('list', 19)
Eye_Temps_Stress = vector('list', 19)

for (i in 1:nrow(Comp_Slopes)){
  Eye_Temps_Control[[i]] = Eye_Temp = c(32.9 + Temp*(2.58 + Comp_Slopes$Control_B[i]))
  Eye_Temps_Stress[[i]] = Eye_Temp = c(32.9 + Temp*2.58*(1.14 + Comp_Slopes$Stress_B[i]))
}

Bound = vector('list', 19)

for (i in 1:19){
  Bound[[i]] = data.frame("Bird.ID" = rep(Comp_Slopes$Bird_ID[i],200),
  "Temp" = rep(Temp, 2),
  "Treatment" = c(rep("Control", 100), rep("Stress", 100)),
  "Eye_Temp" = c(Eye_Temps_Control[[i]], Eye_Temps_Stress[[i]])
  )
}

All = bind_rows(Bound)
All$Group_ID = with(All, paste(Bird.ID, Treatment, sep = "_"))

ggplot(All, aes(x = Temp, y = Eye_Temp, linetype = Treatment, colour = Treatment)) + 
  geom_line(aes(group = Group_ID)) +
  scale_colour_manual(values = c(viridis::viridis(n = 10)[4], "black")) +
  theme(panel.grid.major = element_blank()) +
  xlab("Ambient Temperature (°C)") +
  ylab("Maximum Eye Temperature (°C)") +
  theme_bw() +
  my.theme  

# Perhaps slightly too variable.

Comp_Slopes %>%
  group_by(U_R) %>%
  data_grid(Control_B = seq_range(Control_B, n = 101)) %>%
  slice(rep(row_number(),6)) %>% 
  mutate("Locale" = c(rep("Guelph", 101), rep("Ruthven.Park", 101),
      rep("Brantford", 101), rep("Corwhin", 101),
      rep("Erin", 101), rep("Cambridge", 101))
  ) %>% 
  mutate(Locale = factor(Locale)) %>%
  add_fitted_draws(Model, n = 100) %>%
  ggplot(aes(x = Control_B, y = Stress_B, color = ordered(U_R))) +
  geom_line(aes(y = .value, group = paste(U_R, Locale, .draw)), alpha = .1) +
  geom_point(data = Comp_Slopes) +
  scale_color_brewer(palette = "Dark2", name = "Ecotype") + 
  theme_bw() + my.theme + theme(panel.grid.major = element_blank()) + 
  xlab("Control Slope") + 
  ylab("Stress Slope") + 
  scale_y_continuous( 
    breaks = c(-7.15, -5.15, -3.14, -1.14, 1.14, 3.14, 5.14, 7.14), 
    labels = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
  scale_x_continuous( 
    breaks = c(-7.15, -5.15, -3.14, -1.14, 1.14, 3.14, 5.14, 7.14), 
    labels = c(-6, -4, -2, 0, 2, 4, 6, 8)) +  
  geom_hline(yintercept = -1.14, size = 0.5, colour = "black", linetype = "dashed") +
  annotate("text", x = -1, y = -7.1, label = "Negative Response", size = 5, family = "Noto Sans") + 
  annotate("text", x = 0.95, y = 7, label = "Positive Response", size = 5, family = "Noto Sans")
