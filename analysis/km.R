# survival analysis


# setup -------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
conflict_prefer("filter", "dplyr")


# data format -------------------------------------------------------------

df <-
  read_csv("~/Dropbox (Partners HealthCare)/lactate_detritus/km.csv") %>% 
  mutate(
    time = as.numeric(date - as_date("2021-07-29")), 
    event = if_else(event == "died", 1, 0)
  )


# early deaths ------------------------------------------------------------

df <- 
  df %>% 
  mutate(
    treatment = replace(treatment, time < 8 & disease == "Bleo", "Control"), 
    treatment = factor(treatment, levels = c("Control", "AZD3965", "VB124")), 
    disease = factor(disease, levels = c("Control", "Bleo")), 
    group = interaction(disease, treatment),
  ) %>% 
  filter(disease == "Bleo")


# survival analysis -------------------------------------------------------

f1 <- survfit(Surv(time, event) ~ 1, data = df)
ggsurvplot(
  fit = f1, 
  xlab = "Days", 
  ylab = "Survival Probability"
)

f2 <- survfit(Surv(time, event) ~ (disease + treatment), data = df)
ggsurvplot(
  fit = f2, 
  xlab = "Days", 
  ylab = "Survival Probability"
)


# statistics --------------------------------------------------------------

pairwise_survdiff(Surv(time, event) ~ treatment, data = df)
coxph(Surv(time, event) ~ treatment, data = df)
