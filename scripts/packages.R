# packages.R
# All packages used across project scripts.
# Source this file at the start of a session or from individual scripts.

library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(patchwork)

# Modelling
library(survival)
library(survminer)
library(flexsurv)
library(glmmTMB)
library(emmeans)
library(performance)

# Reporting / tables
library(broom)
library(gtsummary)
library(sjPlot)
library(ggpubr)
