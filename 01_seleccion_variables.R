########################################################
# Trabajo de grado comparación CoxPH vs SRF vs GBS #
########################################################

# Librerias necesarias

library(skimr)
library(tidyverse)
library(dplyr)
library(purrr)
library(ggcorrplot)
library(nnet)
library(truncnorm)
library(MASS)
library(VIM)
library(knitr)
library(tibble)
library(AER)
library(pROC)
library(ResourceSelection)
library(car)
library(AER)
library(survival)
library(stringr)
library(pec)
library(prodlim)
library(purrr)
library(tibble)
library(randomForestSRC)
library(xgboost)
library(ggplot2)
library(forcats)
library(readr)

############################################################
# Carga de Datos

datos <- read_csv("METABRIC_RNA_Mutation.csv")

###########################################################

# 1. Selección de variables elegidas

datos1 <- datos %>% 
  dplyr::select(age_at_diagnosis, type_of_breast_surgery, 
                chemotherapy, neoplasm_histologic_grade, 
                lymph_nodes_examined_positive, her2_status, tumor_size, 
                overall_survival_months, overall_survival )


