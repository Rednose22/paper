#####################################################################
##                        Program set up                           ##
#####################################################################
library(multinma)
library(tidyverse)
library(lazyeval)
library(data.table)
library(kableExtra)
library(rlist) ## save the list as rdata
library(qdapRegex)
library(gridExtra)
library(lemon)

options(mc.cores = parallel::detectCores())

nc <- switch(tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_")), 
             "true" =, "warn" = 2, 
             parallel::detectCores())
options(mc.cores = nc)

setwd("D:\\MPH-Project\\Further work\\paper")
memory.limit(size=40000) # set memory

#####################################################################
##          Redistribute the outcome to increase                   ##
##          treatment effects in  population                ##  
#####################################################################

## checking the missing data
full_ipd <- plaque_psoriasis_ipd %>%
  mutate(   # Check complete cases for covariates of interest
    complete = complete.cases(durnpso, prevsys, bsa, weight, psa, age, male))

# remove the missing data
full_ipd <- filter(full_ipd, complete) %>%
  mutate(cat_pasi = ifelse(pasi75 == 1, "yes", "no"),
         sub_age = ifelse(age>=46, "gt or eq 46", "lt 46"),
         durncat = ifelse(durnpso < 17, "lt 17", "ge 17"))

#####################################################################
##                          Data preparation                       ##
#####################################################################

## Data transformation in IPD
sub_ipd <- full_ipd %>% 
  mutate(# Variable transformations
    bsa = bsa / 100,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    weight = weight / 10,
    age = age / 10,
    durnpso = durnpso / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6, 7) ~ "IL blocker",
                         trtn == 4 ~ "TNFa blocker"),
    sub_durn = ifelse(durnpso < 1.7, "1", "2"),
    male = ifelse(male == "TRUE", 1, 0))

## subgroup result generated from IPD
## five covariate durnpso, prevsys, bsa, weight, psa, age
## include studyc, trt, trtclass, trtn, trtc_long into by vars

sub_eff_func <- function(group_vars, name_) {
  
  sub_ipd %>% 
    group_by(!!!group_vars) %>%
    summarise(age_mean = mean(age),
              age_sd = sd(age),
              bmi_mean = mean(bmi, na.rm = TRUE),
              bmi_sd = sd(bmi, na.rm = TRUE),
              weight_mean = mean(bmi, na.rm = TRUE),
              weight_sd = sd(bmi, na.rm = TRUE),
              durnpso_mean = mean(durnpso),
              durnpso_sd = sd(durnpso),
              pasi75_n = n(),
              pasi75_r = sum(pasi75),
              pasi90_n = n(),
              pasi90_r = sum(pasi90),
              pasi100_n = n(),
              pasi100_r = sum(pasi100), 
              psa = mean(psa), 
              prevsys = mean(prevsys), 
              male = mean(male),
              bsa_mean = mean(bsa), 
              bsa_sd = sd(bsa)) %>%
    mutate(studyc = paste(studyc, name_)) 
  
}

## convert IPD to AgD
ipd_agd <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long), "IPD")

## convert IPD to duration-specific subgroup IPD
agd_sub_durn <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long, sub_durn), "SUB")

## convert IPD to age-specific subgroup IPD
agd_sub_age <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long, sub_age), "SUB")

## combine AGD with study level/treatment level subgroup
sub_agd1 <- bind_rows(ipd_agd,
                      agd_sub_durn, 
                      agd_sub_age)

func.mult <- function(res_var, agd.dt, col_name, value1, value2, ipd) {
  
  if (ipd == "Y") {
      net_org <- combine_network(
        set_ipd(sub_ipd,
              study = studyc, 
              trt = trtc, 
              r = {{res_var}},
              trt_class = trtclass, 
              trt_ref = "PBO"),
              trt_ref = "PBO")
  }
  else {
      filter_criteria1 <- interp(~y %in% x, .values=list(y = as.name(col_name), x = value1))
    
      filter_criteria2 <- interp(~y %in% x, .values=list(y = as.name(col_name), x = value2))
  
      net_org <- combine_network(
        set_ipd(sub_ipd %>% filter_(filter_criteria1),
            study = studyc, 
            trt = trtc, 
            r = {{res_var}},
            trt_class = trtclass, 
            trt_ref = "PBO"),
    
        set_agd_arm(agd.dt %>% filter_(filter_criteria2), 
                study = studyc, 
                trt = trtc, 
                r = pasi75_r, 
                n = pasi75_n,
                trt_class = trtclass,
                trt_ref = "PBO"),
    trt_ref = "PBO"
  )}

  # net_org
  netplot <- plot(net_org, weight_nodes = T, weight_edges = T, show_trt_class = T) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  
  ## Add integration points of to the AgD studies in the network
  
  if (ipd == "Y") {
  ## fixed effect NL-NMR
  FE_net_org <- nma(net_org, 
                    trt_effects = "fixed",
                    link = "probit", 
                    likelihood = "bernoulli2",
                    regression = ~(durnpso + prevsys + bsa + weight + psa + age)*.trt,
                    class_interactions = "common",
                    prior_intercept = normal(scale = 10),
                    prior_trt = normal(scale = 10),
                    prior_reg = normal(scale = 10),
                    init_r = 0.1,
                    QR = TRUE)
  }
  
  else {
    net_org <- add_integration(net_org,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               age = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               n_int = 1000)
    ## fixed effect NL-NMR
    FE_net_org <- nma(net_org, 
                      trt_effects = "fixed",
                      link = "probit", 
                      likelihood = "bernoulli2",
                      regression = ~(durnpso + prevsys + bsa + weight + psa + age)*.trt,
                      class_interactions = "common",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_reg = normal(scale = 10),
                      init_r = 0.1,
                      QR = TRUE)
  }
  
  return(FE_net_org)
  
}

## Using duration of psoriasis as subgroup

#####################################################################
##                        Original IPD+AGD                         ##
#####################################################################

# 4 IPD
FE_ipd <- func.mult(pasi75, ,"studyc" , , ,"Y")

## Choosing Uncover-3 as chosen IPD
## 3 IPD + 1 sub-group IPD aggregated
FE_ml_ipd <- func.mult(pasi75, ipd_agd, "studyc",
                                  c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                                  c("UNCOVER-3 IPD"))

 ## 3 IPD + 1 no-sub-groups IPD
FE_age_ml_sub <- func.mult(pasi75, agd_sub_age, "studyc",
                                  c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                                  c("UNCOVER-3 SUB"))

## save the result
rlist::list.save(FE_ipd, 'D:\\MPH-Project\\Further work\\paper\\FE_ipd.rdata')
rlist::list.save(FE_age_ml_ipd, 'D:\\MPH-Project\\Further work\\paper\\FE_age_ml_ipd.rdata')
rlist::list.save(FE_age_ml_sub, 'D:\\MPH-Project\\Further work\\paper\\FE_age_ml_sub.rdata')



## Using duration of psoriasis as subgroup

#####################################################################
##                        Original IPD+AGD                         ##
#####################################################################

## Choosing Uncover-3 as chosen IPD

## 3 IPD + 1 no-sub-groups IPD
FE_durn_ml_sub <- func.mult(pasi75, agd_sub_durn, "studyc",
                           c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                           c("UNCOVER-3 SUB"))

## save the result
rlist::list.save(FE_durn_ml_sub, 'D:\\MPH-Project\\Further work\\paper\\FE_durn_ml_sub.rdata')

