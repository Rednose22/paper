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

setwd("D:\\MPH-Project\\Further work")
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

## Distribution of duration by study, treatment
full.crs.durn <- full_ipd %>%
  count(studyc, trtc, cat_pasi, durncat)

ggplot(full.crs.durn) +
  geom_bar(aes(x=durncat, y=n, fill = cat_pasi), stat = "identity") +
  facet_wrap(vars(studyc, trtc)) +
  labs(title = "Distribution before simulation by trials, studies")+
  guides(fill=guide_legend(title="pasi75")) +
  xlab("Duration group") +
  ylab("Frequency")


#####################################################################
##                          Data preparation                       ##
#####################################################################
## nest data
full_ipd.muta <- full_ipd %>%
  group_by(studyc) %>%
  nest()

## create model and then predict the f(pasi75)
set.seed(1234)

full_ipd.muta$data <- map(full_ipd.muta$data, function(df){
  a <- glm(pasi75 ~ age*trtc + bmi*trtc + male*trtc + bsa*trtc + weight*trtc + durnpso*trtc, family = "binomial", data = df, na.action = "na.exclude")
  df %>%
    mutate(lp = predict(a, type = "link"))
})

full_ipd.muta$data <- map(full_ipd.muta$data, function(df){
  df %>%
    mutate(
      ## case one: what we have done when the subgroup var is age instead of durnpso
      # lp2 = if_else(durnpso < 17 & trtc %in% c("PBO"),  lp - 0.2, lp),
      lp2 = if_else(durnpso < 17 & trtc %in% c("PBO"),  lp - 0.2, lp),
      lp2 = if_else(durnpso < 17 & !trtc %in% c("PBO"), lp2 + 0.2, lp2),

      ## case two: increase the event rate for durnpso < 17 in placebo group
      # lp2 = if_else(durnpso < 17 & trtc %in% c("PBO"),  lp - 2, lp),
      # lp2 = if_else(durnpso < 17 & !trtc %in% c("PBO"), lp2 + 0.2, lp2))})

      risk = plogis(lp2),
      ## Next line is stochastic
      pasi75_sim = rbinom(n = length(lp2), size = 1, prob = risk))
})

full_ipd <- full_ipd.muta %>%
  unnest(data) %>%
  filter(!is.na(pasi75_sim))

full.crs2 <- full_ipd %>%
  mutate(cat_pasi = ifelse(pasi75_sim == 1, "yes", "no"),
         durncat = ifelse(durnpso < 17, "lt 17", "ge 17")) %>%
  count(studyc, trtc, cat_pasi, durncat)

## histogram by studyc after redistribution
ggplot(full.crs2) +
  geom_bar(aes(x=durncat, y=n, fill = cat_pasi), stat = "identity") +
  facet_wrap(vars(studyc, trtc)) +
  labs(title = "Distribution after simulation by trials, studies") +
  guides(fill=guide_legend(title="pasi75"))


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

## convert IPD to duration-specific subgroup IPD
agd_sub_age <- sub_eff_func(quos(studyc, trtc, trtclass, trtn, trtc_long, sub_age), "SUB")

## combine AGD with study level/treatment level subgroup
sub_agd1 <- bind_rows(ipd_agd,
                      agd_sub_durn, 
                      agd_sub_age)

func.mult <- function(res_var, agd.dt, col_name, value1, value2, mode) {
  
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
  )
  
  # net_org
  # 
  netplot <- plot(net_org, weight_nodes = T, weight_edges = T, show_trt_class = T) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  
  ## Add integration points of to the AgD studies in the network
  
  if (mode == "multi") {
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
  
  else if (mode == "single") {
    net_org <- add_integration(net_org,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               n_int = 1000)
    ## fixed effect NL-NMR
    FE_net_org <- nma(net_org, 
                      trt_effects = "fixed",
                      link = "probit", 
                      likelihood = "bernoulli2",
                      regression = ~(durnpso)*.trt,
                      class_interactions = "common",
                      prior_intercept = normal(scale = 10),
                      prior_trt = normal(scale = 10),
                      prior_reg = normal(scale = 10),
                      init_r = 0.1,
                      QR = TRUE)
  }
  
  return(FE_net_org)
  
}

#####################################################################
##                        Original IPD+AGD                         ##
#####################################################################

# 4 IPD + 5 AgD
FE_net_durn <- func.mult(pasi75, sub_agd, "studyc",
                           c("IXORA", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3"),
                           c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                             "JUNCTURE"), "multi")

## Choosing Uncover-3 as chosen IPD
## 3 IPD + 1 sub-group IPD aggregated + 5 AgD
FE_net_durn_agd.1.1 <- func.mult(pasi75, sub_agd1, "studyc",
                                  c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                    "JUNCTURE", "UNCOVER-3 IPD"), "multi")

 ## 3 IPD + 1 no-sub-groups IPD + 5 AgD
FE_net_durn_agd.1.2 <- func.mult(pasi75, sub_agd1, "studyc",
                                  c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                                  c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                    "JUNCTURE", "UNCOVER-3 SUB"), "multi")

## save the result
rlist::list.save(FE_net_durn, 'D:\\MPH-Project\\Further work\\FE_net_durn.rdata')
rlist::list.save(FE_net_durn_agd.1.1, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.1.1.rdata')
rlist::list.save(FE_net_durn_agd.1.2, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.1.2.rdata')


FE_net_durn_agd.2.1 <- func.mult(pasi75, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "UNCOVER-1 IPD", "UNCOVER-2 IPD", "IXORA IPD"), "multi")

## 3 IPD + 1 no-sub-groups IPD + 5 AgD
FE_net_durn_agd.2.2 <- func.mult(pasi75, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "IXORA SUB", "UNCOVER-1 SUB", "UNCOVER-2 SUB"), "multi")


rlist::list.save(FE_net_durn_agd.2.1, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.2.1.rdata')
rlist::list.save(FE_net_durn_agd.2.2, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.2.2.rdata')

## simulated scenario
## 1 IPD + 3 sub-groups IPD + 5 AgD
FE_net_durn.2 <- func.mult(pasi75_sim, sub_agd, "studyc",
                           c("IXORA", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3"),
                           c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                             "JUNCTURE"), "multi")

FE_net_durn_agd.3.1 <- func.mult(pasi75_sim, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "UNCOVER-1 IPD", "UNCOVER-2 IPD", "IXORA IPD"), "multi")

## 1 IPD + 3 no-sub-groups IPD + 5 AgD
FE_net_durn_agd.3.2 <- func.mult(pasi75_sim, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "IXORA SUB", "UNCOVER-1 SUB", "UNCOVER-2 SUB"), "multi")

rlist::list.save(FE_net_durn.2, 'D:\\MPH-Project\\Further work\\model with higher interaction in plb group\\FE_net_durn.2.rdata')
rlist::list.save(FE_net_durn_agd.3.1, 'D:\\MPH-Project\\Further work\\model with higher interaction in plb group\\FE_net_durn_agd.3.1.rdata')
rlist::list.save(FE_net_durn_agd.3.2, 'D:\\MPH-Project\\Further work\\model with higher interaction in plb group\\FE_net_durn_agd.3.2.rdata')

rlist::list.save(FE_net_durn.2, 'D:\\MPH-Project\\Further work\\FE_net_durn.2.rdata')
rlist::list.save(FE_net_durn_agd.3.1, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.3.1.rdata')
rlist::list.save(FE_net_durn_agd.3.2, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.3.2.rdata')

## Only include durnpso in the NMA

#####################################################################
##                        Original IPD+AGD                         ##
#####################################################################

# 4 IPD + 5 AgD
FE_net_durn.s <- func.mult(pasi75, sub_agd, "studyc",
                         c("IXORA", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3"),
                         c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                           "JUNCTURE"), "single")

## Choosing Uncover-3 as chosen IPD
## 3 IPD + 1 sub-group IPD aggregated + 5 AgD
FE_net_durn_agd.1.1.s <- func.mult(pasi75, sub_agd1, "studyc",
                                 c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "UNCOVER-3 IPD"), "single")

## 3 IPD + 1 no-sub-groups IPD + 5 AgD
FE_net_durn_agd.1.2.s <- func.mult(pasi75, sub_agd1, "studyc",
                                 c("IXORA", "UNCOVER-1", "UNCOVER-2"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "UNCOVER-3 SUB"), "single")

## save the result
rlist::list.save(FE_net_durn.s, 'D:\\MPH-Project\\Further work\\FE_net_durn.s.rdata')
rlist::list.save(FE_net_durn_agd.1.1.s, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.1.1.s.rdata')
rlist::list.save(FE_net_durn_agd.1.2.s, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.1.2.s.rdata')


FE_net_durn_agd.2.1.s <- func.mult(pasi75, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "UNCOVER-1 IPD", "UNCOVER-2 IPD", "IXORA IPD"), "single")

## 3 IPD + 1 no-sub-groups IPD + 5 AgD
FE_net_durn_agd.2.2.s <- func.mult(pasi75, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "IXORA SUB", "UNCOVER-1 SUB", "UNCOVER-2 SUB"), "single")


rlist::list.save(FE_net_durn_agd.2.1.s, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.2.1.s.rdata')
rlist::list.save(FE_net_durn_agd.2.2.s, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.2.2.s.rdata')

## simulated scenario
## 1 IPD + 3 sub-groups IPD + 5 AgD
FE_net_durn.2.s <- func.mult(pasi75_sim, sub_agd, "studyc",
                           c("IXORA", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3"),
                           c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                             "JUNCTURE"), "single")

FE_net_durn_agd.3.1.s <- func.mult(pasi75_sim, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "UNCOVER-1 IPD", "UNCOVER-2 IPD", "IXORA IPD"), "single")

## 1 IPD + 3 no-sub-groups IPD + 5 AgD
FE_net_durn_agd.3.2.s <- func.mult(pasi75_sim, sub_agd1, "studyc",
                                 c("UNCOVER-3"),
                                 c("FIXTURE", "ERASURE", "CLEAR", "FEATURE",
                                   "JUNCTURE", "IXORA SUB", "UNCOVER-1 SUB", "UNCOVER-2 SUB"), "single")

rlist::list.save(FE_net_durn.2.s, 'D:\\MPH-Project\\Further work\\FE_net_durn.2.s.rdata')
rlist::list.save(FE_net_durn_agd.3.1.s, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.3.1.s.rdata')
rlist::list.save(FE_net_durn_agd.3.2.s, 'D:\\MPH-Project\\Further work\\FE_net_durn_agd.3.2.s.rdata')

