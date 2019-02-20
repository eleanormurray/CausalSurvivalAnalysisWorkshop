###########################################################################
# File name: workshop_v5.R
# Programmer: Lucia Petito & Eleanor Murray
# Date: February 20, 2019
# Purpose: Commented code to accompany causal survival workshop
###########################################################################

# Code Section 0 - Data Setup ---------------------------------------

# Load the following packages. If they do not exist on your computer,
# this code will automatically install them into your default library. 
# To set working directory
if(!require(here)) { install.packages("here"); require(here)}
# To use restricted cubic splines
if(!require(rms)) { install.packages("rms"); require(rms)}
# To run Cox PH models
if(!require(survival)) { install.packages("survival"); require(survival)}
# To calculate robust estimates of standard errors
if(!require(sandwich)) { install.packages("sandwich"); require(sandwich)}
if(!require(lmtest)) { install.packages("lmtest"); require(lmtest)}
# To plot results
if(!require(ggplot2)) { install.packages("ggplot2"); require(ggplot2)}
# To plot survival curves
if(!require(survminer)) { install.packages("survminer"); require(survminer)}
# To easily change data from long to wide
if(!require(reshape2)) { install.packages("reshape2"); require(reshape2)}
set_here()

# Code Section 1 - Data Exploration ---------------------------------------

# Load the data from the trial
trial <- read.csv(here("R/trial1.csv"), header=TRUE)

# Print the names of the available variables
names(trial)

# Look at the first 100 observations of the variables 'simid', 'visit', and 'rand'
View(trial[1:100, c('simID', 'visit', 'rand')])

# Total person-time
nrow(trial)

# Sample size
n <- length(unique(trial$simID))
n

# Number of observations at each visit
table(trial$visit)

# Create 'maxVisit' - total amount of time each individual contributed
trial<-trial %>%
  group_by(simID)%>%
  mutate(
    maxVisit = max(visit)
  )

# The variable death is only '1' at end-of-followup
# Create 'deathOverall' - an indicator of whether an individual died at any 
# point during follow-up
trial <- trial %>% 
  group_by(simID)%>%
  mutate(
    deathOverall = max(death)
  )

# Create 'baseline' - data collected at visit 0
baseline <- trial %>%
  dplyr::filter(visit == 0)

#Look at number of individuals in each arm, & deaths by arm
with(baseline, table(rand))
with(baseline, table(rand, deathOverall))
with(baseline, round(100*prop.table(table(rand, deathOverall), 1), 1))


# Code Section 2 - Kaplan Meier -------------------------------------------

# Use kaplan meier to nonparametrically estimate survival in each arm
# Note - this requires 1 observation per person of T (maxVisit - end of follow-up) and
# Delta (an indicator of whether death occurred at end of follow-up)
kmfit <- survfit(Surv(maxVisit, deathOverall) ~ rand, data = baseline)
summary(kmfit)

# Plot the output
ggsurvplot(kmfit, 
           data = baseline,
           conf.int=F,
           legend.labs = c("Placebo", "Treated"),
           ylim = c(0.7, 1),
           surv.scale = 'percent',
           xlab = "Number of Visits",
           title = "Kaplan-Meier Curve showing survival in each trial arm",
           risk.table = TRUE,
           break.time.by=2,
           ggtheme = theme_bw())

png(filename = "ITTkm.png", width = 2*1024, height = 2*1024, units = 'px', res = 72*5)

dev.off()
# Code Section 3 - Conditional Hazard Ratios ------------------------------

# Data processing: create squared time variable [visit2]
trial <- trial%>%
  mutate(
    visit2 = visit*visit
  )

# Calculate the unadjusted hazard ratio from a Cox PH model
cox_fit <- coxph(Surv(maxVisit, deathOverall) ~ rand, data = baseline, method='breslow')
summary(cox_fit)

# Calculate the unadjusted hazard ratio from a pooled logistic regression model
plr_fit <- glm(death ~ visit + visit2 + rand, data = trial, family=binomial())
coeftest(plr_fit, vcov=vcovHC(plr_fit, type="HC1")) # To get robust SE estimates
exp(coef(plr_fit)) # to get Hazard Ratios

# Calculate the baseline covariate-adjusted hazard ratio from a Cox PH model
adj_cox_fit <- coxph(Surv(maxVisit, deathOverall) ~ rand + mi_bin + NIHA_b + HiSerChol_b +
                       HiSerTrigly_b + HiHeart_b + CHF_b + 
                       AP_b + IC_b + DIUR_b + AntiHyp_b + 
                       OralHyp_b + CardioM_b + AnyQQS_b + 
                       AnySTDep_b + FVEB_b + VCD_b, 
                     data = baseline, method='breslow')
summary(adj_cox_fit)

# Calculate the baseline covariate-adjusted hazard ratio from a pooled logistic regression model
adj_plr_fit <- glm(death ~ visit + visit2 + rand + mi_bin + NIHA_b + HiSerChol_b +
                     HiSerTrigly_b + HiHeart_b + CHF_b + 
                     AP_b + IC_b + DIUR_b + AntiHyp_b + 
                     OralHyp_b + CardioM_b + AnyQQS_b + 
                     AnySTDep_b + FVEB_b + VCD_b, data = trial, family=binomial())
coeftest(adj_plr_fit, vcov=vcovHC(adj_plr_fit, type="HC1")) # To get robust SE estimates
exp(coef(adj_plr_fit))



# Code Section 4 - Marginal Effects ---------------------------------------

# Step 1. Data processing: interaction terms
# Create interaction terms between visit, visit2, and randomization
trial <- trial%>%
  mutate(
    randvisit = rand*visit, 
    randvisit2 = rand*visit2
    )

# Step 2. Fit a pooled logistic  regression model with interaction terms between 
# rand and visit & visit1 to allow flexible fitting of baseline hazard
adj_plr_ix_fit <- glm(death ~ visit + visit2 + rand + randvisit + randvisit2 + 
                        mi_bin + NIHA_b + HiSerChol_b +
                        HiSerTrigly_b + HiHeart_b + CHF_b + 
                        AP_b + IC_b + DIUR_b + AntiHyp_b + 
                        OralHyp_b + CardioM_b + AnyQQS_b + 
                        AnySTDep_b + FVEB_b + VCD_b, data = trial, family=binomial)
summary(adj_plr_ix_fit)
exp(coef(adj_plr_ix_fit))

# Step 3. Create simulated data where everyone is treated
# Expand baseline so it contains a visit at each time point for every individual
# where the baseline information has been carried forward at each time
treated <- baseline[rep(1:n,each=15),]
treated$visit <- rep(0:14, times=n) # This recreates the time variable
treated <- treated%>%
  mutate(
    #recreate squared visit term
      visit2 = visit*visit, 
    # Set the treatment assignment to '1' for each individual and
        rand = 1, 
    # recreate the interaction terms
        randvisit = rand*visit,
        randvisit2 = rand*visit2
    ) 

# 'predict' returns predicted "density" of survival at each time
# conditional on covariates
# Turn these into predicted survival density by subtracting from 1
treated$p <- 1 - predict(adj_plr_ix_fit, newdata=treated, type='response')

# We calculate survival by taking the cumulative product by individual
treated <- treated %>%
  dplyr::arrange(simID,visit)%>%
  group_by(simID)%>%
  mutate(
    s = cumprod(p)
  )

# Step 4. Create simulated data where everyone receives placebo
# When simulating data in the placebo arm, only difference from treated is 
# in the randomization assignment, and resulting interaction terms
placebo <- baseline[rep(1:n,each=15),]
placebo$visit <- rep(0:14, times=n) # This recreates the time variable
placebo <- placebo%>%
  mutate(
    #recreate squared visit term
    visit2 = visit*visit, 
    # Set the treatment assignment to '1' for each individual and
    rand = 0, 
    # recreate the interaction terms
    randvisit = rand*visit,
    randvisit2 = rand*visit2
  ) 

# 'predict' returns predicted "density" of survival at each time
# conditional on covariates
# Turn these into predicted survival density by subtracting from 1
placebo$p <- 1 - predict(adj_plr_ix_fit, newdata=placebo, type='response')

# We calculate survival by taking the cumulative product by individual
placebo <- placebo %>%
  dplyr::arrange(simID,visit)%>%
  group_by(simID)%>%
  mutate(
    s = cumprod(p)
  )


# Step 5. Calculate standardized survival at each time
# Create concatenated dataset, only keep s, rand, and visit
both <- dplyr::bind_rows(treated, placebo)
both <- both[,c('s', 'rand', 'visit')]

# Calculate the mean survival at each visit within each treatment arm
results <- both%>%
  group_by(visit, rand)%>%
  dplyr::summarize(mean_survival = mean(s))


# Edit results data frame to reflect that our estimates are for the END of the interval [t, t+1)
# Add a row for each of Placebo and Treated where survival at time 0 is 1.
results <- results%>%
  dplyr::ungroup()%>%
  mutate(
    visit = visit+1
  )

results<-dplyr::bind_rows(c(visit = 0, rand = 0, mean_survival =  1), c(visit = 0, rand = 1, mean_survival =  1), results)

# Add a variable that treats randomization as a factor
results$randf <- factor(results$rand, labels = c("Placebo", "Treated"))

# Step 6. Plot the results
p2 <- ggplot(results, aes(x=visit, y=mean_survival))+
  geom_line(aes(colour=randf)) +
  geom_point(aes(colour=randf))+
  xlab("Number of Visits") +
  scale_x_continuous(limits = c(0, 15), breaks=seq(0,15,2)) +
  ylab("Probability of Survival") +
  ggtitle("Survival Curves Standardized for Baseline Covariate Distribution") +
  labs(colour="Treatment Arm") +
  theme_bw() +
  theme(legend.position="bottom")

png(filename = "ITTsurv.png", width = 2*1060, height = 2*1024, units = 'px', res = 72*5)
p2
dev.off()

# Step 7. Calculate risk difference and hazard ratio at 14 weeks
# Transpose the data so survival in each treatment arm is separate
wideres <- dcast(results, visit ~ randf, value.var = 'mean_survival')
head(wideres)

# Create summary statistics
wideres <- wideres%>%
  mutate(
    RD = (1-Treated) - (1-Placebo),
    logRatio = log(Treated)/log(Placebo),
    CIR = (1-Treated)/ (1-Placebo)
  )
wideres$logRatio[1] <- NA
wideres$cHR <- sapply(0:15, FUN=function(x){mean(wideres$logRatio[wideres$visit <= x], na.rm=T)})

# Print all wide results to fill in table
round(wideres, 3)

# Overall Hazard Ratio
wideres$cHR[wideres$visit==15]
# Risk difference at end of 14 visits
wideres$RD[wideres$visit==15]
# Cumulative incidence ratio at end of 14 visits
wideres$CIR[wideres$visit==15]


# Code Section 5 - Data cleaning for IPW ----------------------------------

# Remove all objects from the R environment EXCEPT the cleaned trial dataframe
rm(adj_cox_fit, adj_plr_fit, adj_plr_ix_fit,baseline, both, cox_fit, kmfit, p2, placebo, plr_fit, results, treated,wideres, n)

# Create a placebo dataset that only contains individuals with placebo
placebo <- trial%>%  
    dplyr::filter(rand == 0)
  

# Check to see how many individuals in your dataset
length(unique(placebo$simID))

# Number of individuals who adhered versus did not adhere at visit 0
table(placebo$adhr_b[placebo$visit==0])

# View first 50 observations of simID, visit, adhr_b, and adhr, and 
# simID, visit, adhr_b, and adhr where adhr_b==1 to understand these variables better
View(placebo[1:50,c("simID", 'visit', 'adhr_b', 'adhr')])
View(placebo[placebo$adhr_b==1,c("simID", 'visit', 'adhr_b', 'adhr')])

#Data set up 
placebo <- placebo%>%
  # Create interaction terms between exposure (adhr_b) and visit, visit2
  # Create censoring variable - indicate when individual deviates from baseline
    mutate(
    adhr0visit = adhr_b*visit,
    adhr0visit2 = adhr_b*visit2,
    adh_change = ifelse(adhr != adhr_b,1,0)
    )%>%
  # Need to recreate maxVisit and deathOverall - slightly more complicated now
  group_by(simID)%>%
  mutate(
    cens_visit = ifelse(is.na((which(adh_change %in% 1)[1])), max(visit),(which(adh_change %in% 1)[1]- 1)),
    cens_new = ifelse(cens_visit == visit, 1, ifelse(cens_visit < visit, NA , 0)),
    
    maxVisit_cens = ifelse(cens_visit == 14, 14, cens_visit ),
    deathOverall_new = ifelse(is.na(cens_new), NA, 
                            ifelse(is.na((which(death %in% 1)[1])), 0, 
                              ifelse((which(death %in% 1)[1]-1)<=  maxVisit_cens, 1, 0)))
  )%>%
  dplyr::ungroup()

View(placebo[1:50,c("simID", 'visit', 'adhr_b', 'adhr', 'cens_new', 'maxVisit_cens','cens_visit', 'adh_change', "death",  "deathOverall_new")])


# Create baseline data
baseline_p <- placebo%>%
  dplyr::filter(visit == 0)

#Look at number of individuals in each arm, & deaths by arm
with(baseline_p, table(adhr_b))
with(baseline_p, table(adhr_b, deathOverall_new))


# Check kaplan-meier fit in each arm
kmfit2 <- survfit(Surv(maxVisit_cens, deathOverall_new) ~ adhr_b, data = baseline_p)
summary(kmfit2)

# Plot the output
p3 <- 
  ggsurvplot(kmfit2, 
           data = baseline_p,
           conf.int=F,
           ylim = c(0.65, 1),
           surv.scale = 'percent',
           legend.labs = c("Nonadherers", "Adherers"),
           xlab = "Number of Visits",
           title = "Kaplan-Meier Curve showing survival by adherence",
           risk.table = TRUE,
           break.time.by=2,
           ggtheme = theme_bw())
p3

png(filename = "PPkm.png", width = 2*1024, height = 2*1024, units = 'px', res = 72*5)
p3
dev.off()

###############FROM HERE -- check p3 graph against STATA & SAS output

# Code Section 6 - Weight Creation ----------------------------------------

# Numerator: Pr(adhr(t)=1|adhr_b, Baseline covariates)
# This model is created in data EXCLUDING the baseline visit

nFit <- glm(adhr ~ visit + visit2 + adhr_b + 
              mi_bin + NIHA_b + HiSerChol_b +
              HiSerTrigly_b + HiHeart_b + CHF_b + 
              AP_b + IC_b + DIUR_b + AntiHyp_b + 
              OralHyp_b + CardioM_b + AnyQQS_b + 
              AnySTDep_b + FVEB_b + VCD_b,
            data=placebo[placebo$visit > 0,],
            family=binomial())

# Create predicted probability at each time point 
# (Pr(adhr(t) = 1 | adhr_b, baseline))
placebo$pnum <- predict(nFit, newdata=placebo, type='response')

# Denominator: Pr(adhr(t)=1|adhr_b, Baseline covariates, Time-varying covariates)
dFit <- glm(adhr ~ visit + visit2 + adhr_b + 
              mi_bin + NIHA_b + HiSerChol_b +
              HiSerTrigly_b + HiHeart_b + CHF_b + 
              AP_b + IC_b + DIUR_b + AntiHyp_b + 
              OralHyp_b + CardioM_b + AnyQQS_b + 
              AnySTDep_b + FVEB_b + VCD_b +
              NIHA + HiSerChol +
              HiSerTrigly + HiHeart + CHF +
              AP + IC + DIUR + AntiHyp +
              OralHyp + CardioM + AnyQQS +
              AnySTDep + FVEB + VCD,
            data=placebo[placebo$visit > 0,],
            family=binomial())
# Create predicted probability at each time point 
# (Pr(adhr(t) = 1 | adhr_b, baseline, time-varying covariates))
placebo$pdenom <-  predict(dFit, newdata=placebo, type='response')

# Sort placebo by simID and visit
placebo <- placebo[order(placebo$simID, placebo$visit),]

# Contribution from adhr(t) = 0 is 1 - p
# Contribution from adhr(t) = 1 is p
placebo$numCont <- with(placebo, adhr*pnum + (1-adhr)*(1-pnum))
placebo$denCont <- with(placebo, adhr*pdenom + (1-adhr)*(1-pdenom))

# Set contribution at baseline visit to 1
placebo$numCont[placebo$visit==0] <- 1; placebo$denCont[placebo$visit==0] <- 1

# Numerator
placebo$k1_0 <- with(placebo, ave(numCont, simID, FUN=cumprod))
# Denominator
placebo$k1_w <- with(placebo, ave(denCont, simID, FUN=cumprod))

# Create both stabilized and unstabilized weights
placebo$stabw <- with(placebo, k1_0 / k1_w)
placebo$unstabw <- with(placebo, 1 / k1_w)

# Check the weights
# Can do this with built-in functions
summary(placebo$stabw); quantile(placebo$stabw, p=c(0.01, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.995))
summary(placebo$unstabw); quantile(placebo$unstabw, p=c(0.01, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.995))

# Or write something for yourself
wt_check_fxn <- function(x, d){
  # x - weight vector
  # d - number of digits to round to
  round(data.frame(N = length(x),
                   Missing = sum(is.na(x)),
                   Distinct = length(unique(x)),
                   Mean = mean(x, na.rm=T),
                   SD = sd(x, na.rm=T),
                   Min = min(x, na.rm=T),
                   Q.01 = quantile(x, p=0.01),
                   Q.25 = quantile(x, p=0.25),
                   Q.50 = quantile(x, p=0.5),
                   Q.75 = quantile(x, p=0.75),
                   Q.90 = quantile(x, p=0.9),
                   Q.95 = quantile(x, p=0.95),
                   Q.99 = quantile(x, p=0.99),
                   Q.995 = quantile(x, p=0.995),
                   Max = max(x, na.rm=T) ), digits=d)
}
wt_check_fxn(placebo$stabw, 3)
wt_check_fxn(placebo$unstabw, 3)

# To prevent high weights from too much influence, truncate
threshold <- quantile(placebo$stabw, 0.99) # Here we chose 99th %tile
placebo$stabw_t <- placebo$stabw
placebo$stabw_t[placebo$stabw > threshold] <- threshold

wt_check_fxn(placebo$stabw_t, 3)

# Code Section 7 - Weighted Conditional Hazard Ratios ---------------------

# Truncated stabilized weights
plrFit_SWT <- glm(death ~ visit + visit2 + adhr_b + 
                    mi_bin + NIHA_b + HiSerChol_b +
                    HiSerTrigly_b + HiHeart_b + CHF_b + 
                    AP_b + IC_b + DIUR_b + AntiHyp_b + 
                    OralHyp_b + CardioM_b + AnyQQS_b + 
                    AnySTDep_b + FVEB_b + VCD_b,
                  data=placebo[placebo$visit <= placebo$maxVisit,],
                  weights = stabw_t,
                  family=quasibinomial())
coeftest(plrFit_SWT, vcov=vcovHC(plrFit_SWT, type="HC1")) # To get robust SE estimates
exp(coef(plrFit_SWT)) # to get Hazard Ratios

# Try changing "weights = stabw_t, " to 
# "weights = stabw_t, " for stabilized fit
# "weights = unstabw, " for unstabilized fit

# Code Section 8 - Weighted Survival Curves -------------------------------

# Step 1. Estimate weighted outcome regression with interactions
plrixFit_USW <- glm(death ~ visit + visit2 + adhr_b + 
                      adhr0visit + adhr0visit2 +
                      mi_bin + NIHA_b + HiSerChol_b +
                      HiSerTrigly_b + HiHeart_b + CHF_b + 
                      AP_b + IC_b + DIUR_b + AntiHyp_b + 
                      OralHyp_b + CardioM_b + AnyQQS_b + 
                      AnySTDep_b + FVEB_b + VCD_b,
                    data=placebo[placebo$visit <= placebo$maxVisit,],
                    weights = stabw_t,
                    family=quasibinomial())

summary(plrixFit_USW)
exp(coef(plrixFit_USW))

# Step 1a. Create dataset with just baseline values 
baseline <- placebo[placebo$visit==0,]

# Step 2. Create simulated data where everyone adheres and doesn't adhere
adherers <- baseline[rep(1:n,each=15),]
adherers$visit <- rep(0:14, times=n) # Recreates the time variable
adherers$visit2 <- adherers$visit * adherers$visit
adherers$adhr_b <- 1
adherers$adhr0visit <- adherers$adhr_b*adherers$visit
adherers$adhr0visit2 <- adherers$adhr_b*adherers$visit2
adherers$p <- 1 - predict(plrixFit_USW, newdata=adherers, type='response')
adherers$s <- ave(adherers$p, adherers$simID, FUN=cumprod)

# Nonadherers
nonadherers <- adherers
nonadherers$adhr_b <- 0
nonadherers$adhr0visit <- nonadherers$adhr_b*nonadherers$visit
nonadherers$adhr0visit2 <- nonadherers$adhr_b*nonadherers$visit2
nonadherers$p <- 1 - predict(plrixFit_USW, newdata=nonadherers, type='response')
nonadherers$s <- ave(nonadherers$p, nonadherers$simID, FUN=cumprod)

# Step 3. Calculate standardized survival at each time
# Create concatenated dataset, only keep s, adhr_b, and visit
both <- rbind(adherers, nonadherers)
both <- both[,c('s', 'adhr_b', 'visit')]

# Calculate the mean survival at each visit within each adherer group
results <- aggregate(s ~ visit + adhr_b, FUN=mean, data=both)

# Edit results data frame to reflect that our estimates are for the END of the interval [t, t+1)
# Add a row for each of Placebo and Treated where survival at time 0 is 1.
results$visit <- results$visit + 1
results <- rbind(c(0, 0, 1), c(0, 1, 1), results)

# Add a variable that treats randomization as a factor
results$adhrf <- factor(results$adhr_b, labels = c("Nonadherers", "Adherers"))

# Plot the results
p4 <- ggplot(results, aes(x=visit, y=s)) + 
  geom_line(aes(colour=adhrf)) +
  geom_point(aes(colour=adhrf))+
  xlab("Number of Visits") +
  scale_x_continuous(limits = c(0, 15), breaks=seq(0,15,2)) +
  ylab("Probability of Survival") +
  ggtitle("Survival Curves Standardized for Baseline Covariate Distribution \nand weighted for time-varying confounders") +
  labs(colour="Adherer Group") +
  theme_bw() +
  theme(legend.position="bottom") 
p4
png(filename = "PPsurv.png", width = 2*1060, height = 2*1024, units = 'px', res = 72*5)
p4
dev.off()

# Step 4. Calculate risk difference and hazard ratio at 14 weeks (EOF)
# Transpose the data so survival in each treatment arm is separate
wideres <- dcast(results, visit ~ adhrf, value.var = 's')
head(wideres)
wideres$RD <- (1 - wideres$Adherers) - (1 - wideres$Nonadherers)
wideres$logRatio <- log(wideres$Adherers) / log(wideres$Nonadherers)
wideres$logRatio[1] <- NA
wideres$cHR <- sapply(0:15, FUN=function(x){mean(wideres$logRatio[wideres$visit <= x], na.rm=T)})

# Overall Hazard Ratio
wideres$cHR[wideres$visit==15]
# Risk difference at 14 visits
wideres$RD[wideres$visit==15]
# Cumulative incidence ratio at end of 14 visits
with(wideres[wideres$visit==15,], (1 - Adherers) / (1 - Nonadherers))


# Extra Code - Bootstrapping Unadjustd Pooled Logistic Regression --------

# Remove all objects from the R environment EXCEPT the cleaned trial dataframe
rm(list=setdiff(ls(), "trial"))

# Data processing: create squared time variable [visit2]
trial$visit2 <- trial$visit*trial$visit

# Calculate the unadjusted hazard ratio from a pooled logistic regression model
plr_fit <- glm(death ~ visit + visit2 + rand, data = trial, family=binomial())
coeftest(plr_fit, vcov=vcovHC(plr_fit, type="HC1")) # To get robust SE estimates
exp(coef(plr_fit)[['rand']]) # to get Hazard Ratio
exp(coefci(plr_fit, parm = 'rand', vcov=vcovHC(plr_fit, type="HC1")))

# Use the ClusterBootstrap package - this is REALLY slow
# library(ClusterBootstrap)
# set.seed(1)
# system.time(bootobj <- clusbootglm(death ~ visit + visit2 + rand, data=trial, 
#                                    clusterid = simID, family=binomial, B=25,
#                                    confint.level = 0.95, n.cores=detectCores()-1))
# exp(coef(plr_fit)[['rand']])
# exp(confint(bootobj, parm = "rand", level=0.95, interval.type = "percentile"))

# Do something more complicated
# Register parallel backend
if(!require(doParallel)) { install.packages("doParallel"); require(doParallel)}
if(!require(foreach)) { install.packages("foreach"); require(foreach)}
registerDoParallel(cores=detectCores()-1)
getDoParWorkers() # Print number of cores
B <- 50
set.seed(123)
# Set up bootstrap samples
ids <- unique(trial$simID)
nind <- length(ids)
bootids <- matrix(ids[sample.int(nind, nind*B, replace=T)], nrow=nind, ncol=B)

# Now bootstrap the results
system.time(bootres <- foreach(i=1:B, .combine=rbind) %dopar% {
  bootdfLoc <- sapply(bootids[,i], function(x) which(trial[,"simID"]==x))
  lengthbid <- sapply(bootids[,i], function(x) length(which(trial[,"simID"]==x))) 
  bootdf <- trial[unlist(bootdfLoc),]
  bootdf$bid <- rep(1:nind, lengthbid)
  
  return(coef(glm(death ~ visit + visit2 + rand, data=bootdf, family=binomial))[['rand']])
})
exp(quantile(bootres, p=c(0.025, 0.975)))


# Now for something even faster!
# Make a list that is B long, where each entry is 
# a matrix indicating how many times each individual occurred in the bth bootstrap resample
if(!require(plyr)) { install.packages("plyr"); require(plyr)}
if(!require(dplyr)) { install.packages("dplyr"); require(dplyr)}
bootidsFreq <- apply(bootids, 2, FUN=function(b){
  temp <- data.frame(simID = b)
  return(plyr::count(temp, vars = 'simID'))
})

system.time(bootres <- foreach(i=1:B, .combine=rbind) %dopar% {
  require(dplyr)
  bootdf <- dplyr::left_join(trial, bootidsFreq[[i]], by = "simID")
  bootdf$freq[is.na(bootdf$freq)] <- 0
  
  return(coef(glm(death ~ visit + visit2 + rand, data=bootdf, family=binomial, weights=freq))[['rand']])
})
exp(quantile(bootres, p=c(0.025, 0.975)))

