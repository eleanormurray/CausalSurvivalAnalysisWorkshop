/************************************************************
* File name: workshop.do
* Programmer: Eleanor Murray
* Date: February 11, 2019
* Purpose: Commented code to go along with survival workshop
************************************************************/

/********************************************/
/* Code Section 0 - Data Setup              */ 
/********************************************/

use "trial1.dta", clear

/*Create some nice labels*/
label define rand_lab 0 "placebo" 1 "treatment"
label values rand rand_lab


/**********************************************/
/* Code Section 1 - Data Exploration*/ 
/********************************************/

/* Load the data from the trial (if you haven't already)*/

use "trial1.dta"

/* Print the names of the available variables*/
describe

/*Look at the first 100 observations of the variables 'simid', 'visit', and 'rand'*/
list simid visit rand in 1/100

/* Total person-time: ie, what are the total number of rows in the dataset?*/
count

/* Sample size: total, i.e. how many individuals are in the dataset.*/
/* To interpret output, look at 'unique values'*/
codebook simid 

/*Number of individuals with at least 10 visits. To interpret, look at 'unique values'*/
codebook simid if visit == 10

/*Sample size: by trial arm*/
codebook simid if rand == 0
codebook simid if rand == 1

/* Create 'maxVisit' - total amount of time each individual contributed*/
by simid, sort: egen maxVisit = max(visit) if visit < .
/* Check by looking at some individuals with 14 or less than 14 visits*/
list simid visit maxVisit in 1/20
list simid visit maxVisit if maxVisit < 14 in 1/200


/* The variable death is only '1' at end-of-followup*/
/* Create 'deathOverall' - an indicator of whether //
  an individual died at any point during follow-up*/
by simid, sort: egen deathOverall = max(death) if visit <.
list simid visit maxVisit deathOverall death in 1/20
list simid visit maxVisit deathOverall death if maxVisit < 14 in 1/200

/*Create a baseline dataset for Kaplan-Meier and Cox models*/

/*In Stata this means set up survival using baseline visits*/
/*Note: our data are formatted so that people who die in the //
  first month have maxVisit = 0. */
/*Stata will drop these events, so we need to add a very small amount //
  of time to these if we want them included*/
gen maxVisit_s = maxVisit + 0.0001
stset maxVisit_s if visit == 0, failure(deathOverall)

/*How  many  individuals  died  overall? */
/*How  many  individuals  died  in  each  treatment  arm?*/
/*What was the cumulative probability of mortality by 14 visits in each treatment arm?*/
tab deathOverall if visit == 0
tab rand deathOverall if visit == 0 


/**********************************************/
/* Code Section 2 - Kaplan Meier*/ 
/********************************************/

/*Make sure you ran stset statment above or this code won't run*/
sts graph, by(rand) risktable ylabel(0.7 (0.1) 1) xlabel(#8) xtitle(Visit) ///
 ytitle(Survival probability) title(Kaplan-Meier Curve showing survival in each trial arm) legend(on)

 sts list, by(rand)


/**********************************************/
/* Code Section 3a - Unadjusted Hazard Ratios*/ 
/********************************************/

/* Data processing: create squared time variable [visit2]*/
gen visit2 = visit*visit

/*Calculate the unadjusted hazard ratio from a Cox PH model*/
stcox rand

/* Calculate the unadjusted hazard ratio from a pooled logistic regression model*/
/*We use the CLUSTER(SIMID) statement to get robust standard errors & valid, but conservative, confidence intervals*/
/*We use the LOGSITIC command to get the exponentiated coefficient, here the Hazard Ratio*/
logistic death visit visit2 rand, cluster(simid)


/**********************************************/
/* Code Section 3b - Conditional Hazard Ratios*/ 
/********************************************/

/* Calculate the baseline covariate-adjusted hazard ratio from a Cox PH model*/
stcox rand mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b
	
/* Calculate the baseline covariate-adjusted hazard ratio from a pooled logistic regression model*/
logistic death visit visit2 rand mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b, cluster(simid)
	
	

/********************************************/
/* Code Section 4 - Marginal Effects*/ 
/********************************************/

/* Step 1. Data processing: interaction terms*/
/* For Stata, the standardization is easier if we don't use separate variables for interactions*/
/* But, we should preserve the dataset so that we can undo the manipulations we'll do below*/
preserve

/* Step 2. Fit a pooled logistic  regression model with interaction terms between rand and visit & visit2 
to allow flexible fitting of baseline hazard, and save the coefficients*/
/*NOTE: We should use bootstraps to get valid, non-conservative confidence intervals*/
logistic death visit rand rand#c.visit c.visit#c.visit rand#c.visit#c.visit mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b, cluster(simid)


/* Step 3. Create simulated data where everyone is treated*/
/* Step 4. Create simulated data where everyone is untreated*/

/* NOTE: in Stata it's easies to do these steps together*/
/* Expand baseline so it contains a visit at each time point for every individual
	where the baseline information has been carried forward at each time*/
drop if visit != 0 
expand 15 if visit ==0 
bysort simid: replace visit = _n -1 

/* Set the treatment assignment to '1' for each individual and recreate  visit and interaction terms*/
/* Set the treatment assignment to '0' for each individual and recreate  visit and interaction terms*/
expand 2, generate(interv) 
replace rand = interv 

/* Calculate the predicted survival at each time and calculate survival from the cumulative product for each individual*/
predict pevent_k, pr 
gen psurv_k = 1-pevent_k 
keep simid visit rand interv psurv_k
sort simid interv visit
gen _t = visit + 1
gen psurv = psurv_k if _t == 1
bysort simid interv: replace psurv = psurv_k*psurv[_t-1] if _t > 1


/* Step 5. Calculate standardized survival at each time*/
/* Create concatenated dataset, only keep s, rand, and visit*/
by interv visit, sort: summarize psurv
keep simid psurv rand interv visit 

/* Calculate the mean survival at each visit within each treatment arm*/
bysort interv visit: egen meanS = mean(psurv)

/* Edit results to reflect that our estimates are for the END of the interval [t, t+1)*/
/* The code below duplicates baseline, sets survival to 100% in the duplicate and adds 1 to the visit count for each other record*/
expand 2 if visit == 0, generate(newvisit)
replace meanS = 1 if newvisit ==1
gen visit2 = 0 if newvisit == 1
replace visit2 = visit + 1 if newvisit ==0

/* Step 6. Plot the results*/
separate meanS, by(interv)


twoway (line meanS0 visit2, sort) (line meanS1 visit2, sort), ylabel(0.75(0.1)1.0) xlabel(#8) ytitle(Survival probability) xtitle(Visit) ///
 title(Survival Curves Standardized for Baseline Covariate Distribution)


/* Step 7. Print risk difference and hazard ratio at 14 weeks*/
bysort rand: egen meanS_14 = mean(psurv) if visit == 14
quietly summarize meanS_14 if(rand==0 & visit ==14)
matrix input observe = (0, `r(mean)')
quietly summarize meanS_14 if(rand==1 & visit == 14)
matrix observe = (observe\1, `r(mean)')
matrix observe = (observe\2, (log(observe[2,2]))/(log(observe[1,2])))
matrix observe = (observe\3, ((1-observe[2,2])/(1-observe[1,2])))
matrix observe = (observe\4, ((1-observe[2,2])-(1-observe[1,2])))
matrix rownames observe = Surv_placebo Surv_treatment cHR CIR risk_diff

matrix list observe
restore


/*********************************************************************************************************/
/**********************************Exercise 3*************************************************************/
/*********************************************************************************************************/

/*****************************************/
/* Code Section 5 - Data cleaning for IPW*/ 
/*****************************************/

use "trial1.dta", clear

/*Create some nice labels*/
label define adhr_lab 0 "Non-adherers" 1 "Adherers"
label values adhr_b adhr_lab

/* For this exercise, we want only the placebo arm. We will drop all individuals with rand = 1*/
drop if rand ==1

/*We don't need the interaction terms for Stata, but we do want to create a censoring variable*/
gen visit2 = visit*visit
gen adh_change = 0
replace adh_change = 1 if adhr != adhr_b
by simid, sort: egen cens_visit = min(cond(adh_change == 1),visit, .)
gen cens_new = 0
replace cens_new = 1 if visit == cens_visit
replace cens_new = . if visit > cens_visit
/*check*/
list simid visit adhr adhr_b adh_change cens_new in 1/100
drop cens_visit adh_change


/* Check to see how many individuals in your dataset*/
codebook simid 

/* Number of individuals who adhered versus did not adhere at visit 0*/
tab adhr_b if visit == 0

/* View first 100 observations of simID, visit, adhr_b, and adhr, and simID, visit, adhr_b, and adhr where adhr_b==1 to understand these variables better*/
list simid visit adhr_b adhr if adhr_b == 1 in 1/30

/* Need to recreate maxVisit and deathOverall*/
/* Create 'maxVisit' - total amount of time each individual contributed WHILE continuously adherent (ie cens_new = 0)*/
by simid, sort: egen maxVisit = max(visit) if cens_new ==0
/* Check by looking at some individuals with 14 or less than 14 visits*/
list simid visit maxVisit cens_new death in 1/100

/* The variable death is only '1' at end-of-followup*/
/* Create 'deathOverall' - an indicator of whether //
  an individual died during follow-up WHILE continuously adherent*/
by simid, sort: egen deathOverall = max(death) if cens_new == 0
list simid visit maxVisit deathOverall death cens_new in 1/100 

/*We don't need a new baseline dataset for Kaplan-Meier and Cox models in Stata*/
/*But we do need to add a small amount to the maxVisit so that deaths in the first month can be counted*/
gen maxVisit_s = maxVisit + 0.0001

/*How  many  individuals  died  overall? */
/*How  many  individuals  died  in  each  treatment  arm?*/
/*What was the cumulative probability of mortality by 14 visits in each treatment arm?*/
tab deathOverall if visit == 0
tab deathOverall rand if visit == 0 

/*Make sure you ran stset statment above or this code won't run*/
preserve
stset maxVisit_s if visit == 0, failure(deathOverall)
sts graph, by(adhr_b) risktable ylabel(0.7 (0.1) 1) xlabel(#8) xtitle(Visit) ///
 ytitle(Survival probability) title(Kaplan-Meier Curve showing survival in placebo arm by adherence) legend(on)
restore


/*****************************************/
/* Code Section 6 - Weight creation      */ 
/*****************************************/

/* Numerator: Pr(adhr(t)=1|adhr_b, Baseline covariates)*/
/* This model is created in data EXCLUDING the baseline visit*/
/* Create predicted probability at each time point */
/* (Pr(adhr(t) = 1 | adhr_b, baseline))*/

logit adhr visit visit2 adhr_b mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b, cluster(simid)
predict pnum, pr
	
/* Denominator: Pr(adhr(t)=1|adhr_b, Baseline covariates, Time-varying covariates)*/
/* Create predicted probability at each time point */
/* (Pr(adhr(t) = 1 | adhr_b, baseline, time-varying covariates))*/

logit adhr visit visit2 adhr_b mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b ///
	niha hiserchol hisertrigly hiheart chf ap ic diur antihyp oralhyp cardiom anyqqs anystdep fveb vcd, cluster(simid)
predict pdenom, pr

/*sort by simID and visit*/
sort simid visit

/*Calculate the weights*/ 
gen numcont = 1 if visit == 0
gen dencont = 1 if visit == 0
replace numcont = adhr*pnum + (1-adhr)*(1-pnum)
replace dencont = adhr*pdenom + (1-adhr)*(1-pdenom)

gen _t = visit + 1
gen k1_0 = 1 if _t == 1
gen k1_w = 1 if _t == 1
replace k1_0 = numcont*k1_0[_n-1] if _t > 1
replace k1_w = dencont*k1_w[_n-1] if _t > 1

gen unstabw = 1.0/k1_w
gen stabw = k1_0/k1_w

quietly summarize stabw, detail
gen stabw_t = stabw
replace stabw_t = r(p99) if stabw_t > r(p99)

summarize unstabw, detail
summarize stabw, detail
summarize stabw_t, detail


/*******************************************************/
/* Code Section 7 - Weighted conditional hazard ratios */ 
/*******************************************************/

/* Calculate the baseline covariate-adjusted hazard ratio from a pooled logistic regression model*/
/*Remember to restrict to only the uncensored person-time*/
/*Also remember to use the weights, as well as robust standard errors*/
/*Note that we are using baseline adherence not adherence over time since we are interested continuously maintaining an adherence level*/
logistic death visit visit2 adhr_b mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b if cens_new == 0 [pweight = stabw_t], cluster(simid)


/*******************************************************/
/* Code Section 8 - Weighted survival curves			*/ 
/*******************************************************/

/*First re-run the outcome model but including interactions between adherence and time*/

/* Calculate the baseline covariate-adjusted hazard ratio from a pooled logistic regression model*/
/*Remember to restrict to only the uncensored person-time*/
/*Also remember to use the weights, as well as robust standard errors*/
logistic death visit visit2 adhr_b adhr_b#c.visit adhr_b#c.visit2 mi_bin niha_b hiserchol_b hisertrigly_b hiheart_b chf_b ap_b ic_b ///
	diur_b antihyp_b oralhyp_b cardiom_b anyqqs_b anystdep_b fveb_b vcd_b if cens_new == 0 [pweight = stabw_t], cluster(simid)


/*Now, we standardize as before. The only difference is we always want adhr_b to be 1 or 0*/

/* Create simulated data where everyone is adherent*/
/* Create simulated data where everyone is non-adherent*/

/* NOTE: in Stata it's easies to do these steps together*/
/* Expand baseline so it contains a visit at each time point for every individual
	where the baseline information has been carried forward at each time*/
drop if visit != 0 
expand 15 if visit ==0 
bysort simid: replace visit = _n -1 

/* Set the adherence to '1' for each individual and recreate  visit and interaction terms*/
/* Set the adherence to '0' for each individual and recreate  visit and interaction terms*/
expand 2, generate(interv) 
replace  adhr_b = interv 

/* Calculate the predicted survival at each time and calculate survival from the cumulative product for each individual*/
predict pevent_k, pr 
gen psurv_k = 1-pevent_k 
keep simid visit rand interv psurv_k
sort simid interv visit
gen _t = visit + 1
gen psurv = psurv_k if _t == 1
bysort simid interv: replace psurv = psurv_k*psurv[_t-1] if _t > 1


/* Step 5. Calculate standardized survival at each time*/
/* Create concatenated dataset, only keep s, rand, and visit*/
by interv visit, sort: summarize psurv
keep simid psurv adhr_b interv visit 

/* Calculate the mean survival at each visit within each treatment arm*/
bysort interv visit: egen meanS = mean(psurv)

/* Edit results to reflect that our estimates are for the END of the interval [t, t+1)*/
/* The code below duplicates baseline, sets survival to 100% in the duplicate and adds 1 to the visit count for each other record*/
expand 2 if visit == 0, generate(newvisit)
replace meanS = 1 if newvisit ==1
gen visit2 = 0 if newvisit == 1
replace visit2 = visit + 1 if newvisit ==0

/* Step 6. Plot the results*/
separate meanS, by(interv)


twoway (line meanS0 visit2, sort) (line meanS1 visit2, sort), ylabel(0.75(0.1)1.0) xlabel(#8) ytitle(Survival probability) xtitle(Visit) ///
 title(Survival Curves Standardized for Baseline Covariate Distribution)


/* Step 7. Print risk difference and hazard ratio at 14 weeks*/
bysort rand: egen meanS_14 = mean(psurv) if visit == 14
quietly summarize meanS_14 if(rand==0 & visit ==14)
matrix input observe = (0, `r(mean)')
quietly summarize meanS_14 if(rand==1 & visit == 14)
matrix observe = (observe\1, `r(mean)')
matrix observe = (observe\2, (log(observe[2,2]))/(log(observe[1,2])))
matrix observe = (observe\3, ((1-observe[2,2])/(1-observe[1,2])))
matrix observe = (observe\4, ((1-observe[2,2])-(1-observe[1,2])))
matrix rownames observe = Surv_placebo Surv_treatment cHR CIR risk_diff

matrix list observe


