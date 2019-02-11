/************************************************************
* File name: workshop.do
* Programmer: Eleanor Murray
* Date: February 11, 2019
* Purpose: Commented code to go along with survival workshop
************************************************************/

/********************************************/
/* Code Section 0 - Data Setup              */ 
/* Read in the data set						*/
/* Tip: if Stata can't find the data do:	*/
/* 1. Close Stata							*/
/* 2. Go to the downloaded workshop folder	*/
/* 3. Double click 'workshop.do'			*/
/* 4. Run 'use' command again				*/
/********************************************/

use "trial1.dta", clear

/*Create some nice labels*/
label define rand_lab 0 "placebo" 1 "treatment"
label values rand rand_lab

/*********************************************************************************************************/
/**********************************Exercise 2*************************************************************/
/*********************************************************************************************************/


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
tab deathOverall rand if visit == 0 


/**********************************************/
/* Code Section 2 - Kaplan Meier*/ 
/********************************************/

/*Make sure you ran stset statment above or this code won't run*/
sts graph, by(rand) risktable ylabel(0.7 (0.1) 1) xlabel(#8) xtitle(Visit) ///
 ytitle(Survival probability) title(Kaplan-Meier Curve showing survival in each trial arm) legend(on)



/**********************************************/
/* Code Section 3 - Conditional Hazard Ratios*/ 
/********************************************/

/* Data processing: create squared time variable [visit2]*/
gen visit2 = visit*visit

/*Calculate the unadjusted hazard ratio from a Cox PH model*/
stcox rand

/* Calculate the unadjusted hazard ratio from a pooled logistic regression model*/
/*We use the CLUSTER(SIMID) statement to get robust standard errors & valid, but conservative, confidence intervals*/
/*We use the LOGSITIC command to get the exponentiated coefficient, here the Hazard Ratio*/
logistic death visit visit2 rand, vce(robust)

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

