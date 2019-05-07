/************************************************************
* File name: workshop_v6.sas
* Programmer: Eleanor Murray
* Date: May 1, 2019
* Purpose: Commented code to go along with survival workshop
************************************************************/

/**********************************************/
/* Code Section 0 - Data Setup               */ 
/* Assign a library to the location of data */
/********************************************/

/*On windows, set user to USERPROFILE; on Unix, set user to HOME*/
%let user = %sysget(USERPROFILE);
/*Set dlgcdir to the location of your extracted file
/*If you unzipped the github repository in your download folder, you don't need to change anything*/
%let rc = %sysfunc(dlgcdir("&user./Downloads/CausalSurvivalWorkshop_2019/SAS"));

libname surv "./";

/* Or, set working directory to your location manually by uncommenting the code below */
/*libname surv "C:\Users\ejmurray\Dropbox\Conferences.Talks\2019\CausalSurvivalWorkshop_2019\SAS";*/

/* Load the data from the trial*/
data trial;
set surv.trial1;
run;


/**********************************************/
/* Code Section 1 - Data Exploration*/ 
/********************************************/

/* Print the names of the available variables*/
proc contents data = trial;
run;


/* Look at the first 100 observations of the variables 'simid', 'visit', and 'rand'*/
proc print data = trial (obs = 100);
var simid visit rand;
run;

/* Total person-time*/
proc means data = trial n;
var simid;
title 'Total person-time';
run;

/* Sample size: total*/
proc freq data=trial nlevels;
	tables simid /noprint;
	title 'Sample size, total';
run;


/*Number of individuals with at least 10 visits*/
proc freq data = trial;
where visit = 10;
tables visit / list nopercent nocum;
title 'Number of people with at least 10 visits';
run;

/*Sample size: by trial arm*/
proc freq data=trial nlevels;
	by rand;
	tables simid /noprint;
	title 'Sample size by randomization arm';
run;

/* Create 'maxVisit' - total amount of time each individual contributed*/
data trial;
set trial;
by simid;
retain vis_count;
if first.simid then vis_count = 1;
else vis_count = vis_count +1;
run;

proc means data = trial max noprint;
by simid;
var vis_count;
output out = maxvisits (keep = simid maxVisit) max = maxVisit;
run;

data trial;
merge trial maxvisits;
by simid;
run;

proc print data = trial (obs = 50);
var simid visit maxVisit vis_count;
title 'Check max visit: all observations';
run;

proc print data = trial (obs = 50);
where death = 1;
var simid visit maxVisit vis_count;
title 'Check max visit: observations with deaths';
run;

/* The variable death is only '1' at end-of-followup
/* Create 'deathOverall' - an indicator of whether an individual died at any point during follow-up*/
proc means data = trial max noprint;
by simid;
var death;
output out = deaths (keep = simid deathOverall) max = deathOverall;
run;

proc print data = deaths(obs = 20);
var simid deathOverall; 
run;

data trial;
merge trial deaths;
by simid;
run;

proc print data = trial (obs = 50);
var simid visit deathOverall;
title 'Check overall death: all observations';
run;

proc print data = trial (obs = 50);
var simid visit deathOverall;
title 'Check overall death: observations with deaths';
run;

proc freq data = trial;
where visit = 0;
tables deathOverall rand*deathOverall /nocol nopercent nocum;
title 'Cumulative mortality, overall and by trial arm';
run;

/* Create 'baseline' - data collected at visit 0*/
data baseline;
set trial;
if visit = 0 then output;
run;


/**********************************************/
/* Code Section 2 - Kaplan Meier*/ 
/********************************************/

/*Create a template for the Kaplan Meier curve with custom title (change at "entrytitle"),
custom y-axis values (change at "viewmin = ", "viewmax = ", "tickvaluelist = "),
and customize legend (change in "layout gridded" at "entry = ")*/
proc template;
      define statgraph Stat.Lifetest.Graphics.ProductLimitSurvival;
         dynamic NStrata xName plotAtRisk plotCensored plotCL plotHW plotEP
            labelCL labelHW labelEP maxTime StratumID classAtRisk plotBand
            plotTest GroupName yMin Transparency SecondTitle TestName pValue GROUPNAME;
      BeginGraph;
            entrytitle "Kaplan-Meier curve";
            if (EXISTS(SECONDTITLE))
            	entrytitle "showing survival in each trial arm ";
       		endif;
			layout lattice / rows=2 columns=1 columndatarange=unionall
                             rowweights=(.85 .15);
            layout overlay / xaxisopts=(shortlabel=XNAME offsetmin=.05
                  linearopts=(viewmax=MAXTIME)) yaxisopts=(label=
                  "Survival Probability" shortlabel="Survival" linearopts=(
                  viewmin=0.6 viewmax=1 tickvaluelist=(.6 .8 1.0)));
                   stepplot y=SURVIVAL x=TIME / group=STRATUM index=STRATUMNUM 
                     name="Survival" rolename=(_tip1=ATRISK _tip2=EVENT) tip=(y
                     x Time _tip1 _tip2);
                   if (PLOTCENSORED)
                     scatterplot y=CENSORED x=TIME / group=STRATUM index=
                     STRATUMNUM markerattrs=(symbol=plus);
                  endif;
                   DiscreteLegend "Survival" / title=GROUPNAME
                     titleattrs=GraphValueText location=inside
                     autoalign=(BOTTOM);
                   dynamic NObs1 NObs2 NObs3 NEvent1 NEvent2;
                  layout gridded / columns=3 border=TRUE autoalign=(TopRight);
                     entry "";      entry "Event";   entry "Total";
                     entry "0: Placebo";     entry NEvent1;   entry NObs1;
                     entry "1: Treatment";     entry NEvent2;   entry NObs2;
                  endlayout;
               endlayout;
                layout overlay / xaxisopts=(display=none);
                  blockplot x=TATRISK block=ATRISK / class=CLASSATRISK
                        repeatedvalues=true display=(label values)
                        valuehalign=start valuefitpolicy=truncate
                        labelposition=left labelattrs=GRAPHVALUETEXT
                        valueattrs=GRAPHDATATEXT(size=7pt)
                        includemissingclass=false;
               endlayout;
            endlayout;
         EndGraph;
      end;
   run;
/* Use kaplan meier to nonparametrically estimate survival in each arm and plot output*/
ods graphics on;
proc lifetest data = baseline plots =(survival(atrisk));
	time maxVisit*deathOverall(0);
	strata rand;
run;
ods graphics off;
/*reset sas Kaplan Meier template*/
proc template;
   delete Stat.Lifetest.Graphics.ProductLimitSurvival;
run;

/**********************************************/
/* Code Section 3 - Conditional Hazard Ratios*/ 
/********************************************/

/* Data processing: create squared time variable [visit2]*/
data trial;
set trial;
	visit2 = visit*visit;
run;

/* Calculate the unadjusted hazard ratio from a Cox PH model*/
proc phreg data = baseline;
	model maxVisit*deathoverall(0) = rand / rl;
	title 'ITT, no adjustment (Cox model)';
run;

/* Calculate the unadjusted hazard ratio from a pooled logistic regression model*/
/*We use the REPEATED statement to get robust standard errors & valid, but conservative, confidence intervals*/
/*We use the ESTIMATE statement to get the exponentiated coefficient, here the Hazard Ratio*/
proc genmod data = trial descending;
	class simid;
	model death = visit visit2 rand /link = logit dist = bin;
	repeated subject = simid / type = ind;
	estimate "Hazard ratio for unadjusted ITT" rand 1 / exp;
	title 'ITT, no adjustment (Pooled logistic)';
run;
	
/* Calculate the baseline covariate-adjusted hazard ratio from a Cox PH model*/
proc phreg data = baseline;
	model maxVisit*deathoverall(0) = rand 
		MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b/ rl;
	title 'ITT, baseline adjustment (Cox model)';
run;

/* Calculate the baseline covariate-adjusted hazard ratio from a pooled logistic regression model*/
proc genmod data = trial descending;
	class simid;
	model death = visit visit2 rand 
		MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b/ link = logit dist = bin;
	repeated subject = simid / type = ind;
	estimate "Hazard ratio for adjusted ITT" rand 1 / exp;
	title 'ITT, baseline adjustment (Pooled logistic)';
run;


/********************************************/
/* Code Section 4 - Marginal Effects*/ 
/********************************************/

/* Step 1. Data processing: interaction terms
/* Create interaction terms between visit, visit2, and randomization*/
data trial;
set trial;
	randvisit = rand*visit;
	randvisit2 = rand*visit2;
run;

/* Step 2. Fit a pooled logistic  regression model with interaction terms between rand and visit & visit2 
to allow flexible fitting of baseline hazard, and save the coefficients*/
/*NOTE: in order to output ParameterEstimates, we can no longer use the REPEATED statement
/*The confidence intervals from this model will not be valid
/*Instead, we should use bootstraps to get valid, non-conservative confidence intervals*/

proc genmod data = trial descending;
	class simid;
	model death = visit visit2 rand randvisit randvisit2
		MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b/ link = logit dist = bin;
	title 'ITT, baseline adjustment (Pooled logistic)';
	ods output ParameterEstimates = adj_plr_ix_fit;
run;

/* Step 3. Create simulated data where everyone is treated
/* Expand baseline so it contains a visit at each time point for every individual where the baseline information has been carried forward at each time
/* Set the treatment assignment to '1' for each individual and recreate  visit and interaction terms
/* Calculate the predicted survival at each time and calculate survival from the cumulative product for each individual*/

/*remove nuisance parameter*/
data adj_plr_ix_fit;
set adj_plr_ix_fit;
if PARAMETER = "Scale" then delete;
run;
/*save parameters and variable names*/
proc sql noprint;
	select ESTIMATE FORMAT =16.12 INTO: IBC_ESTIMATE separated by ' ' from adj_plr_ix_fit;
quit;
proc sql noprint;
	select PARAMETER INTO: MODEL separated by ' ' from adj_plr_ix_fit;
quit;
proc means sum noprint data = adj_plr_ix_fit;	
	var DF;
	output out = nobs (drop = _type_ _freq_ where=(_stat_ ="N"));
run;
proc sql noprint;	
	select DF into:NVAR separated by ' ' from nobs;	quit;


data treated (keep = s ci rand visit);
set trial;
where visit = 0;
	array var{&nvar} &model;
	array coef{&nvar} (&ibc_estimate);
	intercept = 1;
	s=1;

	rand = 1;

	do visit = 0 to 14;
		visit2 = visit*visit;
		xbeta = 0;
        
        randvisit = rand*visit;
        randvisit2 = rand*visit2;

		do i = 1 to dim(var);
			xbeta = xbeta + coef[i] *var[i];
		end;
    	p = 1/(1+exp(-xbeta));
		s = s*(1-p);
		ci = 1-s;
		output;
	end;
run;


/* Step 4. Create simulated data where everyone receives placebo
/* When simulating data in the placebo arm, only difference from treated is in the randomization assignment, and resulting interaction terms*/
data placebo (keep = s ci rand visit);
set trial;
where visit = 0;
	array var{&nvar} &model;
	array coef{&nvar} (&ibc_estimate);
	intercept = 1;
	s=1;

	rand = 0;

	do visit = 0 to 14;
		visit2 = visit*visit;
		xbeta = 0;
        
        randvisit = rand*visit;
        randvisit2 = rand*visit2;
        
		do i = 1 to dim(var);
			xbeta = xbeta + coef[i] *var[i];
		end;
    	p = 1/(1+exp(-xbeta));
		s = s*(1-p);
		ci = 1-s;
		output;
	end;
run;

/* Step 5. Calculate standardized survival at each time
/* Create concatenated dataset, only keep s, rand, and visit*/
data both;
set treated placebo;
by rand;
run;

/* Calculate the mean survival at each visit within each treatment arm*/
proc means data = both mean ;
	class visit rand;
	types visit*rand;
    var s;
    output out = means (drop = _type_ _freq_) mean(s) = s;
run;

/* Edit results data frame to reflect that our estimates are for the END of the interval [t, t+1)*/

proc transpose data=means out = wideres prefix = Surv_;
var s;
id rand;
by visit;
run;

data wideres;
set wideres;
	visit = visit+1;
	merge_id = 1;
    rd = (1-Surv_1) -(1- Surv_0);
	ratio = log(surv_1)/log(surv_0);
	CIR = (1-surv_1)/(1-surv_0);
run;

data wideres;
set wideres;
total + ratio;
cHR = total / _n_;
drop total ratio;
run;
 
data zero;
	visit = 0;
	merge_id = 1;
	risk_1 = 0;
	risk_0 = 0;
	surv_0 = 1;
	surv_1 = 1;
	rd = 0;
	cHR = 1;
	CIR = 1;
run;
data results;
merge wideres zero;
by visit;
drop merge_id _NAME_;
run;

/*Print output*/

/* Step 6. Plot the results*/
goptions reset=all;
axis1 length=5.9in label=(angle=90 height=3 "Survival Curves Standardized for Baseline Covariate Distribution") width=2 order=(0.6 to 1.1 by .1) 
     major=(height=2 width=5) minor=none value=(height=2.5) offset=(2,2);
axis2 length=5.9in label=(height=3 justify=center "Time (visits)") 
     width=2 order=(0 to 15 by 1) major=(height=1 width=5) 
     minor=none value=(height=2.5) offset=(2,2);
symbol1 c=black h=12 l=1 w=5 v=none i=line mode=include;
symbol2 c=black h=12 l=4 w=5 v=none i=line mode=include;
legend1 label=none value=(color=black h=1.5 "Treatment" "Placebo") 
     down=4 position=(bottom left inside) 
     shape=symbol(6,3) offset=(0,1) mode=protect;

proc gplot data=results; plot (surv_1 surv_0 )*visit/overlay vaxis=axis1 haxis=axis2
     legend=legend1 noframe;
run; 
quit; 


/* Step 7. Print risk difference and hazard ratio at 14 weeks*/
proc print data = results round;
var visit surv_0 surv_1 rd;
title 'standardized survival and risk difference estimates at the end of follow-up';
run;
proc print data = results round;
where visit = 15;
var cHR CIR rd;
title 'standardized HR, cumulative incidence ratio, and cumulative incidence difference up to visit 15';
run;


/*********************************************************************************************************/
/**********************************Exercise 3*************************************************************/
/*********************************************************************************************************/

/*****************************************/
/* Code Section 5 - Data cleaning for IPW*/ 
/*****************************************/

/* Create interaction terms between exposure (rand) and visit, visit2
/* Create censoring variable - indicate when individual deviates from baseline*/
data trial_full;
set surv.trial1;
by simid;
	retain cens_new;
	if first.simid then cens_new = 0;
	if adhr = 0 then cens_new = 1;
	visit2 = visit*visit;   
run;

/* Number of individuals who adhered in each arm at visit 0*/
proc sort data = trial_full;
by rand;
run;
proc freq data = trial_full;
	by rand;
	where visit = 0;
	tables adhr_b;
	title 'Sample size by adherence & randomization arm';
run;
/* Check to see how many individuals in your dataset*/
proc freq data=trial_full nlevels;
	tables simid /noprint;
	title 'Sample size by randomization arm';
run;

proc freq data = trial_full;
tables adhr*cens_new;
run;

/* Need to recreate maxVisit */
proc sort data = trial_full; by simid;run;
data trial_full;
set trial_full;
by simid;
retain vis_count;
if first.simid then vis_count = 1;
else do;
	if cens_new = 0 then vis_count = vis_count + 1;
	else vis_count = vis_count;
end;
run;

proc means data = trial_full max noprint;
by simid;
var vis_count;
output out = maxvisits (keep = simid maxVisit_cens) max = maxVisit_cens;
run;

data trial_full;
merge trial_full maxvisits;
by simid;
run;

proc means data = trial_full max noprint;
where visit le maxVisit_cens;
by simid;
var death;
output out = deaths (keep = simid deathOverall) max = deathOverall;
run;

data trial_full;
merge trial_full deaths;
by simid;
run;

/* Create baseline data*/
data baseline;
set trial_full;
where visit = 0;
run;
proc freq data = baseline;
tables deathoverall /missing;
run;


/* Check kaplan-meier fit in each arm*/
proc template;
      define statgraph Stat.Lifetest.Graphics.ProductLimitSurvival;
         dynamic NStrata xName plotAtRisk plotCensored plotCL plotHW plotEP
            labelCL labelHW labelEP maxTime StratumID classAtRisk plotBand
            plotTest GroupName yMin Transparency SecondTitle TestName pValue;
      BeginGraph;
            entrytitle "Kaplan-Meier curve";
            if (EXISTS(SECONDTITLE))
            	entrytitle "showing survival in each trial arm among continuous adherers, unadjusted ";
       		endif;
			layout lattice / rows=2 columns=1 columndatarange=unionall
                             rowweights=(.85 .15);
            layout overlay / xaxisopts=(shortlabel=XNAME offsetmin=.05
                  linearopts=(viewmax=MAXTIME)) yaxisopts=(label=
                  "Survival Probability" shortlabel="Survival" linearopts=(
                  viewmin=0.6 viewmax=1 tickvaluelist=(.6 .8 1.0)));
                   stepplot y=SURVIVAL x=TIME / group=STRATUM index=STRATUMNUM
                     name="Survival" rolename=(_tip1=ATRISK _tip2=EVENT) tip=(y
                     x Time _tip1 _tip2);
                   if (PLOTCENSORED)
                     scatterplot y=CENSORED x=TIME / group=STRATUM index=
                     STRATUMNUM markerattrs=(symbol=plus);
                  endif;
                   DiscreteLegend "Survival" / title=GROUPNAME
                     titleattrs=GraphValueText location=inside
                     autoalign=(BOTTOM);
                   dynamic NObs1 NObs2 NObs3 NEvent1 NEvent2;
                  layout gridded / columns=3 border=TRUE autoalign=(TopRight);
                     entry "";      entry "Event";   entry "Total";
                     entry "0: Placebo";     entry NEvent1;   entry NObs1;
                     entry "1: Treatment";     entry NEvent2;   entry NObs2;
                  endlayout;
               endlayout;
                layout overlay / xaxisopts=(display=none);
                  blockplot x=TATRISK block=ATRISK / class=CLASSATRISK
                        repeatedvalues=true display=(label values)
                        valuehalign=start valuefitpolicy=truncate
                        labelposition=left labelattrs=GRAPHVALUETEXT
                        valueattrs=GRAPHDATATEXT(size=7pt)
                        includemissingclass=false;
               endlayout;
            endlayout;
         EndGraph;
      end;
   run;
/* Use kaplan meier to nonparametrically estimate survival in each arm and plot output*/
ods graphics on;
proc lifetest data = baseline plots =(survival(atrisk));
	time maxVisit_cens*deathOverall(0);
	strata rand;
run;
ods graphics off;
/*reset sas Kaplan Meier template*/
proc template;
   delete Stat.Lifetest.Graphics.ProductLimitSurvival;
run;

/********************************/
/* Code Section 6 - Weight Creation*/ 
/********************************/

/* Numerator: Pr(adhr(t)=1|adhr_b, Baseline covariates)
/* This model is created in data EXCLUDING the baseline visit
/* Create predicted probability at each time point 
/* (Pr(adhr(t) = 1 | adhr_b, baseline))*/
proc logistic data= trial_full (where =(visit >0 & rand = 1 )) descending;
	model adhr = visit visit2 adhr_b
		MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b
		;
	output out = nFit_1 (keep=simid visit pnum_1) p = pnum_1;
	run;

proc logistic data= trial_full (where =(visit >0 & rand = 0 )) descending;
	model adhr = visit visit2 adhr_b
		MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b
		;
	output out = nFit_0 (keep=simid visit pnum_0) p = pnum_0;
	run;

/* Denominator: Pr(adhr(t)=1|adhr_b, Baseline covariates, Time-varying covariates)
/* Create predicted probability at each time point 
/* (Pr(adhr(t) = 1 | adhr_b, baseline, time-varying covariates))*/

proc logistic data= trial_full (where =(visit >0 & rand = 1)) descending;
		model adhr =  visit visit2 adhr_b
        MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b
        
       	NIHA	HiSerChol  HiSerTrigly HiHeart  	CHF 	
		AP  	IC 	DIUR	AntiHyp  	OralHyp
	  	CardioM  	AnyQQS 	AnySTDep	FVEB  	VCD 
		;
		output out = dFIT_1 (keep = simid visit pdenom_1)  p = pdenom_1;
	run;
	
proc logistic data= trial_full (where =(visit >0 & rand = 0)) descending;
		model adhr =  visit visit2 adhr_b
        MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b
        
       	NIHA	HiSerChol  HiSerTrigly HiHeart  	CHF 	
		AP  	IC 	DIUR	AntiHyp  	OralHyp
	  	CardioM  	AnyQQS 	AnySTDep	FVEB  	VCD 
		;
		output out = dFIT_0 (keep = simid visit pdenom_0)  p = pdenom_0;
	run;

/* Sort by simID and visit*/
proc sort data = trial_full;
by simid visit;
run;
proc sort data = nFit_1;
by simid visit;
run;
proc sort data = dFit_1;
by simid visit;
run;
proc sort data = nFit_0;
by simid visit;
run;
proc sort data = dFit_0;
by simid visit;
run;

/*Calculate the weights*/ 
data trial_full;
merge trial_full nFit_1 dFit_1  nFit_0 dFit_0;
by simid visit;

if first.simid then do;
	k1_0 = 1.0;
	k1_w = 1.0;
end;

retain k1_0 k1_w;

else do;
	if pnum_1 = . then pnum_1 = 1;
	if pnum_0 = . then pnum_0 = 1;
	if pdenom_1 = . then pdenom_1 = 1;
	if pdenom_0 = . then pdenom_0 = 1;

	if rand = 1 then do;
		numCont = adhr*pnum_1 +(1-adhr)*(1-pnum_1);
		denCont = adhr*pdenom_1 +(1-adhr)*(1-pdenom_1);
	end;
	else if rand = 0 then do;
		numCont = adhr*pnum_0 +(1-adhr)*(1-pnum_0);
		denCont = adhr*pdenom_0 +(1-adhr)*(1-pdenom_0);
	end;
		k1_0 = k1_0*numCont;
		k1_w = k1_w*denCont;
end;
unstabw = 1.0/k1_w;
stabw = k1_0/k1_w;
run;

/*Create truncated weights*/
proc means data = trial_full p99 noprint;
var stabw;
output out=pctl (keep=p99) p99=p99;
run;
data temp;
set pctl;
call symput ('cutoff', p99);
run;
data trial_full;
set trial_full;
stabw_t = stabw;
if stabw>%sysevalf(&cutoff) then do;
	stabw_t = %sysevalf(&cutoff);
end;
run;

/*Check weights*/
proc means data = trial_full n mean min max p99 std median p25 p75 nmiss;
var unstabw stabw stabw_t;
title 'Distribution of weights';
run;

/*****************************************************/
/*Code Section 7 - Weighted Conditional Hazard Ratios*/
/*****************************************************/

/*Estimate weighted hazard ratio*/
proc genmod data= trial_full (where = (visit < maxVisit_cens))  descending;
	class simid;
		model death =  visit visit2 rand
        MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b /link = logit dist = bin ;
		
		weight stabw_t;
		/*To get robust SE estimates*/
		repeated subject = simid / type = ind;
		/*To get Hazard Ratio*/
		estimate 'PPE' rand 1/exp;	
		ods output GEEFitCriteria = fit;
		title 'Per-protocol effect, with IPW';
	run;

/* Try changing "weights = stabw_t, " to */
/* "weights = stabw_t, " for stabilized fit*/
/* "weights = unstabw, " for unstabilized fit*/


/*********************************************/
/* Code Section 8 - Weighted Survival Curves */
/*********************************************/

/*Create interaction terms with randomization and time*/
data trial_full;
set trial_full;
randvisit = rand*visit;
randvisit2 = rand*visit2;
run;

/* Step 1. Estimate weighted outcome regression with interactions and store parameter estimates*/
proc genmod data= trial_full (where = (visit < maxVisit_cens))  descending;
	class simid;
		model death =  visit visit2 rand randvisit randvisit2
        MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b /link = logit dist = bin ;

		weight stabw_t;		
		ods output ParameterEstimates  = plrFit_ix_SWT;

	
		title 'Per-protocol effect, with IPW';
	run;

data plrFit_ix_SWT;
set plrFit_ix_SWT;
	if PARAMETER = "Scale" then delete;
run;
proc sql noprint;
	select ESTIMATE FORMAT =16.12 INTO: IBC_ESTIMATE separated by ' ' from plrFit_ix_SWT;
quit;
proc sql noprint;
	select PARAMETER INTO: MODEL separated by ' ' from plrFit_ix_SWT;
quit;
proc means sum noprint data = plrFit_ix_SWT;	
	var DF;
	output out = nobs (drop = _type_ _freq_ where=(_stat_ ="N"));
run;
proc sql noprint;	
	select DF into:NVAR separated by ' ' from nobs;	quit;

/* Step 2. Create simulated datasets where everyone is always adherent and never adherent*/

data placebo (keep = s ci rand visit);
set trial_full;
where visit = 0;
	array var{&nvar} &model;
	array coef{&nvar} (&ibc_estimate);
	intercept = 1;
	s=1;

	rand = 0;

	do visit = 0 to 14;
		visit2 = visit*visit;
		xbeta = 0;
        
		randvisit = rand*visit;
		randvisit2 = rand*visit2;
		
		do i = 1 to dim(var);
			xbeta = xbeta + coef[i] *var[i];
		end;
    	p = 1/(1+exp(-xbeta));
		s = s*(1-p);
		ci = 1-s;
		output;
	end;
run;
data treated (keep = s ci rand visit);
set trial_full;
where visit = 0;
	array var{&nvar} &model;
	array coef{&nvar} (&ibc_estimate);
	intercept = 1;
	s=1;

	rand = 1 ;

	do visit = 0 to 14;
		visit2 = visit*visit;
		xbeta = 0;
        
		randvisit = rand*visit;
		randvisit2 = rand*visit2;

		do i = 1 to dim(var);
			xbeta = xbeta + coef[i] *var[i];
		end;
    	p = 1/(1+exp(-xbeta));
		s = s*(1-p);
		ci = 1-s;
		output;
	end;
run;

/* Step 5. Calculate standardized survival at each time
/* Create concatenated dataset, only keep s, rand, and visit*/
data both;
set placebo treated;
by rand;
run;

/* Calculate the mean survival at each visit within each treatment arm*/
proc means data = both mean ;
	class visit rand;
	types visit*rand;
    var s;
    output out = means (drop = _type_ _freq_) mean(s) = s;
run;

proc transpose data=means out = wideres prefix = Surv_;
var s;
id rand;
by visit;
run;

data wideres;
set wideres;
	visit = visit+1;
	merge_id = 1;
    rd = (1-Surv_1) -(1- Surv_0);
	ratio = log(surv_1)/log(surv_0);
	CIR = (1-Surv_1)/(1-Surv_0);
run;

data wideres;
set wideres;
total + ratio;
cHR = total / _n_;
drop total ratio;
run;
 
data zero;
	visit = 0;
	merge_id = 1;
	risk_1 = 0;
	risk_0 = 0;
	surv_0 = 1;
	surv_1 = 1;
	rd = 0;
	cHR = 1;
	CIR = 1;
run;
data results;
merge wideres zero;
by visit;
drop merge_id _NAME_;
run;

/*Plot the results*/
goptions reset=all;
axis1 length=5.9in label=(angle=90 height=3 "Survival Curves Standardized for Baseline Covariate Distribution and weighted for time-varying confounders") width=2 order=(0.6 to 1.1 by .1) 
     major=(height=2 width=5) minor=none value=(height=2.5) offset=(2,2);
axis2 length=5.9in label=(height=3 justify=center "Time (visits)") 
     width=2 order=(0 to 15 by 1) major=(height=1 width=5) 
     minor=none value=(height=2.5) offset=(2,2);
symbol1 c=black h=12 l=1 w=5 v=none i=line mode=include;
symbol2 c=black h=12 l=4 w=5 v=none i=line mode=include;
legend1 label=none value=(color=black h=1.5 "Always adhere" "Never adhere") 
     down=4 position=(bottom left inside) 
     shape=symbol(6,3) offset=(0,1) mode=protect;

proc gplot data=results; plot (surv_1 surv_0 )*visit/overlay vaxis=axis1 haxis=axis2
     legend=legend1 noframe;
run; 
quit; 

/* Step 4. Print risk difference and hazard ratio at 14 weeks*/
proc print data = results round;
var visit surv_0 surv_1 rd;
title 'standardized survival and risk difference estimates at the end of follow-up';
run;
proc print data = results round;
where visit = 15;
var cHR RD CIR ;
title 'IPW Standardized HR, CIR, and risk difference up to visit 15';
run;
