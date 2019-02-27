/************************************************************
* File name: workshop_v4.sas
* Programmer: Eleanor Murray
* Date: October 31, 2018
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
/*libname surv "<yourpathhere>\CausalSurvivalWorkshop_2019\SAS";*/

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

/* Create a placebo dataset that only contains individuals assigned to placebo
/* Create interaction terms between exposure (adhr_b) and visit, visit2
/* Create censoring variable - indicate when individual deviates from baseline*/
data placebo;
set surv.trial1;
where rand = 0;
by simid;
	retain cens_new;
	if first.simid then cens_new = 0;
	if adhr ne adhr_b then cens_new = 1;

	visit2 = visit*visit;
	adhr0visit = adhr_b*visit;
	adhr0visit2 = adhr_b*visit2;
    
run;

proc print data = placebo (obs = 20);
var simid visit adhr adhr_b cens_new;
run;

/* Check to see how many individuals in your dataset*/
proc freq data=placebo nlevels;
	tables simid /noprint;
	title 'Sample size by randomization arm';
run;

/* Number of individuals who adhered versus did not adhere at visit 0*/
proc sort data = placebo;
by adhr_b;
run;
proc freq data = placebo nlevels;
	by adhr_b;
	tables simid / noprint;
	title 'Sample size by randomization arm';
run;

/* View first 30 observations of simID, visit, adhr_b, and adhr, and simID, visit, adhr_b, and adhr where adhr_b==1 to understand these variables better*/
proc print data = placebo (obs = 30);
var simid visit adhr_b adhr;
title 'check data';
run;

proc print data = placebo (obs = 30);
where adhr_b = 1;
var simid visit adhr_b adhr;
title 'check data';
run;

/* Need to recreate maxVisit and deathOverall*/
proc sort data = placebo; by simid;run;
data placebo;
set placebo;
by simid;
retain vis_count;
if first.simid then vis_count = 1;
else do;
	if cens_new = 0 then vis_count = vis_count + 1;
	else vis_count = vis_count;
end;

run;

proc means data = placebo max noprint;
by simid;
var vis_count;
output out = maxvisits (keep = simid maxVisit) max = maxVisit;
run;

data placebo;
merge placebo maxvisits;
by simid;
run;

proc means data = placebo max noprint;
where visit le maxVisit;
by simid;
var death;
output out = deaths (keep = simid deathOverall) max = deathOverall;
run;

data placebo;
merge placebo deaths;
by simid;
run;

/* Create baseline data*/
data baseline;
set placebo;
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
            	entrytitle "showing survival in placebo arm by adherence ";
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
                     entry "0: NonAdherers";     entry NEvent1;   entry NObs1;
                     entry "1: Adherers";     entry NEvent2;   entry NObs2;
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
	strata adhr_b;
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
proc logistic data=  placebo (where =(visit >0 )) descending;
	model adhr = visit visit2 adhr_b
		MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b
		;
	output out = nFit(keep=simid visit pnum) p = pnum;
	run;

/* Denominator: Pr(adhr(t)=1|adhr_b, Baseline covariates, Time-varying covariates)
/* Create predicted probability at each time point 
/* (Pr(adhr(t) = 1 | adhr_b, baseline, time-varying covariates))*/

proc logistic data= placebo (where =(visit >0 )) descending;
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
		output out = dFIT (keep = simid visit pdenom)  p = pdenom;
	run;
	
/* Sort by simID and visit*/
proc sort data = placebo;
by simid visit;
run;
proc sort data = nFit;
by simid visit;
run;
proc sort data = dFit;
by simid visit;
run;

/*Calculate the weights*/ 
data placebo;
merge placebo nFit dFit;
by simid visit;

if first.simid then do;

	k1_0 = 1.0;
	k1_w = 1.0;
end;

retain k1_0 k1_w;

else do;
		numCont = adhr*pnum +(1-adhr)*(1-pnum);
		denCont = adhr*pdenom +(1-adhr)*(1-pdenom);
		k1_0 = k1_0*numCont;
		k1_w = k1_w*denCont;

end;
unstabw = 1.0/k1_w;
stabw = k1_0/k1_w;
run;

/*Create truncated weights*/
proc means data = placebo p99 noprint;
var stabw;
output out=pctl (keep=p99) p99=p99;
run;
data temp;
set pctl;
call symput ('cutoff', p99);
run;
data placebo;
set placebo;
stabw_t = stabw;
if stabw>%sysevalf(&cutoff) then do;
	stabw_t = %sysevalf(&cutoff);
end;
run;

/*Check weights*/
proc means data = placebo n mean min max p99 std median p25 p75 nmiss;
var unstabw stabw stabw_t;
title 'Distribution of weights';
run;

/*****************************************************/
/*Code Section 7 - Weighted Conditional Hazard Ratios*/
/*****************************************************/

/*Estimate weighted hazard ratio*/
proc genmod data= placebo (where = (visit le maxVisit))  descending;
	class simid;
		model death =  visit visit2 adhr_b
        MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b /link = logit dist = bin ;
		
		weight stabw_t;
		/*To get robust SE estimates*/
		repeated subject = simid / type = ind;
		/*To get Hazard Ratio*/
		estimate 'PPE' adhr_b 1/exp;	
		title 'Placebo arm adherence comparison, with IPW';
	run;

/* Try changing "weights = stabw_t, " to */
/* "weights = stabw_t, " for stabilized fit*/
/* "weights = unstabw, " for unstabilized fit*/

/*********************************************/
/* Code Section 8 - Weighted Survival Curves */
/*********************************************/

/* Step 1. Estimate weighted outcome regression with interactions and store parameter estimates*/

proc genmod data= placebo (where = (visit le maxVisit)) descending ;
	class simid;
		model death =  visit visit2 adhr_b adhr0visit adhr0visit2
        MI_bin NIHA_b HiSerChol_b 	
		HiSerTrigly_b HiHeart_b CHF_b
		AP_b IC_b DIUR_b AntiHyp_b 
		OralHyp_b CardioM_b AnyQQS_b 
		AnySTDep_b FVEB_b VCD_b /link = logit dist = bin ;
		
		weight stabw_t;
		ods output ParameterEstimates = aplrixFit_USW;
	run;

data aplrixFit_USW;
set aplrixFit_USW;
	if PARAMETER = "Scale" then delete;
run;
proc sql noprint;
	select ESTIMATE FORMAT =16.12 INTO: IBC_ESTIMATE separated by ' ' from aplrixFit_USW;
quit;
proc sql noprint;
	select PARAMETER INTO: MODEL separated by ' ' from aplrixFit_USW;
quit;
proc means sum noprint data = aplrixFit_USW;	
	var DF;
	output out = nobs (drop = _type_ _freq_ where=(_stat_ ="N"));
run;
proc sql noprint;	
	select DF into:NVAR separated by ' ' from nobs;	quit;

/* Step 2. Create simulated datasets where everyone is always adherent and never adherent*/

data adherers (keep = s ci adhr_b visit);
set placebo;
where visit = 0;
	array var{&nvar} &model;
	array coef{&nvar} (&ibc_estimate);
	intercept = 1;
	s=1;

	adhr_b = 1;

	do visit = 0 to 14;
		visit2 = visit*visit;
		xbeta = 0;
        
        adhr0visit = adhr_b*visit;
        adhr0visit2 = adhr_b*visit2;

		do i = 1 to dim(var);
			xbeta = xbeta + coef[i] *var[i];
		end;
    	p = 1/(1+exp(-xbeta));
		s = s*(1-p);
		ci = 1-s;
		output;
	end;
run;
data nonadherers (keep = s ci adhr_b visit);
set placebo;
where visit = 0;
	array var{&nvar} &model;
	array coef{&nvar} (&ibc_estimate);
	intercept = 1;
	s=1;

	adhr_b = 0;

	do visit = 0 to 14;
		visit2 = visit*visit;
		xbeta = 0;
        
        adhr0visit = adhr_b*visit;
        adhr0visit2 = adhr_b*visit2;

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
set adherers nonadherers;
by adhr_b;
run;

/* Calculate the mean survival at each visit within each treatment arm*/
proc means data = both mean ;
	class visit adhr_b;
	types visit*adhr_b;
    var s;
    output out = means (drop = _type_ _freq_) mean(s) = s;
run;


proc transpose data=means out = wideres prefix = Surv_;
var s;
id adhr_b;
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
var cHR CIR rd;
title 'IPW Standardized HR, CIR, and risk difference up to visit 15';
run;

/*********************************************************************************************************/
/**********************************Extra code: Exercise 2 bootstraps**************************************/
/*********************************************************************************************************/
/* Load the data from the trial and create functional form of time*/
data trial;
set surv.trial1;
run;
data trial;
set trial;
	visit2 = visit*visit;
run;

%macro ITT_boots(inset = , nboot = );
/*make sure data is correctly sorted*/
proc sort data=&inset out=onesample;
	by simid visit;
run;
/*Generate a count of unique individuals*/
data onesample ;
  set onesample end = _end_ ;
  by simid;
 retain _id ;
  if _n_ = 1 then _id = 0;
  if first.simid then do;
  	_id = _id + 1 ;	
  end;
  if _end_ then do ;
     call symput("nids",trim(left(_id)));
  end;
run;
/*create a temporary dataset with a every ID copied &nboot times*/	
data ids ;
   do bsample = 1 to &nboot;
       do _id = 1 to &nids ;
           output ;
       end;
   end;
run;
/*Sample with replacement from list of IDs to create bootstrap samples
/*Temporary dataset _idsamples will contain a list of IDs and the number of times each occurs in each bootstrap sample*/
proc surveyselect data= ids 
         method = urs
         n= &nids
         seed = 1232  
         out = _idsamples (keep = bsample _id  numberhits  ) 
         outall  noprint  ;       
      strata bsample ;
run;
/*Create master output file for storing the results of each bootstrap run*/
data means_all;
bsample = .;
run;

/*Loop through main analysis (bsample = 0) and all bootstrap samples*/
%do bsample = 0 %to &nboot;
	/*create dataset from bootstrap ID lists*/ 
	proc sort data = onesample; 
		by _id;
	run;
	data bootsample;
	merge onesample _idsamples (where = (bsample = &bsample));
		by _id;
	run;
	proc sort data = bootsample; 
		by simid visit;
	run;
	/*For original sample, everyone should show up one time only. We create numberhits = 1 to ensure this*/
	%if &bsample = 0 %then %do;
		data bootsample;
		set bootsample;
			numberhits = 1;
		run;
	%end;

	/*Run pooled logistic regression model. We no longer need the REPEATED statement, but now we need to save the parameter estiamtes*/
	proc genmod data = bootsample descending;
	class simid;
		model death = visit visit2 rand /link = logit dist = bin;
		freq numberhits;
		ods output ParameterEstimates = unadj_plr_fit;
	run;
	proc contents data = unadj_plr_fit;
	run;
	/*Store output from each loop*/
	data unadj_plr_fit_ITT;
		set unadj_plr_fit;
		where PARAMETER = "rand";
		bsample = &bsample;
		keep bsample PARAMETER ESTIMATE;
	run;
	data means_all;
		set means_all unadj_plr_fit_ITT;
		by bsample;
	run;
%end;

/*Exponentiate parameter estimate for RAND*/
data means_all;
set means_all;
	HazardRatio = exp(ESTIMATE);
run;
/*Calculate 2.5th and 97.5th percentiles from bootstraps*/
proc univariate data = means_all (where = (bsample >0));
	by parameter;
	var HazardRatio;
	output out = pctls (keep = parameter HR_2_5 HR_97_5) pctlpre = HR_ pctlpts = 2.5, 97.5;
run;
data sample0;
	set means_all (where = (bsample = 0));
	keep parameter hazardratio;
run;
data final;
	merge sample0 pctls;
	by parameter;
	LowerBound = HR_2_5;
	UpperBound = HR_97_5;
run;
proc print data = final;
	var HazardRatio LowerBound UpperBound;
	title "Unadjusted Intention-to-Treat hazard ratio with 95% CI from &nboot bootstrap samples";
run;
%mend;
%ITT_boots(inset = trial, nboot = 10);
