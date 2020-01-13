libname sasdata "/folders/myfolders/assignment/";

data bchc;
	set sasdata.bchc;
	where mortality and ((Race_ethnicity = 'Black' and place ^= 'San Antonio, TX') or 
		(Race_ethnicity = 'White' and place ^= 'San Antonio, TX') or 
		(Race_ethnicity = 'Hispanic' and place ^= 'San Antonio, TX') or Race_ethnicity = 'Asian/PI');
/* 	where mortality; */
run;

ods select Histogram;
proc univariate data=bchc;
	histogram mortality;
run;

/* PROC SORT DATA=bchc; */
/* 	BY Place Race_Ethnicity Year; */
/* RUN; */
/*  */
/* Analyzing means of mortality */
/* PROC MEANS DATA=bchc NOPRINT; */
/* 	VAR Mortality; */
/* 	BY Place Race_Ethnicity; */
/* 	OUTPUT OUT=MEANS MEAN=Mortality; */
/* RUN; */
/*  */
/* San Antonio has very high mortality rates, should we test for outlier/ommit it? */
/* proc sgplot data=MEANS; */
/*    scatter x=place y=mortality; */
/* run; */

/* Boxplot suggest 2 outliers happening in every race - need to check if it's the same city 
in that case it would be reasonable to get rid of these records */
/* proc sgplot data=bchc; */
/*     vbox mortality/ category=Race_Ethnicity; */
/* run; */
/*  */
/* ods output SolutionR=EBLUP; */
/* PROC MIXED DATA=bchc METHOD=TYPE3; */
/* 	CLASS Place; */
/* 	MODEL Mortality= /SOLUTION outp=PRED; */
/* 	RANDOM PLACE  /SOLUTION; */
/* 	LSMEANS Race_Ethnicity; */
/* RUN; */

/* Annova 2 ways model, Race - Fixed effect, Place - Random effect */
ods output SolutionR=EBLUP;
PROC MIXED DATA=bchc METHOD=TYPE3;
	CLASS Place Race_Ethnicity;
	MODEL Mortality=Race_Ethnicity /SOLUTION ddfm=sat outp=PRED;
	RANDOM PLACE  /SOLUTION;
	LSMEANS Race_Ethnicity / CL;
RUN;

/* Check normality of residuals */
ods listing gpath='/folders/myfolders/assignment/';
ods graphics / imagename="residuals_probplot" imagefmt=png;
PROC UNIVARIATE DATA=PRED NORMAL;
	VAR RESID;
	PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
RUN;

/* Check normality of random effect */
ods listing gpath='/folders/myfolders/assignment/';
ods graphics / imagename="eblups_probplot" imagefmt=png;
PROC UNIVARIATE DATA=EBLUP NORMAL;
	VAR Estimate;
	PROBPLOT Estimate /NORMAL(MU=est SIGMA=est);
RUN;

proc glm data=PRED;
	class Race_Ethnicity;
	model Resid = Race_Ethnicity;
	means Race_Ethnicity / hovtest=BARTLETT;
run;

%RUNSCUC(data=PRED,var=RESID,alpha=0.05);

%RUNSCUC(data=EBLUP,var=Estimate,alpha=0.05);

/* In case of fixed effects hidden in residuals - investigate by city */
/* proc sgplot data=pred; */
/*    scatter x=place y=resid; */
/* run; */
/*  */
/* In case of random effects included in cities */
/* proc sgplot data=eblup; */
/*    scatter x=place y=estimate; */
/* run; */
/*  */
/* Testing randomness of EBLUPS */
/* %RUNSCUC(data=EBLUP,var=Estimate,alpha=0.05); */
/*  */
/* ods select Histogram; */
/* proc univariate data=pred; */
/* 	histogram resid/normal; */
/* run; */
/*  */
/* verify the ANOVA assumption that the residuals are independently distributed */
/* proc sort data=PRED; */
/* 	by Place Race_Ethnicity; */
/* run; */
/*  */
/* %RUNSCUC(data=PRED,var=RESID,alpha=0.05); */
/*  */
/* Are mortality values in different cities comparable? */
/* PROC SORT DATA=bchc; */
/* 	BY Place; */
/* RUN; */
/*  */
/* %RUNSCUC(data=bchc,var=mortality,alpha=0.05); */
/*  */
/* ods select Histogram; */
/* proc univariate data=pred; */
/* 	histogram resid/normal; */
/* run; */

data black;
	set bchc;
	where Race_Ethnicity = 'Black';
run;

data white;
	set bchc;
	where Race_Ethnicity = 'White';
run;

data hispanic;
	set bchc;
	where Race_Ethnicity = 'Hispanic';
run;

data asian;
	set bchc;
	where Race_Ethnicity = 'Asian/PI';
run;

/* Compare destributions between races */
ods select Histogram;
proc univariate data=black;
	histogram mortality/normal;
run;

ods select Histogram;
proc univariate data=white;
	histogram mortality/normal;
run;

ods select Histogram;
proc univariate data=hispanic;
	histogram mortality/normal;
run;

ods select Histogram;
proc univariate data=asian;
	histogram mortality/normal;
run;

/* compare distributions between years */
ods select Histogram;
proc univariate data=bchc;
	where year = 2014;
	histogram mortality/normal;
run;

ods select Histogram;
proc univariate data=bchc;
	where year = 2015;
	histogram mortality/normal;
run;

ods select Histogram;
proc univariate data=bchc;
	where year = 2016;
	histogram mortality/normal;
run;


proc iml;
	sigma_g = 1209.16;
	sigma_e = 471.79;
	
	icc = (sigma_g**2)/((sigma_g**2) + (sigma_e**2));
	print(icc);
quit;
