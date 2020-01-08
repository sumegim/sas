libname sasdata "/folders/myfolders/assignment/";

data bchc;
	set sasdata.bchc;
run;

PROC SORT DATA=bchc;
	BY Place Race_Ethnicity Year;
RUN;

/* Analyzing means of mortality */
PROC MEANS DATA=bchc NOPRINT;
	VAR Mortality;
	BY Place Race_Ethnicity;
	OUTPUT OUT=MEANS MEAN=Mortality;
RUN;

/* Annova 2 ways model, Race - Fixed effect, Place - Random effect */
ods output SolutionR=EBLUP;
PROC MIXED DATA=MEANS METHOD=TYPE3;
	CLASS Place Race_Ethnicity;
	MODEL Mortality=Race_Ethnicity/SOLUTION outp=PRED;
	RANDOM PLACE/SOLUTION;
	LSMEANS Race_Ethnicity;
RUN;

/* Testing randomness of EBLUPS */
%RUNSCUC(data=EBLUP,var=Estimate,alpha=0.05);

/* verify the ANOVA assumption that the residuals are independently distributed */
proc sort data=PRED;
	by Place Race_Ethnicity;
run;

%RUNSCUC(data=PRED,var=RESID,alpha=0.05);

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
