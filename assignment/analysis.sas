libname sasdata "/folders/myfolders/assignment/";
%include "/folders/myfolders/assignment/macro_runstest.sas";

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
/* ods output SolutionR=EBLUP; */
/* PROC MIXED DATA=MEANS METHOD=TYPE3; */
/* 	CLASS Place Race_Ethnicity; */
/* 	MODEL Mortality=Race_Ethnicity/SOLUTION outp=PRED; */
/* 	RANDOM PLACE/SOLUTION; */
/* 	LSMEANS Race_Ethnicity; */
/* RUN; */



proc means data=bchc mean var std nway;
class Place Race_Ethnicity;
var Mortality;
run;

proc sgplot data=bchc;
vline Place / group=Race_Ethnicity stat=mean response=Mortality markers;
run;

proc glm data=bchc plots=diagnostics;
class Place Race_Ethnicity;
model Mortality=Race_Ethnicity Place Race_Ethnicity * Place;
random Place;
lsmeans Race_Ethnicity * Place / diff slice=Place;
store out=interact;
run;
quit;

proc plm restore=interact plots = all;
slice Race_Ethnicity * Place / sliceby=Race_Ethnicity adjust=tukey;
effectplot interaction(sliceby=Race_Ethnicity) / clm;
run;


PROC UNIVARIATE DATA=PRED NORMAL;
    VAR RESID;
    PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
RUN;

proc glm data=PRED;
    class Race_Ethnicity;
    model Resid = Race_Ethnicity;
    means Race_Ethnicity / hovtest=BARTLETT;
run;


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
