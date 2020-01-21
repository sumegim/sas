/* Q 1.1 */
/* What is the total number of children considered in this study? How many of these children */
/* belong to a mother in the control group (TRT = 0), and how many received treatment-M */
/* (TRT = 1) or treatment-C (TRT = 2)? */

LIBNAME DATA "/folders/myfolders/Datasets";

DATA WEEK1;
 SET DATA.IVF;
 where PER = 4; 		/*we only want one observation for each child*/
 drop IMP PER AGE; 		/*this is just for cleanup*/
RUN;

proc freq data=WEEK1;	/*getting the frequencies helps us answer how many kids in each treatment group*/
	tables TRT;			/*I only want the table that gets me the freq by tratment*/
run; 


/* Q 1.1 */
/* In this exercise you will analyze the age of mothers (AGEM) in the IVF dataset. For this */
/* question you may assume that the mother’s age is normally distributed. */
/* (a) Compute the mean and variance of AGEM. Use these to compute a 95% confidence interval */
/* for the average age of a mother, and a 95% prediction interval for a single new observation */
/* of AGEM. */
/* (b) How many mothers were 40 years old or older when they became pregnant? */
/* (c) Compute a 95% confidence interval for the variance of AGEM. */

/*Compute the mean and variance of AGEM.*/
proc means data=WEEK1 n mean var clm;
	var AGEM;
run;

/* Use these to compute a 95% confidence interval for the average age of a mother */
proc univariate data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
run;
/*  */
/* PROC IML; */
/* use WEEK1; */
/* read all var{agem}; */
/* close WEEK1; */
/*  */
/* alpha=0.05; */
/*  */
/* Ybar=mean(agem); */
/* s=var(agem); */
/* n=nrow(agem); */
/* qT=quantile('t',alpha/2,n-1); */
/* UCL=Ybar-qT*sqrt(s/n); */
/* LCL=Ybar+qT*sqrt(s/n); */
/*  */
/* A=Ybar||LCL||UCL; */
/*  */
/* print(Ybar||LCL||UCL); */
/*  */
/* create DATA from A[colname={'mean' 'LCL' 'UCL'}];  */
/* append from A;        */
/* close DATA; */
/* quit; */


/*95% prediction interval for a single new observation*/
proc reg data=WEEK1;
     model AGEM= / cli alpha=0.05;
run;

/* PROC IML; */
/* use WEEK1; */
/* read all var{agem}; */
/* close WEEK1; */
/*  */
/* alpha=0.05; */
/* Ybar=mean(agem); */
/* s=var(agem); */
/* n=nrow(agem); */
/* qT=quantile('t',alpha/2,n-1); */
/* UPL=Ybar-qT*sqrt((n+1)*s/n); */
/* LPL=Ybar+qT*sqrt((n+1)*s/n); */
/* A=Ybar||LPL||UPL; */
/*  */
/* create DATA from A[colname={'mean' 'LPL' 'UPL'}];  */
/* append from A;        */
/* close DATA; */
/* quit; */

/* (b) How many mothers were 40 years old or older when they became pregnant?  */
DATA WEEK1_oldmoms;
 SET WEEK1;
 AGEMB= (AGEM>=40);
RUN;

proc freq data=WEEK1_oldmoms;
	tables AGEMB;
run; 

/* (c) Compute a 95% confidence interval for the variance of AGEM. */
proc univariate data=WEEK1 cibasic;
	var AGEM;
	histogram AGEM/ normal;
run;

/* q 1.3 */
/*In this excercise you will analyze the birth weight (BW)*/
/* compute the first and third quartile of the birth weight together with a corresponding */
/* 95% confidence interval. What is the interquartile range of the birth weight? */
/* proc univariate data=WEEK1 cipctlnormal alpha=0.05; */
/* 	var BW; */
/* run; */
/* Cant use this because the birth weight data is not normally distributed */


*W1-2 CI percentile;
proc univariate data=WEEK1 cipctldf;
	var BW;
run;

/*1.3 b */
/* What is the percentage of children whose birth weight differs at most one interquartile range from the median? */


/* 1.3 c */
/* The birth weight data is not normally distributed. Draw a histogram of the Box-Cox */
/* transformed birth weight, use λ ∈ {−2, −1/2, 0, 1/2, 2}. What is the best value for λ */
/* such that the data becomes approximately normal? */


*W1-2 BOXCOX transformation;
PROC UNIVARIATE DATA=WEEK1 NORMALTEST;
VAR BW;
HISTOGRAM BW /NORMAL;
RUN;
 
DATA WEEK1BOXCOX;
 	SET WEEK1;
 	BWMINUS2= (-1/2)*((BW)**-2-1);
 	BWMINUS1= (-1)*((BW)**-1-1);
 	BWMINUS12= (-2)*((BW)**-(0.5)-1);
 	BW0= log(BW);
 	BW13= (3)*((BW)**(1/3)-1);
 	BW12= (2)*((BW)**(1/2)-1);
 	BW2= (0.5)*((BW)**(2)-1);
RUN;

proc univariate data=WEEK1BOXCOX noprint;
  	histogram BWMINUS2 /normal;
 	histogram BWMINUS1 /normal;
 	histogram BWMINUS12 /normal;
 	histogram BW0 /normal;
 	histogram BW13 /normal;
 	histogram BW12 /normal;
	histogram BW2 /normal;
run;

/* Looks like BW2 (lambda = 2) is the best for normality transformation */

PROC UNIVARIATE DATA=WEEK1BOXCOX NORMALTEST;
VAR BW2;
HISTOGRAM BW2 /NORMAL;
RUN;


/* 1.3 d Using your answer from (c), compute a 95% prediction interval for a single new observation. */
proc reg data=WEEK1BOXCOX;
     model BW2= / cli alpha=0.05;
run;
/* You have to transorm the results back manually afterwards !!!!!!!!!! */

/* 1.3 e */
/* What is the birth weight of the heaviest baby in this dataset, is it a boy or a girl? */

proc freq data=week1;
	tables=BW SEX;
run;
	















