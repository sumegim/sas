libname sasdata "/folders/myfolders/assignment/";

data bchc;
	set sasdata.bchc;
run;


/* San Antonio has very high mortality rates, should we test for outlier/ommit it? */
proc sgplot data=bchc;
   scatter x=place y=mortality;
run;

/* Boxplot suggest 2 outliers happening in every race - need to check if it's the same city 
in that case it would be reasonable to get rid of these records */
proc sgplot data=bchc;
    vbox mortality/ category=Race_Ethnicity;
run;

data race;
	set bchc;
/* 	where Race_Ethnicity = 'White'; */
/* 	where Race_Ethnicity = 'Black'; */
/* 	where Race_Ethnicity = 'Hispanic'; */
	where Race_Ethnicity = 'Asian/PI';
run;

/* use to remove previous runs' outlier  */
/* data race; */
/* 	set race; */
/* 	where mortality ^= 169.2; */
/* run; */

proc means data=race mean std n;
	var mortality;
	output out=L10_sumstat mean=mean median=median std=std n=n;
run;
data L10;
	set race;
	if _n_=1 then set L10_sumstat;
	drop _TYPE_ _FREQ_;
run;

/* here we get rid of outlier */
/* data L10; */
/* 	set L10; */
/* 	where mortality ^= 200; */
/* run; */

data DOORNBOS;
	SET L10;
	U = (mortality-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),	1-CDF('T',DOORNBOS_Y, N-2))*N, 1);
RUN;

PROC SORT DATA =DOORNBOS;
BY DOORNBOS_Y;
RUN;

PROC PRINT DATA=DOORNBOS;
RUN;

proc summary data=race q1 q3;
	var mortality;
	output out=race_summary q1=q1 q3=q3;
run;

proc iml;
	use race;
	read all var{mortality place race_ethnicity year};
	close race;
	
	use race_summary;
	read all var{q1 q3};
	close race_summary;
	
	IQR = q3 - q1;
	lower = q1 - 1.5*IQR;
	upper = q3 + 1.5*IQR;
	
	total = nrow(bw);
	is_outlier = loc(mortality<lower | mortality>upper);
	outliers_m = mortality[is_outlier];
	outliers_p = place[is_outlier];
	outliers_y = year[is_outlier];
/* 	outliers = outliers_m||outliers_p; */
	
	print(total);
	print(outliers_m);
	print(outliers_p);
	print(outliers_y);
	
/* 	no_outliers = loc(bw>=lower & bw<=upper); */
/* 	output=bw[no_outliers]; */
/* 	 */
/* 	create BW_NO_OUTL from output[colname={"bw"}]; */
/* 	append from output; */
/* 	close BW_NO_OUTL; */
quit;

/* Indeed, San Antonio, TX is an outlier for all races */