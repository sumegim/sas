libname sasdata "/folders/myfolders/assignment/";

data bchc;
	set sasdata.bchc;
/* 	where (Race_ethnicity = 'Black' and place ^= 'San Antonio, TX') or  */
/* 		(Race_ethnicity = 'White' and place ^= 'San Antonio, TX') or  */
/* 		(Race_ethnicity = 'Hispanic' and place ^= 'San Antonio, TX') or Race_ethnicity = 'Asian/PI'; */
run;

data race;
	set bchc;
	where Race_ethnicity = 'Black';
run;

PROC UNIVARIATE DATA=race NORMAL;
	VAR mortality;
	PROBPLOT mortality /NORMAL(MU=est SIGMA=est);
RUN;