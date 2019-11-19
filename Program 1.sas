
LIBNAME SAS "/folders/myfolders";

/*Question 1 */ 


DATA WEEK1;
 SET SAS.IVF;
 WHERE PER=4;
 AGEMB= (AGEM<30);
 *KEEP AGEM AGEMB;
 DROP IMP PER AGE
RUN;

*Export .csv file to specified folder;
PROC export data = WEEK1
outfile=  "/folders/myfolders/week1.csv"
dbms= csv
replace;
run;

*Import .csv file from specified folder;
PROC IMPORT OUT= WORK.week12
  DATAFILE= "/folders/myfolders/week1.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
RUN;

*Q1.1
*count the number of people int he control group, in treatmen 1 and treatment 2;
proc freq data=WEEK1;
  table TRT;
run;

*Q1.2;
*a;
PROC IML;
use WEEK1;
read all var{agem};
close WEEK1;

alpha=0.05;

Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
qT=quantile('t',alpha/2,n-1);
UCL=Ybar-qT*sqrt(s/n);
LCL=Ybar+qT*sqrt(s/n);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);

print("mean and frequency:");
print(Ybar||s);
print("LCL and UCL:");
print(LCL||UCL);
print("LPL and UPL:");
print(LPL||UPL);

quit;

*b;
proc freq data=WEEK1;
  table AGEM;
run;



















