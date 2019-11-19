LIBNAME SAS "/folders/myfolders";

/*Question 2 */ 


DATA WEEK1;
 SET SAS.IVF;
 WHERE PER=4;
 AGEMB= (AGEM<30);
 *KEEP AGEM AGEMB;
 DROP IMP PER AGE
RUN;

proc ttest data=WEEK1;
	class FIS;
	var AGEM;
run;

proc glm data=WEEK1;
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=BARTLETT;
run;

proc glm data=WEEK1;
	class FIS;
	model AGEM = FIS;
	means FIS / hovtest=Levene;
run;

