     /* ------------------- */
    /*      Prototypen     */
   /* ------------------- */

#if defined(__STDC__)

/* Speicherverwaltung, Numerical Recipies */
void 	nrerror		(char[]);
int 	*ivector	(long,long);
void 	free_ivector	(int*, long, long);
double 	*dvector	(long,long);
void 	free_dvector	(double*,long,long);
double 	**dmatrix	(long, long, long, long);
void 	free_dmatrix	(double**, long, long, long, long);
void 	upcase		(char*);
void 	lowcase		(char*);

void 	init_quick_075	(void);
double quick_075	(double);
double schicht_power	(double*, double*, double*);

/* nl-Optimierung nach modifiziertem Levenberg-Marquardt */
int 	f_deriv		(int(*)(), double*, double*, int, double*, int*, int, double**);
double 	dpythag		(double,double);
int 	dsvdcmp		(double**, int, int, double*, double**);
void 	dsvbksb		(double**, double*, double**, int, int, double[]);
void 	dsvdinv		(double**, double*, double**, int, double**);
int 	dsvd_solver	(double**, double*, int, int);

void 	covsrt		(double**, int, int*, int);
void 	calc_chisqr	(double*, double*, double*, int, double*, int*, int, double*, double*, double*);
int 	mrqmin		(int(*)(), double*, double*, double*, double*, int, double*, int*, int, int, double*, double*, double**, double**, double*, double*);
int 	magfit		(int(*)(), double*, double*, double*, double*, int*, double*, int*, int*, double*, double*, double*, int*, double*, double*, int*);
int 	mag_solver	(double*, double*, double*, int*, int*, double*, int*, double*, double*, double*, double*);

/* nl-Gleichungssolver, Newton mit variabler Schrittweite */
double 	calc_fmin	(double*, int, int(*)(), double*, double*, double*, double*);
int 	fdjac		(int, double*, double*, double**, int(*)(), double*, double*);
int 	lnsrch		(int, double*, double, double*, double*, double*, double*, double, int*,double*, double*, double*, double*);
int 	newt		(double*, int, int*, double*, double*, double*, double*);

/* Kennlinien */
int 	mag_doppel	(double*, double*, int, double*);
int 	fast_doppel	(double*, double*, int, double*);
int 	resist_doppel	(double*, double*, int, double*);
int 	fast_doppel_nl	(double*, double*, double*, int,int, double*);


/* Interface zu IDL */
int 	run_triple_fit	(int, void**);
int 	idl2mag_err	(int, char*, int);
int 	idl2mag 	(int, void**);


#else

void 	nrerror		();
int 	*ivector	();
void 	free_ivector	();
double 	*dvector	();
void 	free_dvector	();
double 	**dmatrix	();
void 	free_dmatrix	();
void 	upcase		();
void 	lowcase		();

void 	init_quick_075	();
double quick_075	();
double schicht_power	();

int 	f_deriv		();
double 	dpythag		();
int 	dsvdcmp		();
void 	dsvbksb		();
void 	dsvdinv		();
int 	dsvd_solver	();

void 	covsrt		();
void 	calc_chisqr	();
int 	mrqmin		();
int 	magfit		();
int 	mag_solver	();

double 	calc_fmin	();
int 	fdjac		();
int 	lnsrch		();
int 	newt		();

int 	mag_doppel	();
int 	fast_doppel	();
int 	resist_doppel	();
int 	fast_doppel_nl	();


int 	run_triple_fit	();
int 	idl2mag_err	();
int 	idl2mag 	();

#endif	/* #if defined(__STDC__) */
