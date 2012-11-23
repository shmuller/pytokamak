#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#define INT32		int
#define IDL_SHORT	short int
#define IDL_USHORT	unsigned short int
#define IDL_INT		int
#define IDL_FLOAT	float
#define IDL_DOUBLE	double


/* struct idl_string_struct
	{
	IDL_USHORT	slen;
	IDL_USHORT	stype;
	char		*s;
	}; */

typedef struct 
	{
	IDL_USHORT	slen;
	IDL_USHORT	stype;
	char		*s;
	} idl_string_struct ;

#define IDL_STR(x) (((idl_string_struct*)x)->s)
#define IDL_STR_LEN(x) ((int)((idl_string_struct*)x)->slen)


/* physikal. Konstanten */
static double	    c_e    = 1.6022e-19;
static double	    c_eps0 = 8.8542e-12;
static double	    c_me   = 9.1095e-31;
static double	    c_mp   = 1.6726e-27;
static double	    c_quot = 41.3124;  		/* = 0.15^2 * c_mp / c_me */
static double	    c_pi   = 3.1415927;
static double	    c_hpi  = 0.5 * 3.1415927;
static double	    c_dpi  = 3.1415927 + 3.1415927;



/* schnelle Berechnung von ^0.75 durch look-up-table */
static int	n_look_up = 400;
static int 	n_look_max = 398;           /* = n_look_up -2 */
static double	delta_look_up = -0.1;
static double	inv_delta_look_up = -10.;   /* 1/ delta_look_up */
static double	look_up[1000];		    /* eigentlich: look_up[n_look_up] */
static double	look_up_alpha[1000], look_up_beta[1000];
static double	*ptr_alpha, *ptr_beta;

/* look-up-table fuer exp-Funktion */
static int	exp_lt_n = 20002, exp_lt_ninits=0;
static double   exp_lt_min = -10.0, exp_lt_max = 10.0;
static double	exp_lt_delta = 0.001, exp_lt_inv = 1000.0;
static double	exp_lt_alpha[20000], exp_lt_beta[20000];
static double 	*exp_lt_ptr_a, *exp_lt_ptr_b;



/* Funktionen zur Berechnung der Kennlinien, globale Definitionen */
int	mag_doppel();
int	fast_doppel();


#define NR_END 1
#define FREE_ARG char*

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}


/* Numerical Recipes standard error handler */
void nrerror(error_text)
char error_text[];
{
	fprintf(stderr,"\n\nrun-time error: %s\n", error_text);
	fprintf(stderr,"\t --> now exiting to system...\n\n\n");
	exit(1);
}



/* allocate a int vector with subscript range v[nl..nh] */
int *ivector(nl, nh)
long	nl;
long	nh;
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}


/* free a double vector allocated with ivector() */
void free_ivector(v, nl, nh)
int	*v;
long	nl;
long	nh;
{
	free((FREE_ARG) (v+nl-NR_END));
}


/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(nl, nh)
long	nl;
long	nh;
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}


/* free a double vector allocated with dvector() */
void free_dvector(v, nl, nh)
double	*v;
long	nl;
long	nh;
{
	free((FREE_ARG) (v+nl-NR_END));
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)
				((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}



/* free a double matrix allocated by dmatrix() */
void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}



/*	string conversion to upper and to lower case */
int upcase( str )
char	*str;
   {
   int	i, len;
   char *cptr;
   len = strlen(str);
   cptr = str;
   for (i=0; i<len; i++, cptr++) *cptr=toupper(*cptr);
   }	/* upcase */
   
int lowcase( str )
char	*str;
   {
   int	i, len;
   char *cptr;
   len = strlen(str);
   cptr = str;
   for (i=0; i<len; i++, cptr++) *cptr=tolower(*cptr);
   }	/* lowcase */

/*	--------------------------------------------------------

			init_exp_lt
			=========

	Initialisierung der look-up Tabelle fuer die lineare
	Interpolation bei der Berechnung von exp(x)
 	-------------------------------------------------------- 

	M. Weinlich, 09.06.95

 	-------------------------------------------------------- */

int init_exp_lt()


    {
	double	*ptr1, *ptr2;
	double 	x1, y1, x_div, temp;
	int	i;

	/* Initialisierung nur einmal */
	if (exp_lt_ninits > 0 ) return;
	exp_lt_ninits++;

	/* Initialisierung der globalen Pointer */
	exp_lt_ptr_a = exp_lt_alpha;
	exp_lt_ptr_b = exp_lt_beta;

	/* Initialisierung der Arrays */
	x1   = exp_lt_min;
	y1   = exp(x1);
	ptr1 = exp_lt_alpha;
	ptr2 = exp_lt_beta;
	x_div = exp_lt_min * exp_lt_inv;

	for (i=0 ; i<exp_lt_n; i++  ) {

	    /* Funktionswert bei naechstem x */
	    temp = exp(x1 += exp_lt_delta) - y1;

	    *ptr1++ = y1 - temp * (x_div++)  ;
	    *ptr2++ = temp * exp_lt_inv;

	    /* Inkrementieren der Variablen */
	    y1 += temp; 

	    }

    }		/* init_exp_lt */





/*	--------------------------------------------------------

			exp_lt
			=========

	Naeherung fuer exp(x) im Intervall [exp_lt_min,exp_lt_max]


 	-------------------------------------------------------- 

	M. Weinlich, 09.06.95

 	-------------------------------------------------------- */

double exp_lt(x)

    double	x;

    {
    int	  	index;


	
    if ( x <= exp_lt_min ) {
	if ( x <= -35. ) return 0.0;
	return exp(x);
	}

    if ( x >= exp_lt_max ) return exp(x);

    index = (int) ( (x-exp_lt_min) * exp_lt_inv);

    /* Kontrolle der Berechnung */
    /* printf( "exp_lt:  f(x=%g) = %g ( soll = %g) \n", 
		x, *(exp_lt_alpha+index) + *(exp_lt_beta+index)*x, exp(x)); */

    /* lineare Interpolation fuer y */
    return ( *(exp_lt_alpha+index) + *(exp_lt_beta+index)*x ) ;
	
    }		/* exp_lt */

	




/*	-------------------------------------------------	*/
/*	-------------------------------------------------	*/
/*								*/
/*								*/
/*		Berechnung einer Sondenkennlinie		*/
/*		nach der Fitfunktion "mag_doppel"		*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*	-------------------------------------------------	*/





/*	--------------------------------------------------------

			init_quick_075
			=========

	Initialisierung der look-up Tabelle fuer die lineare
	Interpolation bei der Berechnung von x^0.75

 	-------------------------------------------------------- 

	M. Weinlich, 19.09.93

 	-------------------------------------------------------- */

int init_quick_075()


    {
	double	*ptr1, *ptr2;
	double 	x, y1, y2, temp;
	int	i;


	/* Initialisierung der globalen Pointer */
	ptr_alpha = look_up_alpha;
	ptr_beta  = look_up_beta;

	/* Initialisierung der Arrays */
	x    = 0.0;
	y1   = 0.0;
	ptr1 = ptr_alpha;
	ptr2 = ptr_beta;

	for (i=0 ; i<n_look_up; i++  ) {

	    /* Funktionswert bei naechstem x */
	    x += delta_look_up;

	    /* Schichtdicke = const * (phi)^0.75 (Child-Langmuir) */
	    /* y2 = pow( x, 0.75);  */

	    /* Schichtdicke fuer heisse monoenerget. Ionen (M. Weinich) */
	    temp = sqrt( 1. - x-x );
	    temp = ( temp + 2. ) * sqrt( temp - 1. ) - y1;

	    *ptr1++ = y1 - temp * (double)i;
	    *ptr2++ = temp * inv_delta_look_up;

	    /* Inkrementieren der Variablen */
	    y1 += temp;

	    }

    }		/* init_quick_075 */





/*	--------------------------------------------------------

			quick_075
			=========

	Naeherung fuer x^0.75 im Intervall [0.0, 100.0]

	Es wird eine Tabelle angelegt, die x-Werte sind um delta_look_up
	voneinander getrennt, die y-Werte werden nur bei Bedarf berechnet.
	Der eigentliche Funktionswert wird linear zwischen zwei benachbarten 
	und exakt berechneten Werten interpoliert.

       ein Zeit-Vergleich ergibt, dass die umstaendlich erscheinende Version
       ( Nachschauen in Tabelle und lineare Interpolation zwischen zwei bereits
       abgespeicherten Werten ) um fast 30% schneller ist als die direkte
       Berechnung der jeweiligen Potenz.
       Verwendete Parameter: etwa 70000 Aufrufe, n_look_up = 500
    

 	-------------------------------------------------------- 

	M. Weinlich, 19.09.93

 	-------------------------------------------------------- */

double quick_075(x)

    double	x;

    {
    int	  	index;


    /* keine negativen Argumente erlaubt */
    if ( x >= 0. ) return 0.0;

    /* die Werte sind im Abstand delta_x = delta_look_up tabelliert */
    /* index  ausserhalb des erlaubten Bereichs? */
    index = IMIN( (int) ( x * inv_delta_look_up), n_look_max);

    /* Kontrolle der Berechnung */
    /* printf( "quick_075:  f(x=%g) = %g \n", 
		x, *(ptr_alpha+index) + *(ptr_beta+index)*x); */

    /* lineare Interpolation fuer y */
    return ( *(ptr_alpha+index) + *(ptr_beta+index)*x ) ;
	
    }		/* quick_075 */

	


/*	--------------------------------------------------------

			schicht_power
			=============

	Es wird die Abhaengigkeit der Schichtdicke von der 
	anliegenden Sonden-/Wandspannung exakt berechnet.
	Exakt natuerlich nur im Rahmen des vorliegenden Modells.
	Eine Naeherung ist wegen der Abhaengigkeit von dem
	Spannungsabfall in der magnetischen Schicht und dem
	elektrischen Feld an der Grenze zur Debye-Schicht 
	momentan nicht realisiert, da eine Interpolation in 
	zwei Parametern keine wesentliche Zeitersparnis mehr
	mit sich bringt.    

 	-------------------------------------------------------- 

	M. Weinlich, 23.11.94

 	-------------------------------------------------------- */

double schicht_power( phi_sonde, phi_mag, eps )

    double	phi_sonde;
    double	phi_mag;
    double	eps;

    {
    double 	temp ;
    double      eps_square, eta;
    int	  	index;


    /* keine negativen Argumente unter Wurzel erlaubt */
    temp = phi_mag - phi_sonde;
    temp = temp+temp;

    if ( temp < 1. ) return 0.0;

    /* berechne eps^2 */
    eps_square = eps*eps;

    /* berechne eta */
    eta = sqrt( sqrt(1.+temp) - 1. + eps_square );


    /* Potentialabhaengiger Schichtfaktor */
    temp = ( eta*eta + 3.*(1.-eps_square) ) * eta + ( 2.*eps_square-3.)*eps; 
 
    /* Schichtdicke kann nicht negativ werden */
    if ( temp < 0 ) temp = 0.0;

    /* das war's */
    return ( temp ) ;
	
    }		/* schicht_power */

	




/* 	-------------------------------------------------------- 

 			mag_doppel 
			==========

 	berechnet einen Punkt aus Kennlinie einer beidseitig  
 	nichtsaettigenden Doppelsonde 

 	params(0)	Elektronendichte 	log(n_e) 
 	params(1)	Elektronentemperatur 	log(T_e) 
 	params(2)	Korrekturfaktor		log(alpha_1) 
 	params(3)	Flaechenverhaeltnis 	log(beta) 
 	params(4)	Unterschied in V_plas. 	delta_V 
 	params(5)	Korrrekturfaktor 	log(alpha_2) 
 	params(6)	Winkel an Sonde		cos(psi1) 
 	params(7)	Winkel an Gegenelektr.	cos(psi2) 

 	params(10)	Sondenbreite		a_width 
 	params(11)	Sondenhoehe		a_height 
 	params(12)	Massezahl des Fuellg.	a_mass 
	params(13)      Kernladungszahl d. FG.  a_z
	params(14)	Betrag des lok. B's	a_b_loc

        params(15)      Potetialdifferenz z. W. d_phi_w

	params(16) 	Verhaeltnis T_e / T_e	a_tau


	--> Aenderung auf gamma = 3 (Adiabatenkoeffizient)
	    wg. Riemann

 	-------------------------------------------------------- 

	M. Weinlich, 18.01.94

 	-------------------------------------------------------- */

int mag_doppel(vpr, strom, n_vpr, params)

    double 	*vpr;
    double 	*strom;
    int 	n_vpr;
    double 	*params;

    {

    /* Local variables */
    double 	n_e, t_e, beta, alpha1, alpha2, dv;
    double 	a_height, a_width, a_mass, a_proj_area ;
    double 	a_z, a_b_loc, a_tau;
    double	psi1,psi2,cos1,sin1,tan1, cos2, sin2, tan2, tan1p, tan2p;
    double	delta_phi, phi_sonde, d_phi_mag1, d_phi_mag2, 
			d_phi_w, dv_by_te, last_phi ;
    double 	phi1, d_phi, dw1, dw2, zeta1, zeta2, rho_s;
    double 	eps1, eps2;
    double	j_isat, I_isat, j_quot, log_j_quot, cquot, log_nenner, 
	     		ns_add1, ns_add2, ld_by_l, c_i_bohm;
    double	sf1, sf2, h_add1, h_add2, h_diff_phis;
    double	ld, hl1, hl2, hb1, hb2 , delta1, delta2, lfac, bfac;
    double	temp;
    double	zeta;

    int 	i, ind;
    int		n_changes=0;
    int 	verbose;

    /* here we are ... */
    /* printf("Beginn mag_doppel\n"); */

    /* Ueberpruefung auf unsinnige Fit-Parameter:
 	0.002  <  t_e            <  400
        1.e7   <  n_e sqrt(t_e)  <  1.e28
        0.001  <  beta           <  60
        0.000  <  alpha          <  100
	0.000  <  cos(psi)       <  1.000
    */

    verbose = ( 0 == 1 );

    if ( params[1] < -6. ) {
        if (verbose) printf( "Fehler t_e : %g", exp(params[1]) ); 
	/* keine Aenderung in n_e als Folge */
	temp      = params[0] + 0.5*params[1];
        params[1] = -6.;
	params[0] = temp - 0.5*params[1];
        if (verbose) printf( " --> %g \n", exp(params[1]));
	n_changes++;
        }
    else if ( params[1] > 6. ) {
        if (verbose) printf( "Fehler t_e : %g", exp(params[1]) );
	/* keine Aenderung in n_e als Folge */
	temp      = params[0] + 0.5*params[1];
        params[1] = 6.;
	params[0] = temp - 0.5*params[1];
        if (verbose) printf( " --> %g \n", exp(params[1]));
	n_changes++;
        }

    if ( params[0] < 15. ) {
        if (verbose) printf( "Fehler n_e : %g", exp(params[0]-0.5*params[1]) );
	params[0] = 15.;
        if (verbose) printf( " --> %g \n", exp(params[0]-0.5*params[1])); 
	n_changes++;
        }
    else if ( params[0] > 65. ) {
        if (verbose) printf( "Fehler n_e : %g", exp(params[0]-0.5*params[1]) );
	params[0] = 65.;
        if (verbose) printf( " --> %g \n", exp(params[0]-0.5*params[1]));
	n_changes++;
        }

    if ( params[3] < -7. ) {
	if (verbose) printf( "Fehler in beta: %g", exp(params[3]));
	params[3] = -7.;
        if (verbose) printf( " --> %g \n", exp(params[3]));
	n_changes++;
	}
    else if ( params[3] > 4. ) {
	if (verbose) printf( "Fehler in beta: %g", exp(params[3]));
	params[3] = 4.;
        if (verbose) printf( " --> %g \n", exp(params[3]));
	n_changes++;
	}

    if ( params[2] < -50. ) {
	if (verbose) printf( "Fehler in alpha1: %g", exp(params[2]));
	params[2] = -50.;
        if (verbose) printf( " --> %g \n", exp(params[2]));
	n_changes++;
	}
    else if ( params[2] > 5. ) {
	if (verbose) printf( "Fehler in alpha1: %g", exp(params[2]));
	params[2] = 5.;
        if (verbose) printf( " --> %g \n", exp(params[2]));
	n_changes++;
	}

    if ( params[5] < -50. ) {
	if (verbose) printf( "Fehler in alpha2: %g", exp(params[5]));
	params[5] = -50.;
        if (verbose) printf( " --> %g \n", exp(params[5]));
	n_changes++;
	}
    else if ( params[5] > 5. ) {
	if (verbose) printf( "Fehler in alpha2: %g", exp(params[5]));
	params[5] = 5.;
        if (verbose) printf( " --> %g \n", exp(params[5]));
	n_changes++;
	}

    if ( params[6] < 0.001 ) {
        if (verbose) printf( "Fehler in cos1 : %g", params[6] );
        params[6] = 0.001;
        if (verbose) printf(" --> %g \n", params[6]);
	n_changes++;
        }
    else if ( params[6] > 0.999 ) {
        if (verbose) printf( "Fehler in cos1 : %g", params[6] );
        params[6] = 0.999;
        if (verbose) printf(" --> %g \n", params[6]);
	n_changes++;
        }

    if ( params[7] < 0.001 ) {
        if (verbose) printf( "Fehler in cos2 : %g", params[7] );
        params[7] = 0.001;
        if (verbose) printf(" --> %g \n", params[7]);
	n_changes++;
        }
    else if ( params[7] > 0.999 ) {
        if (verbose) printf( "Fehler in cos2 : %g", params[7] );
        params[7] = 0.999;
        if (verbose) printf(" --> %g \n", params[7]);
	n_changes++;
        }

/*    if ( n_changes == 1 ) {
	printf( "Achtung: insgesamt %d Parameter musste geaendert werden!\n",
		n_changes);
	}
    else if ( n_changes > 1 ) {
	printf( "Achtung: insgesamt %d Parameter mussten geaendert werden!\n",
		n_changes);
	if ( n_changes > 3 ) return (-21);
	}
*/

    if ( n_changes > 4 ) {
	/* printf( "Abbruch, da insgesamt %d Parameter ", n_changes); */
	/* printf( "geaendert werden mussten  !\n"); */
	return (-21); 
	}



   /* sind die neuen Parameter Z und B_loc sinnvoll belegt? */

   /* in den meiseten Faellen Wasserstoff als Fuellgas */
   if ( params[13] < 1. ) {
      params[13] = 1. ;
      }

   /* Standard: B_t = 2.0T */
   if (( params[14] < 0.6 ) & ( params[14] > -0.6 )) {
      params[14] = 2.;
      }
 
    /* Ionen- und Elektronentemperatur sind postiv */
    if ( params[16] < 0.1 ) {
	params[16] = 1.;
	}

    /* Zuordnung der Fit-Variablen */
    n_e    = exp(params[0]-0.5*params[1]);
    t_e    = exp(params[1]);
    alpha1 = exp(params[2]);
    beta   = exp(params[3]);
    dv     = params[4];
    alpha2 = exp(params[5]);

    cos1   = params[6];
    cos2   = params[7];
    
    /* wichtige Sonden-Dimensionen */
    a_width  = params[10]; 
    a_height = params[11];
    a_mass   = params[12];
    a_z      = params[13];


    /* Betrag des lokalen Magnetfelds am Ort der Sonde */
    a_b_loc  = fabs(params[14]); 

    /* T_i / T_e */
    a_tau    = params[16];

    sin1 = sqrt( 1. - cos1*cos1 );
    sin2 = sqrt( 1. - cos2*cos2 );
    tan1 = sin1/cos1;
    tan2 = sin2/cos2;

    a_proj_area = a_width*a_height*cos1;


    /* Berechnung der Ionensaettigungsstromdichte 
       ==========================================

       Sekundaerelektronenemissionskoeffizient g = 0.4

       c_i_Bohm = sqrt(  ( Z*T_e + 3 T_i ) / m_i )  wobei 3 = (n+2)/n, n=1
       # c_e_Bohm = sqrt(  (3 T_e + T_i/Z ) / m_e )

       j_isat = n_e * c_i_Bohm * e

       # j_esat = 1/4 n_e * c_e_Bohm * e ( 1-g)

       # j_quot = j_esat / j_isat; 
       # j_quot = (1-g)/4 * sqrt( (a_mass*c_mp) /c_me ); 
 

       neu: 23.11.94 (mnw)
       --> im Doppelsondenbild kommen keine beschleunigten Elektronen zur
	   Elektrode, sondern immer nur thermische Plasmaelektronen.   

       c_e_th = sqrt( (8 T_e) / (pi m_e) )
       j_e_th = 1/4 * e * n_e * c_e_th * (1-g)

       j_quot = j_e_th / j_i_sat =		 mit T_i = T_e
              = 1/4 (1-g) sqrt(8/pi) sqrt(m_p/m_e) sqrt( a_mass/(a_z+3) ) =
	      = 6.8376426 * sqrt( a_mass/(a_z+3) )

    */

    c_i_bohm = sqrt( (3.*a_tau+a_z)* c_e * t_e/ (a_mass*c_mp) );
    j_isat = n_e * c_e * c_i_bohm;
    I_isat = a_proj_area * j_isat;

    j_quot = 6.8376426 * sqrt( a_mass/(a_z+3.*a_tau) ) ;
    log_j_quot = log(j_quot);

    /* Berechne phi_wall: */
    dv_by_te = dv / t_e;
    d_phi_w = log( (1.0+beta) / (1.0+beta*exp(dv_by_te)) ) - log_j_quot;
    params[15] = d_phi_w;	/* Abspeichern fuer spaetere Verwendung */


    /* Debye-Laenge vor Schichten */
    ld = sqrt(t_e * c_eps0 / (n_e * c_e)) ;

    /* Einfluss der Schichten --> Reduzierung in l_debye u.ae. */
    temp  = sqrt( 1. + 3.*a_tau/a_z );
    zeta1 = temp * cos1;
    zeta2 = temp * cos2;

    /* Ionengyroradius bei Schallgeschwindigkeit */
    rho_s = c_i_bohm * ( a_mass * c_mp ) / ( c_e * a_b_loc ) ;

    /* an beiden Sonden Potentialdifferenz                             */
    /* V_wall - V_el.sheath = V_wall - V_plasma + ( V_mag1 - V_plasma) */
    d_phi_mag1 = log(zeta1)	; 	/* i.d.R. < 0.0 */
    d_phi_mag2 = log(zeta2)	; 	/* i.d.R. < 0.0 */

    /* Abschaetzung fuer das elektrische Feld am Eintritt in Debye-Schicht */
    eps1 = - log(zeta1) * ld / ( sqrt(zeta1+zeta1) * rho_s );
    eps2 = - log(zeta2) * ld / ( sqrt(zeta2+zeta2) * rho_s );

    /* Vorfaktoren der Schichtdicken */
    sf1 = ld * sqrt(2./ zeta1) * alpha1 / 3.;
    sf2 = ld * sqrt(2./ zeta2) * alpha2 / 3.;

    /* Potentialdifferenz in el.stat. Schicht an umgebender Wand (V=0)  */
    dw1 = schicht_power( d_phi_w,          d_phi_mag1, eps1 );
    dw2 = schicht_power( d_phi_w+dv_by_te, d_phi_mag2, eps2 );

    /* Hilfsgroessen, muessen nur einmal berechnet werden */
    h_add1 = sf1*dw1;
    h_add2 = sf2*dw2;

    lfac = 1. / ( a_height * sqrt(cos1*beta/cos2) );
    bfac = 1. / ( a_width * sqrt(cos1*beta/cos2) );

    tan1p = tan1 + 2.;
    tan2p = tan2 + 2.;

    /* ab jetzt Berechnung fuer alle Spannungswerte */
    for ( ind=0; ind < n_vpr; ind++ ) 
        {

        /* Spannung an Doppelsonde: phi_sonde */
        phi_sonde = vpr[ind] / t_e;
        /* Verschiebung von phi2 gegen phi1: phi2 = phi1 + d_phi */
        d_phi = dv_by_te - phi_sonde;
        /* Startwert fuer phi1 = V_sonde_1 - V_wall */
        if (phi_sonde < d_phi_w) 
	    phi1 = phi_sonde;
        else 
	    phi1 = d_phi_w;

        /* Hilfsgroessen, muessen nur einmal berechnet werden */
        /* h_diff_phis = d_phi_mag2 - d_phi; */
        /* nenner      = (cquot * exp(d_phi) + 1.) * j_quot ; */
        log_nenner  = log((1.0 + beta*exp(d_phi))) + log_j_quot; 


        /*   	------ eigentliche Iteration ------ 	*/

        i = 0;
        last_phi = phi1 + phi1;

        while( ( fabs(phi1/last_phi-1.) > .001 ) && (i < 20) ) 
	    {

	    /* Schleifenzaehler */
	    last_phi = phi1;
	    ++i;

	    /* Nichtsaettigungskoeffizienten muessen immer groesser als -1 */
	    delta1 = sf1*( schicht_power(phi1,d_phi_mag1,eps1) -dw1);
	    delta2 = sf2*( schicht_power(phi1+d_phi,d_phi_mag2,eps2) -dw2);
	    hl1 = delta1 / a_height;
	    hb1 = delta1 / a_width;
	    hl2 = delta2 * lfac;
	    hb2 = delta2 * bfac;
	    ns_add1 = 1. + hb1+hb1 + (hb1+hb1+1.)*hl1*tan1p ;
	    ns_add2 = 1. + hb2+hb2 + (hb2+hb2+1.)*hl2*tan2p ;
	    if (ns_add1 < 0.) ns_add1 = 0.;
	    if (ns_add2 < 0.) ns_add2 = 0.;

	    /* naechste Naeherung fuer phi1 */
	    temp = ns_add1 + beta*ns_add2;
	    if (temp <= 0. ) break;
	    phi1 = log(temp) - log_nenner;

            }  	/* while ( ... ) */


        /* Berechnung des Strom-Vektors */
	delta1 = sf1*( schicht_power(phi1,d_phi_mag1,eps1) - dw1 );
	hl1 = delta1 / a_height;
	hb1 = delta1 / a_width;
	ns_add1 = 1. + hb1+hb1 + (hb1+hb1+1.)*hl1*(2.+tan1) ;
        if (ns_add1 < 0.) ns_add1 = 0.;
        strom[ind] = -I_isat * (ns_add1 - j_quot * exp(phi1));

        }	/* for ind=0, n_vpr-1 */

    return 0 ;

    }		 /* mag_doppel */






/* 	-------------------------------------------------------- 

 			fast_doppel 
			==========

 	berechnet einen Punkt aus Kennlinie einer beidseitig  
 	nichtsaettigenden Doppelsonde 

  	--> In mag_doppel wird V_1 - V_pl un V_wall - V_pl getrennt
	    bestimmt. Sind diese beiden Potentialdifferenzen bekannt,
	    so kann die effektiv wirksame Schichtdicke d_1 - d_wall
	    berechnet werden, die ihrerseits wiederum die effektive
	    projizierte Sondenflaeche bestimmt.

	--> Die Bestimmung von V_1 - V_pl muss leider iterativ erfolgen,
	    was sehr viel Zeit kostet ( etwa 30% der Gesamtzeit in 
	    mag_doppel )

	--> Um diesen iterativen Aufwand zu vermeiden, wird in fast_doppel
	    d_1 - d_wall durch d(V_Sonde - V_floating) angenaehert. Nun
	    ist die projizierte Flaeche ohne groesseren Aufwand berechenbar.

 	params(0)	Elektronendichte 	log(n_e) 
 	params(1)	Elektronentemperatur 	log(T_e) 
 	params(2)	Korrekturfaktor		log(alpha_1) 
 	params(3)	Flaechenverhaeltnis 	log(beta) 
 	params(4)	Unterschied in V_plas. 	delta_V 
 	params(5)	Korrrekturfaktor 	log(alpha_2) 
 	params(6)	Winkel an Sonde		cos(psi1) 
 	params(7)	Winkel an Gegenelektr.	cos(psi2) 

 	params(10)	Sondenbreite		a_width 
 	params(11)	Sondenhoehe		a_height 
 	params(12)	Massezahl des Fuellg.	a_mass 

	params(13)	Kernladungszahl	d. F.g.	a_z
	params(14)	Toroidalfeld		B_t

        params(15)       Potetialdifferenz z. W. d_phi_w

	params(16)	T_i/T_e			a_tau

	--> Erweiterung, um auch Sonden bei senkrechtem Einfall behandeln
	    zu k"onnen (z.B. Mitelebenenmanipulator, Rohde). Es muss lediglich
	    die projizierte L"ange bzw. H"ohe angegeben werden, da die 
	    Sondenbreite a_width = params(13) bereits als senkrecht zum
	    Magnetfeld vorausgesetzt wird.

	params(17)	Sondenhoehe		p_height

	--> Aenderung auf gamma = 3 (Adiabatenkoeffizient)
	    wg. Riemann

 	-------------------------------------------------------- 

	M. Weinlich, 19.03.94

 	-------------------------------------------------------- */

int fast_doppel(vpr, strom, n_vpr, params)

    double 	*vpr;
    double 	*strom;
    int 	n_vpr;
    double 	*params;

    {

    /* Local variables */
    double 	n_e, t_e, inv_beta, alpha1, alpha2, dv;
    double 	a_height, a_width, a_mass;
    double	p_height;
    double	a_sonde_proj, a_sonde_perp;
    double	delta_phi ;
    double	j_isat, I_isat, ns_add1, ns_add2;
    double	sf1, sf2, lfac, bfac;
    double	ld, delta1, delta2;
    double	temp;

    double	*ptr_vpr, *ptr_strom;
    int 	i, ind;
    int		n_changes=0;

    static double	cos1=0., cos2=0., tan1p=1.e38, tan2p=1.e38;


    /* Ueberpruefung auf unsinnige Fit-Parameter:
 	0.002  <  t_e            <  400
        1.e7   <  n_e sqrt(t_e)  <  1.e28
        0.001  <  beta           <  60
        0.000  <  alpha          <  100
	0.000  <  cos(psi)       <  1.000
    */

    if ( params[1] < -6. ) {
        /* printf( "Fehler t_e : %g", exp(params[1]) ); */
	/* keine Aenderung in n_e als Folge */
	temp      = params[0] + 0.5*params[1];
        params[1] = -6.;
	params[0] = temp - 0.5*params[1];
        /* printf( " --> %g \n", exp(params[1])); */
	n_changes++;
        }
    else if ( params[1] > 6. ) {
        /* printf( "Fehler t_e : %g", exp(params[1]) ); */
	/* keine Aenderung in n_e als Folge */
	temp      = params[0] + 0.5*params[1];
        params[1] = 6.;
	params[0] = temp - 0.5*params[1];
        /* printf( " --> %g \n", exp(params[1])); */
	n_changes++;
        }

    if ( params[0] < 15. ) {
        /* printf( "Fehler n_e : %g", exp(params[0]-0.5*params[1]) ); */
	params[0] = 15.;
        /* printf( " --> %g \n", exp(params[0]-0.5*params[1])); */
	n_changes++;
        }
    else if ( params[0] > 65. ) {
        /* printf( "Fehler n_e : %g", exp(params[0]-0.5*params[1]) ); */
	params[0] = 65.;
        /* printf( " --> %g \n", exp(params[0]-0.5*params[1])); */
	n_changes++;
        }

    if ( params[3] < -7. ) {
	/* printf( "Fehler in beta: %g", exp(params[3])); */
	params[3] = -7.;
        /* printf( " --> %g \n", exp(params[3])); */
	n_changes++;
	}
    else if ( params[3] > 4. ) {
	/* printf( "Fehler in beta: %g", exp(params[3])); */
	params[3] = 4.;
        /* printf( " --> %g \n", exp(params[3])); */
	n_changes++;
	}

    if ( params[2] < -50. ) {
	/* printf( "Fehler in alpha1: %g", exp(params[2])); */
	params[2] = -50.;
        /* printf( " --> %g \n", exp(params[2])); */
	n_changes++;
	}
    else if ( params[2] > 5. ) {
	/* printf( "Fehler in alpha1: %g", exp(params[2])); */
	params[2] = 5.;
        /* printf( " --> %g \n", exp(params[2])); */
	n_changes++;
	}

    if ( params[5] < -50. ) {
	/* printf( "Fehler in alpha2: %g", exp(params[5])); */
	params[5] = -50.;
        /* printf( " --> %g \n", exp(params[5])); */
	n_changes++;
	}
    else if ( params[5] > 5. ) {
	/* printf( "Fehler in alpha1: %g", exp(params[5]));  */
	params[5] = 5.;
        /* printf( " --> %g \n", exp(params[5]));  */
	n_changes++;
	}

    if ( params[6] < 0.001 ) {
        /* printf( "Fehler in cos1 : %g", params[6] ); */
        params[6] = 0.001;
        /* printf(" --> %g \n", params[6]); */
	n_changes++;
        }
    else if ( params[6] > 0.999 ) {
        /* printf( "Fehler in cos1 : %g", params[6] ); */
        params[6] = 0.999;
        /* printf(" --> %g \n", params[6]); */
	n_changes++;
        }

    if ( params[7] < 0.001 ) {
        /* printf( "Fehler in cos2 : %g", params[7] ); */
        params[7] = 0.001;
        /* printf(" --> %g \n", params[7]); */
	n_changes++;
        }
    else if ( params[7] > 0.999 ) {
        /* printf( "Fehler in cos2 : %g", params[7] ); */
        params[7] = 0.999;
        /* printf(" --> %g \n", params[7]); */
	n_changes++;
        }

    if ( n_changes > 3 ) {
	/* printf( "Abbruch, da insgesamt %d Parameter ", n_changes); 
	printf( "geaendert werden mussten  !\n");  */
	return (-21);
	}

    /* Zuordnung der Fit-Variablen */
    n_e    = exp(params[0]-0.5*params[1]);
    t_e    = exp(params[1]);
    alpha1 = exp(params[2]);
    inv_beta = exp(-params[3]) ;
    dv     = params[4];
    alpha2 = exp(params[5]);

    /* Winkelfunktionen muessen nur neu berechnet werden, wenn sich
       auch die Winkel geaendert haben */
    if (cos1 != params[6]) {
        cos1   = params[6];
        tan1p  = sqrt(1./(cos1*cos1) -1.) + 2.;
	}
    if (cos2 != params[7]) {
    	cos2   = params[7];
    	tan2p  = sqrt(1./(cos2*cos2) -1.) + 2.;
	}


    /* wichtige Sonden-Dimensionen */
    a_width  = params[10]; 
    a_height = params[11];
    a_mass   = params[12];

    p_height = params[17];

    /* daraus abgeleitete Fl"achen einer projezierten Sonde bzw. der Sonde */
    a_sonde_proj = a_width * a_height * cos1;
    a_sonde_perp = a_width * p_height;

    /* Berechnung der Ionensaettigungsstromdichte 
       ==========================================

       j_isat = n_e c_s e 
       j_isat = n_e * sqrt(t_e * c_e * 2.0 / (a_mass * c_mp)) * c_e; 
 
       neue Ergebnisse awc 01.09.93 
       c_i_Bohm = sqrt(  ( Z*T_e + 3 T_i ) / m_i )  wobei 3 = (n+2)/n, n=1
       c_e_Bohm = sqrt(  (3 T_e + T_i/Z ) / m_e )
       Sekundaerelektronenemissionskoeffizient g = 0.4
       j_isat = n_e * c_i_Bohm * e
       j_esat = 1/4 n_e * c_e_Bohm * e ( 1-g)
       j_quot = 0.15 * sqrt( (a_mass*c_mp) /c_me );  
    */

    j_isat = n_e * c_e * sqrt( 4.* c_e * t_e/ (a_mass*c_mp) );


    /* L_debye */
    ld = sqrt(t_e * c_eps0 / (n_e * c_e)) ;

   
    /* ---------------------------------------------------------------------
       wenn a_sonde_proj > 0.0 dann wurde eine geometrische Sondenfl"ache
       und ein Projektionswinkel angegeben. Dies hei"st, da"s f"ur eine
       flush mounted Sonde eine Kennlinie unter Ber"ucksichtigung der
       Schichteffekte berechnet werden soll.
       --------------------------------------------------------------------- */

    if ( a_sonde_proj > 0.0 )  
	{

        /* potentialunabhaengige Vorfaktoren vor Schichtdicke */
        sf1 = ld / 3.0 * alpha1 / sqrt(cos1);	
        sf2 = ld / 3.0 * alpha2 / sqrt(cos2);
        temp = sqrt( cos2 * inv_beta / cos1 );
        lfac = temp / a_height;
        bfac = temp / a_width;

        /* ab jetzt Berechnung fuer alle Spannungswerte */
        I_isat = a_sonde_proj * j_isat;

        for ( ind=0, ptr_vpr=vpr, ptr_strom=strom; ind < n_vpr; ind++ ) 
            {

            /* Spannung an Doppelsonde: phi_sonde */
            delta_phi = (dv - *ptr_vpr++) / t_e;

	    /* Nichtsaettigungskoeffizienten muessen immer groesser als -1 */
            if (delta_phi < 0) {
		ns_add1 = 1.0;
	        delta2  = sf2 * quick_075( delta_phi );
	        ns_add2 = ( (delta2+delta2)*bfac +1.) * 
				( 1. + delta2*tan2p*lfac );
		}
	    else {
	    	delta1  = sf1 * quick_075( -delta_phi );
	    	ns_add1 = ( (delta1+delta1)/a_width + 1.) * 
				( 1. + delta1*tan1p/a_height );
		ns_add2 = 1.0;
		}

            /* Berechnung des Strom-Vektors */
	    temp = exp( delta_phi );
	    /* temp = exp_lt( delta_phi ); */ 

            *ptr_strom++ =  I_isat * 
			(ns_add2 - ns_add1 * temp) / ( inv_beta + temp );

            }	/* for ind=0, n_vpr-1 */

        }	/* if (a_sonde_proj > 0.0) .....    */


    /* ---------------------------------------------------------------------
       wenn a_sonde_perp > 0.0 dann wurde bereits eine effektive Sondenfl"ache
       senkrecht zur Feldrichtung eingegeben. Es wird davon ausgegangen, da"s
       es sich hierbei um eine freistehende Sonde handelt, die von beiden
       Seiten Ionen und Elektronen aufsammeln kann. 

       Falls a_sonde_perp > 0.0 und a_sonde_proj > 0.0 dann wird eine 
       Kombination beider Effekte erwartet, d.h. der Gesamtstrom zur Sonde
       setzt sich additiv aus beiden Teilen zusammen.
       --------------------------------------------------------------------- */

    if (a_sonde_perp > 0.0)
	{

        /* Initialisierung, falls n"otig */
        if (a_sonde_proj <= 0.) {
            for ( ind=0, ptr_strom=strom; ind<n_vpr; ind++, *ptr_strom++ =0.0); 
	    }

        /* Kennlinie einer Doppelsonde mit Potentialverschiebung und
	   Fl"achenverh"altnis beta */
        I_isat = 2.* a_sonde_perp * j_isat;
        bfac = alpha2 * ld / a_width; 


        for ( ind=0, ptr_vpr=vpr, ptr_strom=strom; ind < n_vpr; ind++ ) 
	    {

            /* Spannung an Doppelsonde: phi_sonde */
            delta_phi = (dv - *ptr_vpr++) / t_e;
	    temp = exp( delta_phi );

            /* Nichts"attigungsverhalten im Elektronenast? */
	    /* ns_add2 = bfac * quick_075(delta_phi) ; */
            /* Berechnung des Strom-Vektors */
            *ptr_strom++ +=  I_isat * 
		(1. + bfac * quick_075(delta_phi) - temp) / ( inv_beta + temp );

	    }


	}	/* if (a_sonde_perp > 0.0) .....    */


    return 0 ;

    }		 /* fast_doppel */




/* 	-------------------------------------------------------- 

 			f_deriv 

	partielle Ableitung der Funktion an einem Punkt nach 
 	den Parametern 

	In y wird der Funktionswert zurueckgegeben, da er als
	Nebenprodukt der Ableitungs-Berechnung sowieso anfaellt.

 	-------------------------------------------------------- 

	M. Weinlich, 30.08.93

 	-------------------------------------------------------- */


int f_deriv( func, x, yfit, n_x, pars, do_var, n_pars, deriv)

    int		(*func)();
    double 	*x;
    double 	*yfit;
    int 	n_x;
    double 	*pars;
    int		*do_var;
    int 	n_pars;
    double 	**deriv;

    {

    /* Local variables */
    static int 		i, j, rc;
    static double 	save, d_par, temp;
    static double	*y2=NULL;
    static double	min=1.e30;
    static int		n_y2=-1;


    /* Function Body */

    /* Gebe Speicher frei, wenn n_x < 0 */
    if (n_x<=0) {
	printf("f_deriv: Gebe Speicherbereich frei: y2 = %X, Laenge = %d\n",
			y2, n_x+n_x);
	free_dvector(y2,0,-n_x-n_x-1);
        y2 = NULL;
	return 0;
	}

    /* allociere Speicherplatz nur 1x */
    if (y2==NULL) {
        n_y2 = n_x+n_x;
        y2 = dvector(0,n_y2-1);
        printf("f_deriv: Lege Speicherbereich an: y2 = %X, Laenge = %d\n",
					y2, n_x+n_x);
        }

    /* Gebe Speicher frei, wenn Aenderung in n_x  */
    if ( n_x > n_y2 ) {
	printf("f_deriv: Aendere Speicherbereich: y2 = %X, Laenge = %d\n",
			y2, n_x+n_x);
	free_dvector(y2,0,n_y2-1);
        y2 = dvector(0,n_x+n_x-1);
        n_y2 = n_x+n_x;
	}


    /* berechne Funktionswert */
    if ((*func)(x, yfit, n_x, pars) < 0) return (-110);


    /* versuche fuer jeden Paramter partielle Ableitung zu bilden */
    for (j = 0; j < n_pars; ++j) 
	{

	/* einige Parameter sollen festgehalten werden */
	if ( do_var[j] )  {

	    d_par  = (pars[j] == 0. ) ? 1.e-9 : pars[j] *1.e-6;
	    pars[j]  += d_par; 
	    if ( (*func)(x, y2, n_x, pars) < 0) return -111;
	    pars[j]  -= d_par; 

	    /* partielle Ableitung an jedem Punkt */
	    d_par = 1./d_par;
	    for( i=0; i< n_x; i++) deriv[j][i] = (*y2++ - *yfit++) * d_par;	
	    y2   = y2 - n_x;
	    yfit = yfit - n_x;
	
	    }		/* if do_var[j] */

    	}		/* for j=0,n_pars-1 */

    return 0;

    } /* f_deriv */





/*	-------------------------------------------------	*/
/*	-------------------------------------------------	*/
/*								*/
/*								*/
/*		singular value decomposition			*/
/*		zur Loesung eines linearen 			*/
/*		Gleichungssystems				*/
/*								*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*	-------------------------------------------------	*/








/*	-------------------------------------------------	*/
/*								*/
/*			dpythag					*/
/*			======= 				*/
/*								*/
/*	berechnet Betrag (a^2 + b^2)^0.5 unter spezieller 	*/
/*	Beruecksichtigung eventueller Bereichsverletzungen	*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 07.01.93					*/
/*	entnommen aus Numerical Recipies			*/
/*								*/
/*	-------------------------------------------------	*/


double dpythag(a, b)
double a;
double b;

{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if ((absa < 1.e30) && (absb < 1.e30)) return sqrt(absa*absa+absb*absb);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}		/* Ende von dphythag */




/*	-------------------------------------------------	*/
/*								*/
/*			dsvdcmp					*/
/*			======= 				*/
/*								*/
/*	berechnet singular value decomposition einer vor- 	*/
/*	gegebenen Matrix					*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 07.01.93					*/
/*	entnommen aus Numerical Recipies			*/
/*								*/
/*	-------------------------------------------------	*/

int dsvdcmp(a, m, n, w, v)
double 	**a;
int	m;
int	n;
double	w[];
double	**v;

{
	int flag,i,its,j,jj,k,l,nm;
	/* double anorm,c,f,g,h,s,scale,x,y,z,*rv1; */
	double anorm,c,f,g,h,s,scale,x,y,z;

	static double	*rv1=NULL;
	static int	n_rv1=0;

	/* rv1=dvector(1,n); */
	if ( rv1 == NULL ){
	    /* printf("dsvdcmp: Belege Speicherplatz: n = %d\n",n); */
	    rv1=dvector(1,n);
	    n_rv1 = n;
	    }
        if ( n > n_rv1 ) {
	    printf("dsvdcmp: Aendere Speicherplatz: n = %d --> %d\n",
			n_rv1, n);
	    free_dvector(rv1,1,n_rv1);
	    rv1=dvector(1,n);
	    n_rv1 = n;
	    }

	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}

			/* substantieller Fehler, SVD konvergiert nicht */
			if (its == 30) {
			   fprintf(stderr,"\n\n\ndsvcmp: ");
			   fprintf(stderr,"no convergence in 30 iterations");
			   fprintf(stderr," ... sorry!\n\n");
			   fprintf(stderr,"\t--> returning without acttion");
			   fprintf(stderr,"\n\n\n");
			   return (-42);
			   }

			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

	/* free_dvector(rv1,1,n); */

	return (0);

}	/* Ende von dsvdcmp	*/




/*	-------------------------------------------------	*/
/*								*/
/*			dsvbksb					*/
/*			======= 				*/
/*								*/
/*	Rueckwaerts-Substitution fuer Loesung eines lin.	*/
/*	Gleichungssystems mittels singular value decomp.	*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 07.01.93					*/
/*	entnommen aus Numerical Recipies			*/
/*								*/
/*	-------------------------------------------------	*/


void dsvbksb(u, w, v, m, n, b,x)
double 	**u;
double 	w[];
double	**v;
int	m,n;
double	b[];
double	x[];

	{
	int 	jj,j,i;
	/* double 	s,*tmp; */
	double 	s;

	static double 	*tmp=NULL;
	static int	n_tmp;

	if (tmp == NULL ) {
	    /* printf("dsvbksb: Belege Speicherplatz: n = %d \n", n); */
	    tmp=dvector(1,n);
	    n_tmp = n;
	    }
	if ( n > n_tmp ) {
	    printf("dsvbksb: Aendere Speicherplatz: n = %d --> %d\n",
			n_tmp, n);
	    free_dvector(tmp,1,n_tmp);
	    tmp=dvector(1,n);
	    n_tmp = n;
	    }


	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
			}
		tmp[j]=s;
		}

	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
		}

	/* free_dvector(tmp,1,n); */

	}		/* Ende von dsvbksb	*/




/*	-------------------------------------------------	*/
/*								*/
/*			dsvdinv					*/
/*			======= 				*/
/*								*/
/*	inverse Matrix mittels singular value decomposition	*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 07.01.93					*/
/*	entnommen aus Numerical Recipies			*/
/*								*/
/*	-------------------------------------------------	*/


void dsvdinv(u, w, v, n, inv )
double 	**u;
double 	w[];
double	**v;
int	n;
double 	**inv;

{
	int jj,j,i;
	double temp;


	for (i=1; i<=n; i++) {
	    for ( j=1; j<=n; j++ ) {

		temp = 0.0;
		for (jj=1; jj<=n; jj++) {
		    if (w[jj]>0.) temp += v[i][jj]*u[j][jj]/w[jj];
		    }
		inv[i][j] = temp;

		}	/* for i=1, n */
	    }		/* for j=1, n */


}		/* Ende von dsvdinv	*/




/*	-------------------------------------------------	
								
			dsvd_solver				
			=========== 				
								
	Loest lineares Gleichungssystem Ax=b ueber singular	
	value decomposition in double precision.		
								
	Eingabe:	Matrix A				
			Vektor b				
								
								
	Rueckgabe:	A^(-1) in Matrix A 			
			x in Vektor b				
								
	Funktionsreturn:	0 	alles o.k.		
				<0	Fehler in svd		
								
	--> um Zeit zu sparen, werden die internen Felder nur 
	    einmal angelegt und bei folgenden Aufrufen wieder-
	    verwendet. Sollte sich die Dimension n aendern,
	    muss dsvd_solver zuvor mit n<0 aufgerufen werden,
	    um die internen Felder zu loeschen.						
		
	-------------------------------------------------	
								
	M. Weinlich, 07.01.93					
	entnommen aus Numerical Recipies			
								
	-------------------------------------------------	*/


int dsvd_solver(a, b, n)
double 	**a;
double 	*b;
int	n;

{
	int 	rc;
	int	i, j;

        double  det, a00, a01, a02, a10, a11, a12, a20, a21, a22, b0, b1, b2;

	double	wmin, wmax;

	static int	last_n=0;

	/* Hilfsmatrizen fuer svd */
	static double	**u=NULL;
	static double	**v=NULL;
	static double	**ha=NULL;
	static double	*w=NULL;
	static double	*x=NULL;
	static double	*hb=NULL;

	/* printf("Beginn dsvd_solver\n"); */
/* printf("s"); */

	/* allociere nur 1x den benoetigten Speicherplatz */
 	if ( w == NULL ) {
	    printf("dsvd_solver: Lege Speicherplatz an, n = %d\n", n);
	    u  = dmatrix( 1, n, 1, n);
	    v  = dmatrix( 1, n, 1, n);
	    ha = dmatrix( 1, n, 1, n);
	    w  = dvector( 1, n);
	    x  = dvector( 1, n);
	    hb = dvector( 1, n);
	    last_n = n;
	    }

	/* Eingabe von n < 0 : Loesche interne Variablen */
	if ( n <= 0 ) {
	    printf("dsvd_solver: Loesche Speicherplatz, n = %d\n", n);
	    free_dmatrix( u, 1, -n, 1, -n);
	    free_dmatrix( v, 1, -n, 1, -n);
	    free_dmatrix( ha, 1, -n, 1, -n);
	    free_dvector( w, 1, -n);
	    free_dvector( x, 1, -n);
	    free_dvector( hb, 1, -n);
	    u  = NULL;
	    v  = NULL;
	    ha = NULL;
	    w  = NULL;
	    x  = NULL;
	    hb = NULL;
	    last_n = 0;
	    return 0;
	    }

	/* Aenderung in geforderter Dimension : Loesche interne Variablen */
	if ( n > last_n ) {
	    printf("dsvd_solver: Aendere Speicherplatz, n = %d --> n = %d\n", 
			last_n, n);
	    free_dmatrix( u, 1, last_n, 1, last_n);
	    free_dmatrix( v, 1, last_n, 1, last_n);
	    free_dmatrix( ha, 1, last_n, 1, last_n);
	    free_dvector( w, 1, last_n);
	    free_dvector( x, 1, last_n);
	    free_dvector( hb, 1, last_n);
	    u  = dmatrix( 1, n, 1, n);
	    v  = dmatrix( 1, n, 1, n);
	    ha = dmatrix( 1, n, 1, n);
	    w  = dvector( 1, n);
	    x  = dvector( 1, n);
	    hb = dvector( 1, n);
	    last_n = n;
	    }


	/* Sonderbehandlung: 2x2 und 3x3 Matrizen koennen noch
	                     von Hand invertiert werden! */
        if ( n == 2 ) {
	    /* printf("\n\n --> Sonderbehandlung n=2 in dsvd_solver \n\n"); */
	    a00 = a[0][0];
	    a01 = a[0][1];
	    a10 = a[1][0];
	    a11 = a[1][1];
	    det = a00*a11 - a01*a10;
	    if ( fabs(det) > 1.e-30 ) {
		/* inverse Matrix */
		det = 1. / det;
		a[0][0] = a11 * det;
		a[0][1] *= (-det);
		a[1][0] *= (-det);
		a[1][1] = a00 * det;
		/* Loesungsvektor */
		b0 = b[0];
		b1 = b[1];
		b[0] = a[0][0]*b0 + a[0][1]*b1;
		b[1] = a[1][0]*b0 + a[1][1]*b1;
		/* das war's */
		return 0;
		}
	    }
        else if ( n == 3 ) {
	    /* printf("\n\n --> Sonderbehandlung n=3 in dsvd_solver \n\n"); */
	    a00 = a[0][0];
	    a01 = a[0][1];
	    a02 = a[0][2];
	    a10 = a[1][0];
	    a11 = a[1][1];
	    a12 = a[1][2];
	    a20 = a[2][0];
	    a21 = a[2][1];
	    a22 = a[2][2];
	    det = a00*a11*a22 + a01*a12*a20 + a02*a10*a21 -
			a02*a11*a20 - a01*a10*a22 - a00*a12*a21;
	    if ( fabs(det) > 1.e-30 ) {
		/* inverse Matrix */
		det = 1. / det;
		a[0][0] = det * (a11*a22 - a12*a21);
		a[0][1] = det * (-a01*a22 + a02*a21);
		a[0][2] = det * (a01*a12 - a02*a11);
		a[1][0] = det * (-a10*a22 + a12*a20);
		a[1][1] = det * (a00*a22 - a02*a20);
		a[1][2] = det * (-a00*a12 + a02*a10);
		a[2][0] = det * (a10*a21 - a11*a20);
		a[2][1] = det * (-a00*a21 + a01*a20);
		a[2][2] = det * (a00*a11 - a01*a10);
		/* Loesungsvektor */
		b0 = b[0];
		b1 = b[1];
		b2 = b[2];
		b[0] = a[0][0]*b0 + a[0][1]*b1 + a[0][2]*b2;
		b[1] = a[1][0]*b0 + a[1][1]*b1 + a[1][2]*b2;
		b[2] = a[2][0]*b0 + a[2][1]*b1 + a[2][2]*b2;
		/* das war's */
		return 0;
		}
	    }


/*	if (n<4) {
	    printf( "\n\n");
	    printf( "dsvd_solver: invertiere singulaere Matrix \n");
	    printf( "             det = %g\n\n", det);
	    }
*/

	/*	Achtung:	externe Felder haben Indizes
				von 0..n-1, waehrend die neu
				definierten internen Felder 
				der svd-Routine Indizes von
				1..n haben!!!			*/

	/* kopiere a nach u */
	for (i=1; i<=n; i++) {
	    for (j=1; j<=n; j++) u[i][j] = a[i-1][j-1];
	    }

	/* kopiere b nach hb */
	for (i=1; i<=n; i++) hb[i] = b[i-1];

	/* singular value decomposition von A */
	rc = dsvdcmp( u, n, n, w, v);
	if ( rc < 0 ) return (rc);


	/* was haben wir denn fuer Ergebnisse erhalten */
	wmax = 0.0;
	for (i=1; i<=n; i++) wmax = DMAX( wmax, w[i]);
	wmin = wmax * 1.e-12;
	for (i=1; i<=n; i++) if ( w[i] <= wmin ) w[i] = wmin; 
	
	/* Loesung des Gleichungssystems */
	dsvbksb( u, w, v, n, n, hb, x ); 

	/* Berechnung der inversen Matrix */
	dsvdinv( u, w, v, n, ha );	


	/* Rueckschreiben der Ergebnisse - Loesungsvektor x --> b */
	for (i=1; i<=n; i++) b[i-1] = x[i]; 

	/* Rueckschreiben der Ergebnisse - Inverse Matrix --> A */
	/* kopiere a nach u */
	for (i=1; i<=n; i++) {
	    for (j=1; j<=n; j++) a[i-1][j-1] = ha[i][j];
	    }	


	return (0);

}		/* Ende von dsvd_solver	*/




/*	-------------------------------------------------	*/
/*	-------------------------------------------------	*/
/*								*/
/*								*/
/*		Levenberg-Marquard Algorithmus			*/
/*		fuer Nichtlineare Fits mit Neben-		*/
/*		bedingungen					*/
/*								*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*	-------------------------------------------------	*/



/*	-------------------------------------------------	*/
/*								*/
/*			covsrt					*/
/*			====== 					*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 10.01.93					*/
/*	entnommen aus Numerical Recipies			*/
/*								*/
/*	-------------------------------------------------	*/

void covsrt(covar, ma, ia, mfit)

double	**covar;
int	ia[], ma;
int	mfit;

{
	int i,j,k;
	float swap;

	for (i=mfit;i<ma;i++)
		for (j=0;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit-1;

	for (j=ma-1;j>=0;j--) {
		if (ia[j]) {
			for (i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
			}
		}

}		/* covsrt */






/*	-------------------------------------------------	*/
/*								*/
/*			calc_chisqr				*/
/*			===========				*/
/*								*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 10.01.93					*/
/*								*/
/*	-------------------------------------------------	*/

void calc_chisqr(  y, yfit, sig, ndata, pars, do_var, n_pars, 
			nb_values, nb_const, chisq)

double	*y, *yfit, *sig;
int	ndata;
double	*pars;
int	*do_var;
int	n_pars;
double	*nb_values, *nb_const;
double	*chisq;

{
	int 	i,j;

	/* berechen chi^2, mit Sonderbehandlung fuer kleine n's */
	if ( ndata == 2 ) {
	   *chisq = DSQR( (y[0]-yfit[0])/sig[0] ) + 
			DSQR( (y[1]-yfit[1])/sig[1] );
	   }
	else if (ndata == 3) {
	   *chisq = DSQR( (y[0]-yfit[0])/sig[0] ) + 
			DSQR( (y[1]-yfit[1])/sig[1] ) +
			DSQR( (y[2]-yfit[2])/sig[2] );
	   }
 	else {
	   *chisq=0.0;
	    for (i=0;i<ndata;i++) 
		*chisq += DSQR( (y[i]-yfit[i])/sig[i] );
	   }

	/* Nebenbedingungen */
	for (j=0; j<n_pars; j++) {
	    if ( (do_var[j]) && (nb_const[j]>0.) )
		*chisq += nb_const[j]*DSQR( pars[j]-nb_values[j] );
	    }



}		/* calc_chisqr */







/*	-------------------------------------------------	*/
/*								*/
/*			mrqmin					*/
/*			====== 					*/
/*								*/
/*	Hauptroutine fuer Levenberg-Marquardt			*/
/*								*/
/*	Funktionsreturn:	0	alles o.k.		*/
/*				<0	Fehler, wahrsl. in svd	*/
/*								*/
/*								*/
/*	-------------------------------------------------	*/
/*								*/
/*	M. Weinlich, 10.01.93					*/
/*	entnommen aus Numerical Recipies			*/
/*								*/
/*	-------------------------------------------------	*/


int mrqmin( func, x, y, sig, yfit, ndata, a, ia, ma, mfit, nb_values, nb_const, 
		covar, alpha, chisq, alamda)

int	(*func)();
double	x[], y[], sig[];
double  *yfit;
int	ndata;
double	a[];
int 	ia[], ma, mfit;
double	*nb_values, *nb_const;
double	**covar, **alpha;
double	*chisq;
double	*alamda;

	{

	int 		i, j, k, l, m;
	int		rc;
	double		ochisq, temp, dy, wt;

	static double 	*atry=NULL,
			*beta=NULL,
			*da=NULL,
			*oneda=NULL;
	static double 	**dyda=NULL;
	static int	last_ma, last_ndata;



        /* printf("Beginn mrqmin\n"); */

	/* stelle benoetigte Datenfelder bereit */
	if ( beta == NULL ) {
                /* printf( "mrqmin: Lege Speicher an, ma = %d\n",ma); */
		atry  = dvector(0,ma-1);
		beta  = dvector(0,ma-1);
		da    = dvector(0,ma-1);
	        dyda  = dmatrix(0,ma-1,0,2*ndata-1);
		oneda = dvector(0,ma-1);
		last_ma = ma;
		last_ndata = ndata+ndata;
		}

	if (( ma > last_ma ) || ( ndata > last_ndata )) {
                printf( "mrqmin: Aendere Speicher, ma = %d --> %d\n",
					last_ma,ma);
                printf( "                          nd = %d --> %d\n",
					last_ndata,ndata);
		free_dvector(oneda,0,last_ma-1);
		free_dvector(da,0,last_ma-1);
		free_dvector(beta,0,last_ma-1);
		free_dvector(atry,0,last_ma-1);
	        free_dmatrix(dyda,0,last_ma-1,0,last_ndata-1);
		atry  = dvector(0,ma-1);
		beta  = dvector(0,ma-1);
		da    = dvector(0,ma-1);
	        dyda  = dmatrix(0,ma-1,0,2*ndata-1);
		oneda = dvector(0,ma-1);
		last_ma = ma;
		last_ndata = ndata+ndata;
		}

	/* Initalisierung, stelle benoetigte Datenfelder bereit */
	if (*alamda < 0.0) {
                /* printf( "*alamda = %g --> Initialisierung\n",*alamda); */
		*alamda=1.e-4;
		}


	/* Ende der Iterationen --> gebe Speicher frei */
	if (*alamda == 0.0) {
                /* printf( "*alamda = %g --> Abgesang in mrqmin\n",*alamda); */
		covsrt(covar,ma,ia,mfit);
		*alamda = -1.;
		return (0);
		}


	/* berechne aktuellen Fitwert fuer y 
	   sowie alle partiellen Ableitungen */
	rc = f_deriv(func, x, yfit, ndata, a, ia, ma, dyda);
	if ( rc < 0 ) return (rc);

	/* aktuelles chi_sqr */
	calc_chisqr(y,yfit,sig,ndata,a,ia,ma,nb_values,nb_const,&ochisq);


	/* Initialisierung der alpha und beta Matrizen */
	for (j=0;j<mfit;j++) {
	    for (k=0;k<=j;k++) alpha[j][k]=0.0;
	    beta[j]=0.0;
	    }


	/* Berechne alpha, beta */
	for (i=0;i<ndata;i++) {
		dy=y[i]-yfit[i];
		for (j=0,l=0;l<ma;l++) {
			if (ia[l]) {
				wt=dyda[l][i]/(sig[i]*sig[i]);
				for (k=0,m=0;m<=l;m++)
				    if (ia[m]) alpha[j][k++] += wt*dyda[m][i];
				beta[j++] += dy*wt;
				}
			}
		}

	/* fuelle symmetrischen Teil von alpha auf */
	for (j=1;j<mfit;j++)
		for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];

	/* Nebenbedingungen */
	for (j=0,k=0 ; j<ma; j++) {
	    if  (ia[j]) {
		if (nb_const[j]>0.) {
		    beta[k] += nb_const[j]*(nb_values[j]-a[j]);
		    alpha[k][k] += nb_const[j];
		    }
		k++;
		}
	    }


	/* Schleife, um zu akt. alpha und beta Werten Min in chisq zu finden */
	do {

	    for (j=0; j<mfit; j++) {
	        for (k=0; k<mfit; k++) covar[j][k] = alpha[j][k];
	        covar[j][j]=alpha[j][j]*(1.0+(*alamda));
	        oneda[j] = beta[j];
	        }		/* for j=0, mfit-1 */

	    /* Loesung der Gleichung covar * da = oneda */
	    rc = dsvd_solver(covar,oneda,mfit);

	    /* Fehler in svd aufgetreten? wenn ja, dann hat weitere Berechnung
		keinen Sinn mehr --> geordneter Rueckzug aus mrqmin, immer
		moeglich wenn alamda = 0 */
	    if ( rc < 0 ) return rc;

	    /* probiere neue Parameter aus */
	    for (j=0,l=0;l<ma;l++) {
	        atry[l]=a[l];
	        if (ia[l]) atry[l]+=oneda[j++];
	        }
            rc  = (*func)(x, yfit, ndata, atry);
	    if ( rc < 0 ) return(rc);
	    calc_chisqr(y,yfit,sig,ndata,atry,ia,ma,nb_values,nb_const,chisq);

	    /* falls Fit keine Verbesserung bringt */
	    *alamda *= 10.0;
	    if (*alamda > 1.e25) return (-4711);

	    } while ( (*chisq) > ochisq  );


	/* bis jetzt ist alles gut gegangen, neuer Parametersatz 
	   passt besser zu Daten. --> Vorbereitung fuer naechste Iteration */
	if (*alamda > 1.e-20) *alamda *= 0.01;
	for (l=0; l<ma; l++) a[l]=atry[l];

        /* printf("Ende mrqmin\n"); */

	/* alles o.k. */
	return 0;

	}		/* mrqmin */








/* 	-------------------------------------------------------- 

			magfit 
			======

 	Verwaltungsroutine fuer nichtlinearen Fit nach Levenberg-
	Marquardt an Sondenkennlinie entsprechend der Fitfunktion
	mag_doppel

  	-------------------------------------------------------- 

	M. Weinlich, 11.01.93

 	-------------------------------------------------------- */


int magfit( func, x, y, sig, yfit, n_x, pars, do_var, 
		n_pars, nb_values, nb_const,
		chi_sqr, iter_max, eps_abs, eps_rel)

    int		(*func)();
    double 	*x;
    double 	*y;
    double 	*sig;
    double 	*yfit;
    int 	*n_x;
    double 	*pars;
    int 	*do_var;
    int 	*n_pars;
    double	*nb_values, *nb_const;
    double 	*chi_sqr ;
    int 	*iter_max;
    double 	*eps_abs;
    double 	*eps_rel;

    {

    /* Local variables */
    int 	do_continue, iteration, best_estimate, 
		bad_estimate, too_much, do_nb, nb_angepasst;
    int 	i, j, k, ind1, ind2, rc;
    int 	err;
    double 	fac1, fac2;

    int 	iter, mfit;
    double 	chi_1;
    double 	alamda, last_lamda, chi_old, chi_oldold;
    double	*dbl_ptr;
    double	temp;

    static double	**alpha=NULL;
    static double	*curr_nb=NULL;
    static int		*save_do_var=NULL;
    static int		last_npars;

    double	chi_est, mean, min;
    int		pos;

    static double	**sigma=NULL;
    static int		old_n_sigma=-1;


    /* Initialisierung des Hilfsfelds */
    if ( sigma == NULL ) {
        printf("Matrix sigma angelegt mit %dx%d Elementen\n",*n_pars,*n_pars);
        sigma = dmatrix(0, *n_pars-1, 0, *n_pars-1);
	old_n_sigma = *n_pars;
	}
    else if ( *n_pars > old_n_sigma ) {
        printf( "sigam geloescht mit %d Elementen\n", old_n_sigma );
        free_dmatrix( sigma, 0, old_n_sigma-1, 0, old_n_sigma-1);
        sigma = dmatrix(0, *n_pars-1, 0, *n_pars-1);
	old_n_sigma = *n_pars;
        }


    /* Function Body */
/*
    printf( "Startwerte in magfit: \n\n");
    printf( "func = %d \n ", func );
    printf( "n_x = %d \n", *n_x );
    printf( "x = ");
    for (i=0; i<*n_x; i++) printf( "%g  ", x[i] );
    printf("\n");
    printf( "y = ");
    for (i=0; i<*n_x; i++) printf( "%g  ", y[i] );
    printf("\n");
    printf( "sig = ");
    for (i=0; i<*n_x; i++) printf( "%g  ", sig[i] );
    printf("\n");
    printf( "yfit = ");
    for (i=0; i<*n_x; i++) printf( "%g  ", yfit[i] );
    printf("\n");
    printf( "n_pars = %d \n", *n_pars );
    printf( "pars = ");
    for (i=0; i<*n_pars; i++) printf( "%g  ", pars[i] );
    printf("\n");
    printf( "do_var = ");
    for (i=0; i<*n_pars; i++) printf( "%d  ", do_var[i] );
    printf("\n");
    printf( "nb_values = ");
    for (i=0; i<*n_pars; i++) printf( "%g  ", nb_values[i] );
    printf("\n");
    printf( "nb_const = ");
    for (i=0; i<*n_pars; i++) printf( "%g  ", nb_const[i] );
    printf("\n");
    printf("\n");
    printf( " chi_sqr = %g \n", *chi_sqr );
    printf( " iter_max = %d\n", *iter_max );
    printf( " eps_abs = %g \n", *eps_abs );
    printf( " eps_rel = %g \n", *eps_rel);
    printf("\n\n");
*/


    /* Stelle benoetigte Hilfsfelder bereit, wenn noch kein Speicher allok. */
    if ( alpha == NULL ) {
        printf("magfit: stelle Speicher bereit, n_pars = %d\n", *n_pars);
        alpha       = dmatrix( 0, *n_pars-1, 0, *n_pars-1);
        curr_nb     = dvector(0, *n_pars-1);
        save_do_var = ivector(0, *n_pars-1);
	last_npars  = *n_pars;
        }

    /* Passe Speicher bei Aenderung der Anforderungen an */
    if ( *n_pars > last_npars ) {
        printf("magfit: aendere Speicher, n_pars = %d --> \n", 
			last_npars, *n_pars);
        free_dmatrix( alpha, 0, last_npars-1, 0, last_npars-1);
        free_dvector( curr_nb, 0, last_npars-1);
        free_ivector( save_do_var, 0, last_npars-1);
        alpha       = dmatrix( 0, *n_pars-1, 0, *n_pars-1);
        curr_nb     = dvector(0, *n_pars-1);
        save_do_var = ivector(0, *n_pars-1);
	last_npars  = *n_pars;
        }

    /* Initialisierung einiger Flags */
    iteration     = 1;   	/* 1 = TRUE */
    best_estimate = 0;   	/* 0 = FALSE */
    bad_estimate  = 0;   	/* 0 = FALSE */
    nb_angepasst  = 0;   	/* 0 = FALSE */
    too_much   	  = 0;   	/* 0 = FALSE */
    iter          = 0;
    chi_old       = 1.e20;
    chi_oldold    = 1.e22;


    /* Initialisierung der look-up-Tabelle */
    init_quick_075();
    /* init_exp_lt(); */

    /* in allen positiven Parametern logarithmischer Fit */
    if ( pars[0] <= 0.0 ) printf( "Fehler in pars[0] %g \n", pars[0]);
    if ( pars[1] <= 0.0 ) printf( "Fehler in pars[1] %g \n", pars[1]);
    if ( pars[2] <= 0.0 ) printf( "Fehler in pars[2] %g \n", pars[2]);
    if ( pars[3] <= 0.0 ) printf( "Fehler in pars[3] %g \n", pars[3]);
    if ( pars[5] <= 0.0 ) printf( "Fehler in pars[5] %g \n", pars[4]);
    pars[0] = log(pars[0]*sqrt(pars[1]));
    pars[1] = log(pars[1]);
    pars[2] = log(pars[2]);
    pars[3] = log(pars[3]);
    pars[5] = log(pars[5]);
    if (nb_const[0] > 0.0) nb_values[0] = log(nb_values[0]);
    if (nb_const[1] > 0.0) nb_values[1] = log(nb_values[1]);
    if (nb_const[2] > 0.0) nb_values[2] = log(nb_values[2]);
    if (nb_const[3] > 0.0) nb_values[3] = log(nb_values[3]);
    if (nb_const[5] > 0.0) nb_values[5] = log(nb_values[5]);

    /* Fit in cos(psi) statt in psi */
    pars[6] = cos(pars[6]);
    pars[7] = cos(pars[7]);
    if (nb_const[6] > 0.0) nb_values[6] = cos(nb_values[6]);
    if (nb_const[7] > 0.0) nb_values[7] = cos(nb_values[7]);


    /* Fit mit oder ohne Nebenbedingungen */
    for (j=0, do_nb=0;j<*n_pars; j++ ) if (nb_const[j] > 0.0 ) do_nb++;

    /* wenn Nebenbedingungen (do_nb>0):
	Nebenbedingungen wirken auf chisq: chisq = chisq + nb_const*chi_est
	d.h. nb_const gibt den Bruchteil an, den eine bestimmte Abweichung 
	vom eigentliche chisq ausmachen soll. Da das eigentliche chisq aber
	erst nach dem Fit bekannt ist, muss eine Schaetzung (chi_est) vorge-
	nommen werden. Dazu wird davon ausgegangen, das chi_sq im wesentlichen
	von Fluktuationen bestimmt wird, die etwa 30% des Ionensaettigungstroms
	ausmachen koennen. Diese Fluktuationsamplitude wird fuer alle Spannungen
	der Sonde als konstant angesehen. */
    chi_est = 0.0;
    if (do_nb) {
	/* Bestimme min(strom) */
	min = 0.0;
	for (i=0;i<*n_x;i++) if (y[i]<min) min = y[i];
	/* Bis wohin geht in etwa Ionenstrom? */
	pos = *n_x-1;
	while ( y[pos] > 0.1*min ) pos--;
 	/* Berechne mittleren Ionensaettigungsstrom */
	mean = 0.;
	for (i=0;i<=pos;i++) mean += y[i];
	mean /= (double)(pos+1);
	/* Schaetzwert fuer chisq */
	temp = 0.1 * mean;
	chi_est = (*n_x) * temp*temp;
	/* Info, wo liegt der Schaetzwert? */
	printf("Fit mit Nebenbedingungen, Schaetzwert fuer chi^2 = %g\n",
		chi_est);
	}

    /* aktuelle Gewichtsfaktoren, gelten fuer Fit mit und ohne Nebenbed. */
    for (j=0;j<*n_pars;j++) curr_nb[j] = nb_const[j]*chi_est;

    /* aktuelle Anzahl der Fitparameter */
    for (j=0,mfit=0;j<*n_pars;j++) if (do_var[j]) mfit++;


    /* Initialisierung in mrqmin, erzeuge Vektoren, Matrizen, etc. */
    alamda = -1.0;


    /* Iteration, solange bis keine wesentliche Verbesserung in chi^2 */
    while (!too_much)
	{

	/* berechne naechsten Iterationsschritt */
        rc = mrqmin(func, x,y,sig,yfit,*n_x,pars,do_var,*n_pars, mfit, 
			nb_values, curr_nb, sigma,alpha,chi_sqr,&alamda);

	/* rc = -4711 ist Returncode fuer zu grosses alamda, d.h.
	   Schaetzwert laesst sich aus dem einen oder anderen Grund
	   nicht mehr verbessern */
	bad_estimate = ( rc == (-4711) );
	if (bad_estimate) break;

	/* sonstige, noch nicht abgefangenen Fehler */
	if (rc < 0 ) break;

	/* Fit wirklich gut ? */
	best_estimate = ( ( (chi_old) <= chi_oldold ) && 
			    ( ( 1.-(chi_old)/chi_oldold < (*eps_rel) )	|| 
				  ( chi_oldold-(chi_old) < (*eps_abs) ) ));
	best_estimate = best_estimate &&
			( ( (*chi_sqr) <= chi_old ) && 
			    ( ( 1.-(*chi_sqr)/chi_old < (*eps_rel) )	|| 
				  ( chi_old-(*chi_sqr) < (*eps_abs) ) ));
	if (best_estimate) break;

	/* Guete des Fits : der Fit darf nicht unendlich lange dauern */
	too_much = ( iter >= *iter_max );

	/* Vorbereitung fuer naechsten Iterationsschritt */
	iter++;
	chi_oldold = chi_old;
	chi_old = *chi_sqr;

	}	/* while (!too_much) */



    /* nicht nur im Falle einer normalen Beendigung der Schleife muss noch 
       aufgeraeumt werden */
    alamda = 0.0;
    temp = mrqmin(func, x,y,sig,yfit,*n_x, pars,do_var,*n_pars, mfit,
		nb_values, curr_nb, sigma, alpha, chi_sqr, &alamda);


    /* 	Abgesang */
    
    /* "ich lebe noch Zeil" muss nur abgeschlossen werden, wenn sie auch */
    /* wirklich angelegt wurde */
    /* if ( iter >= 10 ) fprintf ( stderr, "\n" ); */

    /* Fit ok, welches Chi? nur interessant fuer Nebenbedingungen, da in
       diesem Fall zu Beginn ein Chi abgeschaetzt werden musste */
    if ( (best_estimate) && (do_nb) ){
	fprintf ( stderr, "Fit nach Levenberg-Marquardt erfolgreich beendet " );
	fprintf ( stderr, "( chi^2 = %g ) \n", *chi_sqr );
	}

    /* sind Fehler aufgetreten, wenn ja, welche? */
    if (rc == -21 ) 
	fprintf(stderr,"Zuviele Parameter ausserhalb ihres Wertebereichs\n");
    else if (rc == -42) 
	fprintf(stderr, "Iteration abgebrochen --- Fehler in svd\n" );
    else if (bad_estimate) 
	fprintf(stderr, "Iteration abgebrochen --- Signal 'bad estimate'\n" );
    else if (too_much) 
	fprintf(stderr, "Iteration nicht konvergiert\n" );
    else if (rc != 0 ) 
	fprintf(stderr, "Fehler waehrend Iteration, Fehlercode = %d\n",rc);

    if (bad_estimate) rc = -33;
    if (too_much) rc = -34;

    /* logarithmischer Fit in allen positiven Parametern */
    pars[0] = exp(pars[0]-0.5*pars[1]);
    pars[1] = exp(pars[1]);
    pars[2] = exp(pars[2]);
    pars[3] = exp(pars[3]);
    pars[5] = exp(pars[5]);

    pars[6] = acos(pars[6]);
    pars[7] = acos(pars[7]);

    /* 	als Info ueber Fit gebe erreichte Werte zurueck */
    *iter_max = iter ;
    *eps_rel  = 1.- (*chi_sqr) / chi_old ;
    *eps_abs  = chi_old - (*chi_sqr);


/*    if ( bad_estimate ) {
        printf("magfit: ne = %g, te = %g, d_V = %g, beta = %g\n",
			pars[0], pars[1], pars[4], pars[3] );
        printf("        alpha1 = %g, alpha2 = %g, psi1 = %g, psi2 = %g\n",
			pars[2], pars[5], pars[6], pars[7]);
        printf("x = ");
	for ( i=0; i<*n_x; i++) printf(" %g ", x[i] );
	printf("\n");
        printf("y = ");
	for ( i=0; i<*n_x; i++) printf(" %g ", y[i] );
	printf("\n");
        printf("yfit = ");
	for ( i=0; i<*n_x; i++) printf(" %g ", yfit[i] );
	printf("\n");
        }
*/

    return (rc);

    } /* magfit */











/* 	-------------------------------------------------------- 

		Fuer die Auswertung von triple oder tetra
		Sonden folgt nun die Loesung eines nicht-
		linearen Gleichungssystem, das die voll-
		staendige Kennlinienform beruecksichtigt

 	-------------------------------------------------------- */












/* 	-------------------------------------------------------- 

 			fast_doppel_nl 
			==============

	abgespeckte Version von fast_doppel aus mag_fit.c um ein
	2x2 oder 3x3 nichtlineares Gleichungssystem zu loesen.

 	berechnet einen Punkt aus Kennlinie einer beidseitig  
 	nichtsaettigenden Doppelsonde 

	Modifikationen: 	in params werden wie in fast_doppel alle
				relevanten Parameter uebergeben, sie werden
				jedoch nur einmal ausgelesen, genau dann wenn
				init = TRUE = 1

				die Prameter, nach denen das nl-System
				geloest werden soll, sind in "vars" abgelegt:
				n=2:	n_e*sqrt(T_e), T_e
				n=3:	n_e*sqrt(T_e), T_e, beta


	vpr		Vektor [n_vars], Spannungswerte
	strom		Vektor [n_vars], Stromwerte --> Rueckgabe
	vars		Vektor [n_vars], aktuelle Variablen zur Loesung des
					 nl-Systems

 	params(0)	Elektronendichte 	log(n_e) 
 	params(1)	Elektronentemperatur 	log(T_e) 
 	params(2)	Korrekturfaktor		log(alpha_1) 
 	params(3)	Flaechenverhaeltnis 	log(beta) 
 	params(4)	Unterschied in V_plas. 	delta_V 
 	params(5)	Korrrekturfaktor 	log(alpha_2) 
 	params(6)	Winkel an Sonde		cos(psi1) 
 	params(7)	Winkel an Gegenelektr.	cos(psi2) 

 	params(10)	Sondenbreite		a_width 
 	params(11)	Sondenhoehe		a_height 
 	params(12)	Massezahl des Fuellg.	a_mass 

	params(13)	Kernladungszahl	d. F.g.	a_z
	params(14)	Toroidalfeld		B_t


	params(16)	T_i/T_e			a_tau

	--> Erweiterung, um auch Sonden bei senkrechtem Einfall behandeln
	    zu k"onnen (z.B. Mitelebenenmanipulator, Rohde). Es muss lediglich
	    die projizierte L"ange bzw. H"ohe angegeben werden, da die 
	    Sondenbreite a_width = params(13) bereits als senkrecht zum
	    Magnetfeld vorausgesetzt wird.

	params(17)	Sondenhoehe		p_height

	--> Aenderung auf gamma = 3 (Adiabatenkoeffizient)
	    wg. Riemann

 	-------------------------------------------------------- 

	M. Weinlich, 12.06.95

 	-------------------------------------------------------- */

int fast_doppel_nl(vpr, strom, vars, n_vars, init, params )

    double 	*vpr;
    double 	*strom;
    double 	*vars;
    int		n_vars;
    int		init;
    double	*params;

    {

    /* Local variables */
    double 	n_e, t_e;
    double	delta_phi ;
    double	j_isat, I_isat, ns_add1, ns_add2;
    double	sf1, sf2, lfac, bfac;
    double	ld, delta1, delta2;
    double	temp;

    double	*ptr_vpr, *ptr_strom;
    int 	i, ind;
    int		n_changes=0;

    static double	inv_beta, alpha1, alpha2, dv;
    static double 	a_height, a_width, p_height, a_mass, a_tau;
    static double	a_sonde_proj, a_sonde_perp;
    static double	cos1, cos2, tan1p, tan2p;

/* printf("f"); */



    /* mehr als drei freie Parameter sind nicht erlaubt */
    if ( n_vars > 3 ) {
	fprintf( stderr, "Achtung: in fast_doppel_nl sind max. ");
	fprintf( stderr, "3 freie Prameter erlaubt! \n");
	fprintf( stderr, "         aktuelle Anzahl n_vars = %d \n", n_vars);
        fprintf( stderr, "\n    ==> Abbruch \n\n");
        exit(-1);
	}


    /* Auslesen von params in static-Variablen nur zur Inititalisierung */
    if (init) {

        /* wichtige Sonden-Dimensionen */
        a_width  = params[10]; 
        a_height = params[11];
        a_mass   = params[12];
        p_height = params[17];

	params[3] = DMAX( params[3], 0.);
	params[3] = DMIN( params[3], 54.);
        inv_beta = 1./params[3] ;

        params[2] = DMAX( params[2], 0.);
	params[2] = DMIN( params[2], 100.);
        alpha1 = params[2];

        params[5] = DMAX( params[5], 0.);
	params[5] = DMIN( params[5], 100.);
        alpha2 = params[5];

	params[6] = DMAX( params[6], 0.001);
	params[6] = DMIN( params[6], 0.999);
        cos1   = params[6];
        tan1p  = sqrt(1./(cos1*cos1) -1.) + 2.;

	params[7] = DMAX( params[7], 0.001);
	params[7] = DMIN( params[7], 0.999);
    	cos2   = params[7];
    	tan2p  = sqrt(1./(cos2*cos2) -1.) + 2.;

        dv     = params[4];

        /* daraus abgeleitete Flaeche einer (projezierten) Sonde  */
        a_sonde_proj = a_width * a_height * cos1;
        a_sonde_perp = a_width * p_height;

        }




    /* Ueberpruefung auf unsinnige Fit-Parameter:
 	0.002  <  t_e            <  400
        1.e7   <  n_e sqrt(t_e)  <  1.e28
        0.001  <  beta           <  60
        0.000  <  alpha          <  100
	0.000  <  cos(psi)       <  1.000
    */

    if ( vars[1] < -6. ) {
        /* printf( "Fehler t_e : %g", exp(vars[1]) ); */
	/* keine Aenderung in n_e als Folge */
	temp      = vars[0] + 0.5*vars[1];
        vars[1] = -6.;
	vars[0] = temp - 0.5*vars[1];
        /* printf( " --> %g \n", exp(vars[1])); */
	n_changes++;
        }
    else if ( vars[1] > 6. ) {
        /* printf( "Fehler t_e : %g", exp(vars[1]) ); */
	/* keine Aenderung in n_e als Folge */
	temp      = vars[0] + 0.5*vars[1];
        vars[1] = 6.;
	vars[0] = temp - 0.5*vars[1];
        /* printf( " --> %g \n", exp(vars[1])); */
	n_changes++;
        }

    if ( vars[0] < 15. ) {
        /* printf( "Fehler n_e : %g", exp(params[0]-0.5*params[1]) ); */
	params[0] = 15.;
        /* printf( " --> %g \n", exp(params[0]-0.5*params[1])); */
	n_changes++;
        }
    else if ( vars[0] > 65. ) {
        /* printf( "Fehler n_e : %g", exp(params[0]-0.5*params[1]) ); */
	params[0] = 65.;
        /* printf( " --> %g \n", exp(params[0]-0.5*params[1])); */
	n_changes++;
        }
    
    /* Zuordnung der Fit-Variablen */
    n_e    = exp(vars[0]-0.5*vars[1]);
    t_e    = exp(vars[1]);

    if ( n_vars == 3 ) {

        if ( vars[2] < -7. ) {
	    /* printf( "Fehler in beta: %g", exp(vars[2])); */
	    vars[2] = -7.;
            /* printf( " --> %g \n", exp(vars[2])); */
	    n_changes++;
	    }
        else if ( vars[2] > 4. ) {
	    /* printf( "Fehler in beta: %g", exp(vars[2])); */
	    vars[2] = 4.;
            /* printf( " --> %g \n", exp(vars[2])); */
	    n_changes++;
	    }

        inv_beta = exp(-vars[2]) ;

 	}



    /* ==> ab hier folgt der identische Teil zu fast_doppel, lediglich
           die Initialisierung bzw. Belegung der Variablen unterscheidet
	   sich geringfuegig */


    /* Berechnung der Ionensaettigungsstromdichte 
       ==========================================

       j_isat = n_e c_s e 
       j_isat = n_e * sqrt(t_e * c_e * 2.0 / (a_mass * c_mp)) * c_e; 
 
       neue Ergebnisse awc 01.09.93 
       c_i_Bohm = sqrt(  ( Z*T_e + 3 T_i ) / m_i )  wobei 3 = (n+2)/n, n=1
       c_e_Bohm = sqrt(  (3 T_e + T_i/Z ) / m_e )
       Sekundaerelektronenemissionskoeffizient g = 0.4
       j_isat = n_e * c_i_Bohm * e
       j_esat = 1/4 n_e * c_e_Bohm * e ( 1-g)
       j_quot = 0.15 * sqrt( (a_mass*c_mp) /c_me );  
    */

    j_isat = n_e * c_e * sqrt( 4.* c_e * t_e/ (a_mass*c_mp) );


    /* L_debye */
    ld = sqrt(t_e * c_eps0 / (n_e * c_e)) ;

   
    /* ---------------------------------------------------------------------
       wenn a_sonde_proj > 0.0 dann wurde eine geometrische Sondenfl"ache
       und ein Projektionswinkel angegeben. Dies hei"st, da"s f"ur eine
       flush mounted Sonde eine Kennlinie unter Ber"ucksichtigung der
       Schichteffekte berechnet werden soll.
       --------------------------------------------------------------------- */

    if ( a_sonde_proj > 0.0 )  
	{

        /* potentialunabhaengige Vorfaktoren vor Schichtdicke */
        sf1 = ld / 3.0 * alpha1 / sqrt(cos1);	
        sf2 = ld / 3.0 * alpha2 / sqrt(cos2);
        temp = sqrt( cos2 * inv_beta / cos1 );
        lfac = temp / a_height;
        bfac = temp / a_width;

        /* ab jetzt Berechnung fuer alle Spannungswerte */
        I_isat = a_sonde_proj * j_isat;

        for ( ind=0, ptr_vpr=vpr, ptr_strom=strom; ind < n_vars; ind++ ) 
            {

            /* Spannung an Doppelsonde: phi_sonde */
            delta_phi = (dv - *ptr_vpr++) / t_e;

	    /* Nichtsaettigungskoeffizienten muessen immer groesser als -1 */
            if (delta_phi < 0) {
		ns_add1 = 1.0;
	        delta2  = sf2 * quick_075( delta_phi );
	        ns_add2 = ( (delta2+delta2)*bfac +1.) * 
				( 1. + delta2*tan2p*lfac );
		}
	    else {
	    	delta1  = sf1 * quick_075( -delta_phi );
	    	ns_add1 = ( (delta1+delta1)/a_width + 1.) * 
				( 1. + delta1*tan1p/a_height );
		ns_add2 = 1.0;
		}

            /* Berechnung des Strom-Vektors */
	    temp = exp( delta_phi );
	    /* temp = exp_lt( delta_phi ); */ 

            *ptr_strom++ =  I_isat * 
			(ns_add2 - ns_add1 * temp) / ( inv_beta + temp );

            }	/* for ind=0, n_vpr-1 */

        }	/* if (a_sonde_proj > 0.0) .....    */


    /* ---------------------------------------------------------------------
       wenn a_sonde_perp > 0.0 dann wurde bereits eine effektive Sondenfl"ache
       senkrecht zur Feldrichtung eingegeben. Es wird davon ausgegangen, da"s
       es sich hierbei um eine freistehende Sonde handelt, die von beiden
       Seiten Ionen und Elektronen aufsammeln kann. 

       Falls a_sonde_perp > 0.0 und a_sonde_proj > 0.0 dann wird eine 
       Kombination beider Effekte erwartet, d.h. der Gesamtstrom zur Sonde
       setzt sich additiv aus beiden Teilen zusammen.
       --------------------------------------------------------------------- */

    if (a_sonde_perp > 0.0)
	{

        /* Initialisierung, falls n"otig */
        if (a_sonde_proj <= 0.) {
            for ( ind=0, ptr_strom=strom; ind<n_vars; ind++, *ptr_strom++ =0.0);
	    }

        /* Kennlinie einer Doppelsonde mit Potentialverschiebung und
	   Fl"achenverh"altnis beta */
        I_isat = 2.* a_sonde_perp * j_isat;
        bfac = alpha2 * ld / a_width; 


        for ( ind=0, ptr_vpr=vpr, ptr_strom=strom; ind < n_vars; ind++ ) 
	    {

            /* Spannung an Doppelsonde: phi_sonde */
            delta_phi = (dv - *ptr_vpr++) / t_e;
	    temp = exp( delta_phi );

            /* Nichts"attigungsverhalten im Elektronenast? */
	    /* ns_add2 = bfac * quick_075(delta_phi) ; */
            /* Berechnung des Strom-Vektors */
            *ptr_strom++ +=  I_isat * 
		(1. + bfac * quick_075(delta_phi) - temp) / ( inv_beta + temp );

	    }


	}	/* if (a_sonde_perp > 0.0) .....    */


    return 0 ;

    }		 /* fast_doppel_nl */







double fmin(x, n, vecfunc, vpr, ipr, icurr, init_fast_doppel_nl, params )

    double	*x;
    int		n;
    int 	(*vecfunc)();
    double	*vpr, *ipr, *icurr;
    int		init_fast_doppel_nl;
    double	*params;

    {
	int 	i;
	double 	sum;

/* printf(":"); */

	(vecfunc)(vpr,icurr, x, n, init_fast_doppel_nl, params);

	for (sum=0.0,i=0;i<n;i++) sum += DSQR(icurr[i]-ipr[i]);

	return 0.5*sum;

        } 	/* fmin */



int fdjac(n, x, icurr, df, vecfunc, vpr, ipr, params )

    int		n;
    double	*x, *icurr, **df;
    int 	(*vecfunc)();
    double	*vpr, *ipr, *params;


    {
	int 	i,j, k;
	double 	h,temp;

	static double	*f=NULL;

	static double	eps=1.e-6;

/* printf("."); */

	/* Standard-Initialisierung, n<=3 in allen Faellen */
	if ( f==NULL ) f=dvector(0,2);


	for (j=0;j<n;j++) {

		temp=x[j];
		h=eps*fabs(temp);

		/* die folgenden drei Zeilen dienen dazu, die Auswirkungen
		   von Rundungsfehlern zu minimieren */
		if (h == 0.0) h=eps;
		x[j]=temp+h;
		h=x[j]-temp;

		/* Funktionsauswertung bei kleinen Abweichungen in x*/
		(*vecfunc)(vpr, f, x, n, 0, params);

		/* partielle Ableitungen */
		for (i=0;i<n;i++) df[i][j]=(f[i]-icurr[i])/h;

		/* Rueckschreiben der Parameter */
		x[j]=temp;

		}


	/* keine Freigabe des Speicherplatzes, erfolgt am PRG-Ende automatisch
	free_dvector(f,1,n); */

    }	/* fdjac */





int lnsrch(n, xold, fold, g, p, x, f, stpmax, check, 
		vpr, ipr, icurr, params)

    int	n;
    double	xold[], fold, g[], p[], x[], *f, stpmax;
    int		*check;
    double	*vpr, *ipr, *icurr, *params;

    {
	int 	i, n_try ;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	static int	n_try_max = 200;
	static double 	tolx=1.e-10;

/* printf("x"); */

	*check=0;
	for (sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=0;i<n;i++)
		slope += g[i]*p[i];

	test=0.0;
	for (i=0;i<n;i++) 
		test = DMAX( test, fabs(p[i])/DMAX(fabs(xold[i]),1.0) );
	alamin=tolx/test;
	alam=1.0;

	/* Endlosschleife, Funktion wird durch Returns innerhalb verlassen */
	for (n_try=0; n_try < n_try_max; n_try++) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];

/* printf("^"); */
		*f=fmin(x,n,fast_doppel_nl, vpr, ipr, icurr, 0, params);

		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
			*check=1;
			return 0;
		} else if (*f <= fold + 1.e-4*alam*slope) return 0;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/
				    	(alam-alam2);
				b = (-alam2*rhs1/(alam*alam)+
				     	alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) 
				    tmplam = -slope/(b+b);
				else {
				    disc=b*b-3.0*a*slope;
				    if (disc<0.0) {
				        /* fprintf(stderr, "Roundoff problem in lnsrch.\n"); */
					return -21;
					}
				    else 
					tmplam=(-b+sqrt(disc))/(3.0*a);
				    }
				tmplam = DMAX( tmplam, 0.5*alam );
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=DMAX(tmplam,0.1*alam);
		}

	/* wird nie erreicht, da Endlosschleife */

        /* um Geschwindigkeit zu gewinnen, wird nach n_try_max Versuchen
	   die Iteration abgebrochen und zum Nebenminimum erklaert */
	*check = 1;
	return 0;

	}	/* lnsrch */








/*	--------------------------------------------------------------------
	
				newt
				====

	global konvergentes Newton-Verfahren zur Loesung eines 
	nichtlinearen Gleichungssystem, entnommen aus Numerical
	Recipes in C, angepasst an spezielles Problem der Kennlinien.


	x 	Parametervektor[n], 2 oder 3dim
	n	Anzahl der Parameter
	check	FALSE (=0), alles ok
		TRUE  (=1), lediglich lokales Minimum gefunden
	vpr	Sondenspannung[n]
	ipr	Sondenstrom[2]
	params	Funktionsparameter[20], vergl. mag_fit

	----------------------------------------------------------------    */


int newt( x, n, check, vpr, ipr, params, f )

    double	*x;
    int		n, *check;
    double	*vpr, *ipr,*params, *f;


    {
	int 	i,its,j;

	double 	d, den, fold, stpmax, sum, temp, test; 


	static int	*indx=NULL;
	static double  	**fjac=NULL, *g=NULL, *p=NULL, *xold=NULL,
			*icurr=NULL, **sfjac; 

	int	init_fast_doppel_nl, rc;

	static int	maxits=100, maxsteps=100;
	static double	tolf=1.e-8, tolmin=1.e-10, tolx=1.e-10;


/* printf("Y"); */

	/* Lege benoetigte Datenfelder an, es gilt immer n<=3 */
	if ( indx == NULL ) {
	    indx = ivector(0,2);
	    fjac = dmatrix(0,2,0,2);
	    sfjac = dmatrix(0,2,0,2);
	    g = dvector(0,2);
	    p = dvector(0,2);
	    xold = dvector(0,2);
	    icurr = dvector(0,2);

	    /* Initialisierung  der loop-up-tables */
	    init_quick_075();

	    }



	/* Fit erfolgt nicht direkt in n_e, T_e und beta */
    	x[0] = log(x[0]*sqrt(x[1]));		/* n_e sqrt(T_e) > 0.0 */
	x[1] = log( fabs(x[1]) ); 		/* T_e > 0.0 */
	if ( n>2 ) x[2] = log( fabs(x[2]) );	/* beta > 0.0 */


	/* Inititalisiere Parameter zur Funktionsberechnung,
	   berechen Funktionswerte "icurr" zu Startparametern */
	init_fast_doppel_nl = 1	; 	/* TRUE */
	*f=fmin(x, n, fast_doppel_nl, vpr, ipr, icurr, 
				init_fast_doppel_nl, params);
	init_fast_doppel_nl = 0	; 	/* FALSE */

	/* maximale Abweichung der Funktionswerte von Sollwerten */
	test=0.0;
	for (i=0;i<n;i++) test = DMAX( test, fabs(icurr[i]-ipr[i]) );
	if (test<0.01*tolf) {
	    /* Variablen muessen wieder retransformiert werden */
    	    x[0] = exp(x[0]-0.5*x[1]);
    	    x[1] = exp(x[1]);
            if (n>2) x[2] = exp(x[2]) ;
	    return 0;
	    }

	/* Schrittweite fuer lineare Suche */
	for (sum=0.0,i=0;i<n;i++) sum += DSQR(x[i]);
	stpmax= maxsteps*DMAX(sqrt(sum),(float)n);

    	/* Beginn der Iteration */
	for (its=0;its<maxits;its++) {

		/* berechne partielle Ableitungen */
		fdjac(n,x,icurr,fjac,fast_doppel_nl, vpr, ipr, params);
		for (i=0;i<n;i++) {
		    for (sum=0.0,j=0;j<n;j++) 
			sum += fjac[j][i]*(icurr[j]-ipr[j]);
		    g[i]=sum;
		    }

		/* Speichere Ist-Zustand zu Vergleichszwecken */
		for (i=0;i<n;i++) xold[i]=x[i];
		fold=*f;


		for (i=0;i<n;i++) p[i] = -(icurr[i]-ipr[i]);

		/* Loese lineares Gleichungssystem: J * dx = -F */
		if ( (rc = dsvd_solver( fjac, p, n)) != 0 ) break;

		/* lineare Suche entlang Newton-Richtung */
		if ( (rc = lnsrch(n, xold, fold, g, p, x, f, stpmax, check, 
				vpr, ipr, icurr, params)) != 0 ) break;

		/* Konvergenzkriterien erfuellt? */
		test=0.0;
		for (i=0;i<n;i++)
			test = DMAX( test, fabs(icurr[i]-ipr[i]));
		if (test < tolf) {
			*check=0;
			rc = 0;
			break;
			}

		/* Gradienten gleich Null */
		if (*check) {
			test=0.0;
			den=DMAX(*f,0.5*n);
			for (i=0;i<n;i++) 
			     test = DMAX( test, 
				       fabs(g[i])*DMAX(fabs(x[i]),1.0)/den );
			*check=(test < tolmin ? 1 : 0);
		    	rc =  (130 + *check);
			break;
			}

		/* keine nennenswerte Aenderung in x mehr moeglich */
		test=0.0;
		for (i=0;i<n;i++) 
		    test = DMAX( test, 
				(fabs(x[i]-xold[i]))/DMAX(fabs(x[i]),1.0) );
		if (test < tolx) {
		    rc = 100;
		    break;
		    }

		}


	/* keine Konvergenz nach endlicher Anzahl von Iterationsschritten? */
	if ( its >= maxits ) {
	   /* printf("newt: failt to converge, its > maxits = %d\n", maxits); */
	   rc = -1;
	   }

	/* Variablen muessen wieder retransformiert werden */
    	x[0] = exp(x[0]-0.5*x[1]);
    	x[1] = exp(x[1]);
        if (n>2) x[2] = exp(x[2]) ;
	
	/* Rueckgabe = aktueller Statuswert */
	return rc;


	}	/* newt */




/* 	-------------------------------------------------------- 

			mag_solver 
			==========

 	Solver fuer nichtlineares Gleichungssystem, um Daten von 
 	Mehrfachsonden, z.B. triple- oder tetra-probes zu bearbeiten. 
 	Die Anzahl der Gleichungen (=Wertepaare Strom und Spannung)  
 	sowie der freien Parameter (n_e, T_e und ggf. beta) muss 
 	uebereinstimmen.


  	-------------------------------------------------------- 

	M. Weinlich, 19.06.95

 	-------------------------------------------------------- */

int mag_solver( data_x, data_y, cos_psi, n_points, n_data, pars, n_pars, 
                ne_result, te_result, beta_result, fmin )


    double	*data_x, *data_y, *cos_psi, *fmin;
    double 	*ne_result, *te_result, *beta_result;
    int		*n_points, *n_data;
    double	*pars;
    int		*n_pars;
    
    {


    double	*ptr_x, *ptr_y, *ptr_ne, *ptr_te, *ptr_beta, *ptr_f, *ptr_cp;
    double	vars[6];	/* Variable, in denen nl. Gl.s. geloest wird */

    int		rc, rc_counter, i, j;
    int		check;
    double	save_ne, save_te, save_beta;



/*
    printf( "Uebergebene Werte in mag_solver: \n\n");
    printf( "n_points = %d \n", *n_points );
    printf( "n_data = %d \n", *n_data );
    printf( "n_pars = %d \n", *n_pars );
    printf( "x = ");
    for (i=0; i<10; i++) printf( "%g  ", data_x[i] );
    printf("\n");
    printf( "y = ");
    for (i=0; i<10; i++) printf( "%g  ", data_y[i] );
    printf("\n");
    printf( "pars = ");
    for (i=0; i<20; i++) printf( "%g  ", pars[i] );
    printf("\n");
    printf("\n\n");
*/

   /* hier kann Verzweigung in Fit fuer penta-probes erfolgen,
      vorlaeufig nicht implementiert */
   if ( *n_points != *n_pars ) {
	fprintf( stderr, "\n\n\n");
	fprintf( stderr, "run_nl_solver: n_points != n_pars \n");
	fprintf( stderr, "               n_points = %d \n", *n_points);
	fprintf( stderr, "               n_pars = %d \n\n", *n_pars);
	fprintf( stderr, "		 Abbruch ... \n\n\n");
	return -42;
	}

    /* Eigentlicher nl-Solver wird nun n_data-Mal durchgelaufen und 
       fittet in jedem Durchgang n_points Datenpunkte */
    ptr_x = data_x;
    ptr_y = data_y;
    ptr_cp = cos_psi;
    ptr_f = fmin;

    ptr_ne   = ne_result;
    ptr_te   = te_result;
    ptr_beta = beta_result;

    save_ne   = pars[0];
    save_te   = pars[1];
    save_beta = pars[3];

    vars[0] = save_ne; 
    vars[1] = save_te; 
    vars[2] = save_beta;

    for ( i=0; i<*n_data; i++ ) {

        rc_counter = -1;

	/* Winkel muss jeweils aktualisiert werden */	
	vars[6] = *ptr_cp;
	vars[7] = *ptr_cp;

        do {

	    switch (rc_counter) {
		case 0 :	vars[0]=save_ne; vars[1]=save_te; 
				vars[2]=save_beta; break;
		case 1 :	vars[0]=1.e18; vars[1]=10.; vars[2]=4.; break;
		case 2 :	vars[0]=1.e20; vars[1]=10.; vars[2]=4.; break;
		case 3 :	vars[0]=1.e18; vars[1]=30.; vars[2]=4.; break;
		case 4 :	vars[0]=1.e20; vars[1]=30.; vars[2]=4.; break;
		default:	;
		}

	    rc = newt( vars , *n_pars, &check, ptr_x, ptr_y, pars, ptr_f);

	    rc_counter++;

	    }

	    /* rc    = lokale Fehler in newt und Subroutinen
	       check = TRUE, wenn lokales Minimum gefunden
	       fmin  = 0.001 ist etwa Messgenauigkeit */
            while( ((rc != 0) || (check) || (*ptr_f > 1.e-3) ) && 
			(rc_counter < 5) );
	

	/* Hallo, ich lebe noch .... */
        if (i%1000 == 0) printf(".");

        *ptr_ne++   = vars[0];
        *ptr_te++   = vars[1];
        *ptr_beta++ = vars[2];

	if (rc != 0) {
	    vars[0] = save_ne; 
	    vars[1] = save_te; 
	    vars[2] = save_beta;
	    *ptr_f = -*ptr_f;
	    }

	ptr_x = ptr_x + *n_points;
	ptr_y = ptr_y + *n_points;
	ptr_f++;
	ptr_cp++;


	}
 
    printf("\n");

    return 0;

    }		/* mag_solver */









/* 	-------------------------------------------------------- 


		Es folgen die Interface-Routine zu IDL
		======================================

 	-------------------------------------------------------- */








/* 	-------------------------------------------------------- 

			run_triple_fit 
			==============

 	Interface zu call_triple_fit. Wird von IDL aus aufgerufen und 
	liest uebergebene Parameter aus und startet eine leicht vereinfachte
        Version des Fits, um Daten von Mehrfachsonden, z.B. triple- oder
	penta-probes zu bearbeiten.


  	-------------------------------------------------------- 

	M. Weinlich, 01.09.93

 	-------------------------------------------------------- */

int run_triple_fit( argc, argv )

    IDL_INT	argc;
    void	*argv[];

    {

    /* user defined functions */
    extern int	magfit();

    static double	**sigma_pars=NULL;
    static int		old_n_sigma;

    IDL_DOUBLE	*data_x, *data_y, *w, *yfit;
    IDL_DOUBLE 	*ne_result, *te_result, *beta_result;
    IDL_INT	*n_points, *n_data;
    IDL_DOUBLE	*pars;
    IDL_INT	*do_var;
    IDL_INT	*n_pars;
    IDL_DOUBLE	*nb_values, *nb_const;
    IDL_DOUBLE	*chi_sqr;
    IDL_INT	*iter_max;
    IDL_DOUBLE	*eps_abs, *eps_rel;

    double	*ptr_x, *ptr_y, *ptr_ne, *ptr_te, *ptr_beta;

    int		rc, rc_counter, i, j;
    int		use_iter_max;
    double	use_eps_abs, use_eps_rel;
    double	save_ne, save_te;

    void  	*old_fpe_state;


    /*   Auslesen der uebergebenen Argumente   */
    data_x      = (IDL_DOUBLE *) argv[0];
    data_y      = (IDL_DOUBLE *) argv[1];
    n_points    = (IDL_INT *) 	 argv[2];
    n_data      = (IDL_INT *) 	 argv[3];
    pars        = (IDL_DOUBLE *) argv[4];
    do_var      = (IDL_INT *)    argv[5];
    n_pars      = (IDL_INT *) 	 argv[6];
    ne_result   = (IDL_DOUBLE *) argv[7];
    te_result   = (IDL_DOUBLE *) argv[8];
    beta_result = (IDL_DOUBLE *) argv[9];
    chi_sqr     = (IDL_DOUBLE *) argv[10];
    iter_max    = (IDL_INT *) 	 argv[11];
    eps_abs     = (IDL_DOUBLE *) argv[12];
    eps_rel     = (IDL_DOUBLE *) argv[13];



    /* Hilfsfeld */
    if ( sigma_pars == NULL ) {
        sigma_pars = dmatrix(0, *n_pars-1, 0, *n_pars-1);
	old_n_sigma = *n_pars;
	}
    if ( *n_pars > old_n_sigma ) {
        printf( "sigam_par geloescht mit %d Elementen\n", old_n_sigma );
        free_dmatrix( sigma_pars, 0, old_n_sigma-1, 0, old_n_sigma-1);
        sigma_pars = dmatrix(0, *n_pars-1, 0, *n_pars-1);
	old_n_sigma = *n_pars;
        }

    w    = dvector(0, *n_points-1);
    yfit = dvector(0, *n_points-1);

    nb_values = dvector(0, *n_pars-1);
    nb_const = dvector(0, *n_pars-1);

    for ( i=0; i<*n_points; i++) {
	w[i] = 1.;
	yfit[i] = 0.;
 	}
    for ( i=0; i<*n_pars; i++) {
	nb_values[i] = 0.;
	nb_const[i] = 0.;
 	}

    /* Eigentliche Fitroutine wird nun n_data-Mal durchgelaufen und 
       fittet in jedem Durchgang n_points Datenpunkte */
    ptr_x = data_x;
    ptr_y = data_y;
    ptr_ne = ne_result;
    ptr_te = te_result;
    ptr_beta = beta_result;

    save_ne = pars[0];
    save_te = pars[1];

    /* ignoriere floationg point errors, da diese IDL durcheinander bringen */
    /* if ( ( signal( SIGFPE, SIG_IGN )) == SIG_ERR ) 
	fprintf( stderr, "Fehler beim Umbiegen der floating execptions\n"); */

    for ( i=0; i<*n_data; i++ ) {

        rc_counter = 0;

        do {

/* 	printf("rc_counter = %d \n", rc_counter); */

	    switch (rc_counter) {
		case 1 :	pars[0]=1.e18; pars[1] = 10.; break;
		case 2 :	pars[0]=1.e20; pars[1] = 10.; break;
		case 3 :	pars[0]=1.e18; pars[1] = 30.; break;
		case 4 :	pars[0]=1.e20; pars[1] = 30.; break;
		default:	;
		}

            use_iter_max = *iter_max;
 	    use_eps_abs  = *eps_abs;
            use_eps_rel  = *eps_rel;

	    rc = magfit(fast_doppel, ptr_x, ptr_y, w, yfit, n_points, 
			pars, do_var, n_pars, nb_values,  nb_const, sigma_pars, 
			chi_sqr, &use_iter_max, &use_eps_abs, &use_eps_rel);

	    rc_counter++;

	    }
            while( (rc != 0) && (rc_counter < 5) );
	
/*	if (( *chi_sqr > 1.e-3 ) && (rc != 0)) exit(-1); */

        *ptr_ne++ = pars[0];
        *ptr_te++ = pars[1];
        *ptr_beta++ = pars[3];

	if (rc == 0) {
	    save_ne = pars[0];
	    save_te = pars[1];
	    }
 	else {
	    *chi_sqr = -*chi_sqr;
	    pars[0] = save_ne;
	    pars[1] = save_te;
	    }

	ptr_x = ptr_x + *n_points;
	ptr_y = ptr_y + *n_points;

        chi_sqr++;

	}


    free_dvector( w, 0, 20 );
    free_dvector( yfit, 0, 20);
    free_dvector( nb_values, 0, 20);
    free_dvector( nb_const, 0, 20);


    /* Standard-Behandlung der floationg point errors */
    /* signal ( SIGFPE, SIG_DFL ); */
 
    return (rc);

    }		/* run_triple_fit */









/* 	-------------------------------------------------------- 


		Interface zu IDL: 
	 	=================

	prueft ob Routine bekannt,
	prueft ob Anzahl der uebergebenen Parameter ok
	ruft entsprechende Subroutine auf

	Syntax in idl: 
	rc=CALL_EXTERNAL('/u/mnw/bin/@sys/libmag.so','idl2mag','name',parameter)

	name         Name of the c subroutine
	parameter    Parameters for the subroutine, all integers must be Long!
	result        0 : successful call
		     -1 : unknown subroutine
		     -2 : wrong number of parameters
		     -3 : number of parameters < 1
		     otherwise the error code
		     
	-------------------------------------------
	
	original (libddww.so)   19.09.94 by arb
	current version 	28.07.95 by mnw
	
 	-------------------------------------------------------- */






/*	--------------------------------------------------------

			idl2mag_err
			===========
			
	beende Routine idl2mag bei Auftreten eines Fehlers
	
	mnw, 28.07.95
	
	------------------------------------------------------- */

int idl2mag_err( err, routine, nargs )

int 	err;
char	*routine;
int	nargs;
 
    {
    int		err_no_routine=-3;
    int		err_wrong_num_arg=-2;
    int		err_unknown_routine=-1;
    char	*err_log="error in idl2mag: ";
    /* Testphase */
    /* printf("idl2mag_err, err = %d \n", err); */
    
    /* Routine ohne Fehler beendet */
    if (err == 0 ) {
       /* fprintf( stderr, "\n\nidl2mag: routine %s successful completed! \n\n",
       		routine); */
       return err;
       }
       
    /* Fehler in Routine oder bei Aufruf der Routine
       ==> kurze Nachricht an user ueber stderr,
           Rueckgabe des entsprechenden Fehlercodes */
    if ( err == err_no_routine ) {
       fprintf( stderr, "%s no subroutine specified\n", err_log);
       }
    else if ( err == err_wrong_num_arg ) {
       fprintf( stderr, "%s routine %s called with wrong number of arguments (%d)\n",
       		err_log, routine, nargs);
       }
    else if ( err == err_unknown_routine ) {
       fprintf( stderr, "%s called routine %s not supplied\n", err_log, routine);
       }
    else {
       /* fprintf( stderr, "%s error %d in routine %s\n", err_log, err, routine); */
       }
 
    return err;
    
    }	/* idl2mag_err */
    
    
    
    

/*	--------------------------------------------------------

			idl2mag
			=======
			
	eigentliches Interface
	
	mnw, 28.07.95
	
	------------------------------------------------------- */
	
int idl2mag (argc, argv)

int 	argc;
void	*argv[];

    {
    int		nargs;
    char	routine[80];
    char	routine_lowcase[80];
    
    int		err_no_routine=-3;
    int		err_wrong_num_arg=-2;
    int		err_unknown_routine=-1;
    
    int 	i, err, rc;
    
    /* number of arguments */
    nargs = argc-1;
    
    /* Testphase */
/*    printf( "\nidl2mag, uebergebenen Argumente: \n");
    printf( "argc = %d \n", argc);
    printf( "argv = ");
    for (i=0; i<argc; i++) printf( "%s, ", argv[i]);
    printf( "\n\n");
*/

    /* keine Argumente, keine Routine */
    if (argc < 1 ) return idl2mag_err( err_no_routine, routine, nargs );
    
    /* welche Routine wird gewuenscht */
    strcpy( routine, IDL_STR(argv[0]) );
    strcpy( routine_lowcase, IDL_STR(argv[0]) );
    lowcase(routine_lowcase);
    
    /* printf("idl2mag: aufzurufende Routine = %s \n\n", routine); */
    
    if ( strcmp( routine_lowcase, "mag_doppel") == 0 ) {
	/* Routine benoetigt insgesamt 14 Uebergabe-Variablen:
	   x, y, w, yfit, n_x, pars, do_var, n_pars, nb_values,  
	   nb_const, chi_sqr, iter_max, eps_abs, eps_rel
	*/
	if ( nargs != 14 ) {
	    return idl2mag_err( err_wrong_num_arg, routine, nargs );
	    }
	else {
	    rc = magfit( mag_doppel, argv[1], argv[2], argv[3], argv[4], 
	    		 argv[5], argv[6], argv[7], argv[8], argv[9], argv[10],
	    		 argv[11], argv[12], argv[13], argv[14]);
	    return idl2mag_err(rc, routine, nargs);
	    }
        }
    else if ( strcmp( routine_lowcase, "fast_doppel") == 0 ) {
	/* Routine benoetigt insgesamt 14 Uebergabe-Variablen:
	   x, y, w, yfit, n_x, pars, do_var, n_pars, nb_values,  
	   nb_const, chi_sqr, iter_max, eps_abs, eps_rel
	*/
	if ( nargs != 14 ) {
	    return idl2mag_err( err_wrong_num_arg, routine, nargs );
	    }
	else {	    
	    rc = magfit( fast_doppel, argv[1], argv[2], argv[3], argv[4], 
	    		 argv[5], argv[6], argv[7], argv[8], argv[9], argv[10],
	    		 argv[11], argv[12], argv[13], argv[14]);
	    return idl2mag_err(rc, routine, nargs);
	    }
        }
        
      else if ( strcmp( routine_lowcase, "mag_solver") == 0 ) {
	/* Routine benoetigt insgesamt 10 Uebergabe-Variablen:
	   data_x, data_y, cos_psi, n_points, n_data, pars, n_pars, 
           ne_res, te_res, beta_res, fmin 
        */
	if ( nargs != 11 ) {
	    return idl2mag_err( err_wrong_num_arg, routine, nargs );
	    }
	else {	
	    rc = mag_solver( argv[1], argv[2], argv[3], argv[4], 
	    		 argv[5], argv[6], argv[7], argv[8], argv[9], 
			 argv[10], argv[11] );
	    return idl2mag_err(rc, routine, nargs);
	    }
        }      
    else {
        /* printf( "idl2mag: unbekannte Subroutine\n"); */
	return idl2mag_err(err_unknown_routine, routine, nargs);
        }
    
    }	/* idl2mag */
    
    
    
    




/* 	-------------------------------------------------------- 


		kurzes Hauptprogramm, das die Routinen testet
		=============================================

 	-------------------------------------------------------- */





/* 	Dummy Hauptprogramm	*/

int main()
{

	int	i, j, n;
	int	rc;
	double	el_t;


    /* Local variables */
    int 		iter_max, k;
    double 	w[200], **sigma_par, chi, fit[200];
    int			do_var[20];
    int 		n_x, err;
    double 	inv[100], vpr[200], ipr[200], save[100], x[3];
    int 		used, unit;
    double 	fixed[13], strom[200], params[20];
    double 	eps_abs, eps_rel;

    int			n_pars;
    double		*nb_values, *nb_const;

    FILE		*demo_file;

    int 		c_6 = 6;
    char		*f_name = "/afs/ipp/u/mnw/idl/single/demo.dat";
    char		*f_mode = "r";

    float	temp;

    double	fmin;
    int		check;


	/*
		-------------------------
		Test fuer linearen Solver
		-------------------------
	*/


    double	**mat, **mat_inv;
    double	*vec, *solve;


    init_exp_lt();

    goto test_fit; 
    /* goto test_nl; */

    el_t = ((double) clock()) ;
    for (i=0; i<100000; i++) temp = exp(i*0.0001 -5.);
    el_t = ((double) clock()) - el_t;
    printf("Zeitaufwand exp(x): %g \n", el_t*1.e-6);

    el_t = ((double) clock()) ;
    for (i=0; i<100000; i++) temp = exp_lt(i*0.0001 -5.);
    el_t = ((double) clock()) - el_t;
    printf("Zeitaufwand exp_lt(x): %g \n", el_t*1.e-6);



   /* erzeuge Testmatrix */
   n = 3;
   mat = dmatrix(0,n-1,0,n-1);
   mat_inv = dmatrix(0,n-1,0,n-1);
   vec = dvector(0,n-1);
   solve = dvector(0,n-1);

    /* Initialisierung */
    mat[0][0] = 1.;
    mat[0][1] = 5.;
    mat[0][2] = 4.;
    mat[1][0] = 2.;
    mat[1][1] = 7.;
    mat[1][2] = 1.;
    mat[2][0] = 4.;
    mat[2][1] = 2.;
    mat[2][2] = 5.;

    vec[0]  = 3.;
    vec[1]  = -1;
    vec[2]  = -4;


	/* Ausgabe der Startwerte */
	printf("\n\n\n\n\n\n\n\n");
	printf("Loese lineares Geleichungssystem Ax=b ueber singular value decomposition\n\n");

	printf(" A =");
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		printf( "\t %g ", mat[i][j]);
		}
	    printf("\n");
	    }

	printf("\n b =");
	for (i=0; i<n; i++ ) printf( "\t %g", vec[i] );
	printf("\n\n");

    /* loese das System mat * x = vec */
    el_t = ((double) clock()) ;
    for (i=0; i<10000; i++) rc = dsvd_solver( mat, vec, n);
    el_t = ((double) clock()) - el_t;
    el_t *= 1.e-6; 



	/* Loeschen der Felder */
	rc = dsvd_solver( mat, vec, -n );

	printf("Loesung des Systems: \n");
	printf(" Inv =");
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		printf( "\t %g ", mat[i][j]);
		}
	    printf("\n");
	    }

	printf("\n x =");
	for (i=0; i<n; i++ ) printf( "\t %g", vec[i] );
	printf("\n\n");

    printf("\n");
    printf("Zeitaufwand = %g sec\n\n", el_t);

if (n>0)  return 0;

	/*
		-------------------------
		Test fuer Fit-Algorithmus
		-------------------------
	*/


test_fit:


    /* 	Einlesen eines Demo-Files */

    demo_file = fopen( f_name, f_mode);

    /* Fehler beim Oeffnen? */
    if ( ! demo_file ) {
	printf ("Fehler beim Oeffnen des Files %s \n\n", f_name );
	}

    /* Laenge der Datensaetze */
    fscanf( demo_file, "%d", &n_x );

    /* Einlesen der Daten */
    for ( i=0; i<n_x; i++ )
	{
	fscanf ( demo_file, "%f", &temp );
	vpr[i] = (double) temp;
	}
    for ( i=0; i<n_x; i++ )
	{
	fscanf ( demo_file, "%f", &temp );
	strom[i] = (double) temp;
	}

    fclose( demo_file );

/*
    printf( "Eingelesene Laenge des Datensatzes: %d \n\n", n_x );

    printf( "    Sondenspannung        Sondenstrom \n\n" );
    for ( i=0; i<n_x; i++ )
	{
	printf( "\t %f \t %f \n", vpr[i], strom[i]);
	}
    printf( "\n\n");
*/

    /* 	Vorbereitung auf den Fit */

    params[0] = 1e18;			/* n_e */
    params[1] = 10.;			/* t_e */
    params[2] = 1.;
    params[3] = 1.;
    params[4] = 0.;
    params[5] = 1.;

    params[6] = 1.5603243;		/* 89.4 Grad */
    params[7] = 1.563815;		/* 89.6 Grad */

    params[10] = .005;			/* a_width */
    params[11] = .04;			/* a_height */
    params[12] = 2.;			/* a_mass */
    params[13] = 1.;			/* a_z */
    params[14] = 2.;			/* a_bt */

    params[17] = 0.00041887;		/* a_height * cos(psi1) */

    for (i=0; i<6; i++ ) do_var[i]=1;
    for (i=6; i<20; i++ ) do_var[i]=0;
    
    /* kein Fit in alpha1, alpha2, psi1, psi2, verwende projizierte Fl"ache */
/*    do_var[2] = 0.;
    do_var[5] = 0.;
    do_var[6] = 0.;
    do_var[7] = 0.;
    params[11] = -1.; */
 

    /* nur flush mountet probe */
    params[17] = -1.; 

    /* konstante Gewichtung */
    for (i=0; i < n_x; ++i) {
	w[i] = 1.;
    	}

    n_pars = 20;
    sigma_par = dmatrix(0,n_pars-1,0,n_pars-1);

    nb_values = dvector(0, n_pars-1);
    nb_const  = dvector(0, n_pars-1);
    for (i=0; i<n_pars; i++) nb_values[i]=0.0;
    for (i=0; i<n_pars; i++) nb_const[i]=0.0;

/*    do_var[2] = 0; */
/*    nb_values[2] = params[2];
    nb_const[2] = 20.; 
    do_var[6] = 1;
    nb_values[6] = params[6];
    nb_const[6] = 200.;
*/

    /* 	eigentlicher Fit */

    printf ( "und los geht`s mit dem Fit \n");

    el_t = ((double) clock()) ;

    /* Fuer Geschwindigkeitsvergleich 50 x Fit */
    for (i=0; i<50; i++ ) {

      params[0] = 1.e18 * ( 1. + i*0.02);
      params[1] = 10. * (3. - i*0.05);
    /* bis hier Einschub fuer Geschwindigkeitsvergliech */

    eps_rel  = 1.e-4;
    eps_abs  = 0.;
    iter_max = 200;


/*    magfit(mag_doppel, vpr, strom, w, fit, &n_x, params, do_var, &n_pars, 
		nb_values, nb_const, 
		&chi, &iter_max, &eps_abs, &eps_rel);
*/
    magfit(fast_doppel, vpr, strom, w, fit, &n_x, params, do_var, &n_pars, 
		nb_values, nb_const, &chi, &iter_max, &eps_abs, &eps_rel);


    /* Fuer Geschwindigkeitsvergleich 50 x Fit mit vergleichbaren Startwerten */
    }


    el_t = ((double) clock()) - el_t;
    el_t *= 1.e-6;

    /* 	Ausgabe der Ergebnisse */

    printf( "Ergebnisse des Fits: \n");
    printf( "-------------------- \n\n\n");

    printf( "chi     = %g \n", chi );
    printf( "eps_rel = %g \n", eps_rel );
    printf( "eps_abs = %g \n", eps_abs );
    printf( "iter    = %d \n", iter_max );
    printf( "Zeitaufwand: %g sec CPU-Zeit\n", el_t);

    printf( "\n\n \t Parameter \t sigma \n");
    for ( i=0; i<6; i++) {
 	printf( "\t %g \t %g \n", params[i], sigma_par[i][i] );
 	}


    printf("\n\n\n");

    return 0;


    /* nichtlineare Gleichungssolver fuer Kennlinie einer Triple-Sonde */

test_nl:

    x[0] = 1.e19;
    x[1] = 10.;
    x[2] = 4.;

    vpr[0] = -50.;
    vpr[1] = 10.;
    vpr[2] = -10.;

    ipr[0] = -1.;
    ipr[1] = 1.;
    ipr[2] = 0.3;

    n = 2;

    for (i=0; i<20; i++) params[i]=0.;
    params[2] = 1.;
    params[3] = 4.;
    params[4] = 0.;
    params[5] = 1.;
    params[6] = 1./30.;		/* cos(psi1) */
    params[7] = params[6];	/* cos(psi2) */
    params[10] = 0.005;
    params[11] = 0.030;
    params[12] = 2.;
    params[13] = 1.;
    params[14] = 2.;
    params[16] = 1.;

    el_t = ((double) clock()) ;
    for (i=0;i<1;i++) {
	/* printf("i = %d\n",i); */
	x[0] = 1.e18;
	x[1] = 1000.;
	rc = newt( x, n, &check, vpr, ipr, params, &fmin );
	}
    el_t = (((double) clock()) - el_t) * 1.e-6;

    printf( "\n\n nach Aufruf von newt: \n\n");
    printf( "vpr = ");
    for (i=0; i<n; i++) printf( " %g ", vpr[i]);
    printf( "\n");
    printf( "ipr = ");
    for (i=0; i<n; i++) printf( " %g ", ipr[i]);
    printf( "\n");
    printf( "x = ");
    for (i=0; i<n; i++) printf( " %g ", x[i]);
    printf( "\n");
    printf( "fmin = %g \n\n", fmin);


    printf( "rc = %d, check = %d \n\n", rc, check );
    printf( "\n");

    printf( "Zeitaufwand: %g sec CPU-Zeit\n", el_t);

    return 0;


    }		/* main */


