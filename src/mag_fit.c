/* this version is required for IDL > IDL 5.4 */


#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <time.h>

/* #include </afs/ipp/u/mnw/C/marquardt/mag_fit.h> */

#define INT32		int
#define IDL_SHORT	short int
#define IDL_USHORT	unsigned short int
#define IDL_INT		int
#define IDL_FLOAT	float
#define IDL_DOUBLE	double

/* 2006-01-31: hwm */
/* the structure of the IDL strings changed with IDL 5.5 */
typedef struct 
	{
	IDL_INT 	slen;
	IDL_SHORT	stype;
	char		*s;
	} idl_string_struct ;

#define IDL_STR(x) (((idl_string_struct*)x)->s)
#define IDL_STR_LEN(x) ((int)((idl_string_struct*)x)->slen)


/* physikal. Konstanten */
#define	    c_e    1.6022e-19
#define	    c_eps0 8.8542e-12
#define	    c_me   9.1095e-31
#define	    c_mp   1.6726e-27
#define	    c_pi   3.1415927



/* schnelle Berechnung von ^0.75 durch look-up-table */
#define		n_look_up 		8000
#define		delta_look_up 		(-0.04)
#define		inv_delta_look_up 	(1.0/delta_look_up)
static double	*look_up_alpha=NULL, *look_up_beta=NULL;

/* Grenzen fuer Fit und NL-Solver */
#define		GTOLF	1.e-8
#define		GTOLX	1.e-8
#define		GTOLMIN	1.e-6


/* Macros fuer max, min, etc. */

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


/* dynamische Speicherverwaltung */

#define NR_END 1
#define FREE_ARG char*

/* allocate an integer vector with subscript range v[nl..nh] */
int *ivector(nl, nh)
long	nl;
long	nh;
{
	int *v;       
                
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
                             

/* free an integer vector allocated with ivector() */
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
void upcase(str)
char	*str;
   {
   int	i, len;
   len = strlen(str);
   for (i=0; i<len; i++) str[i]=toupper(str[i]);
   }	/* upcase */
   
   
void lowcase(str)
char	*str;
   {
   int	i, len;
   len = strlen(str);
   for (i=0; i<len; i++) str[i]=tolower(str[i]);
   }	/* lowcase */








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
			==============

	Initialisierung der look-up Tabelle fuer die lineare
	Interpolation bei der Berechnung von x^0.75
 
 	-------------------------------------------------------- 

	M. Weinlich, 19.09.93

 	-------------------------------------------------------- */

void init_quick_075()

    {
    
    double  x=0.0, y1=0.0, y2, temp;
    int     i;

    if ( look_up_alpha != NULL ) {
       /* printf( "init_quick_075: keine Initialisierung mehr noetig! \n"); */
       return;
       }    
    
    /* Stelle noetigen Speicherplatz bereit */   
    look_up_alpha = dvector(0, n_look_up);
    look_up_beta  = dvector(0, n_look_up);

    for (i=0 ; i<n_look_up; i++  ) {

        /* Schichtdicke = const * (phi)^0.75 (Child-Langmuir) 
           x = delta_look_up * (float)(i); 
           y2 = pow( x, 0.75);  */

        /* Schichtdicke fuer heisse monoenerget. Ionen (M. Weinich) */
        /* temp = ( temp + 2. ) * sqrt( temp - 1. ) - y1;*/
        /* 
        x += delta_look_up; 
        temp = sqrt( 1. - x-x );
        y2 = ( temp + 2. ) * sqrt( temp - 1. );
        temp = y2 - y1;
        */
        
        x += delta_look_up+delta_look_up; 
        temp = sqrt(1.-x)-1;
        y2 = (temp+3.) * sqrt(temp);
        temp = y2 - y1;
        
        look_up_alpha[i] = y1 - temp * (double)i;
        look_up_beta[i]  = temp * inv_delta_look_up;

        /* Inkrementieren der Variablen */
        y1 = y2;

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
    double	temp;

    /* keine negativen Argumente erlaubt */
    if ( x >= 0. ) return 0.0;

    /* die Werte sind im Abstand delta_x = delta_look_up tabelliert */
    /* index  ausserhalb des erlaubten Bereichs? */
    index = (int)(x * inv_delta_look_up);
    
    /* lineare Interpolation fuer y */
    if (index < n_look_up-1) {
        return ( look_up_alpha[index] + look_up_beta[index]*x );
        }
    else {
	temp = sqrt( 1. - x-x );
	return ( ( temp + 2. ) * sqrt( temp - 1. ));
        }
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

double schicht_power( phi_de_w, eps_de_sqr, eps_de )
 
    double	*phi_de_w;
    double	*eps_de_sqr, *eps_de;

    {
    double 	temp;
    double      eps_sqr;

    /* keine negativen Argumente unter Wurzel erlaubt: phi_de_w > -1   */
    /* Theorie nur fuer negative Sonde definiert, da im einfachen      */
    /* Bild keine Schicht vor positiver Sonde:         phi_de_w >  0   */
    /*                                                                 */
    /* lineare Extrapolation der Steigung der spannungsabhaengigen     */
    /* Schichtdicke als Steigung des Schichtwachstums der magnetischen */
    /* Schicht                                                         */
    /*                                                                 */
    if ( (*phi_de_w) < 0. ) {
        /* konstanter additiver Beitrag zur Schichtdicke, der bei Berechnung */
        /* fuer phi_de_w > 0 weggelassen wird, da er sich weghebt und sonst  */
        /* nur zu zusaetzlicher Rechenzeit fuehrt */
        temp = (6. - (*eps_de_sqr) ) * ((*eps_de)+(*eps_de));
        return ( 12*(*phi_de_w) / (*eps_de)  + temp );
        }

    
    /* el. Feld phi' = eps */
    temp = 4. - (*eps_de_sqr);
    eps_sqr = 4.*sqrt(1.+(*phi_de_w))-temp;

    /* eps_sqr > 0 und schicht_power >0 */
    if (eps_sqr < 0) {
       printf("\n\neps_sqr < 0, das duerfte nie(!) passieren\n\n\n");
       return 0;
       }
        
    /* nur Differenzen werden gebraucht, d.h. konstante Teile */
    /* koennen weggelassen werden  */
    return ((eps_sqr + temp+temp+temp) * sqrt(eps_sqr));
       	
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

	params(16) 	Verhaeltnis T_i / T_e	a_tau
	
	params(17) 	Hoehe der Sonde senkrecht
			zur Oberflaeche		p_height
			


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
    
    static double	n_e=1., log_ne=0;
    static double	t_e=1., log_te=0;
    static double	inv_beta=1., log_invbeta=0.;
    static double 	alpha1=1., log_alpha1=0.;
    static double	alpha2=1., log_alpha2=0.;
    
    static double	cos1=0., cos1_sqr=0., sin1=1., tan1=1.;
    static double	cos2=0., cos2_sqr=0., sin2=1., tan2=1.;
    
    double 	dv;
    double 	a_height, a_width, p_height, a_mass, a_proj_area, a_sonde_perp;
    double 	a_z, a_b_loc, a_tau;
    double	last_phi, save_exp ;
    double	phi_me_w, phi_me_pr, phi_me_de;
    double	phi_me2_de, phi_diff;
    double 	delta_phi, phi_norm, dw1, dw2; 
    double	zeta1, zeta2, zeta1_sqr, zeta2_sqr, rho_s;
    double 	eps1, eps2, eps1_sqr, eps2_sqr, q;
    double	j_isat, I_isat, log_j_quot, log_inv_exp, 
	     		ns_add1, ns_add2,  c_i_bohm, log_add;
    double 	norm_e1, norm_e2, norm_l1, norm_l2;
    double	sf1, sf2;
    double	ld_sqr;
    double	hl1, hl2, hb1, hb2;
    double 	delta1, delta2, delta1_fac, delta2_fac, lfac, bfac;
    double	temp, tt1, tt2, tt3;

    int 	i, ind;
    int		n_changes=0;
    
    int 	verbose=0;

    /* here we are ... */
    if (verbose) printf("Beginn mag_doppel\n"); /* */

    /* Zuordnung der Fit-Variablen */
    /* Ueberpruefung auf unsinnige Fit-Parameter:
 	0.002  <  t_e            <  400
        1.e7   <  n_e sqrt(t_e)  <  1.e28
        0.001  <  beta           <  148
        0.000  <  alpha          <  100
	0.000  <  cos(psi)       <  1.000
    */


    if (log_te != params[1]) {
       log_te = params[1];
        /* Bereichsueberpruefung */
        if (log_te < -6) {
           params[1] = log_te = -6;
           n_changes++;
           }
        else if (log_te > 6.) {
           params[1] = log_te = 6.;
           n_changes++;
           }
       t_e = exp(log_te);
       }
       

    temp = params[0]-0.5*log_te;
    if (log_ne != temp) {
       log_ne = temp;
        /* Bereichsueberpruefung */
        if (log_ne < 15) {
           log_ne = 15;
           params[0] = log_ne+0.5*log_te;
           n_changes++;
           }
        else if (log_ne > 65.) {
           log_ne = 65.;
           params[0] = log_ne+0.5*log_te;
           n_changes++;
           }
       n_e    = exp(log_ne);
       }
       
    if (log_alpha1 != params[2]) {
       log_alpha1 = params[2];
        /* Bereichsueberpruefung */
        if (log_alpha1 < -50) {
           log_alpha1 = -50;
           params[2] = log_alpha1;
           n_changes++;
           }
        else if (log_alpha1 > 5.) {
           log_alpha1 = 5.;
           params[2] = log_alpha1;
           n_changes++;
           }
       alpha1 = exp(log_alpha1);
       }
       
    if (log_invbeta != -params[3]) {
       /* beta     = exp(params[3]); */
       log_invbeta = -params[3];
        /* Bereichsueberpruefung */
        if (log_invbeta < -5) {
           log_invbeta = -5;
           params[3] = -log_invbeta;
           n_changes++;
           }
        else if (log_alpha1 > 5.) {
           log_invbeta = 5.;
           params[3] = -log_invbeta;
           n_changes++;
           }
       inv_beta    = exp(log_invbeta);
       }
       
    dv = params[4];   
    
    if (log_alpha2 != params[5]) {
        log_alpha2 = params[5];
        /* Bereichsueberpruefung */
        if (log_alpha2 < -50) {
           params[5] = log_alpha2 = -50;
           n_changes++;
           }
        else if (log_alpha2 > 5.) {
           params[5] = log_alpha2 = 5.;
           n_changes++;
           }
        alpha2 = exp(log_alpha2);
        }
       
    if (cos1 != params[6]) {
        cos1 = params[6];
        /* Bereichsueberpruefung */
        if (cos1 < 1.e-5) {
           params[6] = cos1 = 1.e-5;
           n_changes++;
           }
        else if (cos1 > 1.) {
           params[6] = cos1 = 1.;
           n_changes++;
           }
        /* abgeleitete Groessen */
        cos1_sqr = cos1*cos1;
        sin1     = sqrt(1. - cos1_sqr);
        tan1     = sin1/cos1;
        }
        
    if (cos2 != params[7]) {    
        cos2 = params[7];
        /* Bereichsueberpruefung */
        if (cos2 < 1.e-6) {
           params[7] = cos2 = 1.e-6;
           n_changes++;
           }
        else if (cos2 > 1.) {
           params[7] = cos2 = 1.;
           n_changes++;
           }
        /* abgeleitete Groessen */
        cos2_sqr = cos2*cos2;
        sin2     = sqrt(1. - cos2_sqr);
        tan2     = sin2/cos2;
        }

    
    /* Abbruch, da keine sinnvollen Input-Parameter? */    
    if (verbose) {
       if ( n_changes == 1 ) {
	  printf( "Achtung: insgesamt %d Parameter musste geaendert werden!\n",
		n_changes);
	  }
       else if ( n_changes > 1 ) {
	  printf( "Achtung: insgesamt %d Parameter mussten geaendert werden!\n",
		n_changes);
	  }
       }

    if ( n_changes > 4 ) {
        if (verbose) {
	   printf( "Abbruch, da insgesamt %d Parameter ", n_changes);
	   printf( "geaendert werden mussten  !\n");
	   }
	return (-21); 
	}


    /* ab hier kann einfach ausgelesen werden, keine Umwandlung mehr noetig */    
        
    /* wichtige Sonden-Dimensionen */
    a_width  = params[10]; 
    a_height = params[11]; 
    p_height = params[17];

    /* projezierte Flaeche einer flush mounted Sonde */
    a_proj_area  = a_width*a_height*cos1;	/* <= 0, wenn nicht definiert */

    /* projezierte Flaeche einer pin probe */
    a_sonde_perp = a_width*p_height*sin1;	/* <= 0, wenn nicht definiert */
    
    /* Massenzahl des Fuellgases, meist Deuterium */
    if ( params[12] < 1. ) a_mass = params[12] = 2. ;
    		else a_mass = params[12];
    
    /* Kernladungszahl, Fuellgas meist Wassertoff-Isotop */
    if ( params[13] < 1.) a_z = params[13]=1.;
    		else a_z  = params[13];

    
    /* Betrag des lokalen Magnetfelds am Ort der Sonde */
    /* Standard: |B_t| = 2.0T */
    if ( fabs(params[14]) < 0.6 ) a_b_loc = params[14] = 2.;
   		else a_b_loc  = fabs(params[14]); 

    /* T_i / T_e,  Ionen- und Elektronentemperatur sind postiv */
    if ( params[16] < 0.05 ) a_tau = params[16] = 1.;
    		else a_tau = params[16];

     



    /* Berechnung der Ionensaettigungsstromdichte 
       ==========================================

       Sekundaerelektronenemissionskoeffizient g = 0.6

       c_i_Bohm = sqrt(  ( Z*T_e + 3 T_i ) / m_i )  wobei 3 = (n+2)/n, n=1
       j_isat = n_e * c_i_Bohm * e

       --> im Doppelsondenbild kommen keine beschleunigten Elektronen zur
	   Elektrode, sondern immer nur thermische Plasmaelektronen.   

       c_e_th = sqrt( (8 T_e) / (pi m_e) )
       j_e_th = 1/4 * e * n_e * c_e_th * (1-g)

       j_quot = j_e_th / j_i_sat =		 mit T_i = T_e
              = 1/4 (1-g) sqrt(8/pi) sqrt(m_p/m_e) sqrt( a_mass/(a_z+3) ) =
	      = 6.8376426 * sqrt( a_mass/(a_z+3) )

    */

    c_i_bohm = sqrt( (3.*a_tau+a_z)* c_e * t_e/ a_mass /c_mp );
    j_isat = n_e * c_e * c_i_bohm;

    /* log_j_quot = log(j_quot); */
    log_j_quot = 1.922443023 + 0.5*log(a_mass/(a_z+3.*a_tau));

    /* Debye-Laenge vor Schichten */
    ld_sqr = t_e * c_eps0 / n_e / c_e ;


    if ( a_proj_area > 0. ) {
    
        I_isat = a_proj_area * j_isat;
        
        /* Einfluss der Schichten --> Reduzierung in l_debye u.ae. */
        /* allgemeinerer Ausdruck fuer zeta = n_i,De / n_i,me */
        q = a_z/(6.*a_tau);
        zeta1_sqr = sqrt(q*q + (q+q+1)*cos1_sqr ) - q;
        zeta2_sqr = sqrt(q*q + (q+q+1)*cos2_sqr ) - q;
        zeta1 = sqrt(zeta1_sqr);
        zeta2 = sqrt(zeta2_sqr);
        
        /* Normierung des el. Feldes:                                            */
        /* alle Potentialgroessen wie phi_me_w sind in [eV/k_B T_e] ausgedrueckt */
        /* damit sie mit der Normierung aus dem paper uebereinstimmen, muessen   */
        /* sie noch durch eine Faktor (s.u.) geteilt werden                      */
        norm_e1 = 0.5 + zeta1_sqr/(4.*q) ;
        norm_e2 = 0.5 + zeta2_sqr/(4.*q) ;
        
        /* Normierung der Laenge                                                 */
        /* alle Laengen werden auf Debyelaenge an der mag. Schicht multipliziert */
        /* mit einer dimensionslosen Groesse normiert. Hier enthaelt die         */
        /* Normierung beides, die Debye-Laenge und die Skalierung, da alle       */
        /* Laengen in [m] ausgedrueckt werden                                    */
        norm_l1 = sqrt( ld_sqr * norm_e1 / zeta1 ) ;
        norm_l2 = sqrt( ld_sqr * norm_e2 / zeta2 ) ;
                
        /* Potentialabfall in magnetischer Schicht, in eV/k_B T_e */
        phi_me_de  = -log(zeta1);
        phi_me2_de = -log(zeta2);
    
        /* Ionengyroradius bei Schallgeschwindigkeit */
        rho_s = c_i_bohm * a_mass * c_mp  /  c_e / a_b_loc  ;
    
        /* Abschaetzung fuer das elektrische Feld am Eintritt in Debye-Schicht */
	/* muss immer positiv sein                                             */
        eps1 = DMAX( (phi_me_de / norm_e1) * (norm_l1 / rho_s), 0.) ;
        eps2 = DMAX( (phi_me2_de / norm_e2) * (norm_l2 / rho_s), 0.) ;

/* printf("old: eps1=%g, eps2=%g \n", eps1, eps2); */
   
        /* neue Abschaetzung, awc 02.01.97 */
/*        temp = -4*(1-cos1)*(1-cos1) + 2.*sin1*sin1*phi_me_de/norm_e1;
        eps1 = pow( sqrt(temp)*norm_l1/rho_s, 1./3. );
        
        temp = -4*(1-cos2)*(1-cos2) + 2.*sin2*sin2*phi_me2_de/norm_e2;
        eps2 = pow( sqrt(temp)*norm_l2/rho_s, 1./3. );
        
printf("new: eps1=%g, eps2=%g \n", eps1, eps2);
*/        		
        eps1_sqr = eps1*eps1;
        eps2_sqr = eps2*eps2;
        
        /* Potentialdifferenz zwischen den beiden Referenz-Plasmapotentialen */
        /* phi_me2 = phi_me1 - delta_phi */
        delta_phi = dv / t_e;
           
        /* Berechne phi_me - phi_wall fuer V=0 */
        log_add    = log_j_quot - log(1.+inv_beta);
        phi_me_w   = log_add + log( inv_beta+exp(delta_phi) )  ;
        params[15] = phi_me_w;      /* Abspeichern fuer spaetere Verwendung */

         
        /* Schichtdicke an der umgebenden Wand, d.h. bei V=0 */
        phi_norm = (phi_me_w - phi_me_de) / norm_e1;
        dw1 = schicht_power( &phi_norm, &eps1_sqr, &eps1 );
        phi_norm = (phi_me_w - delta_phi - phi_me2_de) / norm_e2;
        dw2 = schicht_power( &phi_norm, &eps2_sqr, &eps2 );
    
        /* Vorfaktoren der Schichtdicken, incl. Normierung der Laenge */
        sf1 = norm_l1 / 12. * alpha1 ;
        sf2 = norm_l2 / 12. * alpha2 ;
    
        /* Hilfsgroessen, muessen nur einmal berechnet werden */
        temp = sqrt(cos2*inv_beta/cos1);
        lfac = temp / a_height;
        bfac = temp / a_width;
        
        delta1_fac = tan1 / a_height;
        delta2_fac = tan2 * lfac;

/*            
if (cos1 < 0.017) {
   printf("\n\n");
   printf("cos1=%g, phi_me_de=%g, phi_me_w=%g \n", 
   		cos1, phi_me_de, phi_me_w );
   printf("eps1=%g, eps1_sqr=%g, rho_s=%g, norm_e1=%g, norm_l1=%g \n", 
		eps1, eps1_sqr, rho_s, norm_e1, norm_l1);
   printf("te=%g, ne=%g, beta=%g, delta_phi=%g, alpha1=%g \n", 
   		t_e, n_e, 1./inv_beta, delta_phi, alpha1);
   printf("d1_fac=%g, d2_fac=%g, sf1=%g, sf2=%g \n", 
   		delta1_fac, delta2_fac, sf1, sf2);
   phi_norm = (phi_me_w - phi_me_de) / norm_e1;
   printf("phi_norm=%g, dw1=%g \n",phi_norm, dw1);
   }
*/
         
 tt1 = 0.;
 tt2 = 0.;
 tt3 = 0.;
 
        /* ab jetzt Berechnung fuer alle Spannungswerte */
        /* phi_me passt sich jeweils so der von aussen angelegten Sondenspannung */
        /* an, dass die Doppelsondengleichung erfuellt ist, d.h. phi_me ist      */
        /* ist eine implizite Funktion von phi_sonde-phi_wall                    */
        for ( ind=0; ind < n_vpr; ind++ ) 
            {
    
            /* Spannung an Doppelsonde: phi_sonde */
            /* phi_pr_w = vpr[ind] / t_e;           */
            /* phi_diff = delta_phi - phi_pr_w;     */
            phi_diff = (dv - vpr[ind]) / t_e;
            
            /* Startwert fuer Potentialdifferenz zwischen mag. Schicht und Sonde */
            /* stimmt exakt fuer saettigende Doppelsonde, nur bei sehr hohem     */
            /* Nichtsaettigungsanteil ist noch eine weitere Iteration noetig     */
            if ( phi_diff < -300. ) {
                save_exp    = 0.;
                log_inv_exp = log_invbeta;  	/* log_invbeta = log(inv_beta)*/
                }
            else if ( phi_diff > 300. ) {
                save_exp    = 1.95e130;               /* = exp(300) */
                log_inv_exp = phi_diff;
                }
            else { 
                save_exp    = exp(phi_diff);
                log_inv_exp = log(inv_beta + save_exp);
                }
            phi_me_pr = log_inv_exp + log_add;
                
       
            /*      ------ eigentliche Iteration ------     */
    
/* keine Iteration mehr: phi_me_pr wird lediglich fuer Berechnung der Nicht- */
/* saettung benoetigt. Diese wird als eine Korrektur aufgefasst, so dass     */
/* die Abschaetzung voellig genuegt.                                         */
/* Test an Kennlinie mit psi=89.4 und vpr=-200,200 zeigt, dass die Iteration */
/* nur Aenderungen (in T_e und n_e) um max. 3% bringt. Dafuer lohnt sich der */
/* Aufwand an Rechenzeit aber nicht.                                         */

            i = 0;
/*            do {
*/    
                /* Schleifenzaehler */
                last_phi = phi_me_pr;      /*2006-11-17 hwm:  only for iteration */
                i++;
    
                /* Berechne Schichdicken und deren Differenz zur Umgebung */
                phi_norm = (phi_me_pr - phi_me_de) / norm_e1;
                delta1   = sf1*( schicht_power(&phi_norm,&eps1_sqr,&eps1) - dw1);
                phi_norm = (phi_me_pr - phi_diff - phi_me2_de) / norm_e2;
                delta2   = sf2*( schicht_power(&phi_norm,&eps2_sqr,&eps2) - dw2);
    
    /*
                hl1 = delta1 / a_height;
                hb1 = delta1 / a_width; 
                hl2 = delta2 * lfac;
                hb2 = delta2 * bfac; 
    */          
                /* keine Expansion der Schicht parallel zur Oberflaeche */
                ns_add1 = 1. + delta1 * delta1_fac;
                ns_add2 = 1. + delta2 * delta2_fac; 
                
                /* Expansion entlang Oberflaeche = Expansion senkrecht */
    /*          ns_add1 = 1. + hb1+hb1 + (hb1+hb1+1.)*hl1*(tan1+2) ;
                ns_add2 = 1. + hb2+hb2 + (hb2+hb2+1.)*hl2*(tan2+2) ; 
    */
                /* Expansion entlang Oberflaeche = 0.5*Expansion senkrecht
                   diese Version wird von PIC-Rechnungen favorisiert (axb) */
    /*          ns_add1 = (1. + hb1)*(1+ hl1*(1+tan1)) ;
                ns_add2 = (1. + hb2)*(1+ hl2*(1+tan2)) ; 
    */          
                /* Nichtsaettigung(Sondenflaeche) kann nicht kleiner NULL  */
                /* werden. Sondenflaeche kann auch nie gleich Null werden, */
                /* da ein paar Ionen immer den Weg zur Sonde finden werden */
                ns_add1 = DMAX(ns_add1,0.);
                ns_add2 = DMAX(ns_add2,0.);
    
                /* naechste Naeherung fuer phi1 */
                /* nur fuer starke Nichtsaettigung, sonst genuegt es, wenn in */
                /* Pseudo-Schleife ns_add1,2 berechnet wurden */
/*                if (ns_add1 > 1.3)  phi_me_pr = log_inv_exp + log_j_quot
                    			       - log(inv_beta*ns_add1+ns_add2);
*/
                   
/*                phi_me_pr = log_inv_exp + log_j_quot
                    			       - log(inv_beta*ns_add1+ns_add2);
*/    
 
/* keine Iteration mehr: phi_me_pr wird lediglich fuer Berechnung der Nicht- */
/* saettung benoetigt. Diese wird als eine Korrektru aufgefasst, so dass     */
/* die Abschaetzung voellig genuegt.                                         */
/* Test an Kennlinie mit psi=89.4 und vpr=-200,200 zeigt, dass die Iteration */
/* nur Aenderungen (in T_e und n_e) um max. 3% bringt. Dafuer lohnt sich der */
/* Aufwand an Rechenzeit aber nicht.                                         */
              
/*                } while( ( fabs(phi_me_pr/last_phi-1.) > .05 ) && (i < 10) ) ;
*/   
            /* Berechnung des Strom-Vektors, ns_add1,2 wird aus letzter */
            /* Iteration uebernommen, 5% genauigkeit in phi_pr_me langt */
            strom[ind] = I_isat * 
                           (ns_add2 - ns_add1*save_exp) / (inv_beta + save_exp);
                           
   
/*
if ((cos1 < 0.017) && (ind == 0)){
   phi_norm = (phi_me_pr - phi_me_de) / norm_e1;
   printf("vpr=%g, phi=%g, dp1=%g, delta1=%g, delta2=%g, ns1=%g, ns2=%g\n", 
   		vpr[ind], phi_norm,
   		schicht_power(&phi_norm,&eps1_sqr,&eps1),
   		delta1, delta2, ns_add1, ns_add2);
   }
*/
   
    /*
    if ( (ns_add1 < 1.) && (ns_add2 < 1.)) {
       printf("ns < 1: vpr=%g, strom=%g\n", vpr[ind], strom[ind]);
       printf("delta1=%g, delta2=%g, ns1=%g, ns2=%g\n",
       	delta1, delta2, ns_add1, ns_add2);
       printf("eps1=%g, eps2=%g \n", eps1, eps2);
       }
   */    
/*           if ( (fabs(strom[ind]) < 0.001) && (vpr[ind] > 2*t_e) ) {
               printf("vpr=%g, strom=%g\n", vpr[ind], strom[ind]);
               }
*/
            }       /* for ind=0, n_vpr-1 */
    
        } /* if a_proj_area > 0 */
        
        
    /* ---------------------------------------------------------------------
       wenn a_sonde_perp > 0.0 dann wurde bereits eine effektive Sondenfl"ache
       senkrecht zur Feldrichtung eingegeben. 

       original Kommentar:
       Es wird davon ausgegangen, da"s
       es sich hierbei um eine freistehende Sonde handelt, die von beiden
       Seiten Ionen und Elektronen aufsammeln kann. 
       2012-01-17: hwm
                   Diese Annahme geht in der C routine nicht mehr ein.
                   Die Flaeche ist Hoehe*Breite*sin(Winkel).

       Falls a_sonde_perp > 0.0 und a_sonde_proj > 0.0 dann wird eine 
       Kombination beider Effekte erwartet, d.h. der Gesamtstrom zur Sonde
       setzt sich additiv aus beiden Teilen zusammen.
       --------------------------------------------------------------------- */

    if (a_sonde_perp > 0.0)
	{

        /* Initialisierung, falls n"otig */
        if (a_proj_area <= 0.) {
            for ( ind=0; ind<n_vpr; strom[ind++]=0.0); 
	    }

        /* Kennlinie einer Doppelsonde mit Potentialverschiebung und
	   Fl"achenverh"altnis beta */
        I_isat = a_sonde_perp * j_isat;
        bfac = alpha2 * sqrt(ld_sqr) / a_width; 


        for ( ind=0; ind < n_vpr; ind++ ) 
	    {

            /* Spannung an Doppelsonde: phi_sonde */
            delta_phi = (dv - vpr[ind]) / t_e;

            /* Nichts"attigungsverhalten im Elektronenast? */
            /* Berechnung des Strom-Vektors */
            /* *ptr_strom++ +=  I_isat * 
		(1. + bfac * quick_075(delta_phi) - temp) / (inv_beta+temp); */

	    if ( delta_phi < -300. ) {
                strom[ind] +=  I_isat * (1.+bfac*quick_075(delta_phi))
                                   / inv_beta  ;
		}
	    else if ( delta_phi > 300. ) {
                strom[ind] +=  - I_isat ;
	        }
	    else { 
	        temp = exp( delta_phi );
                strom[ind] +=  I_isat * 
		       (1.+ bfac*quick_075(delta_phi) - temp) / (inv_beta+temp);
		}
	    }


	}	/* if (a_sonde_perp > 0.0) .....    */

    return 0 ;

    }		 /* mag_doppel */






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
    static int 		i, j;
    static double 	d_par, temp;
    static double	*y2=NULL;
    static int		n_y2=0;


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
    /* Gebe Speicher frei, wenn Aenderung in n_x  */
    if ( n_x > n_y2 ) {
	/* printf("f_deriv: Aendere Speicherbereich: y2 = %X, Laenge = %d\n",
			y2, n_x+n_x); */
	if ( n_y2 > 0 ) free_dvector(y2,0,n_y2-1);
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

	    /* die folgenden Zeilen dienen dazu, die Auswirkungen
		   von Rundungsfehlern zu minimieren */
	    temp=pars[j];
	    d_par=1.e-6*fabs(temp);
	    if (d_par == 0.0) d_par=1.e-9;
	    pars[j] = temp+d_par;
	    d_par   = pars[j]-temp;
		
	    if ( (*func)(x, y2, n_x, pars) < 0) return -111;
	    pars[j]  = temp; 
	    		
	    /* partielle Ableitung an jedem Punkt */
	    d_par = 1./d_par;
	    for( i=0; i< n_x; i++) deriv[j][i] = (y2[i] - yfit[i]) * d_par;	
	
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
	double anorm,c,f,g,h,s,scale,x,y,z;

	static double	*rv1=NULL;
	static int	n_rv1=0;

	/* rv1=dvector(1,n); */
        if ( n > n_rv1 ) {
	    if (n_rv1 > 0) free_dvector(rv1,0,n_rv1-1);
	    rv1=dvector(0,n-1);
	    n_rv1 = n;
	    }

	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i < m && i != (n-1) ) {
			for (k=l;k<n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for ( i=(n-1);i>=0;i--) {
		if (i < (n-1)) {
			if (g) {
				for (j=l;j<n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
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
					for (j=0;j<m;j++) {
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
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}

			/* substantieller Fehler, SVD konvergiert nicht */
			if (its == 30) {
			   /*
			   fprintf(stderr,"dsvcmp: ");
			   fprintf(stderr,"no convergence in 30 iterations");
			   fprintf(stderr," ... sorry!\n");
			   fprintf(stderr,"\t--> returning without action\n");
			   */
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
				for (jj=0;jj<n;jj++) {
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
				for (jj=0;jj<m;jj++) {
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


void dsvbksb(u, w, v, m, n, b)
double 	**u;
double 	w[];
double	**v;
int	m,n;
double	b[];

	{
	int 	jj,j,i;
	double 	s;

	static double 	*tmp=NULL;
	static int	n_tmp=0;


	if ( n > n_tmp ) {
	    if (n_tmp > 0) free_dvector(tmp,0,n_tmp-1);
	    tmp=dvector(0,n+n-1);
	    n_tmp = n+n;
	    }


	for (j=0;j<n;j++) {
		if (w[j]) {
			for (s=0.,i=0;i<m;i++) s += u[i][j]*b[i];
			tmp[j] = s*w[j];
			}
		else {
			tmp[j] = 0.;
			}
		}

	for (j=0;j<n;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
		b[j]=s;
		}

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


	for (i=0; i<n; i++) {
	    for ( j=0; j<n; j++ ) {
		for (temp=0.,jj=0; jj<n; jj++) temp += v[i][jj]*u[j][jj]*w[jj];
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
			Dimension n
			calc_inverse	Inverse wird nur berechnet, wenn
					"true" (=0)					
								
	Rueckgabe:	A^(-1) in Matrix A 			
			x in Vektor b				
								
	Funktionsreturn:	0 	alles o.k.		
				<0	Fehler in svd		
								
	--> um Zeit zu sparen, werden die internen Felder nur 
	    einmal angelegt und bei folgenden Aufrufen wieder-
	    verwendet. Ist in einem folgenden Aufruf die Dimension
	    groesser, so werden die Felder automatisch angepasst.
	    Sollen die Felder vor Programm-Ende geloescht werden,
	    so muss beim Funktionsaufruf eine negative Dimension
	    n angegeben werden.
	    
	--> fuer kleine Systeme (2x2 und 3x3) ist die Loesung
	    explizit fest vorgegeben. Ist die Determinante nahe
	    Null, so wird die Loesung an die eigentliche svd-Routine
	    weitergegeben.						
		
	-------------------------------------------------	
								
	M. Weinlich, 29.08.95					
	entnommen aus Numerical Recipies			
								
	-------------------------------------------------	*/


int dsvd_solver(a, b, n, calc_inverse)
double 	**a;
double 	*b;
int	n;
int	calc_inverse;
{
	int 	rc;
	int	i, j;

        double  det, a00, a01, a02, a10, a11, a12, a20, a21, a22, b0, b1, b2;

	double	wmin, wmax, temp;


	/* Hilfsmatrizen fuer svd */
	static int	last_n=0;
	static double	**u=NULL;
	static double	**v=NULL;
	static double	*w=NULL;

   
	/* Sonderbehandlung: 2x2 und 3x3 Matrizen koennen noch
	                     von Hand invertiert werden! */
        if ( n == 2 ) {
	    /* printf("\n\n --> Sonderbehandlung n=2 in dsvd_solver \n\n"); */
	    a00 = a[0][0];
	    a01 = a[0][1];
	    a10 = a[1][0];
	    a11 = a[1][1];
	    det = a00*a11 - a01*a10;
	    /* printf("Sonderbehandlung n=2 in dsvd_solver det = %g\n",det); */
	    if ( fabs(det) > 1.e-60 ) {
		b0 = b[0];
		b1 = b[1];
		/* inverse Matrix */
		det = 1. / det;
		a[0][0] = a11 * det;
		a[0][1] *= (-det);
		a[1][0] *= (-det);
		a[1][1] = a00 * det;
		/* Loesungsvektor */
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
	    if ( fabs(det) > 1.e-60 ) {
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


	/* Aenderung in geforderter Dimension : Loesche interne Variablen */
	if ( n > last_n ) {
	    /* printf("dsvd_solver: Lege Speicherplatz an, n = %d\n", n); */
	    if (last_n > 0) {
	        /* printf("dsvd_solver: Aendere Speicherplatz, n = %d --> n = %d\n", 
			last_n, n); */
	        free_dmatrix( u, 0, last_n-1, 0, last_n-1);
	        free_dmatrix( v, 0, last_n-1, 0, last_n-1);
	        free_dvector( w, 0, last_n-1);
	        }
	    u  = dmatrix( 0, n-1, 0, n-1);
	    v  = dmatrix( 0, n-1, 0, n-1);
	    w  = dvector( 0, n-1);
	    last_n = n;
	    }

	/* Eingabe von n < 0 : Loesche interne Variablen */
	if ( n <= 0 ) {
	    if ( last_n >= 0) return 0;
	    /* wenn wirklich Speicher belegt wurde */
	    free_dvector( w , 0, last_n-1);
	    printf("1x free_dvector\n");
	    printf("dsvd_solver: Loesche Speicherplatz, n = %d\n", n);
	    free_dmatrix( u , 0, last_n-1, 0, last_n-1);
	    printf("1x free_dmatrix\n");
	    free_dmatrix( v , 0, last_n-1, 0, last_n-1);
	    printf("2x free_dmatrix\n");

	    u  = NULL;
	    v  = NULL;
	    w  = NULL;
	    last_n = 0;
	    return 0;
	    }

/*
if (n < 4) {
   printf("keine Standard-Loesung fuer n=%d gefunden, det = %g\n",n,det);
   }
*/

	/* singular value decomposition von A */
	rc = dsvdcmp( a, n, n, w, v);
	if ( rc < 0 ) return (rc);


	/* was haben wir denn fuer Ergebnisse erhalten */
	wmax = 0.0;
	for (i=0; i<n; i++) wmax = DMAX( wmax, w[i]);
	wmin = wmax * 1.e-12;
	
	/* eigentlich wird nicht w (>=0.0) sondern 1/w in den folgenden
	   Routinen benoetigt, dabei muss jedoch 1/0 = 0 gesetzt werden */
	for (i=0; i<n; i++) {
	   temp = DMAX( wmin, w[i]);
	   w[i] = ( (temp > 0.) ? 1./temp : 0.0);
	   }
	   
	/* Loesung des Gleichungssystems */
	dsvbksb( a, w, v, n, n, b); 
        
        /* weiter nur, wenn auch Inverse berechnet werden muss */
        if (calc_inverse != 0) return (0);
        
        printf( "Berechne auch Inverse\n");
        
	/* kopiere a nach u */
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) u[i][j] = a[i][j];
	    }

        /* Berechnung der inversen Matrix */
	dsvdinv( u, w, v, n, a );


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
int	ma, ia[];
int	mfit;

{
	int i,j,k;
	double swap;
	
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


	*chisq=0.0;
	for (i=0;i<ndata;i++) 
		*chisq += DSQR( (y[i]-yfit[i])/sig[i] );
	  

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
	double		ochisq, dy, wt;

	static double 	*atry=NULL,
			*beta=NULL,
			*da=NULL,
			*oneda=NULL;
	static double 	**dyda=NULL;
	static int	last_ma=0, last_ndata=0;



        /* printf("Beginn mrqmin\n"); */
        
	/* stelle benoetigte Datenfelder bereit */
	if (( ma > last_ma ) || ( ndata > last_ndata )) {
                /* printf( "mrqmin: Aendere Speicher, ma = %d --> %d\n",
					last_ma,ma); */
                /* printf( "                          nd = %d --> %d\n",
					last_ndata,ndata); */
		if (last_ma > 0) {
		    free_dvector(oneda,0,last_ma-1);
		    free_dvector(da,0,last_ma-1);
		    free_dvector(beta,0,last_ma-1);
		    free_dvector(atry,0,last_ma-1);
	            free_dmatrix(dyda,0,last_ma-1,0,last_ndata-1);
	            }
		atry  = dvector(0,ma-1);
		beta  = dvector(0,ma-1);
		da    = dvector(0,ma-1);
	        dyda  = dmatrix(0,ma-1,0,ndata+ndata-1);
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

        /* gar kein Fit verlangt? */
        if (mfit == 0) return (0);

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
	        covar[j][j] *= (1.0+(*alamda)); 
	        oneda[j] = beta[j];
	        }		/* for j=0, mfit-1 */

	    /* Loesung der Gleichung covar * da = oneda */
	    rc = dsvd_solver(covar,oneda,mfit, 1);

	    /* Fehler in svd aufgetreten? wenn ja, dann hat weitere Berechnung
		keinen Sinn mehr --> geordneter Rueckzug aus mrqmin, immer
		moeglich wenn alamda = 0 */
	    if ( rc < 0 ) return rc;

	    /* probiere neue Parameter aus */
	    for (j=0,l=0;l<ma;l++) 
	        atry[l] =  (ia[l])?  a[l]+oneda[j++] : a[l];
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
		chi_sqr, iter_max, eps_abs, eps_rel, verbose)

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
    int		*verbose;

    {

    /* Local variables */
    static double	n_e=1., log_ne=0.;
    static double	t_e=1.,	log_te=0.;
    static double	beta=1, log_beta=0.;
    static double	alpha1=1., log_alpha1=0.;
    static double	alpha2=1., log_alpha2=0.;
    static double	psi1=0., cos1=1.;
    static double	psi2=0., cos2=1.;
    
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
    static int		last_npars=0;

    double	chi_est, mean, min;
    int		pos;

    static double	**sigma=NULL;
    static int		old_n_sigma=-1;


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


    /* Stelle benoetigte Hilfsfelder bereit, wenn noch kein Speicher allokiert,
       bzw. passe Speicher bei Aenderung der Anforderungen an */
    if ( *n_pars > last_npars ) {
	if ( last_npars > 0) {
            free_dmatrix( alpha, 0, last_npars-1, 0, last_npars-1);
            free_dvector( curr_nb, 0, last_npars-1);
            free_ivector( save_do_var, 0, last_npars-1);
            free_dmatrix( sigma, 0, last_npars-1, 0, last_npars-1);
            }
        alpha       = dmatrix( 0, *n_pars-1, 0, *n_pars-1);
        curr_nb     = dvector(0, *n_pars-1);
        save_do_var = ivector(0, *n_pars-1);
        sigma 	    = dmatrix(0, *n_pars-1, 0, *n_pars-1);
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

/*
    if ( pars[0] <= 0.0 ) printf( "Fehler in pars[0] %g \n", pars[0]);
    if ( pars[1] <= 0.0 ) printf( "Fehler in pars[1] %g \n", pars[1]);
    if ( pars[2] <= 0.0 ) printf( "Fehler in pars[2] %g \n", pars[2]);
    if ( pars[3] <= 0.0 ) printf( "Fehler in pars[3] %g \n", pars[3]);
    if ( pars[5] <= 0.0 ) printf( "Fehler in pars[5] %g \n", pars[4]);
*/
    /* in allen positiven Parametern logarithmischer Fit */
    if ( pars[0] <= 0.0 ) pars[0]=1.e-100 ;
    if ( pars[1] <= 0.0 ) pars[1]=1.e-100 ;
    if ( pars[2] <= 0.0 ) pars[2]=1.e-100 ;
    if ( pars[3] <= 0.0 ) pars[3]=1.e-100 ;
    if ( pars[5] <= 0.0 ) pars[5]=1.e-100 ;

    if (n_e != pars[0]) {
       n_e = pars[0];
       log_ne = log(n_e);
       }
    if (t_e != pars[1]) {
       t_e = pars[1];
       log_te = log(t_e);
       }
    if (beta != pars[2]) {
       beta = pars[2];
       log_beta = log(beta);
       }
    if (alpha1 != pars[3]) {
       alpha1 = pars[3];
       log_alpha1 = log(alpha1);
       }
    if (alpha2 != pars[5]) {
       alpha2 = pars[5];
       log_alpha2 = log(alpha2);
       }
    if (psi1 != pars[6]) {
       psi1 = pars[6];
       cos1 = cos(psi1);
       }
    if (psi2 != pars[7]) {
       psi2 = pars[7];
       cos2 = cos(psi2);
       }    
       
    pars[0] = log_ne + 0.5*log_te;
    pars[1] = log_te;
    pars[2] = log_beta;
    pars[3] = log_alpha1;
    pars[5] = log_alpha2;
    pars[6] = cos1;
    pars[7] = cos2;


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

	/* gegebenenfalls muessen die Nebenbedingungen angepasst werden */
    	if (nb_const[0] > 0.0) nb_values[0] = log(nb_values[0]);
    	if (nb_const[1] > 0.0) nb_values[1] = log(nb_values[1]);
    	if (nb_const[2] > 0.0) nb_values[2] = log(nb_values[2]);
    	if (nb_const[3] > 0.0) nb_values[3] = log(nb_values[3]);
    	if (nb_const[5] > 0.0) nb_values[5] = log(nb_values[5]);
 
        /* Fit in cos(psi) statt in psi */
        if (nb_const[6] > 0.0) nb_values[6] = cos(nb_values[6]);
        if (nb_const[7] > 0.0) nb_values[7] = cos(nb_values[7]);

    	
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
	if (*verbose ) printf("Fit mit Nebenbedingungen, Schaetzwert fuer chi^2 = %g\n",
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

        /* kein Fit verlangt? */
        if (mfit == 0) break;

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
    
    /* "ich lebe noch Zeile" muss nur abgeschlossen werden, wenn sie auch */
    /* wirklich angelegt wurde */
    /* if ( iter >= 10 ) fprintf ( stderr, "\n" ); */

    /* Fit ok, welches Chi? nur interessant fuer Nebenbedingungen, da in
       diesem Fall zu Beginn ein Chi abgeschaetzt werden musste */
    if ( (best_estimate) && (do_nb) && (*verbose) ){
	fprintf ( stderr, "Fit nach Levenberg-Marquardt erfolgreich beendet " );
	fprintf ( stderr, "( chi^2 = %g ) \n", *chi_sqr );
	}

    /* sind Fehler aufgetreten, wenn ja, welche? */
    if ( *verbose ) {
        if (rc == -21 ) 
	    fprintf(stderr,
			"Zuviele Parameter ausserhalb ihres Wertebereichs\n");
        else if (rc == -42) 
	    fprintf(stderr,"Iteration abgebrochen --- Fehler in svd\n");
        else if (bad_estimate) 
	    fprintf(stderr,"Iteration abgebrochen --- Signal 'bad estimate'\n");
        else if (too_much) 
	    fprintf(stderr,"Iteration nicht konvergiert\n");
        else if (rc != 0 ) 
	    fprintf(stderr, "Fehler waehrend Iteration, Fehlercode = %d\n",rc);
	}

    if (bad_estimate) rc = -33;
    if (too_much) rc = -34;

    /* logarithmischer Fit in allen positiven Parametern */
    pars[0] = exp(pars[0]-0.5*pars[1]);
    if (pars[1] == log_te)     pars[1]=t_e;    else pars[1] = exp(pars[1]);
    if (pars[2] == log_beta)   pars[2]=beta;   else pars[2] = exp(pars[2]);
    if (pars[3] == log_alpha1) pars[3]=alpha1; else pars[3] = exp(pars[3]);
    if (pars[5] == log_alpha2) pars[5]=alpha2; else pars[5] = exp(pars[5]);
    if (pars[6] == cos1)       pars[6]=psi1;   else pars[6] = acos(pars[6]);
    if (pars[7] == cos2)       pars[7]=psi2;   else pars[7] = acos(pars[7]);


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
 	
 	
 	
 	
 	



double calc_fmin(x, n, vecfunc, vpr, ipr, icurr, params )

    double	*x;
    int		n;
    int 	(*vecfunc)();
    double	*vpr, *ipr, *icurr;
    double	*params;

    	{
	int 	i;
	double 	sum;
	
	params[0] = x[0];
	params[1] = x[1];
	if (n == 3) params[3]=x[2];
	
	(vecfunc)(vpr,icurr, n, params);
	for (sum=0.0,i=0;i<n;i++) sum += DSQR(icurr[i]-ipr[i]); 
	return (0.5*sum);
	
        } 	/* calc_fmin */



int fdjac(n, x, icurr, df, vecfunc, vpr, params )

    int		n;
    double	*x, *icurr, **df;
    int 	(*vecfunc)();
    double	*vpr, *params;


    {
	int 	i,j;
	double 	h,temp;

	static double	*f=NULL;

	static double	eps=1.e-5;

  	
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
                params[0] = x[0];
                params[1] = x[1];
                if (n == 3) params[3]=x[2];
                
                (vecfunc)(vpr,f, n, params);

		/* partielle Ableitungen */
		h = 1./h; 
		for (i=0;i<n;i++) df[i][j]=(f[i]-icurr[i])*h; 

		/* Rueckschreiben der Parameter */
		x[j]=temp;

		}

/*
printf("fdjac: df[0][0]=%g, df[1][0]=%g, df[0][1]=%g, df[1][1]=%g\n",
		df[0][0], df[1][0], df[0][1], df[1][1]);
*/

	return 0;

    }	/* fdjac */





int lnsrch(n, xold, fold, g, p, x, f, stpmax, check, 
		vpr, ipr, icurr, params)

    int	n;
    double	xold[], fold, g[], p[], x[], *f, stpmax;
    int		*check;
    double	*vpr, *ipr, *icurr, *params;

    {
	int 	i, n_try ;
	double 	a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	static int	n_try_max = 100;
	/* static double 	tolx=1.e-10; */


	*check=0;
	for (sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
	if (sum > stpmax*stpmax) {
		temp = stpmax/sqrt(sum);
		for (i=0;i<n;i++) p[i] *= temp;
		}
	for (slope=0.0,i=0;i<n;i++)
		slope += g[i]*p[i];

	for (test=0.0, i=0;i<n;i++) 
		test = DMAX( test, fabs(p[i])/DMAX(fabs(xold[i]),1.0) );
		
	/* minimale Schrittweite */
	alamin=GTOLX/test;
	
	/* Beginne immer mit vollem Newton-Schritt, laenger gehts nicht */
	alam=1.0;

	/* Endlosschleife, Funktion wird durch Returns innerhalb verlassen */
	for (n_try=0; n_try < n_try_max; n_try++) {
	
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];

		*f=calc_fmin(x,n,mag_doppel, vpr, ipr, icurr, params);

		/* Abbruchbedingungen: 1) lambda zu klein 
				       2) f_min klein genug */
		if (alam < alamin) {
		    for (i=0;i<n;i++) x[i]=xold[i];
		    /* printf("lnsrch: Abbruch da lambda zu klein, n_try =%d\n", n_try); */
		    *check=1;
		    return 0;
		    } 
		else if (*f <= fold + 1.e-4*alam*slope) {
		    /* printf("lnsrch: gewolltes Ende, da Verringerung in f\n"); */
		    return 0;
		    }

		/* Berechne neues lambda aus kubischer Gleichung.
		   Im ersten Durchgang (alam=1.) muss lediglich
		   quadratische Gleichung geloest werden */
		if (alam == 1.0)
		    tmplam = -slope/(2.0*(*f-fold-slope));
		else {
		    temp = alam-alam2;
		    rhs1 = (*f-fold-alam*slope)/(alam*alam*temp);
		    rhs2 = (f2-fold-alam2*slope)/(alam2*alam2*temp);
		    a = rhs1-rhs2;
		    b = -alam2*rhs1 + alam*rhs2;
		    /* loese kubische Gleichung in lambda */
		    if (a == 0.0) 
		        tmplam = -slope/(b+b);
		    else {
		        temp = a+a+a;
			disc=b*b-temp*slope;
			if (disc<0.0) {
			    /* printf("Roundoff problem in lnsrch, disc=%g, slope=%g, a=%g, b=%g\n",disc,slope, a, b); */
			    return -21; 
			    }
			else 
			    tmplam=(-b+sqrt(disc))/temp;
			}
		    /* tmplam = DMIN( tmplam, 0.5*alam ); */
		    }

		tmplam = DMIN( tmplam, 0.5*alam ); 
		alam2 = alam;
		f2    = *f;
		alam  = DMAX(tmplam,0.1*alam);
		}


        /* um Geschwindigkeit zu gewinnen, wird nach n_try_max Versuchen
	   die Iteration abgebrochen und zum Nebenminimum erklaert */
	*check = 1;
	/* printf("lnsrch: Abbruch nach n_try_max Versuchen\n"); */
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


	static double  	**fjac=NULL, *g=NULL, *p=NULL, *xold=NULL,
			*icurr=NULL; 

	int	rc;

	static int	maxits=100, maxsteps=100;
	static double	tolf=GTOLF;



/*	printf( "Beginn von newt: vpr = %g, %g, ipr = %g, %g \n", 
			vpr[0], vpr[1], ipr[0], ipr[1] );
	printf( "                  x = %g %g\n", x[0], x[1]);
*/

	/* Lege benoetigte Datenfelder an, es gilt immer n<=3 */
	if ( xold == NULL ) {
	    fjac  = dmatrix(0,2,0,2);
	    g     = dvector(0,2);
	    p     = dvector(0,2);
	    xold  = dvector(0,2);
	    icurr = dvector(0,2);
	    /* Initialisierung  der loop-up-tables */
	    /* init_quick_075(); */
	    }


	/* Fit erfolgt nicht direkt in n_e, T_e und beta */
    	/* x[0] = log(fabs(x[0]*sqrt(x[1]))); */	/* n_e sqrt(T_e) > 0.0 */
	x[1] = log( fabs(x[1]) ); 		/* T_e > 0.0 */
    	x[0] = log(fabs(x[0])) + 0.5*x[1];	/* n_e sqrt(T_e) > 0.0 */
	if ( n>2 ) x[2] = log( fabs(x[2]) );	/* beta > 0.0 */

        /* in allen positiven Parametern logarithmischer Fit */
        /* if ( params[0] <= 0.0 ) params[0]=1.e-30 ; */
        /* if ( params[1] <= 0.0 ) params[1]=1.e-30 ; */
        if ( params[2] <= 0.0 ) params[2]=1.e-30 ;
        if ( params[3] <= 0.0 ) params[3]=1.e-30 ;
        if ( params[5] <= 0.0 ) params[5]=1.e-30 ;
    
        /* params[0] = log(params[0]*sqrt(params[1])); */
        /* params[1] = log(params[1]); */
        params[2] = log(params[2]);
        params[3] = log(params[3]);
        params[5] = log(params[5]);
    
        /* Fit in cos(psi) statt in psi */
        params[6] = cos(params[6]);
        params[7] = cos(params[7]);


	/* Inititalisiere Parameter zur Funktionsberechnung,
	   berechen Funktionswerte "icurr" zu Startparametern */
	*f=calc_fmin(x, n, mag_doppel, vpr, ipr, icurr, params);


	/* die minimale Abweichung tolf, ab der die Ergebnisse akzeptiert
	   werden, muss relativ zu ipr festgelegt werden */
	for (tolf=0, i=0; i<n; i++) tolf = DMAX( tolf, fabs(ipr[i]));
	tolf = 1.e-8*tolf;

	/* maximale Abweichung der Funktionswerte von Sollwerten */
	for (test=0.0, i=0;i<n;i++) test = DMAX( test, fabs(icurr[i]-ipr[i]) );
	if (test<0.01*tolf) {
	    /* Variablen muessen wieder retransformiert werden */
    	    x[0] = exp(x[0]-0.5*x[1]);
    	    x[1] = exp(x[1]);
            if (n>2) x[2] = exp(x[2]) ;
            /* logarithmischer Fit in allen positiven Parametern */
            params[0] = exp(params[0]-0.5*params[1]);
            params[1] = exp(params[1]);
            params[2] = exp(params[2]);
            params[3] = exp(params[3]);
            params[5] = exp(params[5]);
        
            params[6] = acos(params[6]);
            params[7] = acos(params[7]);
            /* printf("newt: return ohne Suche, da Optimum bereits erreicht\n"); */
	    return 0;
	    }

	/* Schrittweite fuer lineare Suche */
	for (sum=0.0,i=0;i<n;i++) sum += DSQR(x[i]);
	stpmax= maxsteps*DMAX(sqrt(sum),(float)n);

    	/* Beginn der Iteration */
	for (its=0;its<maxits;its++) {

		/* berechne partielle Ableitungen */
		fdjac(n,x,icurr,fjac,mag_doppel, vpr, params);

		/* 1) negative Abweichung im Strom (p) und */		
		/* 2) Speichere Ist-Zustand (x) zu Vergleichszwecken */
		/* 3) Gradientenrichtung (g) fuer lineare Suche und */
		for (i=0;i<n;i++) {
		    /* p[i] = -(icurr[i]-ipr[i]); */
		    p[i] = ipr[i]-icurr[i];
		    xold[i]=x[i];
		    }
		fold=*f;
		for (i=0;i<n;i++) {
		    for (sum=0.0,j=0;j<n;j++) 
			sum -= fjac[j][i]*p[j];
		    g[i]=sum;
		    }


		/* Loese lineares Gleichungssystem: J * dx = -F */
		if ( (rc = dsvd_solver( fjac, p, n, 1)) != 0 ) break;
		
		/* lineare Suche entlang Newton-Richtung */
		if ( (rc = lnsrch(n, xold, fold, g, p, x, f, stpmax, check, 
				vpr, ipr, icurr, params)) != 0 ) break;

		/* Konvergenzkriterien erfuellt? */
		for (test=0.0, i=0;i<n;i++)
			test = DMAX( test, fabs(icurr[i]-ipr[i]));
		if (test < tolf) {
		    *check=0;
		    rc = 0;
		    /* printf("newt: normales Ende, Konvergenzkriterium erfuellt\n"); */
		    break;
		    }

		/* Gradienten gleich Null */
		if (*check) {
			/* warning */
		    	rc =  130;
			den=DMAX(*f,0.5*n);
			for (test=0.0, i=0;i<n;i++) 
			     test = DMAX( test, 
				       fabs(g[i])*DMAX(fabs(x[i]),1.0) );
			*check=( (test/den) < (GTOLMIN) ? 1 : 0);
			/* error wegen Nebenminimum */
		    	if (*check) rc = -130;
		    	/* printf("newt: Abbruch wegen Null-Gradient\n"); */
			break;
			}

		/* keine nennenswerte Aenderung in x mehr moeglich */
		for (test=0.0, i=0;i<n;i++) 
		    test = DMAX( test, 
				(fabs(x[i]-xold[i]))/DMAX(fabs(x[i]),1.0) );
		if (test < GTOLX) {
		    rc = 100;
		    /* printf("newt: Abbruch da keine Aenderung in x\n"); */
		    break;
		    }

		}		/* for (its=0;...) */


	/* keine Konvergenz nach endlicher Anzahl von Iterationsschritten? */
	if ( its >= maxits ) {
	   /* printf("newt: failt to converge, its > maxits = %d\n", maxits); */
	   rc = -1;
	   }

	/* Variablen muessen wieder retransformiert werden */
    	x[0] = exp(x[0]-0.5*x[1]);
    	x[1] = exp(x[1]);
        if (n>2) x[2] = exp(x[2]) ;
	
        /* logarithmischer Fit in allen positiven Parametern */
        params[0] = exp(params[0]-0.5*params[1]);
        params[1] = exp(params[1]);
        params[2] = exp(params[2]);
        params[3] = exp(params[3]);
        params[5] = exp(params[5]);
    
        params[6] = acos(params[6]);
        params[7] = acos(params[7]);

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

	n_points 	Anzahl der Punkte, die ein Gleichungssystem
			ergeben, Anzahl der Gleichungen, die simultan
			geloest werden muessen
			
	n_data		Anzahl der (Zeit-)Punkte, zu denen jeweils
			ein Gleichungssystem geloest werden soll
			
			
  	-------------------------------------------------------- 

	M. Weinlich, 19.06.95

 	-------------------------------------------------------- */



int mag_solver( data_x, data_y, psi, n_points, n_data, pars, n_pars, 
                ne_result, te_result, beta_result, fmin )


    double	*data_x, *data_y, *psi;
    int		*n_points, *n_data;
    double	*pars;
    int		*n_pars;
    double 	*ne_result, *te_result, *beta_result;
    double	*fmin;
    
    {
    double	vars[6];	/* Variable, in denen nl. Gl.s. geloest wird */

    int		rc, rc_counter, i, j, k, data_ind;
    int		check;
    double	save_ne, save_te, save_beta;
    double	i_isat, i_min, a_probe_proj, ne_start, te_start;
    
    static double	*nb_const=NULL, *nb_values=NULL, 
    			*fit_vec, *w_vec;
    static int 		*do_var, n_max_pars=20;
    double		fit_epsa, fit_epsr;
    int			fit_imax;
    double		fit_ok;
    int			do_fit, do_nls, n_fit=0, n_nls=0, n_nls_failed=0;
	
    static int	mag_infos = 0; 	/* ggf. kein Output von magfit */
    
    double		el_t; /* Zeitmessung */

/*
printf("neue Version mag_solver, mag_doppel etc., n_pars = %d \n", *n_pars);
*/

/*
    printf( "Uebergebene Werte in mag_solver: \n\n");
    printf( "n_points = %d \n", *n_points );
    printf( "n_data = %d \n", *n_data );
    printf( "n_pars = %d \n", *n_pars );
    printf( "x = ");
    for (i=0; i<(*n_data)*(*n_points); i++) printf( "%g  ", data_x[i] );
    printf("\n");
    printf( "y = ");
    for (i=0; i<*(n_data)*(*n_points); i++) printf( "%g  ", data_y[i] );
    printf("\n");
    printf("cos(psi) = ");
    for (i=0; i<*n_data; i++) printf( "%g  ", cos(psi[i]) );
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

    /* Datenfelder, fuer evtl. Fit */
    if ( nb_const == NULL ) {
        nb_const  = dvector(0,n_max_pars-1);
        nb_values = dvector(0,n_max_pars-1);
        do_var    = ivector(0,n_max_pars-1);
        for (i=0; i<n_max_pars; i++) {
            nb_const[i] = nb_values[i] = 0.;
            do_var[i] = 0;
            }
        do_var[0] = do_var[1] = 1;
        fit_vec   = dvector(0,2);
        w_vec     = dvector(0,2);
        w_vec[2] = w_vec[1] = w_vec[0] = 1.;
        }

    /* wenn Fit, in wievielen Parametern? */
    if (*n_points > 2) do_var[3] = 1;
    		else do_var[3] = 0;

    /* Initialisierung und Sicherheitskopie der Startwerte */
    vars[0] = save_ne   = pars[0];
    vars[1] = save_te   = pars[1];
    vars[2] = save_beta = pars[3];


    /* Eigentlicher nl-Solver wird nun n_data-Mal durchgelaufen und 
       fittet in jedem Durchgang n_points Datenpunkte */

   for ( i=0; i<*n_data; i++ ) {

	/* Winkel muss jeweils aktualisiert werden */	
	pars[6] = pars[7] = psi[i];
	/* pars[7] = psi[i]; */
	    
	/* projezierte Flaeche der Sonde */
	a_probe_proj = 0.;
	if ( pars[11] > 1.e-6 ) a_probe_proj += pars[10]*pars[11]*cos(psi[i]);
	if ( pars[17] > 1.e-6 ) a_probe_proj += pars[10]*pars[17] ;

	/* Startpunkt der aktuellen Teilvektoren */    
	data_ind = i*(*n_points);
	
	/* Extremwerte in Strom */
	i_isat = 0.;
	for (j=0; j<*n_pars; j++ ) {
	   i_isat = DMIN( i_isat, data_y[data_ind+j] );
	   }
	    
	/* Startwerte, fuer Deuterium-Plasma */
	te_start = save_te;
	ne_start = fabs(i_isat) / a_probe_proj / 
			sqrt( (2.*c_e*c_e*c_e/c_mp) * save_te );
/*
printf("\n");
printf("%g %g %g\n", pars[10],pars[11],cos(psi[i]));
printf("ne_start=%g, te_start=%g \n", ne_start, te_start);
*/
	
	/* beta ist eigentlich eine stabile Groesse. Damit sich evtl. Fehler
	   nicht zusehr fortpflanzen koennen, sollte jedoch immer mit einer
	   sinnvollen Groesse begonnen werden */
	if ( save_beta > 20 ) save_beta = 20;
	
	/* ist Loesung ueber nl-solver ueberhaupt zu erwarten, oder die liegen
	   die Punkte mit Sicherheit nicht auf einer Kennlinie? */
        if (data_y[data_ind]*data_y[data_ind+1] < 0.0) {
            if ( fabs(data_x[data_ind]) < fabs(data_x[data_ind+1]) ) {
                do_nls = ( data_y[data_ind+1] < 
                          (data_y[data_ind]*data_x[data_ind+1]/data_x[data_ind]) );
                }
            else {
                do_nls = ( data_y[data_ind+1] > 
                          (data_y[data_ind]*data_x[data_ind+1]/data_x[data_ind]) );
                }
	    }
	  
	/* zu grosse beta-Werte (oder zu hohe Nichtsaettigung im 
	   Elektronensaettigungsast lassen keine direkte Loesung zu */
	do_nls = ( do_nls && 
	           (fabs(data_y[data_ind+(*n_points)-1]/data_y[data_ind])<60));
	    	   
	/* default: Fit, wenn nicht noetig, wird Fit explizit abgewaehlt */
	do_fit = (0 == 0);
	         
  	/* erster Loesungs-Versuch */
  	if (do_nls) {
  	    n_nls++;
            vars[0] = ne_start;
            vars[1] = te_start;
            vars[2] = save_beta;
            rc = newt( vars, *n_pars, &check, 
                          data_x+data_ind, data_y+data_ind, pars, fmin+i);
	    /* im Fehlerfall: zweiter Versuch mit Fit */
	    /* Fehler: rc < 0 oder Nebenminimum mit T_e > 1eV */
            do_fit = ( (rc < 0) || ((check) && (vars[1]>1.)) );  
            }
          
/*        if ( do_fit && do_nls ) {
           printf( "i=%d, rc=%d \t %f %f %f \t %f %f %f \n", i, rc, 
           	data_x[data_ind], data_x[data_ind+1], data_x[data_ind+2],  
           	data_y[data_ind], data_y[data_ind+1], data_y[data_ind+2]);
           	}  
*/
           
	/* zweiter Versuch mit Fit */
        if ( do_fit ) { 

            pars[0] = ne_start * sqrt(te_start/10.);	/* n_e */
            pars[1] = 10.;				/* T_e */
            pars[3] = save_beta;			/* beta */
            pars[6] = pars[7] = psi[i];			/* psi1, psi2 */
            fit_imax = 100;
            fit_epsa = 0.;
            fit_epsr = 1.e-4;
      
            rc = magfit( mag_doppel, data_x+data_ind, data_y+data_ind, 
            		w_vec, fit_vec, n_points, pars, do_var, &n_max_pars,
            		nb_values, nb_const, fmin+i, 
			&fit_imax, &fit_epsa, &fit_epsr, &mag_infos);

	    /* Loesung wird in vars erwartet */
	    vars[0] = pars[0];
	    vars[1] = pars[1];
	    vars[2] = pars[3];

	    }


	/* Hallo, ich lebe noch .... */
        if (i%1000 == 999) printf(".");


	/* Nur wenn obige Schleife mit akzeptablem Ergebnis verlassen wurde,
	   dann sollen die Ergebnisse uebernommen werden */
	if (rc < 0) {
	    /* setzte Ergebinsse < 0, als Zeichen fuer Fehler in Auswertung
	       Signatur ist eindeutig, da alle Parameter positiv sind */
            ne_result[i]   = -0.01;
            te_result[i]   = -0.01;
            beta_result[i] = -0.01;
	    fmin[i] *= -1.;
            }
        else
            {
            /* printf("mag_solver: Ergebnis uebernommen\n"); */
            ne_result[i]   = vars[0];
            te_result[i]   = vars[1];
            beta_result[i] = vars[2];
            /* Aenderungen in Startwerten sollen etwas gedaempft werden */
/*            save_te = 0.5*(save_te+vars[1]);*/
/*            save_beta = 0.5*(save_beta+vars[2]);*/
            /* einige Statistik */
            if ((do_fit) && (n_fit%20 == 19)) printf( "+" );
            if (do_fit) n_fit++;
            if (do_fit && do_nls) n_nls_failed++;
	    }  
   
	}	/* for */
 
    if (n_fit+n_nls > 100) {
       printf("\n\nmag_solver: n_fit = %d ,  n_nls = %d,  n_nls_failed = %d\n\n", 
    		n_fit, n_nls, n_nls_failed);
       }
       
    return 0;

    }		/* mag_solver */








/* 	-------------------------------------------------------- 


		Es folgen die Interface-Routine zu IDL
		======================================

 	-------------------------------------------------------- */





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
    
    static int	mag_infos = ( 0 == 1);	/* ggf. keine Infos aus magfit */

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
	    		 argv[11], argv[12], argv[13], argv[14], &mag_infos);
	    return idl2mag_err(rc, routine, nargs);
	    }
        }        
      else if ( strcmp( routine_lowcase, "mag_solver") == 0 ) {
	/* Routine benoetigt insgesamt 11 Uebergabe-Variablen:
	   data_x, data_y, psi, n_points, n_data, pars, n_pars, 
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
    double 		chi;
    int			do_var[20];
    int			mag_infos;
    int 		n_x, err;
    double		x[3], x_nl[3];
    double	*vpr,*ipr,*w, *fit, *strom;

    int 		used, unit;
    double 	fixed[13], params[20];
    double 	eps_abs, eps_rel;

    int			n_pars;
    double		*nb_values, *nb_const;

    FILE		*demo_file;

    int 		c_6 = 6;
   char		*f_name = "./demo.dat"; 
/*    char		*f_name = "/afs/ipp/u/mnw/C/marquardt/lsm_sample.dat";*/
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


    goto test_fit; /* */
/*   goto test_nl ; /* */
/*   goto test_kennlinie; /* */

   /* erzeuge Testmatrix */
   n = 3;
   mat     = dmatrix(0,n-1,0,n-1);
   mat_inv = dmatrix(0,n-1,0,n-1);
   vec     = dvector(0,n-1);
   solve   = dvector(0,n-1);
                            
    /* Ausgabe der Startwerte */
    printf("Manfred's Testprogramm für mag_fit.c\n\n");
/*    printf("\n\n\n\n\n\n\n\n");
    printf("Loese lineares Geleichungssystem Ax=b ueber singular value decomposition\n\n");
*/
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
    vec[1]  = -1.;
    vec[2]  = -4.;


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

    for (i=0;i<n;i++){
       solve[i] = vec[i];
       for (j=0;j<n;j++) mat_inv[i][j] = mat[i][j];
       }
    /* loese das System mat * x = vec */
    rc = dsvd_solver( mat_inv, solve, n, 0);

    printf("Loesung des Systems: \n");
    printf(" Inv =");
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            printf( "\t %g ", mat_inv[i][j]);
            }
        printf("\n");
        }

    printf("\n x =");
    for (i=0; i<n; i++ ) printf( "\t %g", solve[i] );
    printf("\n\n");



    printf(" A Inv = " );
    for (i=0;i<n;i++){
        for (j=0;j<n;j++) {
            temp=0.;
            for (k=0;k<n;k++) temp += mat[i][k]*mat_inv[k][j];
            printf( "  %g  ", temp);
            }
        printf("\n         ");
        }
    printf("\n");

    printf(" Ax   = " );
    for (i=0;i<n;i++){
        temp=0.;
        for (j=0;j<n;j++) temp += mat[i][j]*solve[j];
        printf( "  %g  ", temp);
        }
    printf("\n\n");
    

    /* Loeschen der Felder */
    rc = dsvd_solver( mat, vec, -n, 1 );



    /* auf jeden Fall ist's hier jetzt zu Ende */
    if (n>0)  return 0;





	/*
		-----------------------------------
		Test fuer Berechnung der Kennlinie
		-----------------------------------
	*/


test_kennlinie:

    n_x = 100;
    vpr = dvector(0,n_x-1);
    w = dvector(0,n_x-1);
    fit = dvector(0,n_x-1);
    strom = dvector(0,n_x-1);


    /* 	Vorbelegung der Parameter */

    params[0] = 1e19;			/* n_e */
    params[1] = 10.;			/* t_e */
    params[2] = 1.;			/* alpha_1 */
    params[3] = 4.;			/* beta */
    params[4] = 0.;			/* delta_V */
    params[5] = 1;			/* alpha_2 */

    params[6] = 1.553343034;		/* 89. Grad */
    params[6] = 1.;		
    params[7] = 1.553343034;		/* 89. Grad */

    params[10] = .005;			/* a_width */
    params[11] = .04;			/* a_height */
    params[12] = 2.;			/* a_mass */
    params[13] = 1.;			/* a_z */
    params[14] = 2.;			/* a_bt */

    params[17] = -0.1;     		/* nur flush mountet probe */

    
    /* kein Fit, lediglich Berechnung einer Kennlinie */
    for (i=0; i<20; i++ ) do_var[i]=0;
    
     for (i=0; i < n_x; ++i) {
        /* Vorgabe der Sondenspannung */
        vpr[i] = 0.2113*((float)i-(0.5*(float)n_x))*params[1];
        /* konstante Gewichtung */
	w[i] = 1.;
	/* Dummy fuer Strom */
	strom[i] = 0.;
    	}
    	    
    n_pars = 20;

    nb_values = dvector(0, n_pars-1);
    nb_const  = dvector(0, n_pars-1);
    for (i=0; i<n_pars; i++) nb_values[i]=0.0;
    for (i=0; i<n_pars; i++) nb_const[i]=0.0;

    eps_rel  = 1.e-4;
    eps_abs  = 0.;
    iter_max = 200;
    mag_infos = ( 0 == 0 );

    magfit(mag_doppel, vpr, strom, w, fit, &n_x, params, do_var, &n_pars, 
		nb_values, nb_const, 
		&chi, &iter_max, &eps_abs, &eps_rel, &mag_infos); /* */

/*    magfit(fast_doppel, vpr, strom, w, fit, &n_x, params, do_var, &n_pars, 
		nb_values, nb_const, 
		&chi, &iter_max, &eps_abs, &eps_rel, &mag_infos); /* */


    printf( "\n\n");
    printf( "     Kennlinie: \n");
    printf( "-------------------- \n\n\n");

    printf( "  V_sonde   \t  I_sonde \n");
    for (i=0; i<n_x; i++) printf("%g \t %g \n", vpr[i], fit[i] );
    printf("\n\n");
    
    /* Ende dieses speziellen Tests */
    return 0;





	/*
		-------------------------
		Test fuer Fit-Algorithmus
		-------------------------
	*/


test_fit:


    vpr = dvector(0,199);
    w = dvector(0,199);
    fit = dvector(0,199);
    strom = dvector(0,199);

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
    params[2] = 1.;			/* alpha_1 */
    params[3] = 1.;			/* beta */
    params[4] = 0.;			/* delta_V */
    params[5] = 1;			/* alpha_2 */

    params[6] = 1.5603243;		/* 89.4 Grad */
    params[7] = 1.563815;		/* 89.6 Grad */

    params[10] = .005;			/* a_width */
    params[11] = .04;			/* a_height */
    params[12] = 2.;			/* a_mass */
    params[13] = 1.;			/* a_z */
    params[14] = 2.;			/* a_bt */

    params[16] = 1;			/* a_tau */
    params[17] = 0.00041887;		/* a_height * cos(psi1) */

    for (i=0; i<6; i++ ) do_var[i]=1;
    for (i=6; i<20; i++ ) do_var[i]=0;
    
    /* kein Fit in alpha1, alpha2, psi1, psi2, verwende projizierte Fl"ache */
    do_var[2] = 0; 
/*    do_var[5] = 0; */
    do_var[6] = 0;
    do_var[7] = 0;
/*    params[11] = -1.; */
 

    /* nur flush mountet probe */
    params[17] = -1.; 

    /* konstante Gewichtung */
    for (i=0; i < n_x; ++i) {
	w[i] = 1.;
    	}

    n_pars = 20;

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
/* */
    /* bis hier Einschub fuer Geschwindigkeitsvergleich */

    eps_rel  = 1.e-4;
    eps_abs  = 0.;
    iter_max = 200;
    mag_infos = ( 0 == 0 );

    magfit(mag_doppel, vpr, strom, w, fit, &n_x, params, do_var, &n_pars, 
		nb_values, nb_const, 
		&chi, &iter_max, &eps_abs, &eps_rel, &mag_infos); /* */



    /* Fuer Geschwindigkeitsvergleich 50 x Fit mit vergleichbaren Startwerten */
    } /* */


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

    printf( "\n\n \t Parameter  \n");
    for ( i=0; i<6; i++) {
 	printf( "\t %g  \n", params[i] );
 	}


/*    printf( "\n\n V   \t I \t I_fit \n");
    for (i=0; i<n_x; i++)
        printf( "%8.3f \t %8.3f \t %8.3f \n", vpr[i], strom[i], fit[i]);
*/        
    printf("\n\n\n");

    return 0;


    /* nichtlineare Gleichungssolver fuer Kennlinie einer Triple-Sonde */

test_nl:

    vpr = dvector(0,10);
    ipr = dvector(0,10);

    x[0] = 1.67707e+17;
    x[1] = 10.;
    x[2] = 4.;

    for (i=0; i<3; i++) x_nl[i]=x[i];
    
    vpr[0] = -38.;
    vpr[1] =  5;
    vpr[2] =  30.;

    ipr[0] = -0.5;
    ipr[1] =  0.5;
    ipr[2] =  5.;

    vpr[0] = -70.;
    vpr[1] =  10.;
    vpr[2] =  80.;

    ipr[0] = -0.25;
    ipr[1] =  0.25;
    ipr[2] =  0.65;

    /* Fit/Loesung in n Parametern */
    n = 2;

    /* Auswertung fuer einen einzigen Zeitpunkt */
    n_x = 1;
    
    n_pars = 20;
    for (i=0; i<n_pars; i++) params[i]=0.;
    params[2] = 1.;
    params[3] = 4.;
    params[4] = 0.;
    params[5] = 1.;
    params[6] = 1.5603243;		/* 89.4 Grad */
    params[7] = 1.563815;		/* 89.6 Grad */
    params[10] = 0.005;
    params[11] = 0.030;
    params[12] = 2.;
    params[13] = 1.;
    params[14] = 2.;
    params[16] = 1.;

    el_t = ((double) clock()) ;
/*     for (i=0;i<50;i++) {
	/* printf("i = %d\n",i); */
	/* Loesung ueber Fit */
	params[0] = 6.72188e17;
	params[1] = 5.509299;
        rc = mag_solver( vpr, ipr, &params[6], &n, &n_x, params, &n,
		x, x+1, x+2,  &fmin);   /* */

/*
int mag_solver( data_x, data_y, psi, n_points, n_data, pars, n_pars, 
                ne_result, te_result, beta_result, fmin )
*/
	/* Loesung ueber NL-solver */
	x_nl[0] = 6.72188e19;
	x_nl[1] = 5.509299;
	x_nl[2] = 4.;
	rc = newt( x_nl, n, &check, vpr, ipr, params, &fmin );   /* */
	
/* 	} /* */
 
    el_t = (((double) clock()) - el_t) * 1.e-6;

    printf( "\n\n nach Aufruf von newt: \n\n");
    printf( "vpr = ");
    for (i=0; i<n; i++) printf( " %g ", vpr[i]);
    printf( "\n");
    printf( "ipr = ");
    for (i=0; i<n; i++) printf( " %g ", ipr[i]);
    printf( "\n");
    printf( "x_nl = ");
    for (i=0; i<n; i++) printf( " %g ", x_nl[i]);
    printf( "\n");
    printf( "x    = ");
    for (i=0; i<n; i++) printf( " %g ", x[i]);
    printf( "\n");
    printf( "fmin = %g \n\n", fmin);


    printf( "rc = %d, check = %d \n\n", rc, check );
    printf( "\n");

    printf( "Zeitaufwand: %g sec CPU-Zeit\n", el_t);

    return 0;


    }		/* main */


