#include <stdio.h>
#include <math.h>

#include "mag_fit.h"

/* physikal. Konstanten */
#define	    c_e    1.6022e-19
#define	    c_eps0 8.8542e-12
#define	    c_me   9.1095e-31
#define	    c_mp   1.6726e-27
#define	    c_pi   3.1415927


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


int mag_doppel2(double *vpr, double *strom, int n_vpr, double *params)
{

    /* Local variables */
    
    static double	n_e=1., log_ne=0;
    static double	t_e=1., log_te=0;
    static double	inv_beta=1., log_invbeta=0.;
    static double	alpha2=1., log_alpha2=0.;
    
    static double	cos1=0., cos1_sqr=0., sin1=1.;
    
    double 	dv;
    double 	a_height, a_width, p_height, a_mass, a_proj_area, a_sonde_perp;
    double 	a_z, a_b_loc, a_tau;
    double 	delta_phi; 
    double	j_isat, I_isat, c_i_bohm;
    double	ld_sqr;
    double 	bfac;
    double	temp;

    int 	ind;
    int		n_changes=0;
    
    int 	verbose=0;

    /* here we are ... */
    if (verbose) printf("Beginn mag_doppel2\n");

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
        // Bereichsueberpruefung
        if (log_te < -6) {
            params[1] = log_te = -6;
            n_changes++;
        } else if (log_te > 6.) {
            params[1] = log_te = 6.;
            n_changes++;
        }
        t_e = exp(log_te);
    }   

    temp = params[0]-0.5*log_te;
    if (log_ne != temp) {
        log_ne = temp;
        // Bereichsueberpruefung
        if (log_ne < 15) {
            log_ne = 15;
            params[0] = log_ne+0.5*log_te;
            n_changes++;
        } else if (log_ne > 65.) {
            log_ne = 65.;
            params[0] = log_ne+0.5*log_te;
            n_changes++;
        }
        n_e = exp(log_ne);
    }
       
    if (log_invbeta != -params[3]) {
        log_invbeta = -params[3];
        // Bereichsueberpruefung
        if (log_invbeta < -5) {
            log_invbeta = -5;
            params[3] = -log_invbeta;
            n_changes++;
        } else if (log_invbeta > 5.) {
            log_invbeta = 5.;
            params[3] = -log_invbeta;
            n_changes++;
        }
        inv_beta = exp(log_invbeta);
    }
       
    dv = params[4];   
    
    if (log_alpha2 != params[5]) {
        log_alpha2 = params[5];
        // Bereichsueberpruefung
        if (log_alpha2 < -50) {
            params[5] = log_alpha2 = -50;
            n_changes++;
        } else if (log_alpha2 > 5.) {
            params[5] = log_alpha2 = 5.;
            n_changes++;
        }
        alpha2 = exp(log_alpha2);
    }
       
    if (cos1 != params[6]) {
        cos1 = params[6];
        // Bereichsueberpruefung
        if (cos1 < 1.e-5) {
            params[6] = cos1 = 1.e-5;
            n_changes++;
        } else if (cos1 > 1.) {
            params[6] = cos1 = 1.;
            n_changes++;
        }
        // abgeleitete Groessen
        cos1_sqr = cos1*cos1;
        sin1     = sqrt(1. - cos1_sqr);
    }
        
    
    // Abbruch, da keine sinnvollen Input-Parameter?
    if (verbose) {
        if (n_changes >= 1) {
	        printf( "Achtung: insgesamt %d Parameter mussten geaendert werden!\n",
		    n_changes);
	    }
    }

    if (n_changes > 4) {
        if (verbose) {
	        printf( "Abbruch, da insgesamt %d Parameter ", n_changes);
	        printf( "geaendert werden mussten  !\n");
	    }
	    return (-21); 
	}


    // ab hier kann einfach ausgelesen werden, keine Umwandlung mehr noetig
        
    // wichtige Sonden-Dimensionen
    a_width  = params[10]; 
    a_height = params[11]; 
    p_height = params[17];

    // projezierte Flaeche einer flush mounted Sonde
    a_proj_area  = a_width*a_height*cos1;	// <= 0, wenn nicht definiert

    // projezierte Flaeche einer pin probe
    a_sonde_perp = a_width*p_height*sin1;	// <= 0, wenn nicht definiert
    
    // Massenzahl des Fuellgases, meist Deuterium
    if (params[12] < 1.) params[12] = 2.;
    a_mass = params[12];
    
    // Kernladungszahl, Fuellgas meist Wassertoff-Isotop
    if (params[13] < 1.) params[13] = 1.;
    a_z = params[13];
    
    // Betrag des lokalen Magnetfelds am Ort der Sonde
    // Standard: |B_t| = 2.0T
    if (fabs(params[14]) < 0.6) params[14] = 2.;
   	a_b_loc = fabs(params[14]); 

    // T_i / T_e,  Ionen- und Elektronentemperatur sind postiv
    if (params[16] < 0.05) params[16] = 1.;
    a_tau = params[16];


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


    // Debye-Laenge vor Schichten
    ld_sqr = t_e * c_eps0 / n_e / c_e;
   
        
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

        // Initialisierung, falls n"otig
        if (a_proj_area <= 0.) {
            for (ind=0; ind<n_vpr; strom[ind++] = 0.0);
	    }

        // Kennlinie einer Doppelsonde mit Potentialverschiebung und
	    // Fl"achenverh"altnis beta
        I_isat = a_sonde_perp * j_isat;
        bfac = alpha2 * sqrt(ld_sqr) / a_width; 


        for ( ind=0; ind < n_vpr; ind++ ) 
	    {
            // Spannung an Doppelsonde: phi_sonde
            delta_phi = (dv - vpr[ind]) / t_e;

	        if (delta_phi < -300.) {
                strom[ind] +=  I_isat * 
                    (1. + bfac*quick_075(delta_phi)) / inv_beta;
		    } else if (delta_phi > 300.) {
                strom[ind] += -I_isat;
	        } else { 
	            temp = exp(delta_phi);
                strom[ind] +=  I_isat * 
                    (1. + bfac*quick_075(delta_phi) - temp) / (inv_beta + temp);
		    }
	    }


    }

    return 0;

}		 // mag_doppel


