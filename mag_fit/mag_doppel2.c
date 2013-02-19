#include "mag_fit.h"

/* physikal. Konstanten */
#define	    c_e    1.6022e-19
#define	    c_eps0 8.8542e-12
#define	    c_me   9.1095e-31
#define	    c_mp   1.6726e-27
#define	    c_pi   3.1415927

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))


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


int mag_doppel2(vpr, strom, n_vpr, params)

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


