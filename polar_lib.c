/* put together by Monika Moscibrodzka, 
last update: April 14 2020
All functions below (except jar_calc_pl()) 
are part of ipole - a code by Monika Moscibrodzka and Charles Gammie
(https://github.com/moscibrodzka/ipole)
*/

#include "decs.h"

//Planck function, invariant
double Bnu_inv(double nu, double Thetae)
{
    double x;
    x = HPL * nu / (ME * CL * CL * Thetae);

    if (x < 2.e-3)              /* Taylor expand */
      return ((2. * HPL / (CL * CL)) /
	      (x / 24. * (24. + x * (12. + x * (4. + x)))));
    else
      return ((2. * HPL / (CL * CL)) / (exp(x) - 1.));
}



/* from grtrans */
double g(double Xe){ return 1. - 0.11 * log(1 + 0.035 * Xe);}

double h(double Xe)
{
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
	cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
	0.011 * exp(-Xe / 47.2);
}

double Je(double Xe){ return 0.43793091 * log(1. + 0.00185777 * pow(Xe, 1.50316886));}

double jffunc(double Xe)
{
    double extraterm;
    extraterm =
	(0.011 * exp(-Xe / 47.2) -
	 pow(2., -1. / 3.) / pow(3.,
				 23. / 6.) * M_PI * 1e4 * pow(Xe + 1e-16,
							      -8. / 3.)) *
	(0.5 + 0.5 * tanh((log(Xe) - log(120.)) / 0.1));
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
	cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
	0.011 * exp(-Xe / 47.2) + extraterm;
}

double I_I(double x)
{
    return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) +
		     0.9977 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
								    1. /
								    3.));
}

double I_Q(double x)
{
    return 2.5651 * (1 + 0.93193 * pow(x, -1. / 3.) +
		     0.499873 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
								      1. /
								      3.));
}

double I_V(double x)
{
    return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
	    0.0292545 * pow(x, -0.5) + 2.03773 * pow(x,
						     -1. / 3.)) *
	exp(-1.8899 * pow(x, 1. / 3.));
}

double besselk_asym(int n, double x)
{

    if (n == 0)
	return -log(x / 2.) - 0.5772;

    if (n == 1)
	return 1. / x;

    if (n == 2)
	return 2. / x / x;

    fprintf(stderr,"this cannot happen\n");
    exit(1);
}
/* end from grtrans */


// emissivity, absroptivity, rotativity for thermal distribution function
// emissivities from Dexter 2016
void jar_calc_th(double nu,double Ne,double B,double Thetae,double theta,
		 double *jI, double *jQ, double *jU, double *jV,
		 double *aI, double *aQ, double *aU, double *aV,
		 double *rQ, double *rU, double *rV){

    double nusq,x, Xe, omega0, nuc, Bnuinv, Thetaer, wp2;
  
    nusq = nu * nu;
    Thetaer = 1. / Thetae;
    omega0 = EE * B / ME / CL;
    wp2 = 4. * M_PI * Ne * EE * EE / ME;
    
    /* Faraday rotativities for thermal plasma */
    Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));
   
    *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
	jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
		      6. * Thetae) * sin(theta) * sin(theta);
    
    *rU = 0.0;
    
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	(gsl_sf_bessel_Kn(0,Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2,Thetaer) * cos(theta);

    /* invariant rotativities (in fluid frame, nu measured in fluid frame) */
    *rQ *= nu;
    *rU *= nu;
    *rV *= nu;

    
    /*synchrotron emissivity */
    nuc = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
    x = nu / nuc;

    *jI = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_I(x);	// [g/s^2/cm = ergs/s/cm^3]
    *jQ = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_Q(x);
    *jU = 0.0;		// convention; depends on tetrad
    *jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / S3 / CL / Thetae / Thetae / Thetae * I_V(x);
    
    *jI /= nusq;
    *jQ /= nusq;
    *jU /= nusq;
    *jV /= nusq;
  
    /* invariant synchrotron absorptivity */
    Bnuinv = Bnu_inv(nu, Thetae);   /* Planck function */
    *aI = *jI / Bnuinv;
    *aQ = *jQ / Bnuinv;
    *aU = *jU / Bnuinv;
    *aV = *jV / Bnuinv;

}

// emissivity, absroptivity, rotativity for power-law distribution function
// rQV - Dexter 2016, Jones & Odell 1977, error sign in rQ in both, here corrected
// jIQV, and aIQV are approximate from Pandya et al. 2016 with signs of jQ and jV corrected
void jar_calc_pl(double nu,double Ne,double B,double gmin,double gmax,double p,double theta,
		 double *jI, double *jQ, double *jU, double *jV,
		 double *aI, double *aQ, double *aU, double *aV,
		 double *rQ, double *rU, double *rV){

    double nusq = nu * nu;
    double nuc = EE*B/(2.*M_PI*ME*CL);
    double emcon=Ne*EE*EE*nuc/CL;
    double j_pl,j_pl_I,j_pl_Q,j_pl_U,j_pl_V;

    j_pl = emcon*pow(3,p/2.)*(p-1)*sin(theta)/2./(p+1)/(pow(gmin,1-p)-pow(gmax,1-p))*
	gsl_sf_gamma((3*p-1)/12)*gsl_sf_gamma((3*p+19)/12)*
	pow(nu/nuc/sin(theta),-(p-1)/2);
    
    j_pl_I = j_pl;
    j_pl_Q = j_pl*(p+1.)/(p+7./3.);
    j_pl_U = 0.0;
    j_pl_V = 171./250.*sqrt(p)/tan(theta)*pow(nu/3./nuc/sin(theta),-0.5)*j_pl;
    
    /* invatiant emissivity */
    *jI = j_pl_I/nusq;
    *jQ = j_pl_Q/nusq;
    *jU = j_pl_U/nusq;
    *jV = j_pl_V/nusq;

      
    /* absorption */
    double a_pl,a_pl_I,a_pl_Q,a_pl_U,a_pl_V;
    emcon=Ne*EE*EE/nu/ME/CL;
    
    a_pl=emcon*pow(3,(p+1)/2.)*(p-1)/4./(pow(gmin,1-p)-pow(gmax,1-p))*
	gsl_sf_gamma((3*p+12)/12)*gsl_sf_gamma((3*p+22)/12)*
	pow(nu/nuc/sin(theta),-(p+2)/2);
    a_pl_I=a_pl;
    a_pl_Q=a_pl*0.75*pow(p-1,43/500);
    a_pl_U=0.0;
    a_pl_V=7/4.*pow(71/100*p+22/625,197/500)*
	pow(pow(sin(theta),-48/25)-1,64/125)*
	pow(nu/nuc/sin(theta),-0.5)*a_pl;
    
    /*invariant absorptivity*/
    *aI = a_pl_I*nu;
    *aQ = a_pl_Q*nu;
    *aU = a_pl_U*nu;
    *aV = a_pl_V*nu;
      
      
    /*add rotatativity for power-law*/
    double nub=EE*B/2./M_PI/ME/CL;
    double Cv=2.*(p+2.)/(p+1.);
    double rho1=Ne*EE*EE/ME/CL/nub/sin(theta)*(p-1.)/(pow(gmin,1.-p)-pow(gmax,1.-p));
    double numin=gmin*gmin*nub*sin(theta);//ok
    
    /* invariant rotativity, sign corrected */
    *rQ =  (rho1 * pow(nub*sin(theta)/nu,3) * pow(gmin,2-p) * (1.-pow(numin/nu,p/2.-1.)) *pow(p/2.-1.,-1)) *nu ;
    *rU = 0.0;
    *rV =  (Cv *rho1*pow(nub*sin(theta)/nu,2)*pow(gmin,-p-1)*log(gmin)/tan(theta)) * nu;

}




/* main function to integrate Stokes equations in fluid frame */
/* one step */
void integrate_Stokes(double dlam,double L_unit,
		      double jI,double jQ,double jU,double jV,
		      double aI,double aQ,double aU,double aV,
		      double rQ,double rU,double rV,
		      double SI0,double SQ0,double SU0,double SV0,
		      double *SI,double *SQ,double *SU,double *SV){

    if(INT_SPLIT){
	
	/* the 3 steps scheme for stokes integration*/
	/* need to be corrected for L_unit due to units*/
	
	/* step 1 */
	/* apply the Faraday rotation solution for a half step */
	double x = dlam * 0.5 *HPL*L_unit/(ME*CL*CL);
	
	//additional Stokes parameters
	double SI1, SQ1, SU1, SV1;
	double SI2, SQ2, SU2, SV2;
	
	double rdS = rQ * SQ0 + rU * SU0 + rV * SV0;
	double rho2 = rQ * rQ + rU * rU + rV * rV;
	double rho = sqrt(rho2);
	double c, s, sh;
	c = cos(rho * x);
	s = sin(rho * x);
	sh = sin(0.5 * rho * x);
	if (rho2 > 0) {
	    SI1 = SI0;
	    SQ1 = SQ0 * c + 2 * rQ * rdS / rho2 * sh * sh + (rU * SV0 - rV * SU0) / rho * s;
	    SU1 = SU0 * c + 2 * rU * rdS / rho2 * sh * sh + (rV * SQ0 - rQ * SV0) / rho * s;
	    SV1 = SV0 * c + 2 * rV * rdS / rho2 * sh * sh + (rQ * SU0 - rU * SQ0) / rho * s;
	} else {
	    SI1 = SI0;
	    SQ1 = SQ0;
	    SU1 = SU0;
	    SV1 = SV0;
	}
	/* done rotation solution half step */

	/* step 2 */
	/* apply full absorption/emission step */
	x = dlam *HPL*L_unit/(ME*CL*CL) ;
	double aI2 = aI * aI;
	double aP2 = aQ * aQ + aU * aU + aV * aV;
	double aP = sqrt(aP2);
	double ads0 = aQ * SQ1 + aU * SU1 + aV * SV1;
	double adj = aQ * jQ + aU * jU + aV * jV;
	
	if (aP > SMALL) { /* full analytic solution has trouble if polarized absorptivity is small */
	    double expaIx = exp(-aI * x);
	    double sinhaPx = sinh(aP * x);
	    double coshaPx = cosh(aP * x);
	    
	    SI2 = (SI1 * coshaPx * expaIx
		   - (ads0 / aP) * sinhaPx * expaIx
		   + adj / (aI2 - aP2) * (-1 + (aI * sinhaPx + aP * coshaPx) / aP * expaIx)
		   + aI * jI / (aI2 - aP2) * (1 - (aI * coshaPx + aP * sinhaPx) / aI * expaIx));
	    
	    SQ2 = (SQ1 * expaIx
		   + ads0 * aQ / aP2 * (-1 + coshaPx) * expaIx
		   - aQ / aP * SI1 * sinhaPx * expaIx
		   + jQ * (1 - expaIx) / aI
		   + adj * aQ / (aI * (aI2 - aP2)) * (1 - (1 - aI2 / aP2) * expaIx
						      - aI / aP2 * (aI * coshaPx + aP * sinhaPx) * expaIx)
		   + jI * aQ / (aP * (aI2 - aP2)) * (-aP + (aP * coshaPx + aI * sinhaPx) * expaIx));
	    
	    SU2 = (SU1 * expaIx
		   + ads0 * aU / aP2 * (-1 + coshaPx) * expaIx
		   - aU / aP * SI1 * sinhaPx * expaIx
		   + jU * (1 - expaIx) / aI
		   + adj * aU / (aI * (aI2 - aP2)) *
		   (1 - (1 - aI2 / aP2) * expaIx -
		    aI / aP2 * (aI * coshaPx +
				aP * sinhaPx) * expaIx)
		   + jI * aU / (aP * (aI2 - aP2)) *
		   (-aP + (aP * coshaPx + aI * sinhaPx) * expaIx));
	    
	    SV2 = (SV1 * expaIx
		   + ads0 * aV / aP2 * (-1 + coshaPx) * expaIx
		   - aV / aP * SI1 * sinhaPx * expaIx
		   + jV * (1 - expaIx) / aI
		   + adj * aV / (aI * (aI2 - aP2)) * (1 -
						      (1 - aI2 / aP2) * expaIx -
						      aI / aP2 * (aI * coshaPx +
								  aP * sinhaPx) * expaIx)
		   + jI * aV / (aP * (aI2 - aP2)) *
		   (-aP + (aP * coshaPx + aI * sinhaPx) * expaIx));
	    
	} else { /* this should really be a series expansion in aP */
	    SI2 = SI1 + x * jI;
	    SQ2 = SQ1 + x * jQ;
	    SU2 = SU1 + x * jU;
	    SV2 = SV1 + x * jV;
	}
	/* done absorption/emission full step */

	/* step 3 */
	/* apply second rotation half-step */
	x = dlam * 0.5 * HPL*L_unit/(ME*CL*CL);
	rdS = rQ * SQ2 + rU * SU2 + rV * SV2;
	rho2 = rQ * rQ + rU * rU + rV * rV;
	rho = sqrt(rho2);
	c = cos(rho * x);
	s = sin(rho * x);
	sh = sin(0.5 * rho * x);
	if (rho2 > 0) {
	    *SI = SI2;
	    *SQ = SQ2 * c + 2 * rQ * rdS / rho2 * sh * sh + (rU * SV2 - rV * SU2) / rho * s;
	    *SU = SU2 * c + 2 * rU * rdS / rho2 * sh * sh + (rV * SQ2 - rQ * SV2) / rho * s;
	    *SV = SV2 * c + 2 * rV * rdS / rho2 * sh * sh + (rQ * SU2 - rU * SQ2) / rho * s;
	} else {
	    *SI = SI2;
	    *SQ = SQ2;
	    *SU = SU2;
	    *SV = SV2;
	}
	/* done second rotation half-step */
    }
    

    /* full integration step from Landi Degl’Innocenti & Landi Degl’Innocenti (1985).*/
    /* with some regularizers to hangle small optical depth */
    if(INT_FULL){
    
	int k,l;
	double M1[NDIM][NDIM];
	double M2[NDIM][NDIM];
	double M3[NDIM][NDIM];
	double M4[NDIM][NDIM];
	double O[NDIM][NDIM];
	double P[NDIM][NDIM];
    
	double alpha2=aQ*aQ + aU*aU + aV*aV;
	double rho2  =rQ*rQ + rU*rU + rV*rV;
	double alphadrho = aQ*rQ + aU*rU + aV*rV;
	double sig=copysign(1.,alphadrho);
	double T=2.*sqrt(pow(alpha2 - rho2,2)/4.+ alphadrho*alphadrho ) ;
	double ith=1./(T+SMALL); //so that if a=r=0  M234->0, not NaN
        
	double L1=sqrt(T*0.5+(alpha2-rho2)*0.5)+SMALL; //so that if a=0 alone fac1 != 1/0
	double L2=sqrt(T*0.5-(alpha2-rho2)*0.5)+SMALL; //so that if a=r=0 fac2 != 0 
  
	//M1
	for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++) M1[k][l]=0.0;
	M1[0][0]=1.0;
	M1[1][1]=1.0;
	M1[2][2]=1.0;
	M1[3][3]=1.0;
      
	//M2
	M2[0][0]=0.0;
	M2[0][1]=ith*(L2*aQ - sig*L1*rQ);
	M2[0][2]=ith*(L2*aU - sig*L1*rU);
	M2[0][3]=ith*(L2*aV - sig*L1*rV);
	
	M2[1][0]=ith*(L2*aQ - sig*L1*rQ);
	M2[1][1]=0.0;
	M2[1][2]=ith*(sig*L1*aV + L2*rV);
	M2[1][3]=ith*(-sig*L1*aU - L2*rU);
    
	M2[2][0]=ith*(L2*aU - sig*L1*rU);
	M2[2][1]=ith*(-sig*L1*aV - L2*rV);
	M2[2][2]=0.0;
	M2[2][3]=ith*(sig*L1*aQ + L2*rQ);
	
	M2[3][0]=ith*(L2*aV - sig*L1*rV);
	M2[3][1]=ith*(sig*L1*aU + L2*rU);
	M2[3][2]=ith*(-sig*L1*aQ - L2*rQ);
	M2[3][3]=0.0;
	

	//M3 
	M3[0][0]=0.0;
	M3[0][1]=ith*(L1*aQ + sig*L2*rQ);
	M3[0][2]=ith*(L1*aU + sig*L2*rU);
	M3[0][3]=ith*(L1*aV + sig*L2*rV);
      
	M3[1][0]=ith*(L1*aQ + sig*L2*rQ);
	M3[1][1]=0.0;
	M3[1][2]=ith*(-sig*L2*aV + L1*rV);
	M3[1][3]=ith*(sig*L2*aU - L1*rU);

	M3[2][0]=ith*(L1*aU + sig*L2*rU);
	M3[2][1]=ith*(sig*L2*aV - L1*rV);
	M3[2][2]=0.0;
	M3[2][3]=ith*(-sig*L2*aQ + L1*rQ);
	
	M3[3][0]=ith*(L1*aV + sig*L2*rV);
	M3[3][1]=ith*(-sig*L2*aU + L1*rU);
	M3[3][2]=ith*(sig*L2*aQ - L1*rQ);
	M3[3][3]=0.0;
    
	//M4
	M4[0][0]=2.*ith*(alpha2+rho2)/2.;
	M4[0][1]=2.*ith*(aV*rU - aU*rV);
	M4[0][2]=2.*ith*(aQ*rV - aV*rQ);
	M4[0][3]=2.*ith*(aU*rQ - aQ*rU);

	M4[1][0]=2.*ith*(aU*rV - aV*rU);
	M4[1][1]=2.*ith*(aQ*aQ + rQ*rQ - (alpha2+rho2)/2.);
	M4[1][2]=2.*ith*(aQ*aU + rQ*rU);
	M4[1][3]=2.*ith*(aV*aQ + rV*rQ);
  
	M4[2][0]=2.*ith*(aV*rQ - aQ*rV);
	M4[2][1]=2.*ith*(aQ*aU + rQ*rU);
	M4[2][2]=2.*ith*(aU*aU + rU*rU - (alpha2+rho2)/2.);
	M4[2][3]=2.*ith*(aU*aV + rU*rV);
	
	M4[3][0]=2.*ith*(aQ*rU - aU*rQ);
	M4[3][1]=2.*ith*(aV*aQ + rV*rQ);
	M4[3][2]=2.*ith*(aU*aV + rU*rV);
	M4[3][3]=2.*ith*(aV*aV + rV*rV - (alpha2+rho2)/2.);

	double fac1=1./(aI*aI - L1*L1);
	double fac2=1./(aI*aI + L2*L2);
	double x=dlam*HPL*L_unit/(ME*CL*CL);
	double l1dlam=L1*x;
	double l2dlam=L2*x;
	double EaIdlam=exp(-aI*x);
	double coshl1dlam=cosh(l1dlam);
	double sinhl1dlam=sinh(l1dlam);
	double cosl2dlam=cos(l2dlam);
	double sinl2dlam=sin(l2dlam);

	for(k=0;k<4;k++)for(l=0;l<4;l++){
		
		
		O[k][l] = EaIdlam*( 0.5*(coshl1dlam+cosl2dlam)*M1[k][l]
				    - sinl2dlam*M2[k][l]
				    - sinhl1dlam*M3[k][l]
				    + 0.5*(coshl1dlam-cosl2dlam)*M4[k][l]
		    );
	
     

		P[k][l] = (-L1*fac1*M3[k][l] + 0.5*aI*fac1*(M1[k][l] + M4[k][l])) +
		    (-L2*fac2*M2[k][l] + 0.5*aI*fac2*(M1[k][l] - M4[k][l])) -
		    EaIdlam*(
			( -L1*fac1*M3[k][l] + 0.5*aI*fac1*(M1[k][l] + M4[k][l]) )*coshl1dlam +
			( -L2*fac2*M2[k][l] + 0.5*aI*fac2*(M1[k][l] - M4[k][l]) )*cosl2dlam  +
			( -aI*fac2*M2[k][l] - 0.5*L2*fac2*(M1[k][l] - M4[k][l]) )*sinl2dlam  -
			(  aI*fac1*M3[k][l] - 0.5*L1*fac1*(M1[k][l] + M4[k][l]) )*sinhl1dlam
			);
	  
	    }
    
	*SI = P[0][0]*jI + P[0][1]*jQ + P[0][2]*jU + P[0][3]*jV + O[0][0]*SI0 + O[0][1]*SQ0 + O[0][2]*SU0 + O[0][3]*SV0;
	*SQ = P[1][0]*jI + P[1][1]*jQ + P[1][2]*jU + P[1][3]*jV + O[1][0]*SI0 + O[1][1]*SQ0 + O[1][2]*SU0 + O[1][3]*SV0;
	*SU = P[2][0]*jI + P[2][1]*jQ + P[2][2]*jU + P[2][3]*jV + O[2][0]*SI0 + O[2][1]*SQ0 + O[2][2]*SU0 + O[2][3]*SV0;
	*SV = P[3][0]*jI + P[3][1]*jQ + P[3][2]*jU + P[3][3]*jV + O[3][0]*SI0 + O[3][1]*SQ0 + O[3][2]*SU0 + O[3][3]*SV0;

    }

	
}


/// the rest of the functions are just tools to build tetrad and Lorentz transformations 
void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov)
{

  ucov[0] = Gcov[0][0] * ucon[0]
    + Gcov[0][1] * ucon[1]
    + Gcov[0][2] * ucon[2]
    + Gcov[0][3] * ucon[3];
  ucov[1] = Gcov[1][0] * ucon[0]
    + Gcov[1][1] * ucon[1]
    + Gcov[1][2] * ucon[2]
    + Gcov[1][3] * ucon[3];
  ucov[2] = Gcov[2][0] * ucon[0]
    + Gcov[2][1] * ucon[1]
    + Gcov[2][2] * ucon[2]
    + Gcov[2][3] * ucon[3];
  ucov[3] = Gcov[3][0] * ucon[0]
    + Gcov[3][1] * ucon[1]
    + Gcov[3][2] * ucon[2]
    + Gcov[3][3] * ucon[3];

  return;
}

/** tetrad making routines **/

/* 

   econ/ecov index key:

   Econ[k][l]
   k: index attached to tetrad basis
   index down
   l: index attached to coordinate basis 
   index up

   Ecov[k][l]
   k: index attached to tetrad basis
   index up
   l: index attached to coordinate basis 
   index down

*/

/* 

   make orthonormal basis for plasma frame.

   e^0 along U
   e^2 along b
   e^3 along spatial part of K

*/

void make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM],
			double Bcon[NDIM], double Gcov[NDIM][NDIM],
			double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM])
{
  int k, l;
  void normalize(double *vcon, double Gcov[4][4]);
  void project_out(double *vcona, double *vconb, double Gcov[4][4]);
  double check_handedness(double Econ[NDIM][NDIM],
			  double Gcov[NDIM][NDIM]);

  /* start w/ time component parallel to U */
  void set_Econ_from_trial(double Econ[4], int defdir, double trial[4]);
  set_Econ_from_trial(Econ[0], 0, Ucon);
  normalize(Econ[0], Gcov);

  /*** done w/ basis vector 0 ***/

  /* now use the trial vector in basis vector 3 */
  /* cast a suspicious eye on the trial vector... */
  set_Econ_from_trial(Econ[3], 3, Kcon);

  /* project out econ0 */
  project_out(Econ[3], Econ[0], Gcov);
  normalize(Econ[3], Gcov);

  /*** done w/ basis vector 3 ***/

  /* repeat for x2 unit basis vector */
  set_Econ_from_trial(Econ[2], 2, Bcon);

  /* project out econ0,3 */
  project_out(Econ[2], Econ[0], Gcov);
  project_out(Econ[2], Econ[3], Gcov);
  normalize(Econ[2], Gcov);

  /*** done w/ basis vector 2 ***/

  /* whatever is left is econ1 */
  for (k = 0; k < 4; k++)/* trial vector */
    Econ[1][k] = 1.;
  /* project out econ[0-2] */
  project_out(Econ[1], Econ[0], Gcov);
  project_out(Econ[1], Econ[2], Gcov);
  project_out(Econ[1], Econ[3], Gcov);
  normalize(Econ[1], Gcov);

  /* check handedness */
  double dot = check_handedness(Econ, Gcov);

  if (fabs(fabs(dot) - 1.) > 1.e-10) {
    fprintf(stderr, "that's odd: %g\n", fabs(dot) - 1.);
  }

    /* we expect dot = 1. for right-handed system.  
				     If not, flip Econ[1] to make system right-handed. */
  if (dot < 0.) {
    for (k = 0; k < 4; k++)
      Econ[1][k] *= -1.;
  }

  /*** done w/ basis vector 3 ***/

  /* now make covariant version */
  for (k = 0; k < 4; k++) {

    /* lower coordinate basis index */
    lower(Econ[k], Gcov, Ecov[k]);
  }

  /* then raise tetrad basis index */
  for (l = 0; l < 4; l++) {
    Ecov[0][l] *= -1.;
  }

  /* paranoia: check orthonormality */
  /*
       double sum ;
       int m ;
       fprintf(stderr,"ortho check [plasma]:\n") ;
       for(k=0;k<NDIM;k++)
       for(l=0;l<NDIM;l++) {
       sum = 0. ;
       for(m=0;m<NDIM;m++) {
       sum += Econ[k][m]*Ecov[l][m] ;
       }
       fprintf(stderr,"sum: %d %d %g\n",k,l,sum) ;
       }
       fprintf(stderr,"\n") ;
       for(k=0;k<NDIM;k++)
       for(l=0;l<NDIM;l++) {
       fprintf(stderr,"%d %d %g\n",k,l,Econ[k][l]) ;
       }
       fprintf(stderr,"done ortho check.\n") ;
       fprintf(stderr,"\n") ;
  */

  /* done! */

}

void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
				   double Ecov[NDIM][NDIM],
				   double complex T_tetrad[NDIM][NDIM])
{
  int i, j, k, l;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      T_tetrad[i][j] = 0. + I * 0.;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	      T_tetrad[i][j] +=
		T_coord[k][l] * Ecov[i][k] * Ecov[j][l];

  return;
}

void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
				   double Econ[NDIM][NDIM],
				   double complex T_coord[NDIM][NDIM])
{
  int i, j, k, l;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      T_coord[i][j] = 0. + I * 0.;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	      T_coord[i][j] +=
		T_tetrad[k][l] * Econ[k][i] * Econ[l][j];

  return;
}

void stokes_to_tensor(double fI, double fQ, double fU, double fV,
		      double complex f_tetrad[NDIM][NDIM])
{
  int i, j;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      f_tetrad[i][j] = 0. + I * 0.;

  f_tetrad[1][1] = (fI + fQ + 0. * I);
  f_tetrad[1][2] = (fU - I * fV);
  f_tetrad[2][1] = (fU + I * fV);
  f_tetrad[2][2] = (fI - fQ + 0. * I);

}

void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM],
		      double *fI, double *fQ, double *fU, double *fV)
{

  /*here I divide by two to agree with above */
  *fI = creal(f_tetrad[1][1] + f_tetrad[2][2]) / 2;
  *fQ = creal(f_tetrad[1][1] - f_tetrad[2][2]) / 2;
  *fU = creal(f_tetrad[1][2] + f_tetrad[2][1]) / 2;
  *fV = cimag(f_tetrad[2][1] - f_tetrad[1][2]) / 2;

}
/* input and vectors are contravariant (index up) */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
			  double K_tetrad[NDIM])
{
  int k;

  for (k = 0; k < 4; k++) {
    K_tetrad[k] = Ecov[k][0] * K[0] 
                  + Ecov[k][1] * K[1] 
                  + Ecov[k][2] * K[2] 
      + Ecov[k][3] * K[3];
  }
}

/* input and vectors are contravariant (index up) */
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
			  double K[NDIM])
{
  int l;

  for (l = 0; l < 4; l++) {
    K[l] = Econ[0][l] * K_tetrad[0] 
           + Econ[1][l] * K_tetrad[1] 
           + Econ[2][l] * K_tetrad[2] 
      + Econ[3][l] * K_tetrad[3];
  }

  return;
}



#define SMALL_VECTOR 1.e-30

/*
    Kronecker delta
*/
double delta(int i, int j)
{
if (i == j)
  return (1.);
 else
   return (0.);
}
/* 

    normalize input vector (and overwrite)
    so that |v . v| = 1

*/

void normalize(double *vcon, double Gcov[4][4])
{
int k, l;
double norm;

norm = 0.;
for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
    norm += vcon[k] * vcon[l] * Gcov[k][l];

norm = sqrt(fabs(norm));
for (k = 0; k < 4; k++)
  vcon[k] /= norm;

return;
}

/*

    project out vconb from vcona

    both arguments are index up (contravariant)

    covariant metric is third argument.

    overwrite the first argument on return

*/

void project_out(double *vcona, double *vconb, double Gcov[4][4])
{

double adotb, vconb_sq;
int k, l;

vconb_sq = 0.;
for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
    vconb_sq += vconb[k] * vconb[l] * Gcov[k][l];

adotb = 0.;
for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
    adotb += vcona[k] * vconb[l] * Gcov[k][l];

for (k = 0; k < 4; k++)
  vcona[k] -= vconb[k] * adotb / vconb_sq;

return;
}

/* 

   copy the trial vector into a tetrad basis vector,
   checking to see if it is null, and if it is null
   setting to some default value 

*/
void set_Econ_from_trial(double Econ[4], int defdir, double trial[4])
{
double norm = 0.;
int k;

for (k = 0; k < 4; k++)
  norm += fabs(trial[k]);
for (k = 0; k < 4; k++)/* trial vector */
  if (norm <= SMALL_VECTOR)/* bad trial vector; default to radial direction */
    Econ[k] = delta(k, defdir);
  else
    Econ[k] = trial[k];

return;
}

/* 
    check the handedness of a tetrad basis.

    basis is assumed to be in form e^\mu_{(a)} = Econ[a][mu]

    levi_(ijkl) e0^i e1^j e2^k e3^l will be +1 if spatial
    components are right-handed, -1 if left-handed.

    experience suggests that roundoff produces errors of
    order 10^{-12} in the result.

*/





double check_handedness(double Econ[NDIM][NDIM], double Gcov[NDIM][NDIM])
{
int i, j, k, l;
static int firstc = 1;
void set_levi_civita(double levi_civita[NDIM][NDIM][NDIM][NDIM]);
static double levi_civita[NDIM][NDIM][NDIM][NDIM];

if (firstc) {
firstc = 0;
set_levi_civita(levi_civita);
}

 double g = 1;//gdet_func(Gcov);

/* check handedness */
double dot = 0.;
for (i = 0; i < 4; i++)
  for (j = 0; j < 4; j++)
    for (l = 0; l < 4; l++)
      for (k = 0; k < 4; k++) {
    dot += g * levi_civita[i][j][k][l] *
      Econ[0][i] * Econ[1][j] * Econ[2][k] * Econ[3][l];
}

return (dot);
}

/* the completely antisymmetric symbol; not a tensor
   in the coordinate basis */
void set_levi_civita(double levi_civita[NDIM][NDIM][NDIM][NDIM])
{
  int i, j, k, l, n, do_sort, n_perm, val, n_swap;
  int index[NDIM];

  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
	for (l = 0; l < NDIM; l++) {
	  if (i == j || i == k || i == l || j == k || j == l
	      || k == l) {
	    levi_civita[i][j][k][l] = 0;
	  } else {
	    index[0] = i;
	    index[1] = j;
	    index[2] = k;
	    index[3] = l;
	    do_sort = 1;
	    n_perm = 0;
	    while (do_sort) {
	      n_swap = 0;
	      for (n = 0; n < NDIM - 1; n++) {
		if (index[n] > index[n + 1]) {
		  n_perm++;
		  n_swap++;
		  val = index[n];
		  index[n] = index[n + 1];
		  index[n + 1] = val;
		}
	      }
	      do_sort = n_swap;
	    }
	    levi_civita[i][j][k][l] = (n_perm % 2) ? -1 : 1;
	  }
	}

  /* Test levi-civita : */
  /*                              
       for(i=0;i<NDIM;i++)  
       for(j=0;j<NDIM;j++)  
       for(k=0;k<NDIM;k++)  
       for(l=0;l<NDIM;l++) {                         
       fprintf(stdout,"levi-civita[%d%d%d%d] = %d \n", 
       i,j,k,l,levi_civita[i][j][k][l] );                                                
       fflush(stdout); 
       }       
  */
}

