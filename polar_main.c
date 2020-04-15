/* put together by Monika Moscibrodzka,
 last update: April 14 2020 */

#include "decs.h"

// here define your model for plasma in a slab
void main(int argc, char *argv[]){

  int i,k,l;
  double nu,Thetae,Ne,B,theta,Bx,By,Bz,Vx,Vy,Vz;
  double jI,jQ,jU,jV,aI,aQ,aU,aV,rQ,rU,rV;
  double SI0,SQ0,SU0,SV0,SI,SQ,SU,SV;
  double ucov[NDIM],bcov[NDIM];
  double gcov[NDIM][NDIM],Econ[NDIM][NDIM],Ecov[NDIM][NDIM];
  double complex N_coord[NDIM][NDIM],N_tetrad[NDIM][NDIM];

  /* metric tensor is fixed to Minkowski spacetime and initial N is fixed = do not change*/
  for(k=0;k<4;k++)for(l=0;l<4;l++) gcov[k][l] = 0.;
  gcov[0][0]=-1.; gcov[1][1]=1.; gcov[2][2]=1.; gcov[3][3]=1.;
  /* initial conditions for complex coherency tensor carring stokes parameters */
  for(k=0;k<4;k++)for(l=0;l<4;l++) N_coord[k][l] = 0.0 + I*0.0;

  double freqcgs;//frequency of photons
  double tmin;//global time starting point
  sscanf(argv[1], "%lf", &tmin);
  sscanf(argv[2], "%lf", &freqcgs);  
  
  //these paramters can be changed
  //length unit
  double L_unit=9.0e+14; 
  //integration zone
  double Xmin=1.0; double Xmax=20.0;

  /* initial conditions for K-photon wave vector and X position of photon */
  double K[NDIM]={1*freqcgs * HPL / (ME * CL * CL), 1*freqcgs * HPL / (ME * CL * CL),0.0,0.0};
  double X[NDIM]={tmin,Xmin,0.0,0.0};
 
  double lam = 0.0;
  double dlam = 0.001/K[1];  // fixed step size
  double tauF = 0.0;
  double tauabs = 0.0;

  while(X[1] < Xmax) {
  
      /* advance X, allowed to move only along x axis */
      for(i=0;i<4;i++) X[i] += K[i]*dlam;
      lam += dlam; /* in dimentioless units, counting from zero */
      
      /* define velocity (units c), could be a funtion of coordinate */
      /* here static */
      Vx=0.0;Vy=0.0;Vz=0.0;
      /* total velocity and Lorentz factor in flat space-time */
      double v=sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
      double gammaL=1./sqrt(1.-v*v);
      double ucon[NDIM]={gammaL,gammaL*Vx,gammaL*Vy,gammaL*Vz};
      lower(ucon, gcov, ucov);   

      /* define B field (units Gauss), could be a funtion of coordinate */
      Bx=-sqrt(300.)/X[1]; By=0.0; Bz=sqrt(900.)/X[1];
      double bcon[NDIM]={ Bx*ucov[1] + By*ucov[2] + Bz*ucov[3],
			  (Bx + bcon[0]*ucon[1])/ucon[0],
			  (By + bcon[0]*ucon[2])/ucon[0],
			  (Bz + bcon[0]*ucon[3])/ucon[0]};
      lower(bcon, gcov, bcov);   

      
      /* calculate angle of k and b in the fluid frame */
      B= sqrt(bcon[0]*bcov[0]+ bcon[1]*bcov[1]+ bcon[2]*bcov[2]+ bcon[3]*bcov[3]);
      double kk = fabs(K[0]*ucov[0] + K[1]*ucov[1] + K[2]*ucov[2] + K[3]*ucov[3]);
      double mu = (K[0]*bcov[0] + K[1]*bcov[1] + K[2]*bcov[2] + K[3]*bcov[3]) / (kk * B);
      theta=acos(mu);

      /* define plasma density, could be a function of coordinate */
      Ne=2.e5*pow(X[1],-2);

      /* defince temperature of thermal particles in the plasma, could be a funtion of coordinate */
      Thetae=3.e11*KBOL/ME/CL/CL;
      
      /* photon frequency measured in the fluid frame */
      nu = -(K[0]*ucov[0] + K[1]*ucov[1] + K[2]*ucov[2] + K[3]*ucov[3]) * ME * CL * CL / HPL;

      make_plasma_tetrad(ucon, K, bcon, gcov, Econ, Ecov);
      complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
      tensor_to_stokes(N_tetrad, &SI0, &SQ0, &SU0, &SV0);

      if(THERMAL){
	  jar_calc_th(nu,Ne,B,Thetae,theta,
		      &jI,&jQ,&jU,&jV,
		      &aI,&aQ,&aU,&aV,
		      &rQ,&rU,&rV);
      }
      if(POWERL){
	  double gmin=10.;
	  double gmax=1000.;
	  double p=1.5;
      
	  jar_calc_pl(nu,Ne,B,gmin,gmax,p,theta,
		      &jI,&jQ,&jU,&jV,
		      &aI,&aQ,&aU,&aV,
		      &rQ,&rU,&rV);
	  
      }
      
      /* integrate Stokes */
      integrate_Stokes(dlam,L_unit,
		       jI,jQ,jU,jV,
		       aI,aQ,aU,aV,
		       rQ,rU,rV,
		       SI0,SQ0,SU0,SV0,&SI,&SQ,&SU,&SV);
      
      stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
      complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);

      /* calculate optical and Faraday thickess */
      tauF   += fabs(rV)*dlam*HPL*L_unit/(ME*CL*CL);
      tauabs += (aI)*dlam*HPL*L_unit/(ME*CL*CL);
  
      /* build detector (attachaed to Cartesian coordinates) */
      double Ucam[NDIM]={1.0,0.0,0.0,0.0};
      double trial1[NDIM]={0.0,1.0,0.0,0.0};
      double trial2[NDIM]={0.0,0.0,1.0,0.0};

      make_plasma_tetrad(Ucam, trial1, trial2, gcov, Econ, Ecov);
      complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
      tensor_to_stokes(N_tetrad, &SI0, &SQ0, &SU0, &SV0);

      /* print out Stokes parameters as a function of X[1]=x */
      /*
      fprintf(stdout,"%g  %g %g %g %g   %g %lf \n",
	      X[1],
	      SI0*pow(freqcgs,3),
	      SQ0*pow(freqcgs,3),
	      SU0*pow(freqcgs,3),
	      SV0*pow(freqcgs,3),
	      tauabs,tauF);
      */
      

/* prints fluid variables */
/*
  fprintf(stdout,"plasma: x=%g ne=%g thetae=%g B=%g \n",X[1],Ne,Thetae,B);
*/
      
  }

  /*print out final Stokes parameters at the Xmax */
  fprintf(stdout,"%g %g  %g %g %g %g   %g %lf \n",X[0],X[1],
	  SI0*pow(freqcgs,3),
	  SQ0*pow(freqcgs,3),
	  SU0*pow(freqcgs,3),
	  SV0*pow(freqcgs,3),
	  tauabs,tauF);

  /* END */
  
}
