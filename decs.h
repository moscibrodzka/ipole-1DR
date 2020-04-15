#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>


//one can chose integration scheme
#define INT_FULL  1 //full integration step
#define INT_SPLIT 0 //split: 1/2rot + ea + 1/2rot (from development version of ipole)
//chose eDF
#define THERMAL 1
#define POWERL  0

#define NDIM 4
#define SMALL 1.e-40
#define S2                              (1.41421356237310       ) /* sqrt(2) */
#define S3                              (1.73205080756888       ) /* sqrt(3) */
#define EE                              (4.80320680e-10         ) /* electron charge */
#define CL                              (2.99792458e10          ) /* speed of light */
#define ME                              (9.1093826e-28          ) /* electron mass */
#define MP                              (1.67262171e-24         ) /* proton mass */
#define MN                              (1.67492728e-24         ) /* neutron mass */
#define AMU                             (1.66053886e-24         ) /* atomic mass unit */
#define HPL                             (6.6260693e-27          ) /* Planck constant */
#define HPL_MECL2                       (8.09e-21               ) /* hpl/m_ec^2*/
#define HBAR                            (HPL/(2.*M_PI)          ) /* Planck's consant / 2pi */
#define KBOL                            (1.3806505e-16          ) /* Boltzmann constant */
#define GNEWT                           (6.6742e-8              ) /* Gravitational constant */
#define SIG                             (5.670400e-5            ) /* Stefan-Boltzmann constant */
#define RGAS                            (8.3143e7               )       /* erg K^-1 mole^-1: ideal gas const */
#define EV                              (1.60217653e-12         ) /* electron volt in erg */
#define SIGMA_THOMSON                   (0.665245873e-24        ) /* Thomson cross section in cm^2 */
#define JY                              (1.e-23                 ) /* Jansky (flux/freq. unit) in cgs */


