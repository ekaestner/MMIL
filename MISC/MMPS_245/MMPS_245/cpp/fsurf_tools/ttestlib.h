#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef Linux
#include <sys/types.h>
#endif

/* for t to p interchange from mri_stats.c (afni) */
#define ZERO 0.0
#define ONE  1.0
#define ACU  1.0e-15
#define SAE   -15.0
#define TWO     2.0
#define THREE   3.0
#define FOUR    4.0
#define FIVE    5.0
#define SIX     6.0
#ifndef MAX
#define MAX(a,b) (((a)<(b)) ? (b) : (a))
#endif
#ifndef MIN
#define MIN(a,b) (((a)>(b)) ? (b) : (a))
#endif

/* function prototypes */
double student_p2t( double pp , double dof );
double student_t2p( double tt , double dof );
double lnbeta( double p , double q );
double incbeta_inverse( double alpha , double p , double q , double beta );
double incbeta( double x , double p , double q , double beta );

