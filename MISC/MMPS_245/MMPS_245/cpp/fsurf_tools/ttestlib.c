/*
   functions for calculating p-values for t-stats copied from afni's mri_stats.c 
*/

#include "ttestlib.h"

double student_p2t( double pp , double dof )
/* converts to t value from p value */
{
   double bb , binv , tt ;

   if( pp <= 0.0 || pp >= 0.9999999 || dof < 1.0 ) return 0.0 ;

   bb   = lnbeta( 0.5*dof , 0.5 ) ;
   binv = incbeta_inverse( pp, 0.5*dof , 0.5 , bb ) ;
   tt   = sqrt( dof*(1.0/binv-1.0) ) ;
   return tt ;
}

double student_t2p( double tt , double dof )
/* converts to p value from t value */
{
   double bb , xx , pp ;

   if( tt <= 0.0 || dof < 1.0 ) return 1.0 ;

   bb = lnbeta( 0.5*dof , 0.5 ) ;
   xx = dof/(dof + tt*tt) ;
   pp = incbeta( xx , 0.5*dof , 0.5 , bb ) ;
   return pp ;
}

double lnbeta( double p , double q )
/* compute log of complete beta function with Unix math library's log gamma function */
{
   return (gamma(p) + gamma(q) - gamma(p+q)) ;
}

double incbeta( double x , double p , double q , double beta )
/* copied from afni's mri_stats.c */
/* computes incomplete beta function ratio for arguments x between zero and one,
   p and q positive. log of complete beta function, beta, is assumed to be known
   translated from fortran (appl. statist. (1973), vol.22, no.3, algorithm 63)
   -- goto hell! */
{
   double betain , psq , cx , xx,pp,qq , term,ai , temp , rx ;
   int indx , ns ;

   if( p <= ZERO || q <= ZERO ) return -1.0 ;  /* error! */

   if( x <= ZERO ) return ZERO ;
   if( x >= ONE  ) return ONE ;

   /**  change tail if necessary and determine s **/

   psq = p+q ;
   cx  = ONE-x ;
   if(  p < psq*x ){
      xx   = cx ;
      cx   = x ;
      pp   = q ;
      qq   = p ;
      indx = 1 ;
   } else {
      xx   = x ;
      pp   = p ;
      qq   = q ;
      indx = 0 ;
   }

   term   = ONE ;
   ai     = ONE ;
   betain = ONE ;
   ns     = (int)(qq + cx*psq);

   /** use soper's reduction formulae **/

      rx = xx/cx ;

lab3:
      temp = qq-ai ;
      if(ns == 0) rx = xx ;

lab4:
      term   = term*temp*rx/(pp+ai) ;
      betain = betain+term ;
      temp   = fabs(term) ;
      if(temp <= ACU && temp <= ACU*betain) goto lab5 ;

      ai = ai+ONE ;
      ns = ns-1 ;
      if(ns >= 0) goto lab3 ;
      temp = psq ;
      psq  = psq+ONE ;
      goto lab4 ;

lab5:
      betain = betain*exp(pp*log(xx)+(qq-ONE)*log(cx)-beta)/pp ;
      if(indx) betain=ONE-betain ;

   return betain ;
}


double incbeta_inverse( double alpha , double p , double q , double beta )
/* computes inverse of the incomplete beta function ratio for given positive
   values of the arguments p and q, alpha between zero and one. log of complete
   beta function, beta, is assumed to be known
   translated from fortran (appl. statist. (1977), vol.26, no.1, algorithm 109)
   -- goto hell! */
{
#if 0
   int indx , iex ;
#else
   int indx;
#endif
   double fpu , xinbta , a,pp,qq, r,y,t,s,h,w , acu ,
          yprev,prev,sq , g,adj,tx,xin ;

   fpu = pow(10.0,SAE) ;

   if( p <= ZERO || q <= ZERO || alpha < ZERO || alpha > ONE ) return -1.0 ;

   if( alpha == ZERO ) return ZERO ;
   if( alpha == ONE  ) return ONE ;

   /** change tail if necessary **/

   if( alpha > 0.5 ){
      a    = ONE-alpha ;
      pp   = q ;
      qq   = p ;
      indx = 1 ;
    } else {
      a    = alpha ;
      pp   = p ;
      qq   = q ;
      indx = 0 ;
   }

   /** calculate the initial approximation **/

#if 0
lab2:
#endif
     r = sqrt(-log(a*a)) ;
     y = r - (2.30753 + 0.27061*r) / (ONE+(0.99229+0.04481*r)*r) ;
     if(pp > ONE && qq > ONE) goto lab5 ;

     r = qq+qq ;
     t = ONE/(9.0*qq) ;
     t = r * pow( (ONE-t+y*sqrt(t)) , 3.0 ) ;
     if( t <= ZERO ) goto lab3 ;

     t = (FOUR*pp+r-TWO)/t ;
     if( t <= ONE ) goto lab4 ;

     xinbta = ONE-TWO/(t+ONE) ; goto lab6 ;

lab3:
     xinbta = ONE-exp((log((ONE-a)*qq)+beta)/qq) ; goto lab6 ;

lab4:
     xinbta = exp((log(a*pp)+beta)/pp) ; goto lab6 ;

lab5:
     r = (y*y-THREE)/SIX ;
     s = ONE/(pp+pp-ONE) ;
     t = ONE/(qq+qq-ONE) ;
     h = TWO/(s+t) ;
     w = y*sqrt(h+r)/h-(t-s)*(r+FIVE/SIX-TWO/(THREE*h)) ;
     xinbta = pp/(pp+qq*exp(w+w)) ;

     /** solve for x by a modified newton-raphson method **/

lab6:
    r     = ONE-pp ;
    t     = ONE-qq ;
    yprev = ZERO ;
    sq    = ONE ;
    prev  = ONE ;
    if(xinbta < 0.0001) xinbta = 0.0001 ;
    if(xinbta > 0.9999) xinbta = 0.9999 ;

#if 0
    iex = -5.0 / (pp*pp) - 1.0/(a*a) - 13.0 ; if( iex < SAE ) iex = SAE ;
    acu = pow(10.0,iex) ;
#else
    acu = fpu ;
#endif

lab7:
      y = incbeta( xinbta , pp,qq,beta ) ;
      if( y < ZERO ) return -1.0 ;
      xin = xinbta ;
      y = (y-a)*exp(beta+r*log(xin)+t*log(ONE-xin)) ;
      if(y*yprev <= ZERO) prev = MAX(sq, fpu) ;
      g = ONE ;

lab9:
      adj = g*y ;
      sq = adj*adj ;
      if(sq >= prev) goto lab10 ;
      tx = xinbta-adj ;
      if(tx >= ZERO && tx <= ONE) goto lab11 ;

lab10:
      g = g/THREE ; goto lab9 ;

lab11:
      if(tx == ZERO  || tx == ONE ) goto lab10 ;
      if(prev <= acu || y*y <= acu || fabs(xinbta-tx) < fpu) goto lab12 ;
      xinbta = tx ;
      yprev = y ;
      goto lab7 ;

lab12:
      xinbta = tx ;
      if (indx) xinbta = ONE-xinbta ;
#if 0
   printf("alpha = %g  incbeta = %g\n",alpha, incbeta(xinbta,p,q,beta) );
#endif
      return xinbta ;
}

