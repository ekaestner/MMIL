#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <malloc.h>

float *compute_UA(float *U, float *A, int N, int M, int P);
float normal_v(float *v, int n);
float *compute_UT(float *U, int M, int N);
float *dot_UA(float *U, float *A, int N);

/*******************************************************************/
/* This function multiplies the NxM matrix U and the MxP matrix A. */
/* The resulting NxP matrix is returned in UA.                     */
/*******************************************************************/
float *compute_UA(float *U, float *A, int N, int M, int P)
{
   int 		i, j, k;
   double 	*Ud, *Ad;
   float 	*UA;

   Ud = (double *)calloc(N*P, sizeof(double));
   Ad = (double *)calloc(N*P, sizeof(double));
   UA=(float *)calloc(N*P,sizeof(float));

   for(j=0; j<N; j++)
      for(i=0; i<P; i++)
      {
         for(k=0; k<M; k++)
         {
                UA[j*P+i] += U[j*M+k]*A[k*P+i];
         }
      }
   free(Ad);
   free(Ud);
   return(UA);
}


/*******************************************************************/
/* This function multiplies the transfer matrix of MxN matrix U    */
/* The resulting NxM matrix is returned in UT.                     */
/*******************************************************************/
float *compute_UT(float *U, int M, int N)
{
   int 		i, j, k;
   float 	*UT;

   UT=(float *)calloc(N*M,sizeof(float));

   for(j=0; j<M; j++)
        for(k=0; k<N; k++)
        {
              UT[k*M+j] = U[j*N+k];
        }
   return(UT);
}

/*******************************************************************/
/* This function multiplies the transfer matrix of MxN matrix U    */
/* The resulting NxM matrix is returned in UT.                     */
/*******************************************************************/
float normal_v(float *v, int n)
{
        int i;

        float tmp_norm = 0;
        for(i = 0; i < n; i++)
                tmp_norm = tmp_norm + v[i] * v[i];
        tmp_norm = sqrt(tmp_norm);
        return (tmp_norm);
}

/*******************************************************************/
/* This function multiplies the transfer matrix of MxN matrix U    */
/* Returns the scalar product of U and V vectors.                     */
/*******************************************************************/
float dot_UV(float *U, float *V, int N)
{
	int i;

	float slice_normal = 0;	
	for(i=0; i<N; i++)
		slice_normal = slice_normal + U[i] * V[i];
	return slice_normal;
}
