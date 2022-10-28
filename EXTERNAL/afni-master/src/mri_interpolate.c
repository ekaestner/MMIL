#include "mrilib.h"

/*********************************************************************
  * Wrapper functions to interpolate images in various ways.
  * Built on top of functions in mri_genalign_util (thus the GA_).
  * Could be optimized for various sub-cases, but is that
    really worth the effort? We'll see.
*********************************************************************/

/*---------------------------------------------------------------------------*/

void mri_interp_scalar_to_floats_pointset( MRI_IMAGE *fim,
                                           int npp, float *ip, float *jp, float *kp,
                                           int code, float *vv )
{
   MRI_IMAGE *inim=fim ;

ENTRY("mri_interp_scalar_to_floats_pointset") ;

   if( inim == NULL || npp <= 0    ||
       ip   == NULL || jp  == NULL || kp == NULL || vv == NULL ){
     ERROR_message("NULL inputs to mri_interp_scalar_to_floats_pointset()") ;
     EXRETURN ;
   }

   /* convert to float? */

   if( inim->kind != MRI_float ) inim = mri_to_float(fim) ;

   switch( code ){

     case MRI_NN:
       GA_interp_NN     ( inim , npp,ip,jp,kp , vv ) ;
     break ;

     case MRI_LINEAR:
       GA_interp_linear ( inim , npp,ip,jp,kp , vv ) ;
     break ;

     case MRI_CUBIC:
       GA_interp_cubic  ( inim , npp,ip,jp,kp , vv ) ;
     break ;

     case MRI_QUINTIC:
       GA_interp_quintic( inim , npp,ip,jp,kp , vv ) ;
     break ;

     case MRI_WSINC5:
       GA_interp_wsinc5 ( inim , npp,ip,jp,kp , vv ) ;
     break ;

     default:
       ERROR_message("unknown interp code = %d in mri_interp_to_floats_pointset()",code) ;

   }

   if( inim != fim ) mri_free(inim) ;

   EXRETURN ;
}

/*---------------------------------------------------------------------------*/

void mri_interp_to_sametyp_pointset( MRI_IMAGE *fim,
                                     int npp, float *ip, float *jp, float *kp,
                                     int code, void *vv )
{
   register int ii ;

ENTRY("mri_interp_to_sametyp_pointset") ;

   if( fim == NULL || npp <= 0    ||
       ip  == NULL || jp  == NULL || kp == NULL || vv == NULL ){
     ERROR_message("NULL inputs to mri_interp_to_sametyp_pointset()") ;
     EXRETURN ;
   }

   switch( fim->kind ){

     case MRI_float:
       mri_interp_scalar_to_floats_pointset( fim,
                                             npp, ip,jp,kp, code,(float *)vv ) ;
     break ;

     case MRI_int:{
       int *pp = (int *)vv ;
       float *ff = (float *)malloc(sizeof(float)*npp) ;
       mri_interp_scalar_to_floats_pointset( fim ,
                                             npp , ip,jp,kp , code,ff ) ;
       for( ii=0 ; ii < npp ; ii++ ) pp[ii] = (int)rintf(ff[ii]) ;
       free(ff) ;
     }
     break ;

     case MRI_short:{
       short *pp = (short *)vv ;
       float *ff = (float *)malloc(sizeof(float)*npp) ;
       mri_interp_scalar_to_floats_pointset( fim ,
                                             npp , ip,jp,kp , code,ff ) ;
       for( ii=0 ; ii < npp ; ii++ ) pp[ii] = SHORTIZE(ff[ii]) ;
       free(ff) ;
     }
     break ;

     case MRI_byte:{
       byte *pp = (byte *)vv ;
       float *ff = (float *)malloc(sizeof(float)*npp) ;
       mri_interp_scalar_to_floats_pointset( fim ,
                                             npp , ip,jp,kp , code,ff ) ;
       for( ii=0 ; ii < npp ; ii++ ) pp[ii] = BYTEIZE(ff[ii]) ;
       free(ff) ;
     }
     break ;

     case MRI_rgb:{
       byte *pval = (byte *)vv ;
       MRI_IMARR *ppar ;
       float *rr = (float *)malloc(sizeof(float)*npp) ;
       float *gg = (float *)malloc(sizeof(float)*npp) ;
       float *bb = (float *)malloc(sizeof(float)*npp) ;
       ppar = mri_rgb_to_3float(fim) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,0) ,
                                             npp , ip,jp,kp , code, rr ) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,1) ,
                                             npp , ip,jp,kp , code, gg ) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,2) ,
                                             npp , ip,jp,kp , code, bb ) ;
       for( ii=0 ; ii < npp ; ii++ ){
         pval[3*ii+0] = BYTEIZE(rr[ii]) ;
         pval[3*ii+1] = BYTEIZE(gg[ii]) ;
         pval[3*ii+2] = BYTEIZE(bb[ii]) ;
       }
       free(bb) ; free(gg) ; free(rr) ;
     }
     break ;

     case MRI_rgba:{
       rgba *pval = (rgba *)vv ;
       MRI_IMARR *ppar ;
       float *rr = (float *)malloc(sizeof(float)*npp) ;
       float *gg = (float *)malloc(sizeof(float)*npp) ;
       float *bb = (float *)malloc(sizeof(float)*npp) ;
       float *aa = (float *)malloc(sizeof(float)*npp) ;
       ppar = mri_rgba_to_4float(fim) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,0) ,
                                             npp , ip,jp,kp , code, rr ) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,1) ,
                                             npp , ip,jp,kp , code, gg ) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,2) ,
                                             npp , ip,jp,kp , code, bb ) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,3) ,
                                             npp , ip,jp,kp , code, aa ) ;
       for( ii=0 ; ii < npp ; ii++ ){
         pval[ii].r = BYTEIZE(rr[ii]) ;
         pval[ii].g = BYTEIZE(gg[ii]) ;
         pval[ii].b = BYTEIZE(bb[ii]) ;
         pval[ii].a = BYTEIZE(rr[ii]) ;
       }
       free(aa) ; free(bb) ; free(gg) ; free(rr) ;
     }
     break ;

     case MRI_complex:{
       complex *pval = (complex *)vv ;
       MRI_IMARR *ppar ;
       float *xx = (float *)malloc(sizeof(float)*npp) ;
       float *yy = (float *)malloc(sizeof(float)*npp) ;
       ppar =  mri_complex_to_pair(fim) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,0) ,
                                             npp , ip,jp,kp , code, xx ) ;
       mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,1) ,
                                             npp , ip,jp,kp , code, yy ) ;
       for( ii=0 ; ii < npp ; ii++ ){
         pval[ii].r = xx[ii] ;
         pval[ii].i = yy[ii] ;
       }
       free(yy) ; free(xx) ;
     }
     break ;

     case MRI_fvect:{
       MRI_IMARR *ppar ; int nff, kk,ii,qq ; float **ffar, *pval=(float *)vv ;
       ppar = mri_fvect_to_imarr(fim) ;
       nff  = IMARR_COUNT(ppar) ;
       ffar = (float **)malloc(sizeof(float *)*nff) ;
       for( kk=0 ; kk < nff ; kk++ ){
         ffar[kk] = (float *)malloc(sizeof(float)*npp) ;
         mri_interp_scalar_to_floats_pointset( IMARR_SUBIM(ppar,kk) ,
                                               npp , ip,jp,kp , code, ffar[kk] ) ;
       }
       for( qq=ii=0 ; ii < npp ; ii++ ){
         for( kk=0 ; kk < nff ; kk++ ) pval[qq++] = ffar[kk][ii] ;
       }
       for( kk=0 ; kk < nff ; kk++ ) free(ffar[kk]) ;
       free(ffar) ;
     }
     break ;

     default:
       ERROR_message("unknown image kind = %d in mri_interp_to_sametyp_pointset()",fim->kind) ;
     break ;

   }

   EXRETURN ;
}

/*---------------------------------------------------------------------------*/

MRI_IMAGE * mri_interp_to_floats_block( MRI_IMAGE *fim ,
                                        mat44 ijk_to_ijk , int code ,
                                        int ibot , int itop ,
                                        int jbot , int jtop ,
                                        int kbot , int ktop  )
{
   MRI_IMAGE *outim ; float *outar ;
   int nxout, nyout, nzout, nxyz , ii,jj,kk,qq ;
   float *ip , *jp , *kp ;

ENTRY("mri_interp_to_floats_block") ;

   if( fim == NULL || !ISVALID_MAT44(ijk_to_ijk) ||
       ibot > itop || jbot > jtop || kbot > ktop   ) RETURN(NULL) ;

   nxout = itop-ibot+1 ;
   nyout = jtop-jbot+1 ;
   nzout = ktop-kbot+1 ; nxyz = nxout * nyout * nzout ;

   outim = mri_new_vol( nxout, nyout, nzout , MRI_float ) ;
   outar = MRI_FLOAT_PTR(outim) ;

   ip = (float *)malloc(sizeof(float)*nxyz) ;
   jp = (float *)malloc(sizeof(float)*nxyz) ;
   kp = (float *)malloc(sizeof(float)*nxyz) ;

   for( qq=kk=kbot ; kk <= ktop ; kk++ ){
    for( jj=jbot ; jj <= jtop ; jj++ ){
     for( ii=ibot ; ii <= itop ; ii++,qq++ ){
       MAT44_VEC( ijk_to_ijk , ii,jj,kk , ip[qq],jp[qq],kp[qq] ) ;
   }}}

   mri_interp_scalar_to_floats_pointset( fim, nxyz, ip,jp,kp, code, outar ) ;

   free(kp); free(jp); free(ip) ;

   RETURN(outim) ;
}

/*---------------------------------------------------------------------------*/

MRI_IMAGE * mri_interp_to_vectyp_block( MRI_IMAGE *fim ,
                                        mat44 ijk_to_ijk , int code ,
                                        int ibot , int itop ,
                                        int jbot , int jtop ,
                                        int kbot , int ktop  )
{
   MRI_IMAGE *outim ; float *outar ;
   int nxout, nyout, nzout, nxyz , ii,jj,kk,qq ;
   float *ip , *jp , *kp ;

ENTRY("mri_interp_to_vectyp_block") ;

   if( fim == NULL || !ISVALID_MAT44(ijk_to_ijk) ||
       ibot > itop || jbot > jtop || kbot > ktop   ) RETURN(NULL) ;

   if( ! ISVECTIM(fim) ){
     outim = mri_interp_to_floats_block( fim,ijk_to_ijk,code,ibot,itop,jbot,jtop,kbot,ktop) ;
     RETURN(outim) ;
   }

   nxout = itop-ibot+1 ;
   nyout = jtop-jbot+1 ;
   nzout = ktop-kbot+1 ; nxyz = nxout * nyout * nzout ;

   if( fim->kind == MRI_fvect ){
     outim = mri_new_fvectim( nxout , nyout , nzout , fim->vdim ) ;
   } else {
     outim = mri_new_vol( nxout, nyout, nzout , fim->kind ) ;
   }
   outar = mri_data_pointer(outim) ;

   ip = (float *)malloc(sizeof(float)*nxyz) ;
   jp = (float *)malloc(sizeof(float)*nxyz) ;
   kp = (float *)malloc(sizeof(float)*nxyz) ;

   for( qq=kk=kbot ; kk <= ktop ; kk++ ){
    for( jj=jbot ; jj <= jtop ; jj++ ){
     for( ii=ibot ; ii <= itop ; ii++,qq++ ){
       MAT44_VEC( ijk_to_ijk , ii,jj,kk , ip[qq],jp[qq],kp[qq] ) ;
   }}}

   mri_interp_to_sametyp_pointset( fim, nxyz, ip,jp,kp, code, outar ) ;

   free(kp); free(jp); free(ip) ;

   RETURN(outim) ;
}

/*---------------------------------------------------------------------------*/
/* Interpolate 1 point from each image, and return a 1D float array.
   This function is meant to provide a time-series array from an
   interpolated grid point, for use in the AFNI graph viewer.
   It is rather brute force, but as computers are so fast these days .... */
/*---------------------------------------------------------------------------*/

MRI_IMAGE * imarr_interp_scalar_to_floats_onepoint(
                    MRI_IMARR *imar,
                    float iq , float jq , float kq , int code )
{
   MRI_IMAGE *inim, *fim, *outim ; int nim, kk ;
   float *outar , vv[1] , ip,jp,kp ;

ENTRY("imarr_interp_scalar_to_floats_onepoint") ;

   if( imar == NULL || IMARR_COUNT(imar) < 1 ){
     ERROR_message("NULL inputs to imarr_interp_scalar_to_floats_onepoint()") ;
     RETURN(NULL) ;
   }

   nim = IMARR_COUNT(imar) ;
   outim = mri_new( nim , 1 , MRI_float ) ;
   outar = MRI_FLOAT_PTR(outim) ;
   ip = iq ; jp = jq ; kp = kq ;

   for( kk=0 ; kk < nim ; kk++ ){

     inim = IMARR_SUBIM(imar,kk) ;
     if( inim == NULL ) continue ;

     /* convert to float? */
     if( inim->kind != MRI_float ) inim = mri_to_float(fim) ;

     switch( code ){
       case MRI_NN:
         GA_interp_NN     ( inim , 1,&ip,&jp,&kp , vv ) ;
       break ;
       case MRI_LINEAR:
         GA_interp_linear ( inim , 1,&ip,&jp,&kp , vv ) ;
       break ;
       case MRI_CUBIC:
         GA_interp_cubic  ( inim , 1,&ip,&jp,&kp , vv ) ;
       break ;
       case MRI_QUINTIC:
         GA_interp_quintic( inim , 1,&ip,&jp,&kp , vv ) ;
       break ;
       case MRI_WSINC5:
         GA_interp_wsinc5 ( inim , 1,&ip,&jp,&kp , vv ) ;
       break ;
     }
     outar[kk] = vv[0] ;
     if( inim != fim ) mri_free(inim) ;
   }

   RETURN(outim) ;
}
