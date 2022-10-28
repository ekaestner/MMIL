#include "mrilib.h"

#include "mri_pvmap.c"

/*** Functions for doing 4D image to grayplot conversion ***/

/*--------------------------------------------------------------------------*/

/*----- levels = classes of voxels that are grouped together -----*/

static int lev_num ;
static int lev_siz[256] ;

/*----- ordering of voxel time series vertically, within each level -----*/

#define ORDER_IJK  0
#define ORDER_PV   1
#define ORDER_PEEL 2
#define ORDER_LJ   3

static int ordering = 0 ;

static void grayplot_order_by_pvmap( int yesno ){
  ordering = (yesno) ? ORDER_PV : 0 ;
}
static void grayplot_order_by_peels( int yesno ){
  ordering = (yesno) ? ORDER_PEEL : 0 ;
}
static void grayplot_order_by_ijk( int yesno ){
  ordering = 0 ;
}
static void grayplot_order_by_lj( int yesno ){   /* 05 Feb 2020 */
  ordering = (yesno) ? ORDER_LJ : 0 ;
}

static int grayplot_is_ordered_by_peels( void ){ /* [18 Jun 2019 rickr] */
  return (ordering == ORDER_PEEL);
}
static int grayplot_is_ordered_by_pvmap( void ){
  return (ordering == ORDER_PV);
}
static int grayplot_is_ordered_by_lj( void ){    /* 05 Feb 2020 */
  return (ordering == ORDER_LJ);
}

/*----- how time series are normalized [30 Jul 2018] -----*/

static int norming = 0 ;

#define NORM_RMS    0
#define NORM_NONE   1
#define NORM_MAXABS 2

static void grayplot_norming_RMS(void){
  norming = NORM_RMS ;
}

static void grayplot_norming_none(void){
  norming = NORM_NONE ;
}

static void grayplot_norming_maxabs(void){
  norming = NORM_MAXABS ;
}

/*----- grayplot range for 1 unit of data -----*/

static float grange = 0.0f ;

static void grayplot_set_range( float val ){
  grange = MAX(val,0.0f) ;
}

/*----- grayplot raw_with_bounds range -----*/

static float raw_range_bot =  1.0f ;
static float raw_range_top = -1.0f ;
static int   DO_RAW_RANGE  = 0;

static void grayplot_set_raw_range( float bnd_bot, float bnd_top){
   raw_range_bot = bnd_bot;
   raw_range_top = bnd_top;
   DO_RAW_RANGE  = 1;
}

/*----- percentage? -----*/

static int do_percent = 0 ;
static void grayplot_set_percent(void){ do_percent = 1 ; }

/*----- use old resampling method? -----*/

static int use_old_resam = 0 ;
static void grayplot_set_use_old_resam(void){ use_old_resam = 1 ; }

/*--------------------------------------------------------------------------*/
/* Convert to a vectim, using mmask to pick out the levels
   of different voxels to keep together (e.g., GM, WM, CSF)
*//*------------------------------------------------------------------------*/

static MRI_vectim * THD_dset_grayplot_prep( THD_3dim_dataset *dset ,
                                            byte *mmask ,
                                            int polort , float fwhm )
{
   Aint cmval , nmask , nxyz , nts ;
   byte mm,*tmask ; Aint ii,jj,kk ;
   Aint nvim ; MRI_vectim **vim , *vout ;
   float *tsar , *fit=NULL , fac ;

   /* check inputs for plausibility */

ENTRY("THD_dset_grayplot_prep") ;

   lev_num = 0 ;
   if( !ISVALID_DSET(dset) || mmask == NULL ) RETURN(NULL) ;
   nts = DSET_NVALS(dset) ;
   if( nts < 19 ) RETURN(NULL) ;

   nxyz = DSET_NVOX(dset) ;
   nmask = THD_countmask( nxyz , mmask ) ;
   if( nmask < 19 ) RETURN(NULL) ;

   if( do_percent && polort < 0 ) polort = 0 ;

   /* find voxels in each mmask 'level' */

   nvim  = 0 ; vim = NULL ;
   tmask = (byte *)malloc(sizeof(byte)*nxyz) ;           /* binary mask */

   if( polort >= 0 ){
     fit = (float *)malloc(sizeof(float)*(polort+1)) ; /* setup polorts */
   }

   for( ii=1 ; ii <= 255 ; ii++ ){            /* loop over mmask levels */

     /* create binary mask for this level */
     mm = (byte)ii ;
     memset( tmask , 0 , sizeof(byte)*nxyz ) ;
     for( cmval=jj=0 ; jj < nxyz ; jj++ ){
       if( mmask[jj] == mm ){ cmval++ ; tmask[jj] = 1 ; }
     }

     /* extract voxels at this level into a vectim */
     if( cmval > 9 ){
       vim = (MRI_vectim **)realloc( vim , sizeof(MRI_vectim *)*(nvim+1) ) ;
       vim[nvim] = THD_dset_to_vectim( dset , tmask , 0 ) ;
       if( polort >= 0 ){
#if 1
ININFO_message("  Detrending time series with polort=%d",polort) ;
#endif
         for( jj=0 ; jj < vim[nvim]->nvec ; jj++ ){          /* detrend */
           tsar = VECTIM_PTR( vim[nvim] , jj ) ; fit[0] = 0.0f ;
           THD_generic_detrend_LSQ( nts,tsar , polort , 0,NULL,fit ) ;
           if( do_percent && fit[0] > 0.0f ){
             fac = 100.0f / fit[0] ;
             for( kk=0 ; kk < nts ; kk++ ) tsar[kk] *= fac ;
           }
         }
       }

       switch( norming ){  /* normalize */

         case NORM_RMS:
         default:{
#if 1
ININFO_message("  RMS norming time series") ;
#endif
           for( jj=0 ; jj < vim[nvim]->nvec ; jj++ ){ /* set RMS = 1 */
             tsar = VECTIM_PTR( vim[nvim] , jj ) ;
             THD_normRMS( nts , tsar ) ;
           }
         }
         break ;

         case NORM_NONE: /*nada*/
         break ;

         case NORM_MAXABS:{ /* scale so max(abs(x)) = 1 */
           float mab,val ;
#if 1
ININFO_message("  MAXABS norming time series") ;
#endif
           for( jj=0 ; jj < vim[nvim]->nvec ; jj++ ){ /* set RMS = 1 */
             tsar = VECTIM_PTR( vim[nvim] , jj ) ;
             mab = fabsf(tsar[0]) ;
             for( kk=1 ; kk < nts ; kk++){
               val = fabsf(tsar[kk]) ; if( val > mab ) mab = val ;
             }
             if( mab > 0.0f ){
               val = 1.0f / mab ;
               for( kk=0 ; kk < nts ; kk++ ) tsar[kk] *= val ;
             }
           }
         }
         break ;

       }

       if( fwhm > 0.0f ){ /* spatially blur inside this level */
#if 1
ININFO_message("  Blurring (in mask) volumes fwhm=%.3g",fwhm) ;
#endif
         mri_blur3D_vectim( vim[nvim] , fwhm ) ;
       }

       /* re-order spatially, as ordered */

       switch( ordering ){
         default: break ;  /* == IJK */

         case ORDER_PV:{   /* == by coherence with 1st 2 principal vectors */
           MRI_IMAGE *tim , *pim ;
           Aint      *kim ;
           float     *tar , *par ;
#if 1
ININFO_message("  Computing PV order for mask partition #%d - %d voxels", ii,cmval) ;
#endif
           /* copy data vectors into temporary image,
              along with originating index of each vector */
           tim = mri_new(nts,cmval,MRI_float) ; tar = MRI_FLOAT_PTR(tim) ;
           kim = (Aint *)malloc(sizeof(Aint)*cmval) ;
           for( jj=0 ; jj < cmval ; jj++ ){
             memcpy( tar+jj*nts, VECTIM_PTR(vim[nvim],jj), sizeof(float)*nts ) ;
             kim[jj] = jj ;  /* source index */
           }
           /* make the PVmap */
#if 1
ININFO_message("  Computing mri_vec_to_pvmap") ;
#endif
           pim = mri_vec_to_pvmap(tim) ; par = MRI_FLOAT_PTR(pim) ;
           /* sort so largest are first, keeping track of whence they came */
#if 1
ININFO_message("  Sorting mri_vec_to_pvmap") ;
#endif
           for( jj=0 ; jj < cmval ; jj++ ) par[jj] = -par[jj] ;
#if Aintsize == 64
           qsort_floatint64_t( cmval , par , kim ) ;
#else
           qsort_floatint( cmval , par , kim ) ;
#endif
#if 1
ININFO_message("  Putting sorted vectors back into vectim") ;
#endif

           /* copy from temp image back to vectim, in the right order */
           for( jj=0 ; jj < cmval ; jj++ ){
             memcpy( VECTIM_PTR(vim[nvim],jj), tar+kim[jj]*nts, sizeof(float)*nts ) ;
           }
           mri_free(tim) ; free(kim) ; mri_free(pim) ; /* toss the trash */
         }
         break ;

         case ORDER_PEEL:{ /* == by order of peeling from outside */
           short *depth=NULL ; Aint kk ;
           depth = THD_mask_depth( DSET_NX(dset),DSET_NY(dset),DSET_NZ(dset) ,
                                   tmask , 1 , NULL, 2 ) ;
           if( depth != NULL ){
             Aint    *idepth = (Aint *)calloc(sizeof(Aint),cmval) ;
             Aint       *kim = (Aint *)calloc(sizeof(Aint),cmval) ;
             MRI_IMAGE *tim  = mri_new(nts,cmval,MRI_float) ;
             float     *tar  = MRI_FLOAT_PTR(tim) ;
#if 0
             ININFO_message("  Computing PEEL order for mask partition #%d - %d voxels",
                            ii,cmval) ;
#endif
             for( jj=0 ; jj < cmval ; jj++ ){
               memcpy( tar+jj*nts, VECTIM_PTR(vim[nvim],jj), sizeof(float)*nts ) ;
               kim[jj] = jj ;
             }
             for( jj=kk=0 ; jj < nxyz ; jj++ ){
               if( tmask[jj] ) idepth[kk++] = depth[jj] ;
             }
#if Aintsize == 64
             qsort_intint64_t( cmval , idepth , kim ) ;
#else
             qsort_intint( cmval , idepth , kim ) ;
#endif
             for( jj=0 ; jj < cmval ; jj++ ){
               memcpy( VECTIM_PTR(vim[nvim],jj), tar+kim[jj]*nts, sizeof(float)*nts ) ;
             }
             mri_free(tim) ; free(kim) ; free(idepth) ; free(depth) ;
           }
         }
         break ;

         case ORDER_LJ:{   /* == by Ljung-Box statistic [05 Feb 2020] */
           MRI_IMAGE *tim = mri_new(nts,cmval,MRI_float) , *pim ;
           Aint      *kim = (Aint *)malloc(sizeof(Aint)*cmval) ;
           float     *tar = MRI_FLOAT_PTR(tim), *par ;
#if 0
           ININFO_message("  Computing LJ order for mask partition #%d - %d voxels",
                          ii,cmval) ;
#endif
           /* copy data into temporary image */
           for( jj=0 ; jj < cmval ; jj++ ){
             memcpy( tar+jj*nts, VECTIM_PTR(vim[nvim],jj), sizeof(float)*nts ) ;
             kim[jj] = jj ;  /* source index */
           }
           /* make the PVmap */
           pim = mri_vec_to_ljmap(tim) ; par = MRI_FLOAT_PTR(pim) ;
           /* sort so largest are first, keeping track of whence they came */
           for( jj=0 ; jj < cmval ; jj++ ) par[jj] = -par[jj] ;
#if Aintsize == 64
           qsort_floatint64_t( cmval , par , kim ) ;
#else
           qsort_floatint( cmval , par , kim ) ;
#endif
           /* copy from temp image back to vectim, in the right order */
           for( jj=0 ; jj < cmval ; jj++ ){
             memcpy( VECTIM_PTR(vim[nvim],jj), tar+kim[jj]*nts, sizeof(float)*nts ) ;
           }
           mri_free(tim) ; free(kim) ; mri_free(pim) ; /* toss the trash */
         }
         break ;

       }

       lev_siz[nvim] = cmval ;
       nvim++ ;
     }
   }
   free(tmask) ;
   if( fit != NULL ) free(fit) ;
   DSET_unload(dset) ;

   lev_num = nvim ;

   if( nvim == 0 ) RETURN(NULL) ;

   /* glue multiple level vectims into 1, if needed */

   if( nvim == 1 ){
     vout = vim[0] ;
   } else {
     vout = THD_xyzcat_vectims( nvim , vim ) ;
     for( ii=0 ; ii < nvim ; ii++ ) VECTIM_destroy(vim[ii]) ;
   }

   free(vim) ;
   RETURN(vout) ;
}

/*--------------------------------------------------------------------------*/
/* resample a 1D array to a finer (nout > nin) or coarser (nout < nin) grid */
/*--------------------------------------------------------------------------*/

#undef  P_M1     /* cubic interpolation polynomials */
#undef  P_00
#undef  P_P1
#undef  P_P2
#undef  P_FACTOR
#define P_M1(x)  (-(x)*((x)-1.0f)*((x)-2.0f))
#define P_00(x)  (3.0f*((x)+1.0f)*((x)-1.0f)*((x)-2.0f))
#define P_P1(x)  (-3.0f*(x)*((x)+1.0f)*((x)-2.0f))
#define P_P2(x)  ((x)*((x)+1.0f)*((x)-1.0f))
#define P_FACTOR 0.16666667f  /* 1/6 = final scaling factor */

/* M3(x) = minimum sidelobe 3 term window (has no catchy name) */
#undef  PIF
#define PIF   3.1415927f /* PI in float */
#undef  M3
#define M3(x) (0.4243801f+0.4973406f*cosf(PIF*(x))+0.0782793f*cosf(PIF*(x)*2.0f))

static void resample_1D_float( int nin, float *xin, int nout, float *xout )
{
   int ii , jj , nin1,nout1 , jbot,jtop ;
   float ffac , fjmid , fj ;

ENTRY("resample_1D_float") ;

   if( nin < 2 || xin == NULL || nout < 2 || xout == NULL ) EXRETURN ;

   /* nothing to do? */

   if( nin == nout ){
     memcpy( xout , xin , sizeof(float)*nin) ;
     EXRETURN ;
   }

   ffac  = (nin-1.0f)/(nout-1.0f) ;
   nout1 = nout-1 ;
   nin1  = nin -1 ;

   if( ffac < 1.0f ){ /*----- interpolate to finer grid (nout > nin) -----*/

     xout[0]     = xin[0] ;
     xout[nout1] = xin[nin1] ;

     if( use_old_resam ){  /* OLD = linear interp */

       for( ii=1 ; ii < nout1 ; ii++ ){
         fjmid = ii*ffac ;
         jbot  = (int)fjmid ; fj = (fjmid-jbot) ;
         xout[ii] = (1.0f-fj)*xin[jbot] + fj*xin[jbot+1] ;
       }

     } else { /* NEW = cubic interp [Aug 2018] */

       int jx , jx_m1,jx_00,jx_p1,jx_p2 ; float fx , wt_m1,wt_00,wt_p1,wt_p2 ;
       for( ii=1 ; ii < nout1 ; ii++ ){
         fjmid = ii*ffac ; jx = floorf(fjmid) ; fx = fjmid-jx ;
         jx_m1 = jx-1    ; jx_00 = jx         ; jx_p1 = jx+1  ; jx_p2 = jx+2 ;
         if( jx_m1 < 0    ) jx_m1 = 0   ;
         if( jx_m1 > nin1 ) jx_m1 = nin1 ;
         if( jx_00 > nin1 ) jx_00 = nin1 ;
         if( jx_p1 > nin1 ) jx_p1 = nin1 ;
         if( jx_p2 > nin1 ) jx_p2 = nin1 ;
         wt_m1 = P_M1(fx) ; wt_00 = P_00(fx) ;  /* interpolation weights */
         wt_p1 = P_P1(fx) ; wt_p2 = P_P2(fx) ;
         xout[ii] = (  wt_m1*xin[jx_m1] + wt_00*xin[jx_00]
                     + wt_p1*xin[jx_p1] + wt_p2*xin[jx_p2] ) * P_FACTOR ;
       }

     }

   } else { /*---------- coarser grid (nout < nin) ----------*/

     if( use_old_resam ){  /* OLD = not so hot :( */

       int npm = (int)rintf(ffac) ;
       float *xq = (float *)malloc(sizeof(float)*nin) ;
       float *wt = (float *)malloc(sizeof(float)*(npm+1)) ;
       float sum, wsum , ww ;
       memcpy( xq , xin , sizeof(float)*nin) ;
       for( jj=0 ; jj <= npm ; jj++ ) wt[jj] = (npm-jj+0.5f) ;
       for( ii=0 ; ii < nout ; ii++ ){
         jbot = ii-npm ; if( jbot < 0     ) jbot = 0 ;
         jtop = ii+npm ; if( jtop > nout1 ) jbot = nout1 ;
         sum = wsum = 0.0f ;
         for( jj=jbot ; jj <= jtop ; jj++ ){
           ww = wt[abs(jj-ii)] ; wsum += ww ; sum += ww*xin[jj] ;
         }
         xq[ii] = sum / wsum ;
       }
       free(wt) ;
       xout[0]      = xq[0] ;
       xout[nout-1] = xq[nin-1] ;
       for( ii=1 ; ii < nout1 ; ii++ ){
         fjmid = ii*ffac ;
         jbot  = (int)fjmid ; fj = (fjmid-jbot) ;
         xout[ii] = (1.0f-fj)*xq[jbot] + fj*xq[jbot+1] ;
       }
       free(xq) ;

     } else { /*---------- new coarsening method [Aug 2018] ----------*/

       int npm = (int)rintf(0.555f*ffac+0.333f) , jmid ;
       float sum, wsum , ww ;
       static float *wt=NULL ; static int nwt=0 ;

       if( npm < 1 ) npm = 1 ;

       if( npm != nwt ){                /* recompute taper function, if needed */
         nwt = npm ;
         wt  = (float *)realloc(wt,sizeof(float)*(nwt+2)) ;
#if 1
          sum = (npm+0.111f) ;
         wsum = 1.0f / sum ;
         for( jj=0 ; jj <= npm ; jj++ ){
           ww = (sum-jj)*wsum ; wt[jj] = M3(ww) ;
         }
#else
         for( jj=0 ; jj <= npm ; jj++ ) wt[jj] = (npm-jj+0.111f) ;
#endif
       }

       /* resample to coarser grid, averaging over neighbors */

       for( ii=0 ; ii < nout ; ii++ ){
         jmid = (int)rintf(ii*ffac) ;  /* index mapped from xout back to xin */
         jbot = jmid - npm ; if( jbot < 0    ) jbot = 0 ; /* averaging range */
         jtop = jmid + npm ; if( jtop > nin1 ) jtop = nin1 ;
         sum = wsum = 0.0f ;
         for( jj=jbot ; jj <= jtop ; jj++ ){ /* do the local averaging */
           ww = wt[abs(jj-jmid)] ; wsum += ww ; sum += ww*xin[jj] ;
         }
         xout[ii] = sum / wsum ;
       }

     }
   }

   EXRETURN ;
}

/*--------------------------------------------------------------------------*/
/* Convert a vectim to a 2D grayplot */
/*--------------------------------------------------------------------------*/

static MRI_IMAGE * mri_vectim_to_grayplot( MRI_vectim *imts, int nx, int ny )
{
   int nxx=nx , nyy=ny , ntt,nss , ii,jj , domid=0 ;
   MRI_IMAGE *imout ; byte *outar ;
   MRI_IMAGE *imttt ; byte *tttar , *tar ;
   float zbot,ztop , val , zfac , *qar , *zar=NULL , *yar ;

ENTRY("mri_vectim_to_grayplot") ;

   if( imts == NULL ) RETURN(NULL) ;

   if( nxx < 512 ) nxx = 512 ; else if( nxx > 32768 ) nxx = 32768 ;
   if( nyy < 256 ) nyy = 256 ; else if( nyy > 32768 ) nyy = 32768 ;

   ntt = imts->nvals ;
   nss = imts->nvec ; if( ntt < 19 || nss < 19 ) RETURN(NULL) ;

   /* find min and max of pre-processed data */

   if ( DO_RAW_RANGE ){
      zbot = raw_range_bot;
      ztop = raw_range_top;
      INFO_message("Do raw range: [%f, %f]", zbot, ztop);
   }
   else{
      zbot = 6.66e+33 ; ztop = -zbot ;
      for( jj=0 ; jj < nss ; jj++ ){
         qar = VECTIM_PTR(imts,jj) ;
         for( ii=0 ; ii < ntt ; ii++ ){
            if( qar[ii] < zbot ) zbot = qar[ii] ;
            else if( qar[ii] > ztop ) ztop = qar[ii] ;
         }
      }
   }

   if( zbot >= ztop ) RETURN(NULL) ;
   if( grange <= 0.0f ){
     zfac = 255.4f / (ztop-zbot) ; domid = 0 ;
   } else {
     zfac = 127.7f / grange ; domid = 1 ;
   }

   imttt = mri_new( nxx,nss , MRI_byte ) ;
   tttar = MRI_BYTE_PTR(imttt) ;

   if( nxx != ntt ){
     zar = malloc(sizeof(float)*nxx) ;
   }

   /* process each time series from vectim (each row of output) */

   for( jj=0 ; jj < nss ; jj++ ){
     qar = VECTIM_PTR(imts,jj) ;
     if( zar != NULL ){ /* change length of this time series? */
       resample_1D_float( ntt,qar , nxx,zar ) ; yar = zar ;
     } else {
       yar = qar ;
     }
     tar = tttar + jj*nxx ;
     for( ii=0 ; ii < nxx ; ii++ ){ /* scale to 0..255 range */
       if( domid ){
         val = 127.7 + yar[ii]*zfac ;
       } else {
         val = (yar[ii]-zbot)*zfac ;
       }
       tar[ii] = BYTEIZE(val) ;
     }
   }
   if( zar != NULL ){ free(zar); zar = NULL; }

   if( nss == nyy ){ /* number of rows we have == number of rows we want? */
     RETURN(imttt) ;
   }

   /* convert number of rows we have (nss) to number we want (nyy) */

   imout = mri_new( nxx,nyy , MRI_byte ) ;
   outar = MRI_BYTE_PTR(imout) ;

   yar = (float *)malloc(sizeof(float)*nss) ;
   zar = (float *)malloc(sizeof(float)*nyy) ;

   for( ii=0 ; ii < nxx ; ii++ ){  /* resize each column */
     for( jj=0 ; jj < nss ; jj++ )
       yar[jj] = tttar[ii+jj*nxx] ;
     resample_1D_float( nss,yar , nyy,zar ) ;
     for( jj=0 ; jj < nyy ; jj++ )
       outar[ii+jj*nxx] = BYTEIZE(zar[jj]) ;
   }

   free(zar) ; free(yar) ; mri_free(imttt) ;

   /* draw dashed lines between the levels from mmask */

   if( lev_num > 1 ){
     float yfac ; int kk , rr ; byte qq=1 ;
     for( kk=1 ; kk < lev_num ; kk++ ) lev_siz[kk] += lev_siz[kk-1] ;
     yfac = nyy / (float)nss ;
     for( kk=0 ; kk < lev_num-1 ; kk++ ){
       jj = (int)rintf( yfac * lev_siz[kk] ) ;
       if( jj <= 0 ) jj = 1 ; else if( jj >= nyy-1 ) jj = nyy-2 ;
       for( ii=0 ; ii < nxx ; ii++ ){
         rr = ii%19 ;
              if( rr <   9 ) outar[ii+jj*nxx] = outar[ii+(jj-1)*nxx] = qq ;
         else if( rr == 18 ) qq = 255u-qq ;
       }
     }
   }

   RETURN(imout) ;
}

/*--------------------------------------------------------------------------*/
/* Convert a 3D+time dataset to a grayplot, for fun and profit(?). */
/*--------------------------------------------------------------------------*/

MRI_IMAGE * THD_dset_to_grayplot( THD_3dim_dataset *dset ,
                                  byte *mmask ,
                                  int nxout , int nyout ,
                                  int polort , float fwhm )
{
   MRI_vectim *vim ; MRI_IMAGE *imout ;

ENTRY("THD_dset_to_grayplot") ;

   vim = THD_dset_grayplot_prep( dset , mmask , polort , fwhm ) ;

   if( nxout < 128 ) nxout = 1024 ;
   if( nyout <  64 ) nyout =  512 ;

   imout = mri_vectim_to_grayplot( vim, nxout , nyout ) ;

   VECTIM_destroy(vim) ;

   RETURN(imout) ;
}
