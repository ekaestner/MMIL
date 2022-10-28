/*****************************************************************************
   Major portions of this software are copyrighted by the Medical College
   of Wisconsin, 1994-2000, and are released under the Gnu General Public
   License, Version 2.  See the file README.Copyright for details.
******************************************************************************/
   
#ifndef _MCW_MALLOC_HEADER_
#define _MCW_MALLOC_HEADER_

/*----- 24 Jan 2001: modified slightly to add some comments, and
                     to fit in with the hashtable-ized mcw_malloc.c -----*/

#ifdef  __cplusplus
extern "C" {
#endif

#ifndef ALLOW_MCW_MALLOC
#ifndef DONT_USE_MCW_MALLOC  /* don't redefine to avoid warnings */
#define DONT_USE_MCW_MALLOC  /* old way to mark  mcw_malloc usage/non-usage */
#endif
#else
# undef  DONT_USE_MCW_MALLOC
#endif

/*---------------------------------------------------------------------------*/
#ifdef DONT_USE_MCW_MALLOC

#undef USING_MCW_MALLOC

#define MCW_MALLOC_enabled 0

#define MCW_MALLOC_status  NULL
#define MCW_MALLOC_total   0

#undef  mcw_malloc
#define mcw_malloc  malloc

#undef  mcw_malloc_OK
#define mcw_malloc_OK(p) 1

#undef  mcw_realloc
#define mcw_realloc realloc

#undef  mcw_calloc
#define mcw_calloc  calloc

#undef  mcw_free
#define mcw_free(a,b,c)  free(a)

#undef  mcw_strdup
#define mcw_strdup  strdup

#define mcw_malloc_dump_fp(x) /*nada*/

/***** extern void   mcw_malloc_dump_fp(FILE *fp) ; *****/

/*---------------------------------------------------------------------------*/
#else

#define USING_MCW_MALLOC

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "replaceXt.h"

/* Use mrix for X dependent programs */
#ifdef __BUILDING_QUICKLOOK_PLUGIN__
  #include "IntrinsicQuickLook.h"
#else
  #include <X11/Intrinsic.h>
#endif

#include "machdep.h"

/*-- define macros to replace the source code's use of malloc(), etc. --*/

#undef malloc
#undef realloc
#undef calloc
#undef free

#define malloc(a)     mcw_malloc((a),__FILE__,__LINE__)
#define realloc(a,b)  mcw_realloc((a),(b),__FILE__,__LINE__)
#define calloc(a,b)   mcw_calloc((a),(b),__FILE__,__LINE__)
#define free(a)       mcw_free((a),__FILE__,__LINE__)

#undef  strdup
#define strdup(a)     mcw_strdup((a),__FILE__,__LINE__)

/*-- prototypes for interface functions --*/

extern void   enable_mcw_malloc() ;
extern void * mcw_malloc( size_t , char * , int ) ;
extern void * mcw_realloc( void * , size_t , char * , int ) ;
extern void * mcw_calloc( size_t , size_t , char * , int ) ;
extern void   mcw_free( void * , char * , int ) ;
extern char * mcw_strdup( char *, char * , int ) ; /* 06 May 2015 */

extern char * mcw_malloc_status(const char *,int) ;
extern void   mcw_malloc_dump(void) ;
extern void   mcw_malloc_dump_sort(int opt) ;
extern void   mcw_malloc_dump_fp(FILE *fp) ;  /* 23 Apr 2015 */
extern int    mcw_malloc_enabled(void) ;
extern void   pause_mcw_malloc(void);
extern void   resume_mcw_malloc(void);
extern int    mcw_malloc_paused(void);
extern int    mcw_malloc_OK(void *);      /* 10 Jun 2014 */

extern long long mcw_malloc_total(void) ; /* 01 Feb 2007 */

/*-- how to check if the tracking routines are working --*/

#define MCW_MALLOC_enabled mcw_malloc_enabled()

#define MCW_MALLOC_status  mcw_malloc_status(__FILE__,__LINE__)
#define MCW_MALLOC_total   mcw_malloc_total()

/*-- do the same macro thing for the Xt library functions --*/

#undef XtMalloc
#undef XtRealloc
#undef XtFree
#undef XtCalloc

#define XtMalloc(a)     mcw_XtMalloc((a),__FILE__,__LINE__)
#define XtRealloc(a,b)  mcw_XtRealloc((char *)(a),(b),__FILE__,__LINE__)
#define XtCalloc(a,b)   mcw_XtCalloc((a),(b),__FILE__,__LINE__)
#define XtFree(a)       mcw_XtFree((char *)(a))

extern char * mcw_XtMalloc( RwcCardinal , char * ,  int ) ;
extern char * mcw_XtRealloc( char * , RwcCardinal , char * ,  int ) ;
extern char * mcw_XtCalloc( RwcCardinal , RwcCardinal , char * ,  int ) ;
extern void   mcw_XtFree( char * ) ;

#endif /* DONT_USE_MCW_MALLOC */
/*---------------------------------------------------------------------------*/

/*-- some macros used in various AFNI places --*/

#define myXtFree(xp)  (RwcFree((char *)(xp)) , (xp)=NULL)
#define myXtNew(type) ((type *) RwcCalloc(1,(RwcCardinal) sizeof(type)))
#define myfree(xp)    (free((xp)) , (xp)=NULL)

#ifdef  __cplusplus
}
#endif

#endif /* _MCW_MALLOC_HEADER_ */
