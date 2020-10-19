#ifndef MRI_H
#define MRI_H

#include "const.h"
#include "matrix.h"
/* remove unwanted warnings between hips_basic.h vs. volume_io/basic.h */
#undef ABS
#undef SIGN
#ifdef Darwin
  // The result of not defining __MACTYPES__ is a scuba2 build error 
  // complaining about QT conflicting with MNI over the Point typedef. 
  // Notes from the file <mni installation>/include/volume_io/geom_structs.h:
  /* Th 'Point' typedef is annoying to Mac OS users, since Point has been 
   * a basic type on Macs since the beginning.  Testing __MACTYPES__ should
   * work at least with the OS X codebase, I don't know if it existed in
   * earlier versions of the MacTypes.h header.
   */
  #define __MACTYPES__
#endif
#include "volume_io.h"
#include "box.h"
#include "machine.h"

#define BUFTYPE  unsigned char

#define SAMPLE_NEAREST       0
#define SAMPLE_TRILINEAR     1
#define SAMPLE_SINC          2
#define SAMPLE_CUBIC         3 /*E*/
#define SAMPLE_WEIGHTED      4

#define MRI_UCHAR   0
#define MRI_INT     1
#define MRI_LONG    2
#define MRI_FLOAT   3
#define MRI_SHORT   4
#define MRI_BITMAP  5
#define MRI_TENSOR  6

#define MAX_CMDS 1000

typedef struct
{
  int  x ;
  int  y ;
  int  z ;
  int  dx ;
  int  dy ;
  int  dz ;
} MRI_REGION ;

typedef struct
{
  int           width ;
  int           height ;
  int           depth ;     /* # of slices */
  int           type ;      /* data type for slices below */
  int           imnr0 ;     /* starting image # */
  int           imnr1 ;     /* ending image # */
  int           ptype ;     /* not used */
  float         fov ;
  float         thick ;
  float         ps ;   
  float         location ;  /* not used */
  float         xsize ;     /* size of a voxel in the x direction */ 
  float         ysize ;     /* size of a voxel in the y direction */ 
  float         zsize ;     /* size of a voxel in the z direction */ 
  float         xstart ;    /* start x (in xsize units) */
  float         xend ;      /* end x  (in xsize units) */
  float         ystart ;    /* start y   (in ysize units) */  
  float         yend ;      /* end y (in ysize units) */ 
  float         zstart ;    /* start z */  
  float         zend ;      /* end z */
  float         tr ;        /* time to recovery */
  float         te ;        /* time to echo */
  float         ti ;        /* time to inversion */
  char          fname[STR_LEN] ;

  float         x_r, x_a, x_s; /* these are the RAS distances across the whole volume */
  float         y_r, y_a, y_s; /* in x, y, and z                                      */
  float         z_r, z_a, z_s; /* c_r, c_a, and c_s are the center ras coordinates    */
  float         c_r, c_a, c_s; /* ras_good_flag tells if these coordinates are set    */
  int           ras_good_flag; /* and accurate for the volume                         */

  /*  for bshorts and bfloats */
  int           brightness;
  char          subject_name[STRLEN];
  MATRIX        *register_mat;
  char          path_to_t1[STRLEN];
  char          fname_format[STRLEN];

  /* for gdf volumes */
  char          gdf_image_stem[STRLEN];

/* 
  each slice is an array of rows (mri->height of them) each of which is 
  mri->width long.
*/
  BUFTYPE       ***slices ;
  int           scale ;
  char          transform_fname[STR_LEN] ;
  General_transform transform ;   /* the next two are from this struct */
  Transform         *linear_transform ;
  Transform         *inverse_linear_transform ;
  int           free_transform ;   /* are we responsible for freeing it? */
  int           nframes ;          /* # of concatenated images */

  /* these are used to handle boundary conditions (arrays of indices) */
  int           *xi ;
  int           *yi ;
  int           *zi ;
  int           yinvert ;  /* for converting between MNC and coronal slices */
  MRI_REGION    roi ;
  int           dof ;
  double        mean ;   
  double        flip_angle ;  /* in radians */

  void*         tag_data; /* saved tag data */
  int           tag_data_size; /* size of saved tag data */
  MATRIX *i_to_r__; /* cache */
  MATRIX *r_to_i__;
	char   *cmdlines[MAX_CMDS] ;
	int    ncmds;
} MRI_IMAGE, MRI ;

