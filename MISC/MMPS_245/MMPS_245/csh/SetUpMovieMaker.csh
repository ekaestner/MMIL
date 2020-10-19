#/bin/tcsh -k

#	Imagemagick installation found by: 
#		find montage, make sure it's not the wrong one.
#
#
#
#
#

  #	Imagemagick setup
  set imagemagick_dir = `which montage | sed "s/\/montage//"`;
  if ( ! -e "$imagemagick_dir/montage" ) goto no_montage
  if ( ! -e "$imagemagick_dir/convert" ) goto no_convert
  setenv IMAGEMAGICK_DIR $imagemagick_dir;
  
  # 	setup fs_surf_mgh
  if ( "`which fs_surfmgh.csh | grep /`" == "" ) set path = ( $path ~mmildev/csh/fstools )
  if ( "`which ts_make_movie.csh | grep /`" == "" ) set path = ( $path ~mmildev/csh/timesurfer )
    
  # Check the blank.gif file.
  #set moviemaker_dir = `which ts_make_movie.csh | sed "s/\/ts_make_movie.csh//"`;
  #if ( ! -e "$moviemaker_dir/blank.gif" ) then
  #  echo "  Couldn't find blank.gif where ts_make_movie.csh is installed.";
  #  goto exit;
  #endif
  #setenv MOVIEMAKER_DIR $moviemaker_dir;

    
  #   Setup tksurfer and mri_info
  set tksurfer_prog = tksurfer_offscreen
  if ( "`which tksurfer_offscreen | grep /`" != "" ) then 
    setenv TKSURFER_PROG "tksurfer_offscreen"
  else if ( "`which tksurfer | grep /`" != "" ) then
    setenv TKSURFER_PROG "tksurfer"
  else  
    goto no_tksurfer;
  endif

  if ( "`which mri_info | grep /`" == "") goto no_mri_info;
  
  goto exit;
  
no_montage:
  echo "  Imagemagick's "montage" program not found in your path; SetUpMovieMaker.csh cannot complete.";
  goto exit;
	
no_convert:
  echo "  Imagemagick's "convert" program not found in your path; SetUpMovieMaker.csh cannot complete.";
  goto exit;
	
no_tksurfer:
  echo "  Freesurfer's "tksurfer" program not found in your path; SetUpMovieMaker.csh cannot complete.";
  goto exit;
  
no_mri_info:
  echo "  Freesurfer's "mri_info" program not found in your path; $0 cannot complete.";
  goto exit;

exit:
