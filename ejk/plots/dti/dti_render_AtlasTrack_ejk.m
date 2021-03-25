function dti_render_AtlasTrack_ejk(indir,varargin)
%function dti_render_AtlasTrack(indir,[options])
%
% Purpose: render image of AtlasTrack fiber streamlines
%   using tractoview
%
% Required Input:
%   indir: input directory
%     This may be the AtlasTrack output directory
%       or a directory containing fiber path files (grp format)
%     If AtlasTrack directory, will look for fiber path files
%       in subdirectory specified by pathsdir
%
% Optional Parameters to select input files:
%   'fibers': fiber numbers to display (see DTI_Fiber_Legend.csv)
%     {default = [101:110,115:120,123,133,134,141,142,147:150]}
%   'pathsdir': subdirectory containing fiber path files
%     may be full path, otherwise relative to indir
%     {default = 'fiber_paths'}
%   'fiber_suffix': string attached to fiber file name
%     e.g. fiber_001_<fiber_suffix>
%     {default = 'prob_countatlas_pthresh0.08_minflen12_path.grp'}
%   'fstem_image': file stem of underlay image (e.g. 'T1' or 'FA')
%     if 'T1', will look for T1_resDTI.mgh or T1.mgh in indir
%     if 'FA', will look for FA.mgh in indir
%     {default = 'T1'}
%   'fname_image': full path of underlay image volume
%     if supplied, fstem_image will be ignored
%     input image should be in register with fiber streamlines
%     {default = []}
%   'fname_color': full path of text file with RGB color codes for each fiber
%     if empty, will generate temporary color file with pre-set values
%     {default = []}
%
% Optional Parameters to control rendering:
%   'underlay_flag': whether to have a underlay image or not 
%     {default = 1}
%   'plane': plane for underlay image ('sag', 'hor', or 'cor')
%     {default = 'sag'}
%   'view': viewing orientation ('L','R','A','P','I', or 'S')
%      {default = 'L'}
%   'xrot': rotation (deg) about x axis relative to view
%      {default = 0}
%   'yrot': rotation (deg) about y axis relative to view
%      {default = 0}
%   'zrot': rotation (deg) about z axis relative to view
%      {default = 0}
%   'xflip': [0|1] whether to flip streamlines in x direction
%      {default = 0}
%   'yflip': [0|1] whether to flip streamlines in y direction
%      {default = 0}
%   'zflip': [0|1] whether to flip streamlines in z direction
%      {default = 0}
%   'image_alpha': transparency of underlay image (0 - 1)
%     {default = 0.3}
%   'fiber_alpha': transparency of fiber streamlines (0 - 1)
%     {default = 0.4}
%
% Optional Parameters:
%   'fname_tif': output file name for rendered image
%     may be full path, otherwise relative to current directory
%     if supplied, image saved, and tractoview exits
%     if empty, no image saved, and tractoview remains open
%     {default = []}
%   'forceflag': overwrite existing output files
%     {default = 1}
%
% Created:  04/11/12 by Don Hagler
% Last Mod: 07/23/14 by Don Hagler
%

%% todo: create mgh file if only mgz file exists

%% todo: write messages to log file; use mmil_logstr?

%   'fname_log': output file name for log messages
%     may be full path, otherwise relative to current directory
%     if empty, no log file created
%     {default = []}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters and check for problems
parms = check_input(indir,varargin);

% set parameter values based on input
parms = set_parms(parms);

% find underlay file
parms = find_image(parms);

% find fiber files
parms = find_fibers(parms);

% create color file if not specified
parms = set_colors(parms);

% run tractoview
parms = render_fibers(parms);

% remove temporary files
% cleanup_tmpfiles(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(indir,args)
  parms_filter = {...
    'indir',indir,[],...
...
    'fibers',[101:110,115:120,123,133,134,141,142,147:150],[],...
    'pathsdir','fiber_paths',[],...
    'fiber_suffix','prob_countatlas_pthresh0.08_minflen12_path.grp',[],...
    'fstem_image','T1',{'T1','FA'},...
    'fname_image',[],[],...
    'fname_color',[],[],...
...
    'underlay_flag',true,[false true],...
    'plane','sag',{'sag','hor','cor'},...
    'view','L',{'L','R','A','P','I','S'},...
    'xrot',0,[-360,360],...
    'yrot',0,[-360,360],...
    'zrot',0,[-360,360],...
    'xflip',false,[false true],...
    'yflip',false,[false true],...
    'zflip',true,[false true],...
    'image_alpha',0.3,[0 1],...
    'fiber_alpha',0.4,[0 1],...
...
    'fname_tif',[],[],...
    'fname_log',[],[],...
    'forceflag',true,[false true],...
...
    'fname_legend',[],[],...
    'tmpdir',pwd,[],...
    'cleanup_flag',true,[false true],...
    'fname_csh',[],[],...
  };
  parms = mmil_args2parms(args,parms_filter);

  if mmil_isrelative(parms.pathsdir)
    parms.pathsdir = [parms.indir '/' parms.pathsdir];
    if ~exist(parms.pathsdir,'dir')
      parms.pathsdir = parms.indir;
    end;
  elseif ~exist(parms.pathsdir,'dir')
    error('pathsdir %s not found',parms.pathsdir);
  end;
  if ~exist(parms.indir,'dir')
    error('input dir %s not found',parms.indir);
  end;

  if ~isempty(parms.fname_tif) && exist(parms.fname_tif,'file') && ~parms.forceflag
    fprintf('%s: not overwriting existing file %s...\n',mfilename,parms.fname_tif);
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_parms(parms)
  if ~isempty(parms.fname_log)
    parms.log_flag = 1;
  else
    parms.log_flag = 0;
  end;

  if parms.underlay_flag 
    switch parms.plane
      case 'sag'
        parms.ulay_parms = '-norendhor -norendcor';
      case 'hor'
        parms.ulay_parms = '-norendcor -norendsag';
      case 'cor'
        parms.ulay_parms = '-norendsag -norendhor';
      otherwise
        error('invalid plane');
    end;
  else
    parms.ulay_parms = '-norendsag -norendhor -norendcor';
  end;    

  parms.flip_parms = [];
  if parms.xflip
    parms.flip_parms = [parms.flip_parms ' -xflip'];
  else
    parms.flip_parms = [parms.flip_parms ' -noxflip'];
  end;  
  if parms.yflip
    parms.flip_parms = [parms.flip_parms ' -yflip'];
  else
    parms.flip_parms = [parms.flip_parms ' -noyflip'];
  end;  
  if parms.zflip
    parms.flip_parms = [parms.flip_parms ' -zflip'];
  else
    parms.flip_parms = [parms.flip_parms ' -nozflip'];
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = find_image(parms)
  if isempty(parms.fname_image)
    switch parms.fstem_image
      case 'T1'
        parms.fname_image = sprintf('%s/%s_resDTI.mgh',parms.indir,parms.fstem_image);
        if ~exist(parms.fname_image,'file')
          parms.fname_image = sprintf('%s/%s.mgh',parms.indir,parms.fstem_image);
        end;
      case 'FA'
        parms.fname_image = sprintf('%s/%s.mgh',parms.indir,parms.fstem_image);
      otherwise
        error('invalid fstem_image');
    end;
  end;
  if ~exist(parms.fname_image,'file')
    error('image file %s not found',parms.fname_image);
  end;
 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = find_fibers(parms)
  parms.nfibers = length(parms.fibers);
  parms.fname_fiber_list = cell(1,parms.nfibers);
  for f=1:parms.nfibers
    fnum = parms.fibers(f);
    fname_fiber = sprintf('%s/fiber_%03d_%s',parms.pathsdir,fnum,parms.fiber_suffix);
    if exist(fname_fiber,'file')
      parms.fname_fiber_list{f} = fname_fiber;
    else
      fprintf('%s: WARNING: fiber file %s not found\n',mfilename,fname_fiber);
    end;
  end;
  ind_found = find(~cellfun(@isempty,parms.fname_fiber_list));
  if length(ind_found) < parms.nfibers
    parms.fname_fiber_list = parms.fname_fiber_list(ind_found);
    parms.fibers = parms.fibers(ind_found);
    parms.nfibers = length(parms.fibers);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_colors(parms)
  if isempty(parms.fname_color)
    parms.tmp_color_flag = 1;
    parms.legend = DTI_MMIL_Read_Fiber_Legend(parms.fname_legend);
    i = 1;
    while isempty(parms.fname_color) || ((exist(parms.fname_color,'file')) && i<1000)
      [tpath,tstem] = fileparts(tempname);
      parms.fname_color = sprintf('%s/color_%s.dat',parms.tmpdir,tstem);
      i = i + 1;
    end;

    if isfield(parms.legend,'FiberNumber')
      fnums = [parms.legend.FiberNumber];
    else
      error('legend file %s is missing FiberNumber')
    end;
    if isfield(parms.legend,'ColorRGB')
      colors = {parms.legend.ColorRGB};
    else
      error('legend file %s is missing ColorRGB')
    end;
    if isfield(parms.legend,'ColorDescription')
      colornames = deblank({parms.legend.ColorDescription});
    else
      colornames = [];
    end;

    [tmp,ind_fibers,ind_legend] = intersect(parms.fibers,fnums);
    if length(ind_fibers) < parms.nfibers
      error('not all fibers found in legend');
    end;
    [tmp,ind_sort] = sort(ind_fibers);
    colors = colors(ind_legend(ind_sort));
    if ~isempty(colornames)
      colornames = colornames(ind_legend(ind_sort));
    end;

    fid = fopen(parms.fname_color,'wt');
    if fid==-1
      error('unable to write %s',parms.fname_color);
    end;
    if isempty(colornames)
      for f=1:parms.nfibers
        fprintf(fid,'%-20s\n',sprintf('%d ',colors{f}));
      end;
    else
      for f=1:parms.nfibers
        fprintf(fid,'%-20s # %s\n',sprintf('%d ',colors{f}),colornames{f});
      end;
    end;
    fclose(fid);
  else
    if ~exist(parms.fname_color,'file')
      error('color file %s not found',parms.fname_color);
    end;
    parms.tmp_color_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = render_fibers(parms)
  i = 1;
  while isempty(parms.fname_csh) || ((exist(parms.fname_csh,'file')) && i<1000)
    [tpath,tstem] = fileparts(tempname);
    parms.fname_csh = sprintf('%s/run_tracoview_%s.csh',parms.tmpdir,tstem);
    i = i + 1;
  end;

  fid = fopen(parms.fname_csh,'wt');
  if fid==-1
    error('unable to write %s',parms.fname_csh);
  end;
  fprintf(fid,'#!/bin/csh -f\n');
  fprintf(fid,'tractoview \\\n');
  fprintf(fid,' -fnameimage %s %s \\\n',parms.fname_image,parms.ulay_parms);
  fprintf(fid,' -fnamefibers \\\n');
  for f=1:parms.nfibers
    fprintf(fid,'   %s \\\n',parms.fname_fiber_list{f});
  end;  
  fprintf(fid,' -fnamecolor %s \\\n',parms.fname_color);
  fprintf(fid,' -initview %s \\\n',parms.view);
  fprintf(fid,' %s \\\n',parms.flip_parms);
  fprintf(fid,' -xrot %0.2f -yrot %0.2f -zrot %0.2f \\\n',parms.xrot,parms.yrot,parms.zrot);  
  fprintf(fid,' -image_alpha %0.3f \\\n',parms.image_alpha);
  fprintf(fid,' -fiber_alpha %0.3f \\\n',parms.fiber_alpha);
  if ~isempty(parms.fname_tif)
    fprintf(fid,' -savetif -fname_tif %s\\\n',parms.fname_tif);
  end;
  fprintf(fid,'\n');
  fclose(fid);
  
  % change permissions so group can write and execute
  cmd = sprintf('chmod ug+rwx %s',parms.fname_csh);
  [s,r] = mmil_unix(cmd);
  if s
    fprintf('%s: WARNING: failed to set permissions for %s:\n%s\n',...
      mfilename,parms.fname_csh,r);
  end;

  [s,r] = mmil_unix(['source ' parms.fname_csh]);
  if s
    error('failed to run tractoview with %s:\n%s\n',...
      mfilename,parms.fname_csh,r);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpfiles(parms)
  if parms.cleanup_flag
    delete(parms.fname_csh);
    if parms.tmp_color_flag
      delete(parms.fname_color);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

