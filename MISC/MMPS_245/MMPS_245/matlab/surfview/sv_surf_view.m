function sv_surf_view(subj,fname,varargin)
%function sv_surf_view(subj,fname,[options])
%
% Purpose: load FreeSurfer subject surface and overlay values as colors
%
% Required Input:
%   subj: FreeSurfer subject name
%   fname: name of mgh/mgz file containing surface values
%
% Options for specifying surface values:
%   'fname_imag': name of mgh/mgz file containing imaginary component
%     of complex surface values (e.g. retinotopy)
%     {default = []}
%   'frames': vector of frame numbers from fname to load and display
%     if empty, load all and loop over each
%     {default = []}
%
% Options for specifying FreeSurfer subject and surfaces
%   'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%     subjdir/subj should contain the freesurfer subject directory
%     {default = $SUBJECTS_DIR}
%   'surfname': surface to display
%     {default = inflated}
%   'surfdir': input/output directory for surface mat files
%     for faster loading of surfaces on subsequent uses
%     if empty, will load surfaces from binary files
%     {default = []}
%   'hemi': cortical hemisphere; allowed values: {'lh','rh'}
%     if empty, will infer from fname assuming pattern *-?h.mgh
%     {default = []}
%
% Options for controlling colors:
%   'colorscale': type of colorscale
%     allowed values: 'linear','polar','eccen','category'
%     note: for polar and eccen, fname_imag must be specified
%     {default = 'linear'}
%   'cmap': color map matrix with size [n,3]
%     alternatively, may be any string that is a valid input
%       for colormap function (e.g. 'hot', 'gray', 'jet', etc.)
%     if not specified, will use a default cmap depending on colorscale
%       linear:   'mmil_cmap_blueblackred'
%       polar:    'mmil_cmap_polar'
%       eccen:    'mmil_cmap_eccen'
%       category: 'jet'
%     {default = []}
%   'fmax': value corresponding to end of color scale
%     if empty, will be set to max(vals)
%     {default = []}
%   'fmin': value corresponding to start of color scale
%     if empty, will be set to fmax/10
%     {default = []}
%   'fmid': value corresponding to middle of color scale
%     if empty, will be set to (fmin+fmax)/2
%     {default = []}
%   'phase_offset': offset value added to phase values
%     applies to 'polar' and 'eccen' only
%     if not specified, will be set to 0.5
%       if colorscale = 'polar' and hemi = 'rh'
%     {default = []}
%   'phase_revflag': [0|1] reverse phases
%     applies to 'polar' and 'eccen' only
%     if not specified, will be set to 1
%       if colorscale = 'polar' and hemi = 'rh'
%     {default = []}
%   'curvflag': [0|1] display gray-scale curvature
%     {default = 0}
%   'curvfact': gray-scale contrast of curvature values
%     {default = 0.2}
%
% Options for controlling view
%   'view': string defining view angle
%     allowed values: {'left','right','lat','med','pos',...
%                      'ant','sup','ven','infr','infl'}
%     {default = 'lat'}
%   'zoom': zoom factor
%     {default = 1}
%
% Options for saving images:
%   'tif_flag': [0|1] save image as tif
%     {default = 0}
%   'tif_dpi': resolution (dots per inch) for output tif file
%     {default = 200}
%   'outstem': file stem for output tif file
%     if empty, will use stem of fname
%     {default = []}
%   'outdir': output directory for tif files
%     if empty, will use pwd
%     ignored if outstem is supplied with full path
%     {default = []}
%   'visible_flag': [0|1] display image on screen
%     only applies if tif_flag = 1
%     {default = 0}
%   'pause_dur': duration (sec) of pause between frames
%     only applies if multiple frames and tif_flag = 0
%     {default = 0.5}
%   'forceflag': [0|1] overwrite existing output
%     only applies if tif_flag = 1 and visible_flag = 0
%     {default = 0}
%
% Created:  09/25/12 by Don Hagler
% Last Mod: 11/08/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: log transform for eccen

%% todo: colorbar with original scale values for tick labels

%% todo: option to display both hemispheres together
%%       (concat surfs, curvvals, and vals)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

% check input parameters
parms = check_input(subj,fname,varargin);

% check output files
if parms.tif_flag && ~parms.visible_flag && ~parms.forceflag
  allexist = check_output(parms);
  if allexist, return; end;
end;

% load FreeSurfer surface
surf = load_surf(parms);

% load curvature
if parms.curvflag
  parms = load_curv(parms);
end;

% show surfaces with values
switch parms.colorscale
  case 'linear'
    show_surf_vals_linear(surf,parms);
  case {'polar','eccen'}
    show_surf_vals_complex(surf,parms);
  case 'category'
    show_surf_vals_category(surf,parms);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(subj,fname,options)
  parms = mmil_args2parms(options,{,...
    'subj',subj,[],...
    'fname',fname,[],...
  ... % input values
    'fname_imag',[],[],...
    'frames',[],[1,Inf],...
  ... % surfaces
    'subjdir',[],[],...
    'surfname','inflated',[],...
    'surfdir',[],[],...
    'hemi',[],{'lh','rh'},...
  ... % colors
    'colorscale','linear',{'linear','polar','eccen','category'},...
    'cmap',[],[],...
    'fmax',[],[],...
    'fmin',[],[],...
    'fmid',[],[],...
    'curvflag',false,[false true],...
    'curvfact',0.2,[0,1],...
  ... % view
    'view','lat',{'left','right','lat','med','pos','ant','sup','ven','infl','infr'},...
    'zoom',1,[0.01,100],...
  ... % images
    'tif_flag',false,[false true],...
    'tif_dpi',300,[10,1000],...
    'outstem',[],[],...
    'outdir',[],[],...
    'visible_flag',false,[false true],...
    'pause_dur',0.5,[0,10],...
    'forceflag',false,[false true],...
  ...
    'phase_offset',[],[-1,1],...
    'phase_revflag',[],[false true],...
  ...
    'show_tags',{'cvals','view','zoom'},[],...
    'linear_cvals_tags',{'fmax','fmin','fmid','cmap',...
                         'curvvals','curvfact'},[],...
    'complex_cvals_tags',{'fmax','fmin','fmid','cmap',...
                         'curvvals','curvfact',...
                         'phase_revflag','phase_offset'},[],...
    'category_cvals_tags',{'cmap','curvvals','curvfact'},[],...
  });

  if isempty(parms.subjdir)
    parms.subjdir = getenv('SUBJECTS_DIR');
    if isempty(parms.subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  end;

  parms.fspath = [parms.subjdir '/' parms.subj];
  if ~exist(parms.fspath,'dir')
    error('FreeSurfer dir %s not found',parms.fspath);
  end;

  % check input files
  if ~exist(parms.fname,'file'), error('file %s not found',parms.fname); end;
  if ismember(parms.colorscale,{'polar','eccen'})
    if isempty(parms.fname_imag)
      error('fname_imag not specified');
    end;
    if ~exist(parms.fname_imag,'file')
      error('file %s not found',parms.fname_imag);
    end;
  end;

  % get hemi from fname
  if isempty(parms.hemi)
    n = regexp(parms.fname,'.+(?<hemi>[lr]h).mg[hz]','names');
    if isempty(n)
      error('fname does not match expected pattern, must supply hemi');
    end;
    parms.hemi = n.hemi;
  end;

  % set cmap if not specified
  if isempty(parms.cmap)
    switch parms.colorscale
      case 'linear'
        parms.cmap = 'mmil_cmap_blueblackred';
      case 'polar'
        parms.cmap = 'mmil_cmap_polar';
        if strcmp(parms.hemi,'lh')
          if isempty(parms.phase_revflag)
            parms.phase_revflag = true;
          end;
        else
          if isempty(parms.phase_offset)
            parms.phase_offset = 0.5;
          end;
        end;
      case 'eccen'
        parms.cmap = 'mmil_cmap_eccen';
      case 'category'
        parms.cmap = 'jet';
    end;
  end;

  [M,volsz] = fs_read_header(parms.fname);
  if length(volsz)>3
    parms.nframes = volsz(4);
  else
    parms.nframes = 1;
  end;
  if ~isempty(parms.frames)
    badframes = setdiff(parms.frames,[1:parms.nframes]);
    if ~isempty(badframes)
      error('invalid frame numbers: %s',sprintf('%d ',badframes));
    end;
  else
    parms.frames = [1:parms.nframes];
  end;
  if ~isempty(parms.fname_imag)
    [M,volsz] = fs_read_header(parms.fname);
    if volsz(4) ~= parms.nframes
      error('number of frames in fname_imag (%d) does not match fname (%d)',...
        volsz(4),parms.nframes);
    end;  
  end;

  % if view is 'lat' or 'med', replace with 'left' or 'right' depending on hemi
  switch parms.view
    case 'lat'
      switch parms.hemi
        case 'lh', parms.view = 'left';
        case 'rh', parms.view = 'right';
      end;
    case 'med'
      switch parms.hemi
        case 'lh', parms.view = 'right';
        case 'rh', parms.view = 'left';
      end;
    case 'ven'
      switch parms.hemi
        case 'lh', parms.view = 'infl';
        case 'rh', parms.view = 'infr';
      end;
  end;

  if parms.tif_flag
    if isempty(parms.outstem)
      [tmp,parms.outstem] = fileparts(parms.fname);
    end;
    if mmil_isrelative(parms.outstem)
      if isempty(parms.outdir), parms.outdir = pwd; end;
      parms.outstem = [parms.outdir '/' parms.outstem];
    else
      parms.outdir = fileparts(parms.outstem);
    end;
    mmil_mkdir(parms.outdir);
  else
    parms.visible_flag = true;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function allexist = check_output(parms,frame)
  allexist = 1;  
  switch parms.colorscale
    case {'linear','polar','eccen'}
      for f=1:length(parms.frames)
        frame = parms.frames(f);
        fname_out = set_fname_out(parms,frame);
        if ~exist(fname_out,'file')
          allexist = 0;
          break;
        end;
      end;
    otherwise
      fname_out = set_fname_out(parms);
      if ~exist(fname_out,'file'), allexist = 0; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function surf = load_surf(parms)
  surf = [];
  if isempty(parms.surfdir)
    surf = fs_load_subj(parms.subj,parms.hemi,parms.surfname,0,parms.subjdir);
  else
    fname_surf = sprintf('%s/%s-%s-%s.mat',...
      parms.surfdir,parms.subj,parms.surfname,parms.hemi);
    if ~exist(fname_surf,'file')
      mmil_mkdir(parms.surfdir);
      surf = fs_load_subj(parms.subj,parms.hemi,parms.surfname,0,parms.subjdir);
      save(fname_surf,'surf');
    else
      load(fname_surf);
    end;    
  end;
  tmp_surf = surf;
  surf = [];
  surf.vertices = tmp_surf.vertices;
  surf.faces = tmp_surf.faces;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_curv(parms);
  fname_curv = sprintf('%s/%s/surf/%s.curv',...
    parms.subjdir,parms.subj,parms.hemi);
  if ~exist(fname_curv,'file')
    error('file %s not found',fname_curv);
  end;
  parms.curvvals = fs_read_curv(fname_curv);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_surf_vals_linear(surf,parms)
  for f=1:length(parms.frames)
    frame = parms.frames(f);
    % check if output exists or should be overwritten
    [skip_flag,parms] = check_output_file(parms,frame);
    if skip_flag, continue; end;
    % load values
    vals = squeeze(fs_load_mgh(parms.fname,[],frame));
    % set color values
    args = mmil_parms2args(parms,parms.linear_cvals_tags);
    parms.cvals = sv_linear_cvals(vals,args{:});
    % display surface with colors
    show_surf(surf,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_surf_vals_complex(surf,parms)
  for f=1:length(parms.frames)
    frame = parms.frames(f);
    % check if output exists or should be overwritten
    [skip_flag,parms] = check_output_file(parms,frame);
    if skip_flag, continue; end;
    % load values
    vals = squeeze(fs_load_mgh(parms.fname,[],frame));
    vals_imag = squeeze(fs_load_mgh(parms.fname_imag,[],frame));
    vals = cat(2,vals,vals_imag);
    clear vals_imag;
    % set color values
    args = mmil_parms2args(parms,parms.complex_cvals_tags);
    parms.cvals = sv_complex_cvals(vals,args{:});
    % display surface with colors
    show_surf(surf,parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_surf_vals_category(surf,parms)
  % check if output exists or should be overwritten
  [skip_flag,parms] = check_output_file(parms);
  if skip_flag, return; end;
  % load values
  vals = squeeze(fs_load_mgh(parms.fname,[],parms.frames));
  % set color values
  args = mmil_parms2args(parms,parms.category_cvals_tags);
  parms.cvals = sv_category_cvals(vals,args{:});
  % display surface with colors
  show_surf(surf,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [skip_flag,parms] = check_output_file(parms,frame)
  skip_flag = 0;
  if ~exist('frame','var'), frame = []; end;
  if parms.tif_flag
    parms.fname_out = set_fname_out(parms,frame);
    if exist(parms.fname_out,'file') &&...
       ~parms.forceflag && ~parms.visible_flag
      skip_flag = 1;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_fname_out(parms,frame)
  fname_out = [];
  if ~exist('frame','var'), frame = []; end;
  if isempty(frame)
    fname_out = parms.outstem;
  elseif parms.nframes>=10000
    fname_out = sprintf('%s-%05d',parms.outstem,frame);
  elseif parms.nframes>=1000
    fname_out = sprintf('%s-%04d',parms.outstem,frame);
  elseif parms.nframes>=100
    fname_out = sprintf('%s-%03d',parms.outstem,frame);
  elseif parms.nframes>=10
    fname_out = sprintf('%s-%02d',parms.outstem,frame);
  elseif parms.nframes>1
    fname_out = sprintf('%s-%d',parms.outstem,frame);
  else
    fname_out = parms.outstem;
  end;
  fname_out = [fname_out '.tif'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_surf(surf,parms)
  % create figure
  if ~parms.visible_flag
    figure; clf;
    set(gcf,'visible','off');
  else
    clf;
  end;
  % display surface with colors
  args = mmil_parms2args(parms,parms.show_tags);
  sv_show_surf(surf,args{:});
  axis off;
  % save image
  if parms.tif_flag
    print(gcf,'-dtiff',parms.fname_out,sprintf('-r %d',parms.tif_dpi));
  end;
  if ~parms.visible_flag
    close(gcf);
  else
    pause(parms.pause_dur);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


