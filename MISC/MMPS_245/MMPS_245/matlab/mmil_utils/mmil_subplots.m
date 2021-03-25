function mmil_subplots(fnamelist,varargin)
% function mmil_subplots(fnamelist,[options])
%
% Purpose:
%   Creates output 'TIFF' figures of input images.
%
% Required Parameters:
%   fnamelist: list of filenames of input volume file (mgz or mgh file)
%
% Optional Paramters: 
%  'outdir': output directory
%     {default = pwd}
%  'outstem': output file stem
%     {default = 'subplot'}
%  'border_flag': [0|1] allow borders around each subplot
%     {default = 0}
%  'fig_size': figure size in inches
%     if not specified, will use default
%     {default = []}
%  'tif_dpi': resolution of tif files (dots per inch)
%     {default = 300}
%  'frame': frame number for multi-frame volume
%     {default = 1}
%  'planestrings' : plane of the output image
%     {default = {'HOR','SAG','COR'}}
%  'slice_fract': vector of fractional slice position for each plane
%     if empty, will use 0.5 for each
%     {default = []}
%  'clim': vector of lower and upper bounds for imagesc
%     if empty, will auto-scale
%     {default = []}
%  'cmap': colormap
%     {default = 'gray'); 
%  'visible_flag': [0|1] whether to display images as they are generated
%    {default: 0}
%  'forceflag': [0|1] whether to overwrite existing output file
%    {default: 0}
%    
% Created:  05/29/12 by Vijay Venkatraman
% Prev Mod: 04/18/16 by Don Hagler
% Last Mod: 11/05/17 by Don Hagler
%

%% todo: option for cropping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',pwd,[],...
  'outstem','subplot',[],...
  'border_flag',false,[false,true],...
  'fig_size',[],[],...
  'tif_dpi',300,[10,10000],...
  'frame',1,[],...
  'planestrings',{'HOR','SAG','COR'},[],...
  'slice_fract',[],[0,1],...
  'clim',[],[],...
  'cmap','gray',[],...
  'visible_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'plim_flag',true,[false true],...
  'plim',[1,99],[],...
});

% set the output filename
parms.fname_out = sprintf('%s/%s.tif',parms.outdir,parms.outstem);

if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% check whether fnamelist input is cell, otherwise convert to cell
if ~iscell(fnamelist), fnamelist = {fnamelist}; end;

% set number of files
parms.nfiles = length(fnamelist);
if parms.nfiles == 0
  error('%s: input file list is empty',mfilename);
end;
for f=1:parms.nfiles
  fname = fnamelist{f};
  if ~exist(fname,'file')
    error('file %s not found',fname);
  end;
end;

parms.subplot_sizex = parms.nfiles;
parms.subplot_sizey = length(parms.planestrings); % three planes
parms.subplot_sizez = 0; % used as a counter

parms.planestrings = upper(parms.planestrings);
parms.nplanes = length(parms.planestrings);
if isempty(parms.slice_fract)
  parms.slice_fract = 0.5*ones(parms.nplanes,1);
else
  if length(parms.slice_fract) ~= parms.nplanes
    error('slice_fract vector must match number of planestrings');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

figure;
if ~parms.visible_flag
  set(gcf,'Visible','off');
end;
for n = 1: parms.nfiles
  fname = fnamelist{n};
  [vol,M] = fs_load_mgh(fname,[],parms.frame);
  [vol,M]= fs_reorient(vol,M,'PRS');

  % set intensity scaling limits
  if isempty(parms.clim) && parms.plim_flag
    vals = [];
    for i=1:parms.nplanes
      img = extract_slice(vol,parms.planestrings{i},parms.slice_fract(i));
      vals = cat(1,vals,mmil_colvec(img));
    end;
    clim = prctile(vals,parms.plim);
  else
    clim = parms.clim;
  end;
  if any(isnan(clim)) || length(unique(clim))==1
    clim = [];
  end;

  for i=1:parms.nplanes
    img = extract_slice(vol,parms.planestrings{i},parms.slice_fract(i));
    parms.subplot_sizez = parms.subplot_sizez + 1; % increment for subplot
    if parms.border_flag
      subplot(parms.subplot_sizex,parms.subplot_sizey,parms.subplot_sizez); 
    else
      mmil_subplot_tight(parms.subplot_sizex,parms.subplot_sizey,parms.subplot_sizez); 
    end;
    if ~isempty(clim)
      imagesc(img,clim);
    else
      imagesc(img);
    end;
    colormap(parms.cmap); axis off;
  end;
end;   
save_plot(parms);
if ~parms.visible_flag, close(gcf); end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img = extract_slice(vol,plane,slice_fract)
  dim = size(vol);
  s = choose_slice(dim,plane,slice_fract);
  switch plane
    case 'COR'
      img = rot90(squeeze(vol(s,:,:)));
    case 'SAG'
      img = rot90(squeeze(vol(:,s,:)));
    case 'HOR'
      img = squeeze(vol(:,:,s));
  end;
return;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slice_num = choose_slice(dim,plane,slice_fract)
  switch plane
    case 'COR'
      d = 1;
    case 'SAG'
      d = 2;
    case 'HOR'
      d = 3;
  end;
  slice_num = round(dim(d)*slice_fract);
  slice_num = min(max(slice_num,1),dim(d));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms)
  if ~isempty(parms.fig_size)
    if length(parms.fig_size)==1
      parms.fig_size = [parms.fig_size parms.fig_size];
    end;
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [parms.fig_size(1) parms.fig_size(2)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 parms.fig_size(1) parms.fig_size(2)]);
  end;
  print(gcf,'-dtiff',parms.fname_out,sprintf('-r %d',parms.tif_dpi));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

