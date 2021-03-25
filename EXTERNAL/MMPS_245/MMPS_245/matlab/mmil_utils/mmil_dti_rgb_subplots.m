function mmil_dti_rgb_subplots(fname_FA,fname_V0,varargin)
% function mmil_dti_rgb_subplots(fname_FA,fname_V0,[options])
%
% Purpose:
%   Creates output 'TIFF' figures of input images with rgb colormaps.
%
% Required Parameters:
%   fname_FA: full path name input FA volume file (mgz or mgh file)
%   fname_V0: full path name input V0 volume file (mgz or mgh file)
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
%  'planestrings' : plane of the output image
%     {default = {'HOR','SAG','COR'}}
%  'slice_fract': vector of fractional slice position for each plane
%     if empty, will use 0.5 for each
%     {default = []}
%  'clim': vector of lower and upper bounds for imagesc
%     if empty, will auto-scale
%     {default = []}
%  'visible_flag': [0|1] whether to display images as they are generated
%    {default: 0}
%  'forceflag': [0|1] whether to overwrite existing output file
%    {default: 0}
%    
% Created:  04/18/16 by Don Hagler
% Prev Mod: 06/21/16 by Don Hagler
% Last Mod: 05/09/16 by Don Hagler
%

%% todo: option for cropping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',pwd,[],...
  'outstem','subplot',[],...
  'border_flag',false,[false,true],...
  'fig_size',[],[],...
  'tif_dpi',300,[10,10000],...
  'planestrings',{'HOR','SAG','COR'},[],...
  'slice_fract',[],[0,1],...
  'clim',[],[],...
  'visible_flag',false,[false true],...
  'forceflag',false,[false true],...
...
  'plim_flag',true,[false true],...
  'plim',[1,99],[],...
});

% set the output filename
parms.fname_out = sprintf('%s/%s.tif',parms.outdir,parms.outstem);

if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% check that input files exist
if ~exist(fname_FA,'file')
  error('file %s not found',fname_FA);
end;
if ~exist(fname_V0,'file')
  error('file %s not found',fname_V0);
end;

parms.subplot_sizex = 1;
parms.subplot_sizey = length(parms.planestrings);
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

[vol_FA,M_FA] = fs_load_mgh(fname_FA);
[vol_FA,M_FA]= fs_reorient(vol_FA,M_FA,'PRS');

[vol_V0,M_V0] = fs_load_mgh(fname_V0);
[vol_V0,M_V0]= fs_reorient(vol_V0,M_V0,'PRS');

% set intensity scaling limits
if isempty(parms.clim) && parms.plim_flag
  vals = [];
  for i=1:parms.nplanes
    img = extract_slice(vol_FA,parms.planestrings{i},parms.slice_fract(i));
    vals = cat(1,vals,mmil_colvec(img));
  end;
  clim = prctile(vals,parms.plim);
else
  clim = parms.clim;
end;
if any(isnan(clim)) || length(unique(clim))==1
  clim = [];
end;

img_rgb_cat = [];
ind_plane = [];
for i=1:parms.nplanes
  img_FA = extract_slice(vol_FA,parms.planestrings{i},parms.slice_fract(i));
  img_r = abs(extract_slice(vol_V0,parms.planestrings{i},parms.slice_fract(i),1));
  img_g = abs(extract_slice(vol_V0,parms.planestrings{i},parms.slice_fract(i),2));
  img_b = abs(extract_slice(vol_V0,parms.planestrings{i},parms.slice_fract(i),3));
  % scale FA
  if ~isempty(clim)
    img_FA = (img_FA-clim(1))*(clim(2)-clim(1));
    img_FA = min(max(img_FA,0),1);
  end;
  % combine colors
  img_rgb = cat(3,img_r.*img_FA,img_g.*img_FA,img_b.*img_FA);
  [img_rgb2,cmap] = rgb2ind(img_rgb,256);
  % concatenate images
  img_rgb_cat = cat(1,img_rgb_cat,img_rgb);
  nx = size(img_rgb,1);
  ind_plane = cat(1,ind_plane,i*ones(nx,1));
end;
% convert rgb values to linear color scale
%[img_rgb_cat,cmap] = rgb2ind(img_rgb_cat,256);

for i=1:parms.nplanes
  parms.subplot_sizez = parms.subplot_sizez + 1; % increment for subplot
  if parms.border_flag
    subplot(parms.subplot_sizex,parms.subplot_sizey,parms.subplot_sizez); 
  else
    mmil_subplot_tight(parms.subplot_sizex,parms.subplot_sizey,parms.subplot_sizez); 
  end;
  ind = find(ind_plane==i);
%  img_rgb = img_rgb_cat(ind,:);
  img_rgb = img_rgb_cat(ind,:,:);
  % display rgb image
  imagesc(img_rgb);
  colormap(cmap);
  axis off;
end;
save_plot(parms);
if ~parms.visible_flag, close(gcf); end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img = extract_slice(vol,plane,slice_fract,frame)
  if ~exist('frame','var') || isempty(frame), frame = 1; end;
  dim = size(vol);
  s = choose_slice(dim,plane,slice_fract);
  switch plane
    case 'COR'
      img = rot90(squeeze(vol(s,:,:,frame)));
    case 'SAG'
      img = rot90(squeeze(vol(:,s,:,frame)));
    case 'HOR'
      img = squeeze(vol(:,:,s,frame));
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

