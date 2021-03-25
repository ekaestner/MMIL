function mmil_multislice_subplots(fname,varargin)
% function mmil_multislice_subplots(fname,[options])
%
% Purpose:
%   Creates output 'TIFF' figures of input images
%
% Required Parameters:
%   fname: file name of input volume file (mgz or mgh file)
%
% Optional Paramters: 
%  'outdir': output directory
%     {default = pwd}
%  'outstem': output file stem
%     {default = 'multislice_subplot'}
%  'frame': frame index
%     use first if empty
%     {default = 1}
%  'border_flag': [0|1] allow borders around each subplot
%     {default = 0}
%  'fig_size': figure size in inches
%     if not specified, will use default
%     {default = []}
%  'tif_dpi': resolution of tif files (dots per inch)
%     {default = 300}
%  'planestring' : plane of the output image
%     'HOR', 'SAG', or 'COR'
%     {default = 'HOR'}
%  'slice_fracts': vector of fractional slice positions
%     {default = [0.1:0.1:0.9]}
%  'clim': vector of lower and upper bounds for imagesc
%     if empty, will auto-scale
%     {default = []}
%  'plim_flag': [0|1] calculate clim from percentile values in plim
%     {default = 1}
%  'plim': percentile values (0-100) for intensity scaling range
%     {default = [1,99]}
%  'logtrans_flag': [0|1] calculate percentiles on log-transformed values
%     {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output file
%    {default: 0}
%    
% Created:  04/11/13 by Don Hagler
% Last Mod: 06/21/16 by Don Hagler
%

%% todo: option for cropping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',pwd,[],...
  'outstem','multislice_subplot',[],...
  'frame',1,[1,Inf],...
  'border_flag',false,[false,true],...
  'fig_size',[],[],...
  'tif_dpi',300,[10,10000],...
  'planestring','HOR',{'HOR','SAG','COR'},...
  'slice_fracts',[0.1:0.1:0.9],[0,1],...
  'clim',[],[],...
  'plim_flag',true,[false true],...
  'plim',[1,99],[],...
  'logtrans_flag',false,[false true],...
  'forceflag',false,[false true],...
  ...
  'orient','PRS',[],...
  'subsamp',10,[1,1000],... % save memory
});

% set the output filename
parms.fname_out = sprintf('%s/%s.tif',parms.outdir,parms.outstem);

if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% check whether fname exists
if ~exist(fname,'file')
  error('file %s not found',fname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

figure;
set(gcf,'Visible','off');

% load fname
[vol,M] = fs_load_mgh(fname,[],parms.frame(1));

% reorient volume to standard orientation
[vol,M]= fs_reorient(vol,M,parms.orient);

% set intensity scaling limits
if isempty(parms.clim) && parms.plim_flag
  vals = mmil_colvec(vol);
  % apply log transform
  if parms.logtrans_flag
    vals = log(max(1,vals));
  end;
  % calulate intensity limits from percentile values
  clim = prctile(vals(1:parms.subsamp:end),parms.plim);
  % transform back to original scale
  if parms.logtrans_flag
    clim = exp(clim);
  end;
  if any(isnan(clim)) || length(unique(clim))==1
    parms.clim = [];
  else
    parms.clim = clim;
  end;
end;

% determine number of rows and columns
nslices = length(parms.slice_fracts);
nrows = ceil(sqrt(nslices));
ncols = ceil(nslices/nrows);

k = 1;
for i=1:nslices
  img = extract_slice(vol,parms.planestring,parms.slice_fracts(i));
  if parms.border_flag
    subplot(nrows,ncols,k); 
  else
    mmil_subplot_tight(nrows,ncols,k); 
  end;
  k = k + 1;
  if ~isempty(parms.clim)
    imagesc(img,parms.clim);
  else
    imagesc(img);
  end;
  colormap(gray); axis off;
end;

save_plot(parms);
close(gcf);

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

