function rsi_plot_FOD(FODfit,varargin)
%function rsi_plot_FOD(FODfit,[options])
%
% Purpose: plot FOD overlayed on brain images
%
% Required Parameters:
%   FODfit: structure containing FOD calculations with fields:
%     volFOD : FOD fit parameters      size = [nx,ny,nz,nb]
%     volb0 : average b=0 image        size = [nx,ny,nz]
%     volmask : mask volume            size = [nx,ny,nz]
%     volmask_dilated : dilated mask   size = [nx,ny,nz]
%     M : vox2ras matrix
%     volsz : size of input 4D vol          = [nx,ny,nz,nf]
%     qmat_orig : input q matrix       size = [nf,3]
%     qmat : unit vector q matrix      size = [nf,3]
%     bvals : input b values           size = [nf,1]
%     b0_scalefacts : scale factors calculated from b=0 images
%     i_b0 : indices of b=0 images
%     ADC_long : longitudinal ADC
%     ADC_trans : transverse ADC
%     SH_order : spherical harmonic order
%     lambda : regularization parameter
%     icoverts : icosahedral sphere vertex coordinates
%     icostruct : icosahedral sphere surface struct
%     beta2ico : matrix mapping from FOD to ico3 surface
%     nb : number of parameters (betas)
%     A : RSI forward matrix           size = [nf,nb]
%     Ainv : RSI inverse matrix        size = [nb,nf]
%     B : tensor forward matrix        size = [nf,7]
%     Binv : tensor inverse matrix     size = [7,nf]
%
%  NOTE: may also supply MFfit structure containing volMF
%
% Optional Parameters:
%   'vol': brain volume matrix to be used as overlay
%     if not supplied, will use FODfit.volb0
%     {default = []}
%   'clim': vector of lower and upper bounds for color scale used by imagesc
%     {default = []}
%   'plane': slice plane (1=x, 2=y, 3=z)
%     {default = 3}
%   'slices': vector of slice indices to plot
%     if empty, plot all slices (of 'plane')
%     {default = []}
%   'xrange': vector of low and high index numbers in first (x) dimension
%     {default = []}
%   'yrange': vector of low and high index numbers in second (y) dimension
%     {default = []}
%   'zrange': vector of low and high index numbers in third (z) dimension
%     {default = []}
%   'tif_flag': [0|1] whether to save images as tif files
%     {default = 1}
%   'outdir': output directory
%     if not specified, will write to current directory
%     {default = []}
%   'outstem': output file stem
%     {default = 'FOD'}
%   'tif_dpi': resolution of tif files (dots per inch)
%     {default = 300}
%   'fig_size': figure size in inches
%     if not specified, will use default
%     {default = []}
%   'visible_flag': [0|1] display images on screen
%     ignored if tif_flag = 0 (always visible)
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% see also: rsi_fit_FOD, rsi_fit_MFOD
%
% Created:  05/18/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

% based on code by Nate White

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms( varargin, { ...
  'vol',[],[],...
  'clim',[],[],...
  'slices',[],[],...
  'xrange',[],[],...
  'yrange',[],[],...
  'zrange',[],[],...
  'plane',3,[1:3],...
  'tif_flag',true,[false true],...
  'outdir',pwd,[],...
  'outstem','FOD',[],...
  'tif_dpi',300,[10,10000],...
  'fig_size',[],[],...
  'visible_flag',true,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volmask = FODfit.volmask_dilated;

if isfield(FODfit,'volMF')
  volFOD = FODfit.volMF;
else
  volFOD = FODfit.volFOD;
end;

volsz = FODfit.volsz(1:3);

if ~isempty(parms.xrange)
  if any(~ismember(parms.xrange,[1:volsz(1)]))
    error('xrange must be between 1 and %d',volsz(1));
  end;
end;

if ~isempty(parms.yrange)
  if any(~ismember(parms.yrange,[1:volsz(2)]))
    error('yrange must be between 1 and %d',volsz(2));
  end;
end;

if ~isempty(parms.zrange)
  if any(~ismember(parms.zrange,[1:volsz(3)]))
    error('zrange must be between 1 and %d',volsz(3));
  end;
end;

if isempty(parms.vol);
  vol = FODfit.volb0;
else
  vol = parms.vol;
end

tmp_volsz = size(vol);
if any(volsz ~= tmp_volsz)
  error('mismatch in volume sizes');
end;

switch parms.plane
  case 1
    permvec = [3,2,1];
    xrange = parms.zrange;
    yrange = parms.yrange;
    zrange = parms.xrange;
    planestr = 'x';
  case 2
    permvec = [3,1,2];
    xrange = parms.zrange;
    yrange = parms.xrange;
    zrange = parms.yrange;
    planestr = 'y';
  case 3
    permvec = [1,2,3];
    xrange = parms.xrange;
    yrange = parms.yrange;
    zrange = parms.zrange;
    planestr = 'z';
end;

if parms.plane~=3
  volFOD = permute(volFOD,[permvec,4]);
  vol = permute(vol,permvec);
  volmask = permute(volmask,permvec);
  volsz = volsz(permvec);
end;

nx = volsz(1);
ny = volsz(2);
nz = volsz(3);

if isempty(xrange)
  xrange = [1:nx];
elseif numel(xrange)>=2
  xrange = [min(xrange):max(xrange)];
end;
if isempty(yrange)
  yrange = [1:ny];
elseif numel(yrange)>=2
  yrange = [min(yrange):max(yrange)];
end;
if isempty(zrange)
  zrange = [1:nz];
end;

if isempty(parms.slices)
  parms.slices = zrange;
else
  if any(~ismember(parms.slices,[1:nz]))
    error('slices must be between 1 and %d',nz);
  end;
end;

npf = (FODfit.SH_order+1)*(FODfit.SH_order+2)/2;

[TH,PHI,RHO] = ...
  cart2sph(FODfit.icoverts(:,1),FODfit.icoverts(:,2),FODfit.icoverts(:,3));
FOD.faces = FODfit.icostruct.faces;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.tif_flag
  mmil_mkdir(parms.outdir);
end;

for s=1:length(parms.slices)
  z = parms.slices(s);
  if parms.tif_flag
    fname_tif = sprintf('%s/%s_%s%d.tif',...
      parms.outdir,parms.outstem,planestr,z);
    if exist(fname_tif,'file') && ~parms.visible_flag && ~parms.forceflag
      continue;
    end;
  end;
  slicemask = squeeze(volmask(xrange,yrange,z));
  vox = find(slicemask);
  nvox = length(vox);
  if nvox>0
    % restricted FOD
    h1 = figure; clf; % set(h1,'color',[1 1 1]);
    im = squeeze(vol(xrange,yrange,z));
    if ~isempty(parms.clim)
      imagesc(im,parms.clim);
    else
      imagesc(im);
    end;
    colormap gray; axis image; axis off; hold on
    for k = 1:nvox
      [i j] = ind2sub([length(xrange) length(yrange)],vox(k));
      beta = squeeze(volFOD(xrange(i),yrange(j),z,:));
      icoF = max(0,FODfit.beta2ico*beta(1:npf));
      [val ind] = max(icoF);
      V0 = FODfit.icoverts(ind,:)';
      %% todo: normalize by max across all voxels?
      %% todo: scale by F2?
      [tmpX,tmpY,tmpZ] = sph2cart(TH,PHI,0.5*icoF/max(icoF));
      switch parms.plane
        case 1
          X = tmpZ;
          Y = tmpY;
          Z = tmpX;
        case 2
          X = tmpZ;
          Y = tmpX;
          Z = tmpY;
        case 3
          X = tmpX;
          Y = tmpY;
          Z = tmpZ;              
      end;      
      FOD.vertices = [Y+j,X+i,Z+2];
      p = patch(FOD);
      %% todo: normalize by max across all voxels
      %% todo: scale by F2?
      tcolor = repmat(max(icoF/max(icoF),eps),1,3).*abs(FODfit.icoverts);
      set(p,'FaceVertexCData',tcolor,'LineStyle','none','FaceColor','interp','Facelighting','phong');
    end
    axis ij; axis image;
  end

  % save image
  if parms.tif_flag
    if ~parms.visible_flag, set(gcf,'Visible','off'); end;
    if ~isempty(parms.fig_size)
      if length(parms.fig_size)==1
        parms.fig_size = [parms.fig_size parms.fig_size];
      end;
      set(gcf, 'PaperUnits', 'inches');
      set(gcf, 'PaperSize', [parms.fig_size(1) parms.fig_size(2)]);
      set(gcf, 'PaperPositionMode', 'manual');
      set(gcf, 'PaperPosition', [0 0 parms.fig_size(1) parms.fig_size(2)]);
    end;
    print(gcf,'-dtiff',fname_tif,sprintf('-r %d',parms.tif_dpi));
    if ~parms.visible_flag, close(gcf); end;
  end;  
end

