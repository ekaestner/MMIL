function dti_plot_tensor(DTfit,varargin)
%function dti_plot_tensor(DTfit,[options])
%
% Purpose: plot tensor overlayed on brain images
%
% Required Parameters:
%   DTfit: structure containing tensor calculations with fields:
%     volDT : tensor fit parameters    size = [nx,ny,nz,7]
%     volb0 : average b=0 image        size = [nx,ny,nz]
%     volmask : mask volume            size = [nx,ny,nz]
%     volmask_dilated : dilated mask   size = [nx,ny,nz]
%     volsz : size of input 4D vol          = [nx,ny,nz,nf]
%     qmat : input q matrix            size = [nf,3]
%     bvals : input b values           size = [nf,1]
%     Q : tensor fit forward matrix    size = [nf,7]
%     censor_mat : censored slices per frame       size = [nz,nf]
%     censor_err : fit error per slices per frame  size = [nz,nf] 
%     censor_niter : number of censoring iterations
%
% Optional Parameters:
%   'vol': brain volume matrix to be used as overlay
%     if not supplied, will use DTfit.volb0
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
%     {default = 'DT'}
%   'tif_dpi': resolution of tif files (dots per inch)
%     {default = 300}
%   'fig_size': figure size in inches
%     if not specified, will use default
%     {default = []}
%   'visible_flag': [0|1] display images on screen
%     ignored if tif_flag = 0 (always visible)
%     {default = 1}
%
% see also: dti_fit_tensor
%
% Created:  05/19/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

% based on code from Nate White

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
  'outstem','DT',[],...
  'tif_dpi',300,[10,10000],...
  'fig_size',[],[],...
  'visible_flag',true,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volDT = DTfit.volDT(:,:,:,1:6);
volmask = DTfit.volmask_dilated;

volsz = DTfit.volsz(1:3);

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
  vol = DTfit.volb0;
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
  volDT = permute(volDT,[permvec,4]);
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

icostruct = icosaSphere(3);
icoverts = icostruct.vertices;
nverts = size(icoverts,1);

[TH,PHI,RHO] = ...
  cart2sph(icoverts(:,1),icoverts(:,2),icoverts(:,3));
DT.faces = icostruct.faces;

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
      T = dti_components2tensor(squeeze(volDT(xrange(i),yrange(j),z,:)));
      icoT = zeros(nverts,1);
      for v = 1:nverts
        icoT(v) = icoverts(v,:)*T*icoverts(v,:)';
      end
      %% todo: normalize by max across all voxels
      %% todo: scale by FA? by ADC?
      [tmpX,tmpY,tmpZ] = sph2cart(TH,PHI,0.5*icoT/max(icoT));
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
      DT.vertices = [Y+j,X+i,Z+2];
      p = patch(DT);
      %% todo: normalize by max across all voxels
      %% todo: scale by FA? by ADC?
      tcolor = repmat(max(icoT/max(icoT),eps),1,3).*abs(icoverts);
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

