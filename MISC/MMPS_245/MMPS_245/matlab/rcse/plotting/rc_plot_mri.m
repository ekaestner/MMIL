function rc_plot_mri(subject,slicenum,plane);
%function rc_plot_mri(subject,slicenum,[plane]);
%
% Purpose: plotting images of MRI slices
%
% Required Input:
%   subject: freesufer subject name
%   slicenum: number between 1 and 256
%
% Optional Input:
%   plane: 'cor', 'sag', or 'hor'
%     {default = 'cor'}
%
% Early Mod: 12/18/07 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

%% todo: varargin
%% todo: rootdir, cropping
%% todo: mkdir
%% todo: use fs_load_mgh not fs_read_cor

if ~mmil_check_nargs(nargin,2), return; end;
if length(slicenum)~=1
  error('slicenum must be a single number between 1 and 256');
end;
if slicenum<1 | slicenum>256
  error('slicenum must be between 1 and 256');
end;
if ~exist('plane','var') || isempty(plane), plane = 'cor'; end;
if ~ismember(plane,{'cor','sag','hor'})
  error('plane must be cor, sag, or hor');
end;

if (strcmp(plane,'cor'))
  viewplane = 2;
elseif (strcmp(plane,'sag'))
  viewplane = 1;
elseif (strcmp(plane,'hor'))
  viewplane = 3;
end;

!mkdir -p matfiles

% set variables
crop = 1;
mridir = '.';
mristem = sprintf('mri');

if crop
  if viewplane==1 % sag
    x0 = 1;
    x1 = 256;
    y0 = 26;
    y1 = 80;
    z0 = 101;
    z1 = 155;
  elseif viewplane==2 % cor
    x0 = 88;
    x1 = 162;
    y0 = 1;
    y1 = 256;
    z0 = 92;
    z1 = 166;
  elseif viewplane==3 % hor
    x0 = 88;
    x1 = 152;
    y0 = 26;
    y1 = 90;
    z0 = 1;
    z1 = 256;
  end;
else
  % no cropping
  x0 = 1;
  x1 = 256;
  y0 = 1;
  y1 = 256;
  z0 = 1;
  z1 = 256;
end;

ax_pos = [0.05,0.05,0.9,0.9];

matname=sprintf('matfiles/mri.mat');
if exist(matname,'file')
  load(matname);
else
  SubjectsDir = deblank(getenv('SUBJECTS_DIR'));
  if(isempty(SubjectsDir))
    fprintf('%s: Cannot find SUBJECTS_DIR environment variable\n',mfilename);
    return;
  end
  fname = sprintf('%s/%s/mri/T1',SubjectsDir,subject);
  if ~exist(fname,'file')
    fprintf('%s: file %s not found\n',mfilename,fname);
    return;
  end;
  % load mri file
  fprintf('%s: loading mri volume %s...\n',mfilename,fname);
  mri_vol = fs_read_cor(fname);
  mri_vol = permute(mri_vol, [1 3 2]);
  mri_vol = flipdim(mri_vol,3);
  
  % save to mat file
  save(matname,'mri_vol');
end;

% select slice
if viewplane==1
  img = squeeze(mri_vol(slicenum,y0:y1,z0:z1));
  img = img';
  X_max = y1-y0+1;
  Y_max = z1-z0+1;
elseif viewplane==2
  img = squeeze(mri_vol(x0:x1,slicenum,z0:z1));
  img = img';
  img = fliplr(img);
  X_max = x1-x0+1;
  Y_max = z1-z0+1;
elseif viewplane==3;
  img = squeeze(mri_vol(x0:x1,y0:y1,slicenum));
  img = img';
  img = fliplr(img);
  X_max = x1-x0+1;
  Y_max = y1-y0+1;
end;

% make figure
crange = [40 125];
ax_mri = axes('position',ax_pos);
hold on;
h = imagesc(img,crange);
colormap(ax_mri,gray);
axis xy;
axis([1 X_max 1 Y_max]);
axis off;
hold off;

