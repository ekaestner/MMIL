% Early Mod:03/09/12 by Don Hagler
% Last Mod: 09/09/12 by Don Hagler

%% todo: use mmil_rbreg_vol2vol_jpdf instead of mmil_atlas_jpdfreg_func2struct

dir_atlas = [getenv('MMPS_DIR') '/atlases'];
ContainerRootDir = '/home/mmildev/data/MMILDB/DTI/Containers';
ContainerDir = 'MRIPROC_epd08_070126b_20070126.100512_1';
ftype = 'GECAL';
stype = 'MPR';
reinitflag = 0;

fname_jpdf = sprintf('%s/jpdfs/%s_%s_jpdf.mat',dir_atlas,stype,ftype);
if exist(fname_jpdf,'file'), delete(fname_jpdf); end;

fname_F = sprintf('%s/%s/GEB1CAL1.mgh',ContainerRootDir,ContainerDir);
fname_S = sprintf('%s/%s/MPR1.mgh',ContainerRootDir,ContainerDir);

volmask = [];

[mvolF,M_F]=fs_load_mgh(fname_F,[],1);
[mvolS,M_S]=fs_load_mgh(fname_S);
volS = ctx_mgh2ctx(mvolS,M_S);

% resample F to S
M = eye(4);
M(1,4) = 0; % L/R
M(2,4) = 0; % A/P
M(3,4) = 0; % I/S
volF = ctx_mgh2ctx(mvolF,M*M_F);
volFres = vol_resample_pad(volF,volS,eye(4),1);
viewplane = 'hor';
if 0
  % adjust M_F to get good initial alignment
  r=120;
  c=128;
  s=128;

  switch viewplane
    case 'cor'
      tmpF = squeeze(volFres.imgs(r,:,:));
      tmpS = squeeze(volS.imgs(r,:,:));
    case 'sag'
      tmpF = squeeze(volFres.imgs(:,:,s))';
      tmpS = squeeze(volS.imgs(:,:,s))';
    case 'hor'
      tmpF = squeeze(volFres.imgs(:,c,:));
      tmpS = squeeze(volS.imgs(:,c,:));
  end;
  for i=1:100
    imagesc(tmpF);
    pause;
    imagesc(tmpS);
    pause;
  end;
end;

[M_S_to_F,RegInfo,volmask]=mmil_atlas_jpdfreg_func2struct(volF,volS,...
  'ftype',ftype,'stype',stype,'volmask',volmask,'reinitflag',reinitflag)

volFres = vol_resample_pad(volF,volS,M_S_to_F,1);
fname_out = sprintf('%s/%s/GEB1CAL1_res.mgh',ContainerRootDir,ContainerDir);
mvolFres =volFres.imgs;
fs_save_mgh(mvolFres,fname_out,M_S);

