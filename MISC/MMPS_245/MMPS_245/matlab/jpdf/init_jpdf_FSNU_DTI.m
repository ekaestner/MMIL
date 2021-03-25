% Rcnt Mod: 03/09/12 by Don Hagler
% Last Mod: 06/09/14 by Don Hagler

ContainerRootDir = '/home/mmildev/DTI/proc';
ContainerDir = 'MRIPROC_ma01_060908_20060908.141137_1';
FSRootDir = '/home/mmildev/data/DTI/fsurf';
subjname = 'FSURF_ma01_060908_20060908.141137_1';
DTIScanNum = 3;
infix = 'mc_reg3_ecc_B0uw';

parms = [];
parms.type1 = 'FSNU';
parms.type2 = 'DTI';
parms.jpdfdir = [pwd '/atlases/jpdfs'];
parms.initflag = 1;

resinfo_fname = sprintf('%s/%s/MPR_resinfo.mat',ContainerRootDir,ContainerDir);
if exist(resinfo_fname,'file')
  load(resinfo_fname);
else
  fprintf('%s: WARNING: %s not found\n',mfilename,resinfo_fname);
  M_res_to_mpr = eye(4);
end;

fname_F = sprintf('%s/%s/DTI%d_%s.mgh',ContainerRootDir,ContainerDir,DTIScanNum,infix);
fname_S = sprintf('%s/%s/mri/nu.mgz',FSRootDir,subjname);
maskfname = sprintf('%s/%s/mri/nu-brainmask.mgz',FSRootDir,subjname);

volmask = [];
if exist(maskfname,'file')
  fprintf('%s: loading brain mask volume %s...\n',mfilename,maskfname);
  volmask = ctx_load_mgh(maskfname);
end;

[mvolF,M_F]=fs_load_mgh(fname_F,[],1);
[mvolS,M_S]=fs_load_mgh(fname_S);
volS = ctx_mgh2ctx(mvolS,M_S);
volF = ctx_mgh2ctx(mvolF,M_F);

% resample F to S
M = eye(4);
M(1,4) = -2; % L/R
M(2,4) = 4; % A/P
M(3,4) = 11; % I/S

volF.Mvxl2lph = inv(M)*inv(M_res_to_mpr)*volF.Mvxl2lph;

if 0
  volFres = vol_resample_pad(volF,volS,eye(4),1);

  viewplane = 'sag';
  % adjust M_F to get good initial alignment
  r=128;
  c=128;
  s=100;

  switch viewplane
    case 'sag'
      tmpF = squeeze(volFres.imgs(s,:,:));
      tmpS = squeeze(volS.imgs(s,:,:));
    case 'cor'
      tmpF = squeeze(volFres.imgs(:,:,r))';
      tmpS = squeeze(volS.imgs(:,:,r))';
    case 'hor'
      tmpF = flipud(squeeze(volFres.imgs(:,c,:))');
      tmpS = flipud(squeeze(volS.imgs(:,c,:))');
  end;
  for i=1:100
    imagesc(tmpF);
    pause;
    imagesc(tmpS);
    pause;
  end;
end;

args = mmil_parms2args(parms);
[M_S_to_F,volmask]=mmil_rbreg_vol2vol_jpdf(volS,volF,args{:});

if 0
  volFres2 = vol_resample_pad(volF,volS,M_S_to_F,1);

  viewplane = 'sag';
  r=128;
  c=128;
  s=100;

  switch viewplane
    case 'sag'
      tmpF = squeeze(volFres2.imgs(s,:,:));
      tmpS = squeeze(volS.imgs(s,:,:));
    case 'cor'
      tmpF = squeeze(volFres2.imgs(:,:,r))';
      tmpS = squeeze(volS.imgs(:,:,r))';
    case 'hor'
      tmpF = flipud(squeeze(volFres2.imgs(:,c,:))');
      tmpS = flipud(squeeze(volS.imgs(:,c,:))');
  end;
  for i=1:100
    imagesc(tmpF);
    pause;
    imagesc(tmpS);
    pause;
  end;
end;

if 0
  viewplane = 'sag';
  r=128;
  c=128;
  s=100;

  switch viewplane
    case 'sag'
      tmpM = squeeze(volmask.imgs(s,:,:));
      tmpS = squeeze(volS.imgs(s,:,:));
    case 'cor'
      tmpM = squeeze(volmask.imgs(:,:,r))';
      tmpS = squeeze(volS.imgs(:,:,r))';
    case 'hor'
      tmpM = flipud(squeeze(volmask.imgs(:,c,:))');
      tmpS = flipud(squeeze(volS.imgs(:,c,:))');
  end;
  for i=1:100
    imagesc(tmpM);
    pause;
    imagesc(tmpS);
    pause;
  end;
end;

