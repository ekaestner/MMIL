% Created:  06/05/14 by Don Hagler
% Last Mod: 06/09/14 by Don Hagler

ContainerRootDir = '/space/mdeh4/1/halgdev/projects/dhagler/processing/SLP/proc';
ContainerDir = 'MRIPROC_EH_sleep_s2_20040621.151203.172000_1';

fname1 = 'MPR1_uw_nu.mgz';
fname2 = 'MEDIChi1_uw_nu.mgz';

parms = [];
parms.type1 = 'MPR';
parms.type2 = 'MEDICHI';
parms.jpdfdir = [pwd '/atlases/jpdfs'];
parms.initflag = 1;
parms.forceflag = 1;

fullfname1 = sprintf('%s/%s/%s',ContainerRootDir,ContainerDir,fname1);
fullfname2 = sprintf('%s/%s/%s',ContainerRootDir,ContainerDir,fname2);

[mvol1,M1]=fs_load_mgh(fullfname1);
[mvol2,M2]=fs_load_mgh(fullfname2);
if isempty(mvol1) | isempty(mvol2), return; end;
vol1 = ctx_mgh2ctx(mvol1,M1);

% resample v2 to v1
M = eye(4);
M(1,4) = 0; % L/R
M(2,4) = 0; % A/P
M(3,4) = 0; % I/S
vol2 = ctx_mgh2ctx(mvol2,M*M2);
vol2res = vol_resample_pad(vol2,vol1,eye(4),1);

viewplane = 'cor';
if 0
  % adjust M to get good initial alignment
  r=128;
  c=128;
  s=128;

  switch viewplane
    case 'sag'
      tmp1 = squeeze(vol1.imgs(r,:,:));
      tmp2 = squeeze(vol2res.imgs(r,:,:));
    case 'cor'
      tmp1 = squeeze(vol1.imgs(:,:,s))';
      tmp2 = squeeze(vol2res.imgs(:,:,s))';
    case 'hor'
      tmp1 = flipud(squeeze(vol1.imgs(:,c,:))');
      tmp2 = flipud(squeeze(vol2res.imgs(:,c,:))');
  end;
  for i=1:100
    imagesc(tmp2);
    pause;
    imagesc(tmp1);
    pause;
  end;
end;

args = mmil_parms2args(parms);
[M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_jpdf(vol1,vol2,args{:});


