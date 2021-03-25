function MMIL_Process_PET_Exam(PETrawdir,PETprocdir,MRIdir,forceflag)
%function MMIL_Process_PET_Exam(PETrawdir,PETprocdir,MRIdir,forceflag)
%
% Early Mod:  10/12/09 by Don Hagler
% Last Mod:   04/12/12 by Don Hagler
%

if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

fname_atlas = [getenv('MMPS_DIR') '/atlases/jpdfs/MRI_PET_jpdf.mat'];

PET_files = dir(sprintf('%s/vols*',PETrawdir));
volnames = {};
petregstems = {};
for j=1:length(PET_files)
  volnames{end+1} = PET_files(j).name(6:end-4);
end

for j=1:length(volnames)
  petregstems{end+1} = sprintf('PET_reg_%s',volnames{j});
end

if ~exist(PETprocdir,'dir'), mkdir(PETprocdir); end

for h=1:length(volnames)
  fname_PET_out = sprintf('%s/%s.mgh',PETprocdir,petregstems{h});
  fname_PET_info = sprintf('%s/%s.mat',PETprocdir,petregstems{h});
  fname_PET_info_txt = sprintf('%s/%s.txt',PETprocdir,petregstems{h});
  try
    cinfo = sprintf('%s/ContainerInfo_%s.mat',PETrawdir,volnames{h});
    load(cinfo);
  catch
    fprintf(1,'%s: ERROR: cannot read file %s\n',mfilename,cinfo);
    return
  end

  if ~isfield(ContainerInfo,'PIBflag') % check PIB flag
    fprintf('%s: ERROR: missing PIBflag',mfilename);
    return;
  end;
  if ContainerInfo.PIBflag
    fprintf('%s: WARNING: skipping registration for PIB data',mfilename);
    return;
  end;

  load(fname_atlas);

  if exist(fname_PET_out,'file') && exist(fname_PET_info,'file') &&...
     exist(fname_PET_info_txt,'file') && ~forceflag
    return;
  end;

  PETfname = sprintf('%s/vols_%s.mat',PETrawdir,volnames{h});
  load(PETfname);
  if ~ismember(length(vols),[1,5,6,7])
    fprintf('%s: ERROR: incorrect number of PET volumes (length(vols) = %d ~= (1, 5, 6, or 7))\n',...
      mfilename,length(vols));
    return;
  end
  % todo: if nframes = 33, use last 3 only
  vol_PET = 0;
  M = vols{1}.M;
  % register each frame to first
  fprintf('%s: registering PET frames to each other...\n',mfilename);
  tmp_vol = max(0,double(vols{1}.vol));
  vol_PET = ctx_mgh2ctx(tmp_vol,M);
  volmask = vol_PET;
  volmask.imgs = zeros(size(volmask.imgs));
  volsz = size(volmask.imgs);
  x = round(volsz(1)/4);
  y = round(volsz(2)/4);
  z = round(volsz(3)/4);
  volmask.imgs(x:3*x,y:3*y,z:3*z)=1;
  for j = 2:length(vols)
    vol_tmp = max(0,double(vols{j}.vol));
    vol_tmp = ctx_mgh2ctx(vol_tmp,M);
    [M_v1_to_v2, min_cost] = rbreg_vol2vol_icorr_mask(vol_PET,vol_tmp,volmask,0,4);
    vol_tmp = vol_resample(vol_tmp,vol_PET,M_v1_to_v2,2);
    vol_tmp.imgs(find(isnan(vol_tmp.imgs))) = 0;
    vol_PET.imgs = vol_PET.imgs + vol_tmp.imgs;
  end
  vol_PET.imgs = vol_PET.imgs/length(vols);
  [hc,bv] = hist(vol_PET.imgs(:),1000);

  % Compute PET scalefactor
  sf = ComputeHistScalefactor(hc,bv,PEThc,PETbv,PETbn0);
  vol_PET.imgs = sf*vol_PET.imgs;

  MRIfname = sprintf('%s/mri/nu.mgh',MRIdir);
  if ~exist(MRIfname,'file')
    mgzfname = sprintf('%s/mri/nu.mgz',MRIdir);
    cmd = sprintf('mri_convert %s %s\n',mgzfname,MRIfname);
    system(cmd);
  end
  vol_MRI = ctx_load_mgh(MRIfname);
  [hc,bv] = hist(vol_MRI.imgs(:),1000);

  % Compute MRI scalefactor
  sf = ComputeHistScalefactor(hc,bv,MRIhc,MRIbv,MRIbn0);
  vol_MRI.imgs = sf*vol_MRI.imgs;

  % Register PET to MRI data
  fprintf('%s: registering PET to MRI...\n',mfilename);
  volmask = mmil_dct_brainmask(vol_MRI);
  Mreg0 = eye(4);
  tic
  [M_MRI_to_PET, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol_MRI, vol_PET, volmask, Mreg0, 1, 4, [], [], [], [], jpdf_mask, bins1_mask, bins2_mask);
  toc
  vol_PET_res = vol_resample(vol_PET, vol_MRI, M_MRI_to_PET, 1);
  vol_PET_res.imgs(find(isnan(vol_PET_res.imgs))) = 0;
  %showVol(vol_MRI,vol_PET_res)

  % Rinse and repeat...
  for iter = 1:2
    % Calculate jpdf from registered MRI & PET volumes
    [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=vols_jhist_mask_amd(vol_MRI, vol_PET_res, volmask, 1, eye(4), eye(4));
    jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function

    Mreg0 = M_MRI_to_PET;
    [M_MRI_to_PET, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol_MRI, vol_PET, volmask, Mreg0, 1, 4, [], [], [], [], jpdf_mask, bins1_mask, bins2_mask);
    vol_PET_res = vol_resample(vol_PET, vol_MRI, M_MRI_to_PET, 1);
    vol_PET_res.imgs(find(isnan(vol_PET_res.imgs))) = 0;
  end

  PETinfo.MRIdir = MRIdir;
  PETinfo.PETrawdir = PETrawdir;
  PETinfo.M_MRI_to_PET = M_MRI_to_PET;
  PETinfo.min_cost = min_cost;
  PETinfo.jpdf_mask = jpdf_mask;
  PETinfo.bins1_mask = bins1_mask;
  PETinfo.bins2_mask = bins2_mask;
  ctx_save_mgh(vol_PET_res,fname_PET_out);
  save(fname_PET_info,'PETinfo');
  fid = fopen(fname_PET_info_txt,'w');
  if fid>0
    fprintf(fid,'MRIfname = %s\n',MRIfname);
    fclose(fid);
  end
end
