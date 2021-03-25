function DTI_MMIL_Average_Atlas_Fiber(RootDirs,StudyInfo,f,varargin);
%function DTI_MMIL_Average_Atlas_Fiber(RootDirs,StudyInfo,f,varargin);
%
% Required Input:
%  RootDirs: RootDirs for the ProjID
%  StudyInfo: StudyInfo for the ProjID
%  f: fiber number
%
% Optional Parameters:
%  'outdir': path for storing the DTI atlas output
%    {default = pwd}
%  'tensor_smooth_sigma: smoothing applied to tensors after averaging
%    {default = 5}
%  'countflag: use fiber counts as input rather than fiber masks
%    {default = 1}
%  'first_only_flag: use first eigen vector only
%    {default = 1}
%  'min_tensor_count: for calculating mean after smoothing
%    {default = 0.5}
%  'smf': threshold used for the standardized tensors
%    {default = 1e-6}
%  'forceflag: force overwrite of existing output files
%    {default: 0}
%
% Created:  03/07/11 by Vijay Venkatraman
% Last Mod: 08/07/13 by Don Hagler
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'atlasdir',pwd,[],...
  'tensor_smooth_sigma',5,[],... 
  'countflag',1,[],...
  'first_only_flag',1,[],...
  'min_tensor_count',0.5,[],...
  'smf',1e-6,[],...
  'forceflag',false,[true false],...
... % hidden parameters
  'tensor_components',[1,5,9,2,3,6],[1:9],... % do not change
  'warpdir','WarpToAtlas',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms.num_tensor_components = length(parms.tensor_components);

nstudies = length(StudyInfo);
if ~nstudies, error('StudyInfo is empty'); end;

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: reduce redundancy between this section and next

% average fiber masks or fiber counts
M = [];
avgfiber_dir = sprintf('%s/avgfiber',parms.outdir);
mmil_mkdir(avgfiber_dir);
if parms.countflag % 1 means count and 0 means mask
  fname_count = sprintf('%s/fiber_%03d_count_countatlas.mat',avgfiber_dir,f);
  fname_mean = sprintf('%s/fiber_%03d_mean_countatlas.mat',parms.outdir,f);
  fname_log = sprintf('%s/fiber_%03d_countatlas.log',avgfiber_dir,f);
else
  fname_count = sprintf('%s/fiber_%03d_count_maskatlas.mat',avgfiber_dir,f);
  fname_mean = sprintf('%s/fiber_%03d_mean_maskatlas.mat',parms.outdir,f);
  fname_log = sprintf('%s/fiber_%03d_maskatlas.log',avgfiber_dir,f);
end;
if (~exist(fname_mean,'file') | ...
   ~exist(fname_count,'file') | ...
   ~exist(fname_log,'file')) | parms.forceflag
   nsum = 0;
   volsum = single(0);
   %volsum2 = single(0);
   fprintf('%s: averaging fiber %d...\n',mfilename,f);
   flog = fopen(fname_log,'wt');
   if flog == -1
     error('failed to open %s for writing',fname_mean);
   end;
   fprintf(flog,'fiber %d:\n',f);
   for i= 1:nstudies
     SubjID = StudyInfo(i).SubjID;
     PROCContainerDir = StudyInfo(i).proc_dti;
     PROCContainerPath = sprintf('%s/%s',RootDirs.proc_dti,PROCContainerDir);
     ATLContainerPath = sprintf('%s/%s/%s/fiber_atlas',...
      RootDirs.proc_dti,PROCContainerDir,parms.warpdir);
     n = regexp(char(PROCContainerDir),'DTIPROC_(?<SessID>[^]+)','names');
     if isempty(n)
       fprintf('%s: WARNING: %s is an invalid DTIPROC container name',mfilename,PROCContainerDir);
       continue;
     end;
     if isempty(PROCContainerDir) | ~exist(PROCContainerPath,'dir') |...
       ~exist(ATLContainerPath,'dir')
       continue;
     end;
     if parms.countflag
       % load fiber count
       fname_fiber = sprintf('%s/fiber_%03d_count_atlas.mat',ATLContainerPath,f);
     else
       % load fiber mask
       fname_fiber = sprintf('%s/fiber_%03d_mask_atlas.mat',ATLContainerPath,f);
     end;
     if ~exist(fname_fiber,'file'), continue; end;
     fprintf('%s: loading %s...\n',mfilename,fname_fiber);
     tic
     [vol,tmpM] = mmil_load_sparse(fname_fiber);
     toc
     if ~isempty(vol)
       if ~isempty(tmpM) & isempty(M), M = tmpM; end;
       if ~parms.countflag
         %scale to 1
         maxval = max(vol(:));
         vol = vol/maxval;
       else
         %% todo: normalize to median?
         %% todo: normalize to total number of voxels in fiber?
         %% todo: normalize to total number of streamlines in all voxels in fiber?
       end;
        nsum = nsum + 1;
        volsum = volsum + vol;
        %volsum2 = volsum2 + vol.^2;
     end;    
     fprintf(flog,'total number of subjects: %s \n',int2str(nsum));
     if nsum<=0
       fprintf(flog,'zero subjects in fiber %03d average\n',f);
       fclose(flog);
       error('zero subjects in fiber %03d average',f);
     elseif nsum==1
       fprintf('%s: WARNING: only one subject in fiber %03d average\n',...
         mfilename,f);
       volmean = volsum;
       %volstd = zeros(size(volmean));
     else
       fprintf('%s: %d subjects included in average\n',mfilename,nsum);
       volmean = volsum / (eps+nsum);
       %volstd = sqrt((nsum.*volsum2 - volsum.^2)./(eps+nsum.*(nsum-1)));
     end;
     if parms.countflag
      % normalize average fiber count
      volmean = volmean/(max(volmean(:))+eps);
      if 0
       % apply exponent so that median gets value of 0.5
       vals = volmean(volmean>0);
       medval = median(vals);
       p = log(0.5)/log(medval);
       volmean = volmean.^p;
      end;
     end;
     mmil_save_sparse(volsum,fname_count,M);
     mmil_save_sparse(volmean,fname_mean,M);
   end;
   fclose(flog);
end;
clear vol volmean volsum nsum fname_log fname_mean fname_count fname_fiber;
clear tmp_M M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% average tensors (components done separately to save memory)
M = [];
avgtensor_dir = sprintf('%s/avgtensor',parms.outdir);
mmil_mkdir(avgtensor_dir);
for i=1:parms.num_tensor_components
  t = parms.tensor_components(i);
  if parms.first_only_flag
    fname_count = sprintf('%s/fiber_%03d_count_tensor%dV0_atlas.mat',avgtensor_dir,f,i);
    fname_sum = sprintf('%s/fiber_%03d_sum_tensor%dV0_atlas.mat',avgtensor_dir,f,i);
    fname_log = sprintf('%s/fiber_%03d_tensor%dV0_atlas.log',avgtensor_dir,f,i);
  else
    fname_count = sprintf('%s/fiber_%03d_count_tensor%d_atlas.mat',avgtensor_dir,f,i);
    fname_sum = sprintf('%s/fiber_%03d_sum_tensor%d_atlas.mat',avgtensor_dir,f,i);
    fname_log = sprintf('%s/fiber_%03d_tensor%d_atlas.log',avgtensor_dir,f,i);
  end;
  if (~exist(fname_sum,'file') | ...
     ~exist(fname_count,'file') | ...
     ~exist(fname_log,'file')) | parms.forceflag
    flog = fopen(fname_log,'wt');
    if flog==-1
      error('failed to open %s for writing',fname_log);
    end;
    fprintf(flog,'fiber %d:\n',f);
    nsum = 0;
    voln = single(0);
    volsum = single(0);
    %volsum2 = single(0);
    fprintf('%s: averaging tensor component %d for fiber %d...\n',...
      mfilename,i,f);
    for s = 1:nstudies
      SubjID = StudyInfo(s).SubjID;
      PROCContainerDir = StudyInfo(s).proc_dti;
      PROCContainerPath = sprintf('%s/%s',RootDirs.proc_dti,PROCContainerDir);
      ATLContainerPath = sprintf('%s/%s/%s/fiber_atlas',...
        RootDirs.proc_dti,PROCContainerDir,parms.warpdir);
      n = regexp(char(PROCContainerDir),'DTIPROC_(?<SessID>[^]+)','names');
      if isempty(n)
        fprintf('%s: WARNING: %s is an invalid DTIRPC container name',mfilename,PROCContainerDir);
        continue;
      end;
      if isempty(PROCContainerDir) |~exist(PROCContainerPath,'dir') |...
        ~exist(ATLContainerPath,'dir')
        continue;
      end;
      if parms.countflag
        % load fiber count
        fname_fiber = sprintf('%s/fiber_%03d_count_atlas.mat',ATLContainerPath,f);
      else
        % load fiber mask
        fname_fiber = sprintf('%s/fiber_%03d_mask_atlas.mat',ATLContainerPath,f);
      end;
      if ~exist(fname_fiber,'file'), continue; end;
      fprintf('%s: loading %s...\n',mfilename,fname_fiber);
      tic
      [vol,tmpM] = mmil_load_sparse(fname_fiber);
      toc
      if isempty(vol), continue; end;
      volmask = zeros(size(vol));
      volmask(vol>0) = 1;
      clear vol;

      % load tensors
      if parms.first_only_flag
        fname_tensor = sprintf('%s/fiber_%03d_tensorV0_atlas.mat',ATLContainerPath,f);
      else
        fname_tensor = sprintf('%s/fiber_%03d_tensor_atlas.mat',ATLContainerPath,f);
      end;
      if ~exist(fname_tensor,'file'), continue; end;
      fprintf('%s: loading %s...\n',mfilename,fname_tensor);
      tic
      [tmpvol,tmpM,tmp_sz] = mmil_load_sparse(fname_tensor);
      % check if num frames = 6 or 9
      if tmp_sz(4)==6
        volT = tmpvol(:,:,:,i); % load one component
      elseif tmp_sz(4)==9
        volT = tmpvol(:,:,:,t); % load one component
      end;
      toc
      if isempty(volT), continue; end;
      if ~isempty(tmpM) & isempty(M), M = tmpM; end;
      volT = volmask.*volT; % only take tensor from voxels inside fiber
      nsum = nsum + 1;
      voln = voln + volmask;
      volsum = volsum + volT;
      %volsum2 = volsum2 + volT.^2;
      fprintf(flog,'\t%s\n',SubjID);
    end
    fprintf(flog,'total number of subjects= %d \n',nsum);
    if nsum<=0
      fprintf(flog,'zero subjects in fiber %03d average\n',f);
      fclose(flog);
      error('zero subjects in fiber %03d average',f);
    elseif nsum==1
      fprintf('%s: WARNING: only one subject in fiber %03d tensor average\n',...
        mfilename,f);
    else
      fprintf('%s: %d subjects included in tensor sum\n',mfilename,nsum);
    end;
    mmil_save_sparse(voln,fname_count,M);
    mmil_save_sparse(volsum,fname_sum,M);
    fclose(flog);
  end;
end;
clear voln volsum volmask volT fname_fiber fname_tensor fname_count fname_sum fname_log;
clear tmp_vol tmpM tmpsz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine tensor components into single file (and smooth)
M = [];
if parms.first_only_flag
  fname_tensor = sprintf('%s/fiber_%03d_mean_tensorV0_atlas_sm%0.2f.mat',...
    avgtensor_dir,f,parms.tensor_smooth_sigma);
  fname_tensor_norm = sprintf('%s/fiber_%03d_mean_tensorV0_atlas_sm%0.2f_norm.mat',...
    parms.outdir,f,parms.tensor_smooth_sigma);
else
  fname_tensor = sprintf('%s/fiber_%03d_mean_tensor_atlas_sm%0.2f.mat',...
    avgtensor_dir,f,parms.tensor_smooth_sigma);
  fname_tensor_norm = sprintf('%s/fiber_%03d_mean_tensor_atlas_sm%0.2f_norm.mat',...
    parms.outdir,f,parms.tensor_smooth_sigma);
end;
if ~exist(fname_tensor,'file') | parms.forceflag
  vol_tensor = [];
  fprintf('%s: combining tensor components for fiber %d...\n',...
    mfilename,f);
  for i=1:parms.num_tensor_components
    t = parms.tensor_components(i);
    if parms.first_only_flag
      fname_count = sprintf('%s/fiber_%03d_count_tensor%dV0_atlas.mat',avgtensor_dir,f,i);
      fname_sum = sprintf('%s/fiber_%03d_sum_tensor%dV0_atlas.mat',avgtensor_dir,f,i);
    else
      fname_count = sprintf('%s/fiber_%03d_count_tensor%d_atlas.mat',avgtensor_dir,f,i);
      fname_sum = sprintf('%s/fiber_%03d_sum_tensor%d_atlas.mat',avgtensor_dir,f,i);
    end;
    if ~exist(fname_count,'file')
      error('file %s not found',fname_count);
    end;
    if ~exist(fname_sum,'file')
      error('file %s not found',fname_sum);
    end;
    % load atlas data
    [vol_count,M_count] = mmil_load_sparse(fname_count);
    [vol_sum,tmpM] = mmil_load_sparse(fname_sum);
    if isempty(vol_sum)
      error('file %s contains empty volume',fname_sum);
    end;
    if ~isempty(tmpM) & isempty(M), M = tmpM; end;
    if isempty(vol_tensor)
      vol_tensor = zeros([size(vol_sum),parms.num_tensor_components]);
    end;
    % smooth count and mean
    if parms.tensor_smooth_sigma>0
      ps = sqrt(sum(M(1:3,1:3).^2,1)); % voxel dimensions in mm
      vol_sum = mmil_smooth3d(vol_sum,...
        parms.tensor_smooth_sigma/ps(1),...
        parms.tensor_smooth_sigma/ps(2),...
        parms.tensor_smooth_sigma/ps(3));
      vol_count = mmil_smooth3d(vol_count,...
        parms.tensor_smooth_sigma/ps(1),...
        parms.tensor_smooth_sigma/ps(2),...
        parms.tensor_smooth_sigma/ps(3));
    end;
    vol_mean = vol_sum./(vol_count + eps);
    vol_mean(vol_count<parms.min_tensor_count) = 0;
    % insert into combined tensor
    fprintf('%s: calculated mean for tensor component %d\n',mfilename,t);
    vol_tensor(:,:,:,i) = vol_mean;
  end;
  if isempty(vol_tensor)
    error('vol_tensor is empty');
  end;
  mmil_save_sparse(vol_tensor,fname_tensor,M);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_tensor_norm,'file') | parms.forceflag
  fprintf('%s: normalizing tensors for fiber %d...\n',...
    mfilename,f);
  if ~exist('vol_tensor','var') | isempty(vol_tensor)
    [vol_tensor,M]=mmil_load_sparse(fname_tensor);
  end;
  % normalize tensors by product with first eigen vector
  vol_tensor_std = sum(abs(vol_tensor),4);
  volsz = size(vol_tensor_std);
  voxels = find(vol_tensor_std>parms.smf);
  nvox_ROI = length(voxels);
  fprintf('%s: %d voxels with non-zero tensor\n',mfilename,nvox_ROI);
  for m=1:nvox_ROI
    ind=voxels(m);
    [i,j,k] = ind2sub(volsz,ind);
    vec = squeeze(vol_tensor(i,j,k,:));
    T = dti_components2tensor(vec);
    [U,S] = eig(T);
    lamdas = real(diag(S)');
    [sorted_lamdas,ind] = sort(lamdas,2,'descend');
    v = U(:,ind(1));
    scale_fact = v'*T*v;
    if scale_fact<parms.smf
      T = zeros(3);
    else
      T = T/scale_fact;
    end;
    if any(T>ones(3)), T = zeros(3); end;
    vol_tensor(i,j,k,:) = dti_tensor2components(T,parms.num_tensor_components);
  end;
  mmil_save_sparse(vol_tensor,fname_tensor_norm,M);
end;

return;
