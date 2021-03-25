function [vol,M,qmat,bvals,TEs,errcode]=DTI_MMIL_Load_Data(ContainerPath,varargin)
%function [vol,M,qmat,bvals,TEs,errcode]=DTI_MMIL_Load_Data(ContainerPath,[options])
%
% Usage:
%  [vol,M,qmat,bvals,TEs] = DTI_MMIL_Load_Data(ContainerPath,'key1', value1,...);
%
% Required Input Parameters:
%   ContainerPath: full path of directory containing processed diffusion data (mgh format)
%
% Optional Input Parameters:
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr_regT1'
%     {default = []}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'min_ndirs': require at least this many diffusion directions
%     {default = 6}
%   'min_bval': minimum b value a scan must have to be included
%     {default = 1}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': minimum number of b=0 volumes a scan must have
%     {default = 1}
%
% Output:
%   vol: 4D matrix containing volumes for each diffusion direction
%   M: 4x4 RAS to VOX matrix
%   qmat: 2D matrix of diffusion vectors corresponding to each frame of vol
%   bvals: vector of b-values for each frame of vol
%   TEs: vector of echo time (TE) values for each frame of vol
%
% Created:  01/01/07 by Don Hagler
% Last Mod: 09/10/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output variables
vol = [];
M = [];
qmat = [];
bvals = [];
TEs = [];
errcode = 0;

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'infix',[],[],...
  'revflag',0,[0,1,2],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
  'min_nb0',1,[],...
...
  'fnamestem','DTI',[],...
  'orient_ref','LPI',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ContainerInfo, get DTI scan info, determine valid scans, reference scans
tags = {'snums','min_nb0','min_ndirs','min_bval','flex_flag','revflag'};
args = mmil_parms2args(parms,tags);
[ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
if errcode ~= 0, return; end;

% load data and qmat files
seq_type = [];
for s=SessInfo.snums_DT
  revflag_list = set_revflag_list(ScanInfo(s),parms.revflag);
  for revflag=revflag_list
    fstem = sprintf('%s/%s%d',ContainerPath,parms.fnamestem,s);
    if revflag, fstem = [fstem '_rev']; end;
    if ~isempty(parms.infix), fstem = [fstem '_' parms.infix]; end;

    % if file does not exist, skip it
    fname = [fstem parms.ext];
    if ~exist(fname,'file')
      fprintf('%s: file %s not found, skipping...\n',mfilename,fname);
      continue;
    end;

    % load matrix of diffusion vectors (ndirs x 3)
    fname_qmat = [fstem '_qmat.mat'];
    if exist(fname_qmat,'file')
      tmp = load(fname_qmat);
      tmp_qmat = tmp.qmat;
      nq = size(tmp_qmat,1);
      if isfield(tmp,'bvals')
        tmp_bvals = tmp.bvals;
      else
        tmp_bvals = ScanInfo(s).bval * ones(nq,1);
      end;
    else % may not have done motion correction
      fprintf('%s: WARNING: qmat file %s not found, loading original diffusion vectors\n',...
        mfilename,fname_qmat);
      diffdirs = mmil_getfield(ScanInfo(s),'diffdirs',[]);
      tmp_qmat = dti_load_qmat(ScanInfo(s).DTI_Sequence_Type,...
        ScanInfo(s).nb0,ScanInfo(s).ndiffdirs,...
        'diffdirs',diffdirs,...
        'tensor_fnum',ScanInfo(s).tensor_fnum);
      tmp_bvals = [];
    end;

    % load diffusion data
    fprintf('%s: loading data from %s...\n',mfilename,fname);
    [tmp_vol,tmpM] = fs_load_mgh(fname);
    
    % may need to censor some frames
    fname_censor = [fstem '_censor.txt'];
    if exist(fname_censor,'file') && ~exist(fname_qmat,'file')
    % NOTE: if fname_qmat exists, frames were censored before motion correction
      fname_DT = [fstem '_DT.mat'];
      if exist(fname_DT,'file')
        % if fname_DT exists, vol was censored before eddy current correction
        % and censored qmat is stored in DTfit
        load(fname_DT);
        tmp_qmat = DTfit.qmat;
        tmp_bvals = DTfit.bvals;
      else
        % need to censor vol and qmat
        [tmp_vol,tmp_qmat,tmp_bvals,ind_censored_frames] = ...
          dti_censor_frames(fname_censor,tmp_vol,tmp_qmat,tmp_bvals);
      end;
    end;
    % check that number of frames matches for vol and qmat
    nframes = size(tmp_vol,4);
    nq = size(tmp_qmat,1);
    if nframes~=nq
      fprintf('%s: ERROR: number of frames in volume (%d) does not match number of diffusion directions (%d)\n',...
        mfilename,nframes,nq);
      errcode = 1;
      return;
    end;
    
    % concatenate qmat, bvals, TEs, and vol
    qmat = cat(1,qmat,tmp_qmat);
    if isempty(tmp_bvals)
      tmp_bvals = set_bvals(ScanInfo(s),tmp_qmat);
    end;
    bvals = cat(1,bvals,tmp_bvals);
    TEs = cat(1,TEs,ScanInfo(s).TE*ones(nq,1));
    tmp_vol = single(tmp_vol);
    if isempty(M) ||...
      (~SessInfo.revflag && s==SessInfo.reg_ref_for) ||...
      (SessInfo.revflag  && s==SessInfo.reg_ref_rev)
      M = tmpM;
    end;
    vol = cat(4,vol,tmp_vol);
    clear tmp_vol;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(vol)
  fprintf('%s: WARNING: no valid DTI files found in %s with infix "%s"\n',...
    mfilename,ContainerPath,parms.infix);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revflag_list = set_revflag_list(ScanInfo,revflag)
  revflag_list = [];
  if ismember(ScanInfo.pepolar,[0,1])
    revflag_list = ScanInfo.pepolar;
  elseif revflag==2
    if ismember(ScanInfo.DTI_Sequence_Type,[4,5,6])
      if ScanInfo.pepolar==2 % main scan is "forward"
        revflag_list = 0;
      else                   % main scan is "reverse"
        revflag_list = 1;
      end;
    else
      revflag_list = [0,1];
    end;
  else
    revflag_list = revflag;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bvals = set_bvals(ScanInfo,qmat)
  if ScanInfo.DTI_Sequence_Type==6
    bvals = ScanInfo.bval*sum(qmat.^2,2);
  else
    bvals = ScanInfo.bval*ones(size(qmat,1),1);
  end;
return

