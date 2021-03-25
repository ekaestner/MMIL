function [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,varargin)
%function [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,[options])
%
% Purpose: load DTI scan information, select reference scans
%   for between-scan registration and B0 unwarping, and select valid
%   scans for processing
%
% Required Input:
%   ContainerPath: full path of DTIPROC Container
%
% Optional Parameters:
%   'snums': list of scan numbers to process and use for tensor calculations
%     if empty, use all DTI scans in container
%     {default = []}
%   'revflag': [0|1|2] specify whether to process non-rev or rev data
%     0: process only forward phase-encode polarity data
%     1: process only reverse phase-encode polarity data
%     2: process both forward and reverse data
%     {default = 2}
%   'min_ndirs': minimum number of gradient directions allowed
%     for tensor calculations
%     {default = 6}
%   'min_bval': minimum b-value allowed for tensor calculations
%     {default = 1}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': minimum number of b=0 images required for tensor calculations
%     Note: some scans (e.g. 1-dir diffusion) do have a b=0 image but the
%       report 0 b=0 images in the dicom header
%     {default = 1}
%   'verbose': [0|1] display status messages and warnings
%     {default = 1}
%
% Output:
%   ScanInfo: struct array containing info about each DTI scan
%     from ContainerInfo
%   SessInfo: struct containing summary info about session, including:
%     snums_valid
%     snums_DT
%     snums_for
%     snums_rev
%     B0uw_refs_for
%     B0uw_refs_rev
%     reg_ref_for
%     reg_ref_rev
%     regT1_ref
%
%   errcode: 0 if successful, 1 if error
%
% NOTE: Regardless of snums, revflag, min_ndirs, min_bval, and flex_flag
%  any DTI scans can be used for estimating B0 distortion or as between
%  scan motion reference
%
%
% Created:  02/25/10 by Don Hagler
% Last Mod: 03/09/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScanInfo = [];
SessInfo = [];
diffdirs = [];
errcode = 0;

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'snums',[],[],...
  'revflag',2,[0,1,2],...
  'min_nb0',1,[],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
  'verbose',true,[false true],...
... % hidden
  'fnamestem','DTI',[],...
  'copy_tags',{'VisitID','StudyDate','StudyTime','StudyInstanceUID',...
    'MagneticFieldStrength','Manufacturer',...
    'ManufacturersModelName','MagneticFieldStrength'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ContainerInfo
[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode~=0, return; end;

% get info for DTI scans
ScanInfo = mmil_getfield(ContainerInfo.ScanInfo,parms.fnamestem);
if isempty(ScanInfo), return; end;

% copy fields from ContainerInfo
for t=1:length(parms.copy_tags)
  tag = parms.copy_tags{t};
  SessInfo.(tag) = ContainerInfo.(tag);
end;

% identify forward and reverse scans for B0 unwarp, valid scans
SessInfo.nscans = length(ScanInfo);
if isempty(parms.snums), parms.snums = [1:SessInfo.nscans]; end;
SessInfo.snums_for = [];
SessInfo.snums_rev = [];
SessInfo.snums_valid = [];
SessInfo.snums_DT = [];
DTI_ipp_flag = 0;
for s=1:SessInfo.nscans
  % make sure fields are set to some value
  ScanInfo(s).nb0 = mmil_getfield(ScanInfo(s),'nb0',0);
  ScanInfo(s).bval = mmil_getfield(ScanInfo(s),'bval',0);
  ScanInfo(s).ndiffdirs = mmil_getfield(ScanInfo(s),'ndiffdirs',0);
  ScanInfo(s).pepolar = mmil_getfield(ScanInfo(s),'pepolar',0);  

  % exclude scans that failed to convert
  if isfield(ScanInfo,'valid') && ~ScanInfo(s).valid
    if parms.verbose
      fprintf('%s: WARNING: skipping invalid DTI scan %d in %s\n',...
        mfilename,s,ContainerPath);
    end;
    continue;
  end;

  % exclude scans without pepolar info
  if isempty(ScanInfo(s).pepolar)
    if parms.verbose
      fprintf('%s: WARNING: skipping DTI scan %d in %s: missing pepolar info\n',...
        mfilename,s,ContainerPath);
    end;
    continue;
  end;

  % which scans are valid for processing
  SessInfo.snums_valid = [SessInfo.snums_valid,s];

  % which scans are forward and reverse phase encode polarity
  switch ScanInfo(s).pepolar
    case 0 % forward
      SessInfo.snums_for = [SessInfo.snums_for,s];
    case 1 % reverse
      SessInfo.snums_rev = [SessInfo.snums_rev,s];
    case {2,3}
      % integrated pepolar 2: reverse (first) + forward (remaining)
      % integrated pepolar 3: forward (first) + reverse (remaining)
      SessInfo.snums_for = [SessInfo.snums_for,s];
      SessInfo.snums_rev = [SessInfo.snums_rev,s];
      if ismember(ScanInfo(s).DTI_Sequence_Type,[4,5,6])
        DTI_ipp_flag = ScanInfo(s).pepolar;
      end;
  end;
  
  % if scan is DTI_flex type
  if ScanInfo(s).DTI_Sequence_Type == 6
    flex_flag = 1;
  else
    flex_flag = 0;
  end;

  % which scans are valid for tensor calculations
  if ismember(s,parms.snums) &&...
    (ScanInfo(s).nb0>=parms.min_nb0) &&...
    (ScanInfo(s).bval >= parms.min_bval) &&...
    (ScanInfo(s).ndiffdirs >= parms.min_ndirs) &&...
    ((ismember(ScanInfo(s).pepolar,[0,2]) &&...
     ismember(parms.revflag,[0,2])) ||...
    (ismember(ScanInfo(s).pepolar,[1,3]) &&...
     ismember(parms.revflag,[1,2]))) &&...
    (~flex_flag || (flex_flag && parms.flex_flag))
    SessInfo.snums_DT = [SessInfo.snums_DT,s];
  end;

  if ScanInfo(s).ndiffdirs>=6
    % load qmat
    if ~isfield(ScanInfo,'qmat') || isempty(ScanInfo(s).qmat)
      if isfield(ScanInfo(s), 'diffdirs')
        diffdirs = ScanInfo(s).diffdirs;
      else
        diffdirs = [];
      end;
      if isfield(ScanInfo,'tensor_fnum')
        tensor_fnum = ScanInfo(s).tensor_fnum;
      else
        tensor_fnum = [];
      end;
      % load matrix of diffusion vectors (ndirs x 3)
      ScanInfo(s).qmat = dti_load_qmat(ScanInfo(s).DTI_Sequence_Type,...
        ScanInfo(s).nb0,ScanInfo(s).ndiffdirs,...
        'diffdirs',diffdirs,...
        'tensor_fnum',tensor_fnum);
    end;
  else
    ScanInfo(s).qmat = [];
  end;
  ScanInfo(s).B0uw_ref_for = [];
  ScanInfo(s).B0uw_ref_rev = [];
end;
SessInfo.nscans_for = length(SessInfo.snums_for);
SessInfo.nscans_rev = length(SessInfo.snums_rev);

% choose pairs of forward and reverse scans
%  for B0 unwarp and between scan registration
SessInfo.B0uw_refs_for = [];
SessInfo.B0uw_refs_rev = [];
SessInfo.reg_ref_for = [];
SessInfo.reg_ref_rev = [];
SessInfo.revflag = 0;
if SessInfo.nscans_for>0 & SessInfo.nscans_rev>0
  if parms.revflag==1
    SessInfo.revflag = 1;
  elseif parms.revflag==0
    SessInfo.revflag = 0;
  elseif DTI_ipp_flag==2
    SessInfo.revflag = 0;
  elseif DTI_ipp_flag==3
    SessInfo.revflag = 1;
  elseif (SessInfo.nscans_rev > SessInfo.nscans_for)
    SessInfo.revflag = 1;
  elseif (SessInfo.nscans_for == SessInfo.nscans_rev) &&...
         (SessInfo.snums_for(1)<SessInfo.snums_rev(1))
    SessInfo.revflag = 1;
  else
    SessInfo.revflag = 0;
  end;
  if SessInfo.revflag
    snumsA = SessInfo.snums_for;
    snumsB = SessInfo.snums_rev;
  else
    snumsA = SessInfo.snums_rev;
    snumsB = SessInfo.snums_for;
  end;
  for i=1:length(snumsA)
    sA = snumsA(i);
    tmp_diff = abs(sA - snumsB);
    [tmp,ind] = min(tmp_diff);
    % if tie, choose most similar b value
    % if all same, choose later scan
    ind_close = find(tmp_diff==tmp);
    if length(ind_close)>1
      close_bvals = cell2mat({ScanInfo(snumsB(ind_close)).bval});
      curr_bval = ScanInfo(sA).bval;
      if length(unique(close_bvals))==1
        ind = ind_close(2);
      else
        tmp_diff_bval = abs(close_bvals-curr_bval);
        [tmp,ind_bval] = min(tmp_diff_bval);
        ind_tmp = find(tmp_diff_bval==tmp);
        if length(ind_tmp)>1
          ind_bval = ind_tmp(2);
        end;
        ind = ind_close(ind_bval);
      end;
    end;
    sB = snumsB(ind);    
    if SessInfo.revflag
      SessInfo.B0uw_refs_for = [SessInfo.B0uw_refs_for,sA];
      SessInfo.B0uw_refs_rev = [SessInfo.B0uw_refs_rev,sB];
    else
      SessInfo.B0uw_refs_rev = [SessInfo.B0uw_refs_rev,sA];
      SessInfo.B0uw_refs_for = [SessInfo.B0uw_refs_for,sB];
    end;
  end;
  
  % set for and rev references for between scan registration
  nrefs = length(SessInfo.B0uw_refs_for);
  ind = floor(median(1:nrefs));
  SessInfo.reg_ref_for = SessInfo.B0uw_refs_for(ind);
  SessInfo.reg_ref_rev = SessInfo.B0uw_refs_rev(ind);

  % choose B0uw_ref for each scan 
  for s=1:SessInfo.nscans
    if ismember(ScanInfo(s).pepolar,[1,3])
      B0uw_refs = SessInfo.B0uw_refs_rev;
    else
      B0uw_refs = SessInfo.B0uw_refs_for;
    end;
    tmp = B0uw_refs - s;
    if any(tmp<0) % use most recent reference
      tmp(tmp>0) = Inf; % ignore future scans
      [minval,ind]=min(abs(tmp));
    else % ref scan is here or ahead, find closest
      [minval,ind]=min(tmp);
    end;
    ScanInfo(s).B0uw_ref_for = SessInfo.B0uw_refs_for(ind);
    ScanInfo(s).B0uw_ref_rev = SessInfo.B0uw_refs_rev(ind);
  end;
elseif SessInfo.nscans_for>0
  SessInfo.revflag = 0;
  nrefs = length(SessInfo.snums_for);
  ind = floor(median(1:nrefs));
  SessInfo.reg_ref_for = SessInfo.snums_for(ind);
elseif SessInfo.nscans_rev>0
  SessInfo.revflag = 1;
  nrefs = length(SessInfo.snums_rev);
  ind = floor(median(1:nrefs));
  SessInfo.reg_ref_rev = SessInfo.snums_rev(ind);
end;

if ~SessInfo.revflag
  SessInfo.regT1_ref = SessInfo.reg_ref_for;
else
  SessInfo.regT1_ref = SessInfo.reg_ref_rev;
end;

