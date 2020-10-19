function scaninfo = mmil_convert_DTI_GE(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_DTI_GE(serinfo,i_series,i_scan,outdir,forceflag)
%
% Required Input:
%   serinfo: series info struct created by mmil_classify_dicoms
%     must be a single struct, not an array
%   i_series: series index (within study)
%   i_scan: scan index (used in output name)
%
% Optional Input:
%   outdir: output directory
%     {default = pwd}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   scaninfo: struct containing selected information for scan
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 02/19/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scaninfo = [];
if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set DTI sequence type
[DTI_Sequence_Type,ipp_flag] = set_sequence_type(serinfo);

% gather information
mr_parms = mmil_serinfo_mrparms(serinfo);
[nr,nc,ns] = mmil_serinfo_dims(serinfo);
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);
nreps = serinfo.nb0+serinfo.ndiffdirs;
if nreps==0, nreps = floor(serinfo.nimgs/ns); end;
revflag = strcmp(serinfo.SeriesType,'DTI_rev');
if ~ipp_flag
  if revflag
    pepolar = 1;
  else
    pepolar = 0;
  end;
else
  pepolar = serinfo.pepolar;
end;

% create scaninfo struct
scaninfo = struct(...
  'SeriesIndex',i_series,...
  'readoutFOV',readoutFOV,...
  'phaseFOV',phaseFOV,...
  'PhaseDir',PhaseDir,...
  'NumberOfPhaseEncodingSteps',serinfo.NumberOfPhaseEncodingSteps,...
  'SliceThickness',serinfo.SliceThickness,...
  'TE',sort(unique(serinfo.EchoTimes)),...
  'TR',serinfo.RepetitionTime,...
  'FlipAngle',serinfo.FlipAngle,...
  'ScanType',serinfo.SeriesType,...
  'gradwarpinfo',serinfo.gradwarpinfo,...
  'ReceiveCoilName',serinfo.ReceiveCoilName,...
  'DTI_Sequence_Type',DTI_Sequence_Type,...
  'diffdirs',serinfo.diffdirs,...
  'ndiffdirs',serinfo.ndiffdirs,...
  'nb0',serinfo.nb0,...
  'bval',serinfo.bval,...
  'nreps',nreps,...
  'pepolar',pepolar,...
  'tensor_fnum',serinfo.tensor_fnum,...
  'EchoSpacing',serinfo.EchoSpacing,...
  'AcquisitionRows',serinfo.AcquisitionRows,...
  'AcquisitionColumns',serinfo.AcquisitionColumns,...
  'valid',false);

% convert dicoms to mgz
if ipp_flag
  % combining the two files in nreps to split using mmil_read_dicom_4dvol
  nreps = nreps + 1;
  fname = sprintf('%s/DTI%d.mgz',outdir,i_scan);
  fname_rev = sprintf('%s/DTI%d_rev.mgz',outdir,i_scan);
  if ~exist(fname,'file') || ~exist(fname_rev,'file') || forceflag
    fprintf('%s: converting DTI%d (pepolar = %d)...\n',...
      mfilename,i_scan,pepolar);
    % read dicom 4d volume
    [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,nreps);
    if isempty(vol), return; end;
    % split into forward and reverse
    switch pepolar
      case 2
        Vr = vol(:,:,:,1);     % first frame is "reverse"
        V  = vol(:,:,:,2:end); % subsequent frames are "forward"
      case 3
        V  = vol(:,:,:,1);     % first frame is "forward"
        Vr = vol(:,:,:,2:end); % subsequent frames are "reverse"
      otherwise
        error('invalid pepolar (%d) for ipp',pepolar);
    end;
    Vr = Vr(:,end:-1:1,:,:); % flip "reverse" images in y dimension
    if isempty(V) | isempty(Vr), return; end;
    fs_save_mgh(V,fname,M,mr_parms);
    fs_save_mgh(Vr,fname_rev,M,mr_parms);
  end;
else
  if revflag
    fstem = sprintf('DTI%d_rev',i_scan);
  else
    fstem = sprintf('DTI%d',i_scan);
  end;
  fname = sprintf('%s/%s.mgz',outdir,fstem);
  if ~exist(fname,'file') || forceflag
    fprintf('%s: converting %s...\n',mfilename,fstem);
    % read dicom 4d volume
    [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,nreps);
    if isempty(vol), return; end;
    if isfield(serinfo,'SequenceLabel') && strcmp(upper(serinfo.SequenceLabel),'MMILDTI') && revflag
      % flip rev images in y
      vol = vol(:,end:-1:1,:,:);
    end;
    fs_save_mgh(vol,fname,M,mr_parms);
  end;
end;

scaninfo.valid = true;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DTI_Sequence_Type,ipp_flag] = set_sequence_type(serinfo)
  DTI_Sequence_Type = []; ipp_flag = 0; flex_flag = 0;
  if strcmp(serinfo.SeriesType,'DTI_flex')
    flex_flag = 1;
    if serinfo.pepolar>1, ipp_flag = 1; end;
  elseif strcmp(serinfo.SeriesType,'DTI_ipp')
    ipp_flag = 1;
  end;
  if flex_flag
    DTI_Sequence_Type = 6; % variable tensorX.dat file
  elseif ipp_flag || ~isempty(regexp(upper(serinfo.SequenceLabel),'MMILDTI'))
    DTI_Sequence_Type = 5; % different tensor.dat file for 14x
  else
    DTI_Sequence_Type = 4; % original tensor.dat file for 12x
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

