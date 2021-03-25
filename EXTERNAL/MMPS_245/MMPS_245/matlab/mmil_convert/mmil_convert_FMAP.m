function scaninfo = mmil_convert_FMAP(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_FMAP(serinfo,i_series,i_scan,outdir,forceflag)
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
% Last Mod: 02/19/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scaninfo = [];
if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;
ext = '.mgz';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gather information
mr_parms = mmil_serinfo_mrparms(serinfo);
[nr,nc,ns] = mmil_serinfo_dims(serinfo);
if nr == serinfo.Rows
  nreps = floor(serinfo.nimgs/ns);
  mosaic_flag = 0;
else
  nreps = serinfo.nimgs;
  mosaic_flag = 1;
end;
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);

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
  'nreps',nreps,...
  'valid',false);

% convert dicoms to mgz
if strcmp(serinfo.SeriesType,'fm_TE1_NFS')
  fstem = sprintf('FMAP%d_TE1',i_scan);
elseif strcmp(serinfo.SeriesType,'fm_TE2_NFS')
  fstem = sprintf('FMAP%d_TE2',i_scan);
else
  fstem = sprintf('FMAP%d',i_scan);
end;
fname = sprintf('%s/%s%s',outdir,fstem,ext);
if ~exist(fname,'file') || forceflag
  fprintf('%s: converting %s\n',mfilename,fstem);
  % read dicom 4d volume
  [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,nreps,mosaic_flag);
  if isempty(vol), return; end;
  fs_save_mgh(vol,fname,M,mr_parms);
end;

scaninfo.valid = true;
