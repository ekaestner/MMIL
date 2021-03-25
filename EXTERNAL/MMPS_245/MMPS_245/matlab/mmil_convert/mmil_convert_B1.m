function scaninfo = mmil_convert_B1(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_B1(serinfo,i_series,i_scan,outdir,forceflag)
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
% Created:  09/12/12 by Don Hagler
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
ScanType = serinfo.SeriesType;
mr_parms = mmil_serinfo_mrparms(serinfo);
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);

% create scaninfo struct
scaninfo = struct(...
  'SeriesIndex',i_series,...
  'readoutFOV',readoutFOV,...
  'phaseFOV',phaseFOV,...
  'PhaseDir',PhaseDir,...
  'SliceThickness',serinfo.SliceThickness,...
  'TE',sort(unique(serinfo.EchoTimes)),...
  'TR',serinfo.RepetitionTime,...
  'FlipAngle',serinfo.FlipAngle,...
  'ScanType',ScanType,...
  'gradwarpinfo',serinfo.gradwarpinfo,...
  'ReceiveCoilName',serinfo.ReceiveCoilName,...
  'valid',false);

% convert dicoms to mgh
fstem = sprintf('%s%d',ScanType,i_scan);
fname = sprintf('%s/%s%s',outdir,fstem,ext);
if ~exist(fname,'file') || forceflag
  fprintf('%s: converting %s\n',mfilename,fstem);
  [vol,M] = mmil_read_dicom_vol(serinfo.FileNames); 
  if isempty(vol), return; end;
  if ~isempty(strfind(ScanType,'_phi'))
    vol = vol*2*pi/(2^12);
  end
  fs_save_mgh(vol,fname,M,mr_parms);
end;

scaninfo.valid = true;
