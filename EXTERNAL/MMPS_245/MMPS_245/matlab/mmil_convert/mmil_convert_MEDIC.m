function scaninfo = mmil_convert_MEDIC(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_MEDIC(serinfo,i_series,i_scan,outdir,forceflag)
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
PUREflag = ... % only on GE
  ~isempty(regexp(mmil_getfield(serinfo,'ScanOptions',''),'FILTERED_GEMS'));
EchoTimes = serinfo.EchoTimes;
TE = sort(unique(EchoTimes));

% create scaninfo struct
scaninfo = struct(...
  'SeriesIndex',i_series,...
  'readoutFOV',readoutFOV,...
  'phaseFOV',phaseFOV,...
  'PhaseDir',PhaseDir,...
  'SliceThickness',serinfo.SliceThickness,...
  'TE',TE,...
  'TR',serinfo.RepetitionTime,...
  'TI',serinfo.InversionTime,...
  'FlipAngle',serinfo.FlipAngle,...
  'ScanType',ScanType,...
  'gradwarpinfo',serinfo.gradwarpinfo,...
  'ReceiveCoilName',serinfo.ReceiveCoilName,...
  'nslices',serinfo.nimgs,...
  'PUREflag',PUREflag,...
  'valid',false);

% convert dicoms to mgh
fstem = sprintf('%s%d',ScanType,i_scan);
fname = sprintf('%s/%s%s',outdir,fstem,ext);
fname_R2 = sprintf('%s/%s_R2%s',outdir,fstem,ext);
if ~exist(fname,'file') || ~exist(fname_R2,'file') || forceflag
  fprintf('%s: converting %s...\n',mfilename,fstem);
  wM = [-TE'/1000 ones(length(TE),1)];
  wMi = pinv(wM);
  w1 = wMi(1,:);
  w2 = wMi(2,:);
  volokflag = true;
  vol_sum = single(0);
  vol_sum_R2 = single(0);
  for t=1:length(TE)
    imlist = find(EchoTimes==TE(t));
    if length(imlist)>1
      [vol,M] = mmil_read_dicom_vol(serinfo.FileNames(imlist));
      if isempty(vol), return; end;
      try
        vol_sum = vol_sum+single(vol.^2);
      catch
        volokflag = false;
      end;
      try
        vol_sum_R2 = vol_sum_R2+single(log(vol+eps)*w1(t));
      catch
        volokflag = false;
      end;
      fname_tmp = sprintf('%s/%s_%d%s',outdir,fstem,t,ext);
      fs_save_mgh(vol,fname_tmp,M,mr_parms);
    end
  end
  if volokflag
    vol = double(sqrt(vol_sum/length(TE)));
    if isempty(vol), return; end;
    fs_save_mgh(vol,fname,M,mr_parms);
    vol_sum_R2(find(~isfinite(vol_sum_R2))) = 0;
    vol = double(vol_sum_R2);
    if isempty(vol), return; end;
    fs_save_mgh(vol,fname_R2,M,mr_parms);
  end
end

scaninfo.valid = true;
