function scaninfo = mmil_convert_DTI_Philips(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_DTI_Philips(serinfo,i_series,i_scan,outdir,forceflag)
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

tmpdir = [outdir '/tmp_convert_DTI_Philips'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set DTI sequence type
DTI_Sequence_Type = 14;

% gather information
revflag = strcmp(serinfo.SeriesType,'DTI_rev');
mr_parms = mmil_serinfo_mrparms(serinfo);
[nr,nc,ns] = mmil_serinfo_dims(serinfo);
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);
nreps = serinfo.nb0+serinfo.ndiffdirs;
if nreps==0, nreps = floor(serinfo.nimgs/ns); end;

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
  'pepolar',serinfo.pepolar,...
  'tensor_fnum',serinfo.tensor_fnum,...
  'EchoSpacing',serinfo.EchoSpacing,...
  'AcquisitionRows',serinfo.AcquisitionRows,...
  'AcquisitionColumns',serinfo.AcquisitionColumns,...
  'valid',false);

% reshape fnames to make slices change as inner loop
fnames = serinfo.FileNames;
if nreps>1
  nfiles = length(fnames);
  % check if there is an extra rep that needs to be discarded
  if nfiles == ns*(nreps+1)
    nreps = nreps + 1;
    skip_last_flag = 1;
  else
    skip_last_flag = 0;
  end;
  if nfiles ~= ns*nreps
    fprintf('%s: WARNING: nfiles (%d) does not match nslices*nreps (%dx%d)\n',...
      mfilename,nfiles,ns,nreps);
    return;
  end;
  fnames = reshape(reshape(fnames,[nreps,ns])',[nfiles,1]);
else
  skip_last_flag = 0;
end;

% convert dicoms to mgz
fail_flag = 0;
if revflag
  fstem = sprintf('DTI%d_rev',i_scan);
else
  fstem = sprintf('DTI%d',i_scan);
end;
fname = sprintf('%s/%s.mgz',outdir,fstem);
if ~exist(fname,'file') || forceflag
  mmil_mkdir(tmpdir);
  fprintf('%s: converting %s...\n',mfilename,fstem);
  % read dicom 4d volume
  [vol,M] = mmil_read_dicom_4dvol(fnames,nr,nc,ns,nreps);
  if skip_last_flag>1
    % remove extra synthesized image added to Philips diffusion scans
    vol = vol(:,:,:,1:end-1);
  end;
  if isempty(vol), return; end;
  fname_tmp = sprintf('%s/%s.mgz',tmpdir,fstem);
  fs_save_mgh(vol,fname_tmp,M,mr_parms);
  if revflag
    cmd = sprintf('shiftX -i1 %s -o %s -osx -1',fname_tmp,fname_tmp);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: WARNING: failed to x-shift Philips DTI:\n%s\n',...
        mfilename,result);
      fail_flag = 1;
    end;
  end
  if nr == 320 && nc == 320 && ~fail_flag
    mmil_mkdir(tmpdir);
    paramfile = [getenv('MMPS_PARMS') '/MRIREG/inputParamsScale.txt'];
    cmd = sprintf('resize_volume -ip %s -s %s -t %s -od %s',...
      paramfile,fname_tmp,fname_tmp,tmpdir);
    [status,result] = unix(cmd);
    if status
      fprintf('%s: WARNING: failed to resize DTI volume:\n%s\n',...
        mfilename,result);
      fail_flag = 1;
    end;
    fname_tmp = sprintf('%s/%s_SRescaled.mgz',tmpdir,fstem);
  end;                        
  if ~fail_flag
    cmd = sprintf('mv %s %s',fname_tmp,fname);
    [status,result] = unix(cmd);
    if status, error('cmd %s failed:\n%s',cmd,result); end;
  end;
end;

scaninfo.valid = true;

% remove tmpdir
if exist(tmpdir,'dir')
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result] = unix(cmd);
  if status, error('failed to remove tmpdir %s:\n%s',tmpdir,result); end;
end;

