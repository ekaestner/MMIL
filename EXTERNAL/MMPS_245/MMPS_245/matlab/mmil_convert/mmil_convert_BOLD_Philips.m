function scaninfo = mmil_convert_BOLD_Philips(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_BOLD_Philips(serinfo,i_series,i_scan,outdir,forceflag)
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

tmpdir = [outdir '/tmp_convert_DTI_Philips'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gather information
ScanType = serinfo.SeriesType;
mr_parms = mmil_serinfo_mrparms(serinfo);
[nr,nc,ns] = mmil_serinfo_dims(serinfo);
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);
pepolar = serinfo.pepolar;
nreps = floor(serinfo.nimgs/ns);

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
  'ScanType',ScanType,...
  'gradwarpinfo',serinfo.gradwarpinfo,...
  'ReceiveCoilName',serinfo.ReceiveCoilName,...
  'nreps',nreps,...
  'pepolar',pepolar,...
  'EchoSpacing',serinfo.EchoSpacing,...
  'AcquisitionRows',serinfo.AcquisitionRows,...
  'AcquisitionColumns',serinfo.AcquisitionColumns,...
  'slice_tpattern',serinfo.slice_tpattern,...
  'slice_order',serinfo.slice_order,...
  'valid',false);

if nreps==0
  fprintf('%s: nreps = 0... skipping\n',mfilename);
  return;
end

% reshape fnames to make slices change as inner loop
fnames = serinfo.FileNames;
if nreps>1
  nfiles = length(fnames);
  if nfiles ~= ns*nreps
    fprintf('%s: WARNING: nfiles (%d) does not match nslices*nreps (%dx%d)\n',...
      mfilename,nfiles,ns,nreps);
    return;
  end;
  fnames = reshape(reshape(fnames,[nreps,ns])',[nfiles,1]);
end;

% convert dicoms to mgz
fail_flag = 0;
switch ScanType
  case 'BOLD_SE'
    switch pepolar
      case 0    
        fstem = sprintf('BOLD%d_for',i_scan);
        fname = sprintf('%s/%s.mgz',outdir,fstem);
        if ~exist(fname,'file') || forceflag  
          fprintf('%s: converting %s...\n',mfilename,fstem);
          [vol,M] = mmil_read_dicom_4dvol(fnames,nr,nc,ns,nreps);
          if isempty(vol), return; end;               
          fs_save_mgh(vol,fname,M,mr_parms);
        end;
      case 1                                  
        fstem = sprintf('BOLD%d_rev',i_scan);
        fname = sprintf('%s/%s.mgz',outdir,fstem);
        if ~exist(fname,'file') || forceflag
          mmil_mkdir(tmpdir);
          fprintf('%s: converting %s...\n',mfilename,fstem);
          [vol,M] = mmil_read_dicom_4dvol(fnames,nr,nc,ns,nreps);
          if isempty(vol), return; end;               
          fname_tmp = sprintf('%s/%s.mgz',tmpdir,fstem);
          fs_save_mgh(vol,fname_tmp,M,mr_parms);
          cmd = sprintf('shiftX -i1 %s -o %s -osx 1',fname_tmp,fname_tmp);
          [status,result] = unix(cmd);
          if status
            fprintf('%s: ERROR: failed to x-shift Philips DTI:\n%s\n',...
              mfilename,result);
            fail_flag = 1;
          end;
          if ~fail_flag
            cmd = sprintf('mv %s %s',fname_tmp,fname);
            [status,result] = unix(cmd);
            if status, error('cmd %s failed:\n%s',cmd,result); end;
          end;
       end;
    end;
  case 'BOLD_bu'
    fstem = sprintf('BOLD%d_for',i_scan);
    fname = sprintf('%s/%s.mgz',outdir,fstem);
    if ~exist(fname,'file') || forceflag
      fprintf('%s: converting %s...\n',mfilename,fstem);
      % read dicom 4d volume
      [vol,M] = mmil_read_dicom_4dvol(fnames,nr,nc,ns,nreps); 
      if isempty(vol), return; end;
      fs_save_mgh(vol,fname,M,mr_parms);
    end;
  case 'BOLD_td'
    fstem = sprintf('BOLD%d_rev',i_scan);
    fname = sprintf('%s/%s.mgz',outdir,fstem);
    if ~exist(fname,'file') || forceflag
      fprintf('%s: converting %s...\n',mfilename,fstem);
      % read dicom 4d volume
      [vol,M] = mmil_read_dicom_4dvol(fnames,nr,nc,ns,nreps); 
      if isempty(vol), return; end;
      fs_save_mgh(vol,fname,M,mr_parms);
    end;
end;

scaninfo.valid = true;

% remove tmpdir
if exist(tmpdir,'dir')
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result] = unix(cmd);
  if status, error('failed to remove tmpdir %s:\n%s',tmpdir,result); end;
end;

