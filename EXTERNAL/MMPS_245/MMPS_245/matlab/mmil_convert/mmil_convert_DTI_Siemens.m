function scaninfo = mmil_convert_DTI_Siemens(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_DTI_Siemens(serinfo,i_series,i_scan,outdir,forceflag)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set DTI sequence type
[DTI_Sequence_Type,rrr_flag] = set_sequence_type(serinfo);

% gather information
mr_parms = mmil_serinfo_mrparms(serinfo);
[nr,nc,ns] = mmil_serinfo_dims(serinfo);
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);
revflag = strcmp(serinfo.SeriesType,'DTI_rev');
ndiffdirs = serinfo.ndiffdirs;
nb0 = serinfo.nb0;
scalefact = mmil_getfield(serinfo,'scalefact',1000);
splitflag = 0;
rr_flag = 0;
ismosaic = ~isempty(regexp(mmil_getfield(serinfo,'ImageType',''),'MOSAIC'));
if ismosaic
  mosaic_flag = 1;
  nreps = serinfo.nimgs;
  if ~isempty(regexp(serinfo.SequenceName,'epse2d1_96'))
    ndiffdirs = 0;
    nvols = nreps/nb0;
    if nvols == 2
      % if 2 images, then it is a mosaic field map,
      %  with each image containing a different phase encoding direction
      rr_flag = 1;
    end;
  end;
  nvols = nreps/(nb0+ndiffdirs);
  nreps = nb0+ndiffdirs;
  if nvols > 1 && ~rrr_flag && ~rr_flag, splitflag = 1; end;  
else
  mosaic_flag = 0;
  nimgs = serinfo.nimgs;
  nreps = floor(serinfo.nimgs/ns); 
%  if ndiffdirs~=1
%    fprintf('%s: WARNING: setting ndiffdirs = 1, actual value is %d\n',...
%      mfilename,ndiffdirs);
%    keyboard
%  end;
%  ndiffdirs = 1;
  nvols = 1;
end;        
if nreps==0, nreps = floor(serinfo.nimgs/ns); end;
if rr_flag
  pepolar = 3; % forward then reverse
elseif rrr_flag
  pepolar = serinfo.pepolar;
elseif splitflag
  if revflag
    pepolar = 1;
  else
    pepolar = 0;
  end;
else
  if revflag
    pepolar = 1;
  else
    pepolar = 0;
  end;
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
  'ndiffdirs',ndiffdirs,...
  'nb0',serinfo.nb0,...
  'bval',serinfo.bval,...
  'nreps',nreps,...
  'pepolar',pepolar,...
  'tensor_fnum',serinfo.tensor_fnum,...
  'EchoSpacing',serinfo.EchoSpacing,...
  'AcquisitionRows',serinfo.AcquisitionRows,...
  'AcquisitionColumns',serinfo.AcquisitionColumns,...
  'valid',false);

if splitflag
  scaninfo(2) = scaninfo;
end;

% convert dicoms to mgz
if rr_flag
  fname = sprintf('%s/DTI%d.mgz',outdir,i_scan);
  fname_rev = sprintf('%s/DTI%d_rev.mgz',outdir,i_scan);
  if ~exist(fname,'file') || ~exist(fname_rev,'file') || forceflag
    fprintf('%s: converting DTI%d (rr)...\n',mfilename,i_scan);
    % read dicom 4d volume
    [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,2*nreps,...
      mosaic_flag);
    if isempty(vol), return; end;
    % split into forward and reverse (specific to revro sequence)
    nreps = nb0+ndiffdirs;
    V1 = vol(:,:,:,1);
    V2 = vol(:,:,:,2);
    if isempty(V1) | isempty(V2), return; end;
    fs_save_mgh(V1,fname,M,mr_parms);
    fs_save_mgh(V2,fname_rev,M,mr_parms);
  end;
elseif rrr_flag
  fname = sprintf('%s/DTI%d.mgz',outdir,i_scan);
  fname_rev = sprintf('%s/DTI%d_rev.mgz',outdir,i_scan);
  if ~exist(fname,'file') || ~exist(fname_rev,'file') || forceflag
    fprintf('%s: converting DTI%d (rrr)...\n',mfilename,i_scan);
    % read dicom 4d volume
    [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,2*nreps,...
      mosaic_flag);
    if isempty(vol), return; end;
    % split into forward and reverse (specific to rorevro sequence)
    nreps = nb0+ndiffdirs;
    V1 = vol(:,:,:,1:nreps);
    V2 = vol(:,:,:,nreps+1:end);
    if isempty(V1) | isempty(V2), return; end;
    fs_save_mgh(V1,fname,M,mr_parms);
    fs_save_mgh(V2,fname_rev,M,mr_parms);
  end;
elseif splitflag
  if revflag
    fstem1 = sprintf('DTI%d_rev',i_scan);
    fstem2 = sprintf('DTI%d_rev',i_scan+1);
  else
    fstem1 = sprintf('DTI%d',i_scan);
    fstem2 = sprintf('DTI%d',i_scan+1);
  end;
  fname1 = sprintf('%s/%s.mgz',outdir,fstem1);
  fname2 = sprintf('%s/%s.mgz',outdir,fstem2);
  if ~exist(fname1,'file') || ~exist(fname2,'file') || forceflag
    fprintf('%s: converting %s and %s...\n',mfilename,fstem1,fstem2);
    [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,...
      nvols*nreps,mosaic_flag);         
    if isempty(vol), return; end;
    V1 = vol(:,:,:,1:nreps);
    V2 = vol(:,:,:,nreps+1:end);
    fs_save_mgh(V1,fname1,M,mr_parms);                      
    fs_save_mgh(V2,fname2,M,mr_parms);                                     
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
    [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,...
      nvols*nreps,mosaic_flag);         
    if isempty(vol), return; end;
    fs_save_mgh(vol,fname,M,mr_parms);
  end;
end;

if splitflag
  scaninfo(1).valid = true;
  scaninfo(2).valid = true;
else
  scaninfo.valid = true;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DTI_Sequence_Type,rrr_flag] = set_sequence_type(serinfo)
  DTI_Sequence_Type = []; rrr_flag = [];
  if ~isempty(regexp(serinfo.SeriesDescription,'rorevro'))
    if ~isempty(regexp(serinfo.ManufacturersModelName,'Avanto')) |...
       ~isempty(regexp(serinfo.ManufacturersModelName,'Allegra'))
      DTI_Sequence_Type = 1; % MGH Avanto(1.5T)/Allegra(3T) reversing PE polarity sequence
    elseif ~isempty(regexp(serinfo.ManufacturersModelName,'Symphony'))
      DTI_Sequence_Type = 2; % UCSD Symphony reversing PE polarity sequence
    else
      DTI_Sequence_Type = 3; % some other Siemens scanner, rorevro sequence
    end;
    rrr_flag = 1;
  else
    if ~isempty(regexp(serinfo.ManufacturersModelName,'Allegra'))
      DTI_Sequence_Type = 11; % NYU Allegra (3T)
    else
      DTI_Sequence_Type = 13; % some other Siemens sequence / scanner
    end;
    rrr_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

