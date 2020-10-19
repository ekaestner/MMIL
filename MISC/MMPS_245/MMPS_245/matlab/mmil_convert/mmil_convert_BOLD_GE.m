function scaninfo = mmil_convert_BOLD_GE(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_BOLD_GE(serinfo,i_series,i_scan,outdir,forceflag)
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
% Last Mod: 06/01/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scaninfo = [];
if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gather information
ScanType = serinfo.SeriesType;
mr_parms = mmil_serinfo_mrparms(serinfo);
[nr,nc,ns] = mmil_serinfo_dims(serinfo);
[readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo);
if strcmp(ScanType,'BOLD_bu')
  pepolar = 0;
elseif strcmp(ScanType,'BOLD_td');
  pepolar = 1;
else
  pepolar = serinfo.pepolar;
end;
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

% convert dicoms to mgz
fstem = sprintf('BOLD%d',i_scan);
switch ScanType
  case {'BOLD_SE','BOLD_bu','BOLD_td'}
    if pepolar
      fstem = sprintf('BOLD%d_rev',i_scan);
    else
      fstem = sprintf('BOLD%d_for',i_scan);
    end;
    fname = sprintf('%s/%s.mgz',outdir,fstem);
    if ~exist(fname,'file') || forceflag
      fprintf('%s: converting %s (pepolar = %d)...\n',...
        mfilename,fstem,pepolar);
      [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,nreps); 
      if isempty(vol), return; end;
      fs_save_mgh(vol,fname,M,mr_parms);                             
    end;
  case {'BOLD_pep','BOLD_ipp','BOLD_ape'}
    fname_for = sprintf('%s/%s_for.mgz',outdir,fstem);
    fname_rev = sprintf('%s/%s_rev.mgz',outdir,fstem);
    % pepolar:
    %   0 = all frames are forward (bottom-up)
    %   1 = all frames are reverse (top-down)
    %   2 = first frame is rev, rest are for
    %   3 = first frame is for, rest are rev
    % unless BOLD_ape, in which case:
    %   2 = odd frames are rev, even frames are for
    %   3 = odd frames are for, even frames are rev
    switch pepolar
      case 0
        fname_rev = [];
      case 1
        fname_for = [];
    end;
    if (~isempty(fname_for) && ~exist(fname_for,'file')) ||...
       (~isempty(fname_rev) && ~exist(fname_rev,'file')) ||...
       forceflag
      fprintf('%s: converting %s (pepolar = %d)...\n',...
        mfilename,fstem,pepolar);
      % read dicom 4d volume
      [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,nreps);
      if isempty(vol), return; end;
      switch pepolar
        case 0
          fs_save_mgh(vol,fname_for,M,mr_parms);
        case 1
          switch PhaseDir
            case 'ROW'
              vol = vol(end:-1:1,:,:,:); % flip "reverse" images in x dimension
            case 'COL'
              vol = vol(:,end:-1:1,:,:); % flip "reverse" images in y dimension
            otherwise
              fprintf('%s: WARNING: flipping images and PhaseDir = %s\n',...
                mfilename,PhaseDir);
          end;
          fs_save_mgh(vol,fname_rev,M,mr_parms);
        case {2,3}
          nframes = size(vol,4);
          if strcmp(ScanType,'BOLD_ape')
            switch pepolar
              case 2
                ind_rev = [1:2:nframes]; % odd frames are "reverse"
                ind_for = [2:2:nframes]; % even frames are "forward"
              case 3
                ind_for = [1:2:nframes]; % odd frames are "forward"
                ind_rev = [2:2:nframes]; % even frames are "reverse"
            end;
          else
            switch pepolar
              case 2
                ind_rev = 1;           % first frame is "reverse"
                ind_for = [2:nframes]; % subsequent frames are "forward"
              case 3
                ind_for = 1;           % first frame is "forward"
                ind_rev = [2:nframes]; % subsequent frames are "reverse"
            end;
          end;
          Vfor = vol(:,:,:,ind_for);
          Vrev = vol(:,:,:,ind_rev);
          switch PhaseDir
            case 'ROW'
              Vrev = Vrev(end:-1:1,:,:,:); % flip "reverse" images in x dimension
            case 'COL'
              Vrev = Vrev(:,end:-1:1,:,:); % flip "reverse" images in y dimension
            otherwise
              fprintf('%s: WARNING: flipping images and PhaseDir = %s\n',...
                mfilename,PhaseDir);
          end;
          if isempty(Vrev) | isempty(Vfor), return; end;
          fs_save_mgh(Vfor,fname_for,M,mr_parms);
          fs_save_mgh(Vrev,fname_rev,M,mr_parms);
      end;
    end;
  otherwise
    fname = sprintf('%s/%s.mgz',outdir,fstem);
    if ~exist(fname,'file') | forceflag
      fprintf('%s: converting %s (pepolar = %d)...\n',...
        mfilename,fstem,pepolar);
      % read dicom 4d volume
      [vol,M] = mmil_read_dicom_4dvol(serinfo.FileNames,nr,nc,ns,nreps);
      if isempty(vol), return; end;
      fs_save_mgh(vol,fname,M,mr_parms);
    end;
end;

scaninfo.valid = true;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

