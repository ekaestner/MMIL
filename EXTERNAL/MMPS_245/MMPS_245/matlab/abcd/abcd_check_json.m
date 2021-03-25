function [jinfo,fstem_json,fname_tar,fname_pfile_json] = abcd_check_json(fname_json,indir)
%function [jinfo,fstem_json,fname_tar,fname_pfile_json] = abcd_check_json(fname_json,indir)
%
% Required Input:
%   fname_json: full path of input json file
%
% Optional Input:
%   indir: input directory containing tgz files
%     if not supplied, will be assumed to be directory containing fname_json
%
% Output:
%    fname_tar: full path of corresponding tar/tgz file
%    jinfo: struct containing info from json file
%    fname_pfile_json: full path of json for tgz file containing P file
%      empty if no tgz with P file
%
% Created:  07/18/16 by Don Hagler
% Prev Mod: 06/13/17 by Don Hagler
% Last Mod: 09/08/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_tar = []; jinfo = []; fname_pfile_json = [];
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('indir','var'), indir = []; end;
recon_sequence_list = {'research_ABCD_muxepi','research_ABCD_muxepi2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check json file exists
if ~exist(fname_json,'file')
  error('file %s not found',fname_json);
end;
[fdir,fstem_json,fext] = fileparts(fname_json);
if ~strcmp(fext,'.json')
  error('file %s is not json',fname_json);
end;
if isempty(indir), indir = fdir; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verify that json file is ASCII
%cmd = sprintf('file %s',fname_json);
%[s,r] = unix(cmd);
%if s, error('cmd %s failed:\n%s',mfilename,cmd,r); end;
%if isempty(regexp(r,'ASCII'))
%  fprintf('%s: WARNING: wrong format for json file %s\n',mfilename,fname_json);
%  return;
%end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read json file
try
  jinfo = loadjson(fname_json);
catch me
  fprintf('%s: WARNING: error reading json file %s\n',mfilename,fname_json);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize additional output fields
jinfo.invalid_flag = 0;
jinfo.pfile_flag = 0;
jinfo.multiband_flag = 0;
jinfo.recon_required_flag = 0;
jinfo.missing_kspace_flag = 0;
jinfo.pinfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% quit if this is a P-file json file
if isfield(jinfo,'ReconRequired') &&...
    (~isfield(jinfo,'StudyInstanceUID') || ~isfield(jinfo,'Manufacturer'))
  jinfo.pfile_flag = 1;
  jinfo.invalid_flag = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for required fields
if ~isfield(jinfo,'SeriesNumber') ||...
   ~isfield(jinfo,'SeriesDescription')% ||...
%   ~isfield(jinfo,'SeriesTime')
  fprintf('%s: WARNING: json file %s has missing fields\n',mfilename,fname_json);
  jinfo.invalid_flag = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine if is multi-band acquisition based on SeriesDescription
SeriesDescription = jinfo.SeriesDescription;
if ~isempty(regexp(SeriesDescription,'DTI')) ||...
    ~isempty(regexp(SeriesDescription,'dMRI')) ||...
    ~isempty(regexp(SeriesDescription,'Diffusion_Field_Map')) ||...
    (~isempty(regexp(SeriesDescription,'fMRI')) &&...
      isempty(regexp(SeriesDescription,'fMRI_Field_Map')))
  jinfo.multiband_flag = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% decide if recon_required based on Manufacturer and SequenceName
if isfield(jinfo,'SequenceName')
  SequenceName = jinfo.SequenceName;
elseif isfield(jinfo,'PulseSequenceName')
  SequenceName = jinfo.PulseSequenceName;
else
  SequenceName = 'unknown';
end;
SequenceName = regexprep(SequenceName,'/','_');
if ~isempty(regexp(jinfo.Manufacturer,'GE')) &&...
   (ismember(SequenceName,recon_sequence_list) || jinfo.multiband_flag)
  % if OS level is dv26, no recon required
  %% todo: use new tag in json
  if ~isempty(regexp(jinfo.PatientID,'dv26'))
    jinfo.recon_required_flag = 0;
  else
    jinfo.recon_required_flag = 1;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% identify tar file
if jinfo.recon_required_flag
  % search for P-file tgz
  fpat = sprintf('%s/*SUID_%s_*se%s_*.t*',indir,jinfo.StudyInstanceUID,jinfo.SeriesNumber);
  jlist = dir(fpat);
  % search for tar file with alternate naming convention
  if isempty(jlist)
    fpat = sprintf('%s/subjid%s_*se%s_*.t*',indir,jinfo.PatientID,jinfo.SeriesNumber);
    jlist = dir(fpat);
  end;
  if ~isempty(jlist)
    % choose from possibly multiple files
    nfiles = length(jlist);
    for j=1:nfiles
      fname_tar_tmp = sprintf('%s/%s',indir,jlist(j).name);
      k = regexp(fname_tar_tmp,'subjid(?<SubjID>\w+)_ex\d+_se(?<SeriesNummber>\d+)_(?<StudyDate>\d{8})_(?<SeriesTime>\d+)','names');
      if ~isempty(k) || ~isempty(regexp(fname_tar_tmp,'SUID_'))
        if ~isempty(k) && isfield(jinfo,'SeriesTime')
          if ~strcmp(k.StudyDate,jinfo.StudyDate), continue; end;
          if length(k.SeriesTime>6), k.SeriesTime = k.SeriesTime(1:6); end;
          if str2num(k.SeriesTime) < str2num(jinfo.StudyTime), continue; end;
          if str2num(k.SeriesTime) - str2num(jinfo.StudyTime) > 20000, continue; end;
        end;
        fname_tar = fname_tar_tmp;
        break;
      else
        fprintf('%s: WARNING: tar file %s did NOT match the expected pattern\n',mfilename,fname_tar_tmp);
      end;
    end;
  end;
  % look for matching P-file json file
  if ~isempty(fname_tar)
    [tmp,fstem_pfile_tar,fext] = fileparts(fname_tar);
    % remove extra .tar from .tar.gz (if present)
    if strcmp(fext,'.gz') && ~isempty(regexp(fstem_pfile_tar,'.tar$'))
      [tmp,fstem_pfile_tar,fext] = fileparts(fstem_pfile_tar);
    end;
    % check that expected P-file json exists
    fname_pfile_json = sprintf('%s/%s.json',indir,fstem_pfile_tar);
    if ~exist(fname_pfile_json)
      fname_pfile_json = [];
    else
      % set recon_required_flag to 0 if indicated in json file
      try
        jinfo.pinfo = loadjson(fname_pfile_json);
      catch me
        fprintf('%s: WARNING: error reading P-file json file %s\n',mfilename,fname_pfile_json);
      end;
      if isfield(jinfo.pinfo,'ReconRequired') && str2num(jinfo.pinfo.ReconRequired)==0
        jinfo.recon_required_flag = 0;
      end;
    end;
  else
    fprintf('%s: WARNING: P-file tgz file not found for %s\n',mfilename,fname_json);
    jinfo.missing_kspace_flag = 1;
  end;
end;
if ~jinfo.recon_required_flag
  fname_tar = sprintf('%s/%s.tgz',indir,fstem_json);
  if ~exist(fname_tar,'file')
    fprintf('%s: WARNING: tgz file not found for %s\n',mfilename,fname_json);
    fname_tar = []; 
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

