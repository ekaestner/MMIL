function [SeriesInfo,errcode,msg,invalid_files] = mmil_sort_dicoms(indir,varargin)
%function [SeriesInfo,errcode,msg,invalid_files] = mmil_sort_dicoms(indir,[options])
%
% Required Input:
%   indir: full path of input directory containing dicoms
%     recursive search will be used to find dicoms in subdirectories
%
% Optional Parameters:
%   'tags': dicom header info field names to be extracted
%     first entry will be used to sort into series
%     {default = 'SeriesInstanceUID','SeriesNumber',...
%                'StudyInstanceUID','InstanceNumber','EchoTime'}}
%   'regpat': cell array of regular expression for extracting results
%     one for each entry in tags; names in regpat must match tags
%     {default = {'\(0020,000e\) UI \[(?<SeriesInstanceUID>[\.\d]+)\]',...
%                 '\(0020,0011\) IS \[(?<SeriesNumber>\d+)\]',...
%                 '\(0020,000d\) UI \[(?<StudyInstanceUID>[\.\d]+)\]',...
%                 '\(0020,0013\) IS \[(?<InstanceNumber>\d+)\]',...
%                 '\(0018,0081\) DS \[(?<EchoTime>\d+)\]'},[],...
%   'batch_limit': maximum number of lines per unix call
%     {default = 250}
%
% Output:
%   SeriesInfo: struct array containing FileNames and values for each series
%     'SeriesInstanceUID','SeriesNumber', 'StudyInstanceUID',InstanceNumber',
%     and 'errmsg'
%   errcode: returns 1 if error, 0 if successful
%   msg: error message
%   invalid_files: cell array of non-dicom or unreadable files
%
% Created:  03/14/11 by Don Hagler
% Last Mod: 09/10/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'tags',{'SeriesInstanceUID',...
          'SeriesNumber',...
          'StudyInstanceUID',...
          'InstanceNumber',...
          'EchoTime'},[],...
  'regpats',{'\(0020,000e\) UI \[(?<SeriesInstanceUID>[\.\d]+)\]',...
             '\(0020,0011\) IS \[(?<SeriesNumber>\d+)\]',...
             '\(0020,000d\) UI \[(?<StudyInstanceUID>[\.\d]+)\]',...
             '\(0020,0013\) IS \[(?<InstanceNumber>\d+)\]',...
             '\(0018,0081\) DS \[(?<EchoTime>[\.\d]+)\]'},[],...
  'batch_limit',250,[],...
  ...
  'sep','###########',[],...
  'invalstr','invalid',[],...
  'numeric_tags',{'SeriesNumber','InstanceNumber','EchoTime'},[],...
  'series_tags',{'SeriesInstanceUID','SeriesNumber',...
                 'StudyInstanceUID'},[],...
});

SeriesInfo = [];
errcode = 0;
msg = [];
invalid_files = [];

ntags = length(parms.tags);
if ntags~=length(parms.regpats)
  error('lengths of tags and regpats do not match');
end;

if ~exist(indir,'dir')
  errcode = 1;
  msg = sprintf('input directory %s not found',indir);
  return;
end;
flist = recursive_dir(char(indir));

% remove 'DICOMDIR' (dicom media directory file) from flist (if present)
ok_inds = ~find_cell(regexp(flist,'/DICOMDIR$|^DICOMDIR$'));
flist = flist(ok_inds);

nfiles = length(flist);
if ~nfiles
  errcode = 1;
  msg = sprintf('no files in %s',indir);
  return;
end;

% construct cmd
fprintf('%s: preparing to sort %d files...\n',mfilename,nfiles);
cmd = [];
for f=1:nfiles
  cmd = sprintf('%sdcmdump "%s" %s; echo; echo "%s"\n',...
    cmd,flist{f},sprintf(' --search %s',parms.tags{:}),parms.sep);
end;

fprintf('%s: sorting %d files...\n',mfilename,nfiles);
[errcode,msg] = mmil_unix(cmd,parms.batch_limit);

if errcode
  fprintf('%s: WARNING: error sorting files:\n%s\n',mfilename,msg);
  return;
end;

% split msg into cell array
results = regexp(msg,sprintf('\\n%s\\n',parms.sep),'split');
% parse results to get values for each file
values = cell(ntags,nfiles);
for t=1:ntags
  tag = parms.tags{t};
  regpat = parms.regpats{t};
  n = regexp(results,regpat,'names');
  for f=1:nfiles
    if isempty(n{f})
      if t==1
        values{t,f} = {parms.invalstr};
      else
        values{t,f} = {0};
      end;
    else
      tmp = {n{f}.(tag)};
      if t==1 && numel(tmp)>1 % handle Philips scans with two (or more?) SeriesInstanceUIDs
        tmp = tmp(end);
      end;
      values{t,f} = tmp;
    end;
  end;
end;

% exclude invalid files
IDs = [values{1,:}];
ind_invalid = strcmp(IDs,parms.invalstr);
invalid_files = flist(ind_invalid);
ind_valid = ~strcmp(IDs,parms.invalstr);
values = values(:,ind_valid);
IDs = [values{1,:}];
flist = flist(ind_valid);

% find unique IDs
[unique_IDs,ind_A,ind_B] = unique(IDs);
[sorted_ind_A,ind_sort] = sort(ind_A);
% sort unique IDs so earliest files in flist get lower series and study number
unique_IDs = unique_IDs(ind_sort);

% create output struct array containing ID, file names, and other info
for s=1:length(unique_IDs)
  SeriesInfo(s).(parms.tags{1}) = unique_IDs{s};
  ind_files = find(strcmp(IDs,unique_IDs{s}));
  ind_uniq = [];
  for t=2:ntags
    tag = parms.tags{t};
    if ismember(tag,parms.series_tags)
      vals = values{t,ind_files(1)};
    else
      vals = [values{t,ind_files}];
    end;
    if ismember(tag,parms.numeric_tags)
      vals = str2double(vals);
    elseif ismember(tag,parms.series_tags)
      vals = vals{1};
    end;
    SeriesInfo(s).(tag) = vals;
  end;
  SeriesInfo(s).FileNames = flist(ind_files);  
  SeriesInfo(s).EchoTimes = SeriesInfo(s).EchoTime;
  SeriesInfo(s).EchoTime = SeriesInfo(s).EchoTimes(1);
  SeriesInfo(s).InstanceNumbers = SeriesInfo(s).InstanceNumber;
  SeriesInfo(s).InstanceNumber = SeriesInfo(s).InstanceNumbers(1);
  SeriesInfo(s).errmsg = '';
  SeriesInfo(s) = check_info(SeriesInfo(s),s);
  SeriesInfo(s) = remove_duplicates(SeriesInfo(s),s);
end;

if isempty(SeriesInfo)
  errcode = 1;
  msg = 'empty SeriesInfo';
  fprintf('%s: WARNING: %s\n',mfilename,msg);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function serinfo = check_info(serinfo,s)
  nfiles = length(serinfo.FileNames);
  indir = fileparts(serinfo.FileNames{1});
  % handle double InstanceNumber issue with some Philips images
  if length(serinfo.InstanceNumbers) == 2*nfiles
    if length(unique(serinfo.InstanceNumbers(1:2:end))) == nfiles
       serinfo.InstanceNumbers = serinfo.InstanceNumbers(1:2:end);
    elseif length(unique(serinfo.InstanceNumbers(2:2:end))) == nfiles
       serinfo.InstanceNumbers = serinfo.InstanceNumbers(2:2:end);
    end
  end;
  % general solution to problem of extra EchoTimes or InstanceNumbers
  if nfiles ~= length(serinfo.EchoTimes) ||...
     nfiles ~= length(serinfo.InstanceNumbers)
    fprintf('%s: re-reading dicom headers for series %d (%s)...\n',...
      mfilename,s,indir);
    serinfo.EchoTimes = zeros(1,nfiles);
    serinfo.InstanceNumbers = zeros(1,nfiles);
    for f=1:nfiles
      tmpinfo = dicominfo(serinfo.FileNames{f});
      serinfo.EchoTimes(f) = ...
        mmil_getfield(tmpinfo,'EchoTime',-1);
      serinfo.InstanceNumbers(f) = ...
        mmil_getfield(tmpinfo,'InstanceNumber',-1);
    end;
    serinfo.EchoTime = serinfo.EchoTimes(1);
    serinfo.InstanceNumber = serinfo.InstanceNumbers(1);
  end;
  if any(serinfo.InstanceNumbers==-1)
    k = find(serinfo.InstanceNumbers>0);
    nfiles = length(k);
    if nfiles
      serinfo.FileNames = serinfo.FileNames(k);
      serinfo.EchoTimes = serinfo.EchoTimes(k);
      serinfo.InstanceNumbers = serinfo.InstanceNumbers(k);
      serinfo.EchoTime = serinfo.EchoTimes(1);
      serinfo.InstanceNumber = serinfo.InstanceNumbers(1);
    else
      serinfo.errmsg = sprintf(...
        '%s: WARNING: missing/invalid InstanceNumber tags in series %d (%s). Series will be excluded.',...
        mfilename,s,indir);
      fprintf('%s\n',serinfo.errmsg);
    end;
  end;
  if nfiles>1 && length(unique(serinfo.InstanceNumbers))==1
    serinfo.errmsg = sprintf('%s: WARNING: only one unique InstanceNumber tag (%d) in series %d (%s). Series will be excluded.',...
      mfilename,serinfo.InstanceNumber,s,indir);
    fprintf('%s\n',serinfo.errmsg);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function serinfo = remove_duplicates(serinfo,s)
  if ~isempty(serinfo.errmsg), return; end;
  indir = fileparts(serinfo.FileNames{1});
  uniq_InstanceNumbers = unique(serinfo.InstanceNumbers);
  if length(uniq_InstanceNumbers) ~= length(serinfo.FileNames)
    fprintf('%s: WARNING: duplicate InstanceNumbers in series %d (%s). Will use only unique InstanceNumbers.\n',...
      mfilename,s,indir);
    [uniq_vals,ind_uniq] = unique(serinfo.InstanceNumbers,'first');
    serinfo.FileNames = serinfo.FileNames(ind_uniq);
    serinfo.EchoTimes = serinfo.EchoTimes(ind_uniq);
    serinfo.InstanceNumbers = serinfo.InstanceNumbers(ind_uniq);
    serinfo.EchoTime = serinfo.EchoTimes(1);
    serinfo.InstanceNumber = serinfo.InstanceNumbers(1);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

