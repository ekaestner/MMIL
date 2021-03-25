function [SeriesInfo,errcode,msg] = mmil_copy_dicoms(SeriesInfo,varargin)
%function [SeriesInfo,errcode,msg] = mmil_copy_dicoms(SeriesInfo,[options])
%
% Required Input:
%   SeriesInfo: struct array containing FileNames and values for each series
%     'SeriesInstanceUID','SeriesNumber', 'StudyInstanceUID',InstanceNumbers'
%     all non-dicom or unreadable files will be placed in SeriesInfo(end)
%       with SeriesInstanceUID = 'invalid'
%     see mmil_sort_dicoms
%
% Optional Parameters:
%   'outdir': output directory; subdirectories will be created for each series
%     {default = pwd}
%   'linkflag': [0|1] create symbolic links instead of copying
%     {default = 1}
%   'batch_limit': maximum number of lines per unix call
%     {default = 250}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Output:
%   errcode: returns 1 if error, 0 if successful
%   msg: error message
%
% Created:  04/15/11 by Don Hagler
% Last Mod: 08/22/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir',pwd,[],...
  'linkflag',true,[false true],...
  'batch_limit',250,[],...
  'forceflag',false,[false true],...
});
errcode = 0;
msg = [];

if parms.linkflag
  parms.copy_cmd = 'ln -s';
else
  parms.copy_cmd = 'cp';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StudyInstanceUIDs = unique({SeriesInfo.StudyInstanceUID});
for s=1:length(SeriesInfo)
  SeriesNumber = SeriesInfo(s).SeriesNumber;
  StudyNumber = find(strcmp(SeriesInfo(s).StudyInstanceUID,StudyInstanceUIDs));
  SeriesInfo(s).SeriesDirPath = sprintf('%s/st%03d_ser%04d',...
    parms.outdir,StudyNumber,SeriesNumber);
  if isfield(SeriesInfo(s),'errmsg') && ~isempty(SeriesInfo(s).errmsg)
    continue;
  end;
  mmil_mkdir(SeriesInfo(s).SeriesDirPath);
  cmd = [];
  SeriesInfo(s).OrigFileNames = SeriesInfo(s).FileNames;
  for f=1:length(SeriesInfo(s).FileNames)
    fname = SeriesInfo(s).FileNames{f};
    fname_out = sprintf('%s/im%04d.dcm',...
      SeriesInfo(s).SeriesDirPath,SeriesInfo(s).InstanceNumbers(f));
    % remove existing link/copy if forceflag
    if parms.forceflag
      cmd = sprintf('%s set f = "%s"; if (-e "$f") rm -f "$f"\n',...
        cmd,fname_out);
    end;
    cmd = sprintf('%s set f = "%s"; if (! -e "$f") %s "%s" "$f"\n',...
      cmd,fname_out,parms.copy_cmd,fname);
    % update file name
    SeriesInfo(s).FileNames{f} = fname_out;
  end;
  [errcode,msg] = mmil_unix(cmd,parms.batch_limit);
  if errcode
    if isempty(msg), msg = 'unspecified error copying dicoms'; end;
    return;
  end;
  [sorted_INs,ind_sort] = sort(SeriesInfo(s).InstanceNumbers);
  SeriesInfo(s).InstanceNumbers = sorted_INs;
  SeriesInfo(s).FileNames = SeriesInfo(s).FileNames(ind_sort);
  SeriesInfo(s).EchoTimes = SeriesInfo(s).EchoTimes(ind_sort);
end;

