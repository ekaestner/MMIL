function abcd_summarize_pc(indir,varargin)
%function abcd_summarize_pc(indir,[options])
%
% required input:
%   indir: input directory containing subdirectories with
%          protocol compliance json files
%
% optional input:
%   'auxdir': auxiliary directory containing incoming tar files
%     {default = '/space/syn05/1/data/MMILDB/DAL_ABCD/incoming'}
%   'outdir': output directory
%     {default = indir}
%   'outstem': output file stem
%     {default = 'pcinfo'}
%   'col_order': column names in the desired order
%      for output csv file arranged by series
%      additional columns will be moved to the end
%     {default = {'pGUID','VisitID','EventName','SessionType','SiteName',...
%                        'SeriesDescription','SeriesType','ABCD_Compliant','Completed',...
%                        'AdditionalInfo''}}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%
% Created:  08/19/16 by Jose Teruel
% Last Mod: 01/17/17 by Don Hagler
%

% NOTE: based on json_csv_parser_pc, last mod 08/19/16, created by Jose
% Teruel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input
if ~mmil_check_nargs(nargin,1), return; end;
parms = check_input(indir,varargin);

if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end;
mmil_mkdir(parms.outdir);

% check output, maybe quit depending on forceflag
output_exists_flag = check_output(parms);
if output_exists_flag && ~parms.forceflag
  fprintf('%s: WARNING: not overwriting existing output\n',mfilename);
  return;
end;

% find all json files
jlist = find_json_files(parms);

% read each json file
jinfo = read_json_files(jlist,parms);

% write csv file with one row for each series
write_output(jinfo,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(indir,options)
  parms = mmil_args2parms(options,{...
    'indir',indir,[],...
    ...
    'auxdir','/space/syn05/1/data/MMILDB/DAL_ABCD/incoming',[],...
    'outdir',indir,[],...
    'outstem','pcinfo',[],...
    'col_order',{'pGUID','VisitID','EventName','SessionType','SiteName',...
                        'SeriesDescription','SeriesType','ABCD_Compliant','Completed',...
                        'AdditionalInfo'},[],...
    'forceflag',false,[false true],...
  });
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output_exists_flag = check_output(parms)
  output_exists_flag = true;
  suffix_list = {'series','sessions'};
  for i=1:length(suffix_list)
    suffix = suffix_list{i};
    fname_out = sprintf('%s/%s_%s.csv',parms.outdir,parms.outstem,suffix);
    if ~exist(fname_out,'file')
      output_exists_flag = false;
      break;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jlist = find_json_files(parms)
  fprintf('%s: finding json files...\n',mfilename);
  jlist = [];
  sep = '#####';
  cmd = sprintf('echo "%s"; find %s -type f -name "*.json"; echo "%s"',...
    sep,parms.indir,sep);
  [s,r] = unix(cmd);
  if s, error('cmd %s failed:\n%s',cmd,r); end;
  q = regexp(r,sprintf('%s\\n',sep),'split');
  jlist = regexp(strtrim(q{2}),'\n','split');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jinfo = read_json_files(jlist,parms)
  jinfo = [];
  nfiles = length(jlist);
  fprintf('%s: reading %d json files...\n',mfilename,nfiles);
  k = 1;
  for j=1:length(jlist)
    fname = jlist{j};
    % read json file
    try
      tinfo = loadjson(jlist{j});
    catch me
      fprintf('%s: WARNING: error reading json file %s\n',mfilename,jlist{j});
      continue;
    end;
    tinfo = tinfo.ProtocolCompliance;
    % get site ID and name
    [~,tinfo.VisitID] = fileparts(fileparts(fileparts(fname)));
    tinfo.SiteID = tinfo.VisitID(1:4);
    tinfo.SiteName = get_SiteName(tinfo.SiteID);
    % determine the pGUID
    [tinfo.pGUID,tinfo.EventName,tinfo.SessionType,tinfo.fname_json] = get_pGUID(tinfo,parms);
    % censor certain fields by manufacturer
    tinfo.fname_pc_json = fname;
    % get fieldnames from struct
    tags = fieldnames(tinfo);
    % reorder tags
    tags = reorder_tags(tags,parms.col_order);
    % copy values to struct array
    for t=1:length(tags)
      tag = tags{t};
      jinfo(k).(tag) = tinfo.(tag);
    end;
    k = k + 1;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tags = reorder_tags(tags,tag_order)
  [tmp,i_first,i_order] = intersect(tags,tag_order);
  [tmp,i_sort] = sort(i_order);
  i_first = i_first(i_sort);
  i_last = setdiff([1:length(tags)]',i_first);
  i_new = cat(1,mmil_colvec(i_first),mmil_colvec(i_last));
  tags = tags(i_new);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SiteName = get_SiteName(SiteID)
  switch SiteID
    case 'G010'
      SiteName = 'sri';
    case 'G031'
      SiteName = 'daic';
    case 'G032'
      SiteName = 'libr';
    case 'G054'
      SiteName = 'mssm';      
    case 'G075'
      SiteName = 'umich';
    case 'G087'
      SiteName = 'uwm'; 
    case 'P023'
      SiteName = 'vcu';
    case 'P043'
      SiteName = 'uvm';  
    case 'P064'
      SiteName = 'chla';
    case 'S011'
      SiteName = 'fiu';    
    case 'S012'
      SiteName = 'oahu';
    case 'S013'
      SiteName = 'upmc';
    case 'S014'
      SiteName = 'utah';
    case 'S020'
      SiteName = 'umn';      
    case 'S021'
      SiteName = 'washu';   
    case 'S022'
      SiteName = 'cub';
    case 'S042'
      SiteName = 'ohsu';
    case 'S053'
      SiteName = 'yale';      
    case 'S065'
      SiteName = 'ucla';
    case 'S076'
      SiteName = 'ufl';
    case 'S086'
      SiteName = 'musc';   
    case 'S101' % for PCGC only
      SiteName = 'bch';
    case 'S102' % for PCGC only
      SiteName = 'chop';
    case 'G103' % for PCGC only
      SiteName = 'ucsf';
    otherwise
      SiteName = 'unknown';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pGUID,EventName,SessionType,fname_json] = get_pGUID(tinfo,parms)
  pGUID = []; EventName = []; SessionType = []; fname_json = [];
  [pGUID,EventName,SessionType,fname_json] = get_auxdir_pGUID(tinfo,parms);
  if isempty(pGUID), pGUID = 'invalid'; end;
  if isempty(EventName), EventName = 'invalid'; end;
  if isempty(SessionType), SessionType = 'invalid'; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lookup NDAR ID from incoming json file
function [pGUID,EventName,SessionType,fname_json] = get_auxdir_pGUID(tinfo,parms)
  pGUID = []; EventName = []; SessionType = []; fname_json = [];
  fpat = sprintf('%s/%s/*%s_%s.j*',...
    parms.auxdir,tinfo.SiteName,tinfo.StudyInstanceUID,tinfo.SeriesInstanceUID);
%  fprintf('%s: getting auxdir pGUID for %s...\n',mfilename,tinfo.SeriesInstanceUID);
  flist = dir(fpat);
  if isempty(flist) % NOTE: offline recon creates new SeriesInstanceUID
    fpat = sprintf('%s/%s/*%s*.j*',...
      parms.auxdir,tinfo.SiteName,tinfo.StudyInstanceUID);
    flist = dir(fpat);
  end;
  if ~isempty(flist)
    fname_json = flist(1).name;
    %% NOTE: for data with offline recon, fname_json is first from this study
  else
    fname_json = 'missing';
  end;
%  [SubjID,pGUID,EventName,SessionType] = abcd_get_SubjID(fname_json,tinfo);    
  [SubjID,pGUID,EventName,SessionType] = abcd_get_SubjID(fname_json);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_output(jinfo,parms)
  fname_out = sprintf('%s/%s.csv',parms.outdir,parms.outstem);
  fprintf('%s: writing output to %s...\n',mfilename,fname_out);
  mmil_struct2csv(jinfo,fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

