function abcd_check_kspace_tgz(varargin)
%function abcd_check_kspace_tgz(varargin)
% Purpose: Check if a series has issues of:
%        1. missing tgz
%        2. is kspace
%
% Optional input:
%   'indir': metadata directory (contains summary spreadsheets)
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'indir_jq_cache': jq cache dir
%     {default = '/dev/shm/jq_cache'}
%   'outdir': output directory
%     {default = 'MetaData/DAL_ABCD_QC'}
%   'fname_info': input file containing incoming info
%     {default = []}
%   'fname_kt': input file containing about kspace and missingtgz
%     {default = []}
%   'fname_output': output file name
%     {default = []}
%   'fname_ks_only': output file name of kspace only info
%     {default = []}
%   'instem': instem
%     {default = 'DAL_ABCD_QC'}
%   'outstem': outstem
%     {default = 'DAL_ABCD_QC'}
%   'infix': file suffix of input file
%     containing info about incoming info
%     {default = 'combined_incoming_info'}
%   'infix_kt': file suffix of input file
%     containing info about unique events in kspace_tgz directory
%     {default = 'kspace_tgzissue'}
%   'outfix : file suffix for output
%     {default = 'combined_incoming_info_kspace_tgz'}
%   'outfix_ks_only : file suffix for output file of kspace only info
%     {default = 'combined_incoming_info_kspace_only'}
%   'forceflag': ignore existing results and run
%     {default = 1}
%
% Created:  04/26/17 by Feng Xue
% Last Mod: 07/05/17 by Feng Xue
%

  %check lock files

  fname_lck=sprintf('%s/MetaData/DAL_ABCD_QC/.check_kspace_tgz.lck',getenv('HOME'));
  if exist(fname_lck,'file')
    fprintf('%s\n','lock files exist!.');
    return;
  end
  %Place lock file
  fclose(fopen(fname_lck, 'w'));

  parms=check_input(varargin);

  incoming_info = abcd_load_csv(parms.fname_info);
  kspace_tgz_info = mmil_csv2struct(parms.fname_kt);

  [idx idx]=unique(regexprep({kspace_tgz_info.json_fstem},'[Ss]ession[^_]*_',''),'last');
  kspace_tgz_info=kspace_tgz_info(idx);
  json_ks=regexprep({kspace_tgz_info.json_fstem},'[Ss]ession[^_]*_','');
  json_incoming=regexprep({incoming_info.json},'[Ss]ession[^_]*_|\.json','');
  id_common=find(ismember(json_ks,json_incoming));
  kspace_tgz_common=kspace_tgz_info(id_common);
  kspace_tgz_info=kspace_tgz_info(setdiff(1:length(json_ks),id_common));
  mmil_struct2csv(kspace_tgz_info,parms.fname_output);

  id_kspace=[];
  for i=1:length(kspace_tgz_common)
    [~,targets]=unix(sprintf('readlink %s/jq_tmp_hastgz_*/%s.json',parms.indir_jq_cache,regexprep(kspace_tgz_common(i).json_fstem,'[Ss]ession[^_]*_','*')));
    if length(regexp(targets,'syn07|syn08')) == 0
      id_kspace=[id_kspace i];
    end
  end
  if ~isempty(id_kspace)
    kspace_valid=kspace_tgz_common(id_kspace);
    mmil_struct2csv(kspace_valid,parms.fname_ks_only);
  end
  
  delete(parms.fname_kt);

  %delete lock file
  delete(fname_lck);
return

function parms = check_input(options)
  parms = mmil_args2parms(options,{...
    'indir','MetaData/DAL_ABCD_QC',[],...
    'indir_jq_cache','/dev/shm/jq_cache',[],...
    'outdir','MetaData/DAL_ABCD_QC',[],...
    'fname_info',[],[],...
    'fname_kt',[],[],...
    'fname_output',[],[],...
    'fname_ks_only',[],[],...
    'instem','DAL_ABCD_QC',[],...
    'outstem','DAL_ABCD_QC',[],...
    'infix','combined_incoming_info',[],...
    'infix_kt','kspace_tgzissue',[],...
    'outfix','combined_incoming_info_kspace_tgz',[],...
    'outfix_ks_only','combined_incoming_info_kspace_only',[],...
    'forceflag',true,[false true],...
  });
  if parms.outdir(1) ~= '/', parms.outdir = sprintf('%s/%s',getenv('HOME'),parms.outdir); end;
  if parms.indir(1) ~= '/', parms.indir = parms.outdir; end;
  if isempty(parms.fname_info)
    parms.fname_info = sprintf('%s/%s_%s.csv',...
      parms.indir,parms.instem,parms.infix);
  end;
  if isempty(parms.fname_kt)
    parms.fname_kt = sprintf('%s/%s_%s.csv',...
      parms.indir,parms.instem,parms.infix_kt);
  end;
  if isempty(parms.fname_output)
    parms.fname_output = sprintf('%s/%s_%s.csv',...
      parms.outdir,parms.outstem,parms.outfix);
  end;
  if isempty(parms.fname_ks_only)
    parms.fname_ks_only = sprintf('%s/%s_%s.csv',...
      parms.outdir,parms.outstem,parms.outfix_ks_only);
  end;
  if ~exist(parms.fname_info,'file')
    error('info file %s not found',parms.fname_info);
  end;
  if ~exist(parms.fname_kt,'file')
    error('info file %s not found',parms.fname_kt);
  end;
return;
