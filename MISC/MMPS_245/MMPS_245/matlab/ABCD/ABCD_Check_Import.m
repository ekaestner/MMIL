function ABCD_Check_Import(ProjID,varargin)
%function ABCD_Check_Import(ProjID,[options])
%
% Purpose: generate jobs for event level incoming checkup.
%
% Required Input:
%   ProjID: project ID string
%
% Optional Input: ('key', value pairs)
%   'batchname': name of output batchdir
%     {default = 'check_import'}
%   'rootdir': root input directory containing container type dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%   'fname_info': spreadsheet containing series info
%     if empty, will use incoming_unique_enrolled_info from MetaData
%     {default = []}
%   'outdir': output directory
%     if not  specified, output in MetaData/{ProjID}
%     {default = []}
%   'outstem': output file stem
%     if not  specified, {ProjID}_import
%     {default = []}
%   'infix': file suffix of input file
%            containing info about unique events in incoming directory
%     {default = 'incoming_info'}
%   'outfix : file suffix for output file
%            containing info about status of unpack / recon
%     {default = 'import_info'}
%   'sites' : cell array of sites to check
%     if empty, check all sites
%     {default = []}
%   'forceflag': ignore checkpoints and run on all series
%     {default = 1}
%
% Created:  03/03/17 by Don Hagler
% Last Mod: 05/19/17 by Feng Xue
%

%   'forceflag': overwrite existing output
%     {default = 0}
%  'forceflag',false,[false true],...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'batchname','check_import',[],...
  'rootdir','/space/syn05/1/data/MMILDB',[],...
  'fname_info',[],[],...
  'outdir',[],[],...
  'outstem',[],[],...
  'infix','combined_incoming_info_classified',[],...
  'outfix','import_info',[],...
  'sites',[],[],...
  'forceflag',true,[false true],...
  'tags',{'indir','outdir','fname_projinfo',...
    'instem','outstem','infix','outfix','site','forceflag'},[],...
});
%  'tags',{'incoming_dir','unpack_dir','indir','outdir',...
%    'instem','outstem','infix','outfix','site','forceflag'},[],...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get project information and root dirs
[ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID);
% create output batch directory
if ~isempty(ProjID)
  parms.batchname = [ProjID '_' parms.batchname];
end;
batchdir = [RootDirs.batch '/' parms.batchname];
scriptlistfname = sprintf('%s/scriptlist.txt',batchdir);
if exist(batchdir,'dir')
  cmd = sprintf('rm -rf %s/*\n',batchdir);
  fprintf('cmd = %s',cmd);
  unix(cmd);
end;
mmil_mkdir(batchdir);

fid = fopen(scriptlistfname,'w');
if fid==-1
  error('failed to open scriptlist file %s for writing\n',scriptlistfname);
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parms.incoming_dir = mmil_getfield(RootDirs,'incoming',sprintf('%s/%s/incoming',parms.rootdir,ProjID));
%parms.unpack_dir = mmil_getfield(RootDirs,'unpack',sprintf('%s/%s/unpack',parms.rootdir,ProjID));
%if ~exist(parms.incoming_dir,'dir')
%  error('input directory %s not found',parms.incoming_dir);
%end;
%if ~exist(parms.unpack_dir,'dir')
%  error('unpack directory %s not found',parms.unpack_dir);
%end;
parms.indir = sprintf('%s/MetaData/%s',RootDirs.home,ProjID);
parms.outdir = sprintf('%s/MetaData/%s/import_reports',RootDirs.home,ProjID);
parms.instem = ProjID;
parms.outstem = ProjID;
if isempty(parms.fname_info)
  parms.fname_info = sprintf('%s/%s_%s.csv',parms.indir,parms.instem,parms.infix);
end;
if ~exist(parms.fname_info,'file')
  error('info file %s not found',parms.fname_info);
end;

if isempty(parms.sites)
  incoming_info = abcd_load_csv(parms.fname_info);
  %dlist = dir(sprintf('%s/*',parms.incoming_dir));
  %parms.sites = setdiff({dlist.name},{'.','..'});
  parms.sites = unique({incoming_info.site});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create jobs for each site
j = 1;
for i=1:length(parms.sites)
  parms.site = parms.sites{i};
  jstem = parms.site;
  jstem = jstem(1:min(20,length(jstem)));
  jobID = sprintf('job_%03d_%s',j,jstem);
  jobfname = [batchdir '/' jobID '.m'];
  mmil_write_script(jobfname,'abcd_check_import',...
    {},parms.tags,parms);
  fid = fopen(scriptlistfname,'a');
  fprintf(fid,'%s\n',jobID);
  fclose(fid);
  j=j+1;
end;

fprintf('%%%% Now login to a cluster and run this:\n');
fprintf('    qmatjobs %s\n',parms.batchname);

fprintf('%%%% or run:\n');
fprintf('    bmatjobs %s 12\n',parms.batchname);

