function ABCD_Import_Exams(ProjID,varargin)
%function ABCD_Import_Exams(ProjID,varargin)
%
% Required input:
%   ProjID: project ID  string
%
% Optional parameters:
%   'batchname': name of batch directory in /home/mmilrec/batchdirs
%     {default = 'ABCD_Import_Exams'}
%   'fname_sites': spreadsheet containing site names and IDs
%     {default = '/home/abcddaic/ProjInfo/ABCD/ABCD_Sites.csv'}
%   'phantom_flag': [0|1|2] whether to import phantom scans
%     0: exclude phantom scans (pGUID matches ABCDPhantomXXX)
%     1: allow phantom scans
%     2: only phantom scans
%     {default = 0}
%   'humanphantom_flag': [0|1|2] whether to import human phantom scans
%     0: exclude human phantom scans (pGUID contains PhantomTravelingHuman or PhantomLocalHuman))
%     1: allow human phantom scans
%     2: only human phantom scans
%     {default = 0}
%   'testQA_flag': [0|1|2] whether to import test QA scans
%     0: exclude test QA scans (pGUID contains TestQA)
%     1: allow test QA scans
%     2: only test QA scans
%     {default = 0}
%   'PMU_flag':  [0|1|2] whether to import test PMU files
%     0: exclude PMU files (SeriesDescription contains PMU)
%     1: allow PMU files
%     2: only PMU files
%     {default = 0}
%   'require_subjid_flag': [0|1] whether to skip series without a SubjID in fstem
%     {default = 1}
%   'mb_recon_flag': [-1|0|1|2] whether to perform multi-band recon
%     -1: unpack only. This will work with fast-track data
%     0: no multi-band recon
%     1: unpack and perform multi-band recon if required
%     2: do not unpack unless multi-band recon is required
%     {default = 1}
%   'mb_tmpdir': temporary directory that will contain untar'd P-files
%     {default = '/scratch'}
%   'mb_dicoms_flag': [0|1] whether to write dicoms after multi-band recon
%     0: produce mat files only
%     1: produce mat files and dicoms
%     2: write dicoms only (do not produce jobs if recon has not been run yet)
%     {default = 1}
%   'mb_cleanup_flag': [0|1|2] whether to cleanup temporary mb output
%     0: perform no cleanup
%     1: remove P-file when complete
%     2: remove mb_tmpdir when complete
%     {default = 2}
%   'separate_series_flag': [0|1] create separate output directories
%     for each series
%     {default = 0}
%   'check_dates_flag': [0|1] compare date of unpack to tgz in incoming
%     force import if unpack is older than tgz
%     {default = 0}
%   'check_md5sum_flag': [0|1] compare contents of md5sum files in incoming
%     force import if there is a mismatch with file in unpack
%     {default = 0}
%   'check_md5sum_server_flag': [0|1] compare contents of md5sum and md5sum_server files
%     skip import if there is a mismatch or one does not exist
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}     
%   'force_err_flag': [0|1] overwrite existing output if previous errors
%     {default = 0}
%   'force_unpack_err_flag': [0|1] overwrite existing output if unpack errors
%     {default = 0}
%   'force_recon_err_flag': [0|1] overwrite existing output if recon errors
%     {default = 0}
%   'force_dicoms_err_flag': [0|1] overwrite existing output if previous dicoms errors
%     for multi-band recon only
%     {default = 0}
%   'force_rename_err_flag': [0|1] overwrite existing output if previous rename errors
%     {default = 0}
%   'rootdir': root directory containing project directory
%     with incoming, unpack, and orig dirs
%     unless specified in ProjInfo
%     {default = '/space/syn05/1/data/MMILDB'}
%   'sites': cell array of site names
%     if empty, use all site dirs within incoming dir
%     {default = []}
%
% Created:  09/08/16 by Don Hagler
% Prev Mod: 08/25/17 by Don Hagler
% Last Mod: 09/07/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'batchname','ABCD_Import_Exams',[],...
  'fname_sites','/home/abcddaic/ProjInfo/ABCD/ABCD_Sites.csv',[],...
  'phantom_flag',0,[0:2],...
  'humanphantom_flag',0,[0:2],...
  'testQA_flag',0,[0:2],...
  'PMU_flag',0,[0:2],...
  'require_subjid_flag',true,[false true],...
  'mb_recon_flag',1,[-1:2],...
  'mb_tmpdir','/scratch',[],...
  'mb_dicoms_flag',1,[0:2],...
  'mb_cleanup_flag',2,[0:2],...
  'separate_series_flag',false,[false true],...
  'check_dates_flag',false,[false true],...
  'check_md5sum_flag',false,[false true],...
  'check_md5sum_server_flag',false,[false true],...
  'forceflag',false,[false true],...
  'force_err_flag',false,[false true],...
  'force_unpack_err_flag',false,[false true],...
  'force_recon_err_flag',false,[false true],...
  'force_dicoms_err_flag',false,[false true],...
  'force_rename_err_flag',false,[false true],...
  'rootdir','/space/syn05/1/data/MMILDB',[],...
  'sites',[],[],...
  ...
  'phantom_pat_list',{'ABCDPhantom'},[],...
  'humanphantom_pat_list',{'PhantomTravelingHuman','PhantomLocalHuman'},[],...
  'testQA_pat_list',{'TestQA'},[],...
  ...
  'tags',{'batchname','batch_append_flag','batch_msg_flag',...
         'indir','unpackdir','outdir','fname_sites',...
         'phantom_flag','humanphantom_flag','testQA_flag',...
         'PMU_flag','require_subjid_flag'...
         'mb_recon_flag','mb_tmpdir','mb_dicoms_flag','mb_cleanup_flag',...
         'separate_series_flag',...
         'check_dates_flag','check_md5sum_flag','check_md5sum_server_flag',...
         'forceflag','force_err_flag',...
         'force_unpack_err_flag','force_recon_err_flag','force_dicoms_err_flag',....
         'force_rename_err_flag','logfile',...
         'phantom_pat_list','humanphantom_pat_list','testQA_pat_list'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID);
if ~isempty(mmil_getfield(RootDirs,'incoming'))
  indir = RootDirs.incoming;
else
  indir = sprintf('%s/%s/incoming',parms.rootdir,ProjID);
end;
if ~isempty(mmil_getfield(RootDirs,'unpack'))
  unpackdir = RootDirs.unpack;
else
  unpackdir = sprintf('%s/%s/unpack',parms.rootdir,ProjID);
end;
if ~isempty(mmil_getfield(RootDirs,'orig'))
  outdir = RootDirs.orig;
else
  outdir = sprintf('%s/%s/orig',parms.rootdir,ProjID);
end;

if ~exist(indir,'dir')
  error('input directory %s not found',indir);
end;

parms.batchname = [ProjID '_' parms.batchname];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch parms.mb_recon_flag
  case 0
    suffix = 'nomb';
  case 1
    suffix = 'all';
  case 2
    suffix = 'mb';
end;

parms.logfile = sprintf('%s/logs/%s_abcd_import_%s.log',...
                        getenv('HOME'),ProjID,suffix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get list of sites
if isempty(parms.sites)
  dlist = dir(sprintf('%s/*',indir));
  sites = setdiff({dlist.name},{'.','..'});
else
  sites = parms.sites;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove logfile
if exist(parms.logfile,'file')
  delete(parms.logfile);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create jobs for each site
parms.batch_append_flag = 0;
nsites = length(sites);
for s=1:nsites
  site = sites{s};
  if ~isempty(regexp(site,'\.lock')), continue; end;
  parms.indir = [indir '/' site];
  if ~exist(parms.indir,'dir')
    fprintf('%s: WARNING: site dir %s not found\n',...
      mfilename,parms.indir);
    continue;
  end;
  parms.unpackdir = [unpackdir '/' site];
  parms.outdir = outdir;
  if s==nsites
    parms.batch_msg_flag = 1;
  else
    parms.batch_msg_flag = 0;
  end;
  args = mmil_parms2args(parms,parms.tags);
  abcd_import_all(args{:});
  parms.batch_append_flag = 1;
end;

