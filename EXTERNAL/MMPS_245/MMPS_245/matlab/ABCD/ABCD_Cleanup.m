function ABCD_Cleanup(ProjID,varargin)
%function ABCD_Cleanup(ProjID,[options])
%
% Purpose:
%   remove out-of-date processing output
%     i.e. processing results that are older than unpacking/recon
%
% Usage:
%  ABCD_Cleanup(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%
% Optional Parameters:
%   'rootdir': root directory containing project directory
%     with unpack, orig, pc, raw, proc, proc_dti, proc_bold dirs
%     {default = '/space/syn05/1/data/MMILDB'}
%   'sites': cell array of site names
%     if empty, use all site dirs within incoming dir
%     {default = []}
%   'ContainerTypes': cell array of container types to check
%     {default = {'pc' 'raw' 'proc' 'proc_dti' 'proc_bold'}}
%   'check_only_flag': [0|1] check dates only, do not remove anything
%     {default = 0}
%   'keep_qc_flag': do not remove proc, proc_dti, or proc_bold
%     dirs if they contain raw qcinfo files
%     {default = 0}
%   'verbose': [0|1] display status messages and warnings
%     {default: 1}
%   'logfile': name of file created to log messages
%     if empty, will be /home/{user}/logs/{ProjID}_cleanup.log
%     {default = []}
%
% Created:  12/06/16 by Don Hagler
% Prev Mod: 04/06/17 by Don Hagler
% Last Mod: 10/05/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(ProjID,varargin);

% create log file and display starting message
parms = create_log(parms);

% check containers for each site
check_sites(parms);

% close log file
fclose(parms.flog);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms = mmil_args2parms(options,{...
    'ProjID',ProjID,[],...
  ...
    'rootdir','/space/syn05/1/data/MMILDB',[],...
    'sites',[],[],...
    'ContainerTypes',{'pc' 'raw' 'proc' 'proc_dti' 'proc_bold'},{'pc' 'raw' 'proc' 'proc_dti' 'proc_bold' 'fsurf' 'fsico'},...
    'check_only_flag',false,[false true],...
    'keep_qc_flag',false,[false true],...
    'verbose',true,[false true],...
    'logfile',[],[],...
  });
  parms.ntypes = length(parms.ContainerTypes);

  % set RootDirs
  [ProjInfo,RootDirs] = MMIL_Get_ProjInfo(ProjID);
  dtypes = cat(2,{'unpack','orig'},parms.ContainerTypes);
  for i=1:length(dtypes)
    dname = dtypes{i};
    indir = mmil_getfield(RootDirs,dname);
    if isempty(indir)
      RootDirs.(dname) = sprintf('%s/%s/%s',parms.rootdir,ProjID,dname);
    end;
  end;
  parms.RootDirs = MMIL_Check_RootDirs(RootDirs);
  
  % get list of sites
  if isempty(parms.sites)
    dlist = dir(sprintf('%s/*',parms.RootDirs.unpack));
    parms.sites = setdiff({dlist.name},{'.','..'});
  end;
  parms.nsites = length(parms.sites);

  % default logfile
  if isempty(parms.logfile)
    logdir = sprintf('%s/logs',getenv('HOME'));
    mmil_mkdir(logdir);
    parms.logfile = sprintf('%s/%s_cleanup.log',logdir,ProjID);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = create_log(parms)
  % create logfile
  logdir = fileparts(parms.logfile);
  mmil_mkdir(logdir);
  parms.flog = fopen(parms.logfile,'a');
  if parms.flog==-1
    error('failed to open logfile %s for writing\n',parms.logfile);
  end;
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('%s: starting abcd cleanup for %s at %s...\n',...
    mfilename,parms.RootDirs.orig,datestr(now));
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(parms.flog,'%s: starting abcd cleanup for %s at %s...\n',...
    mfilename,parms.RootDirs.orig,datestr(now));
  fprintf(parms.flog,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_sites(parms)
  for s=1:parms.nsites
    site = parms.sites{s};
    % get list of directories for each site
    sitedir = sprintf('%s/%s',parms.RootDirs.unpack,site);
    dlist = dir(sitedir);
    % exclude '.' and '..'
    ind = find(~ismember({dlist.name},{'.','..'}));
    dlist = dlist(ind);
    ndirs = length(dlist);
    if parms.verbose
      fprintf('%s: %d unpack dirs to review for %s\n',...
        mfilename,ndirs,site);
    end;
    fprintf(parms.flog,'%s: %d unpack dirs to review for %s\n',...
      mfilename,ndirs,site);

    % loop over dirlist
    for d=1:ndirs
      dname = dlist(d).name;
      % find .unpacked .reconned and .renamed files
      fname_unpack = sprintf('%s/%s/%s.unpacked',sitedir,dname,dname);
      fname_recon = sprintf('%s/%s/%s.reconned',sitedir,dname,dname);
      fname_rename = sprintf('%s/%s/%s.renamed',sitedir,dname,dname);
      dinfo_unpack = dir(fname_unpack);
      dinfo_recon = dir(fname_recon);
      % find latest date
      unpack_datenum = max(cat(1,[dinfo_unpack.datenum],[dinfo_recon.datenum]));
      % get orig dir name
      if exist(fname_rename,'file')
        fid = fopen(fname_rename);
        origsubpath = fscanf(fid,'%s');
        fclose(fid);
        origpath = fileparts(origsubpath);
        [origrootpath,origdir] = fileparts(origpath);
        VisitID = origdir;
      else
        continue;
      end;
      % compare dates against dates of containers
      for k=1:parms.ntypes
        ContainerType = parms.ContainerTypes{k};
        check_container(VisitID,ContainerType,unpack_datenum,parms);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_container(VisitID,ContainerType,unpack_datenum,parms)
  if strcmp(ContainerType,'pc')
    ContainerPath = sprintf('%s/%s',parms.RootDirs.pc,VisitID);
  else
    ContainerPath = MMIL_Get_Container(parms.RootDirs,VisitID,ContainerType);
  end;
  if ~isempty(ContainerPath)
    if parms.verbose
      fprintf('%s: checking %s %s...\n',mfilename,VisitID,ContainerType);
    end;
    fprintf(parms.flog,'%s: checking %s %s...\n',mfilename,VisitID,ContainerType);
    % get file info about contents
    cont_info = dir(ContainerPath);
    % exclude '.' and '..' from cont
    ind = find(~ismember({cont_info.name},{'.','..'}));
    cont_info = cont_info(ind);
    % exclude 'txt' files from cont
    ind = find(cellfun(@isempty,regexp({cont_info.name},'\.txt$')));
    cont_info = cont_info(ind);
    % compare dates
    if unpack_datenum > min([cont_info.datenum])
      % check whether manual QC has been done
      if parms.keep_qc_flag
        if ismember(ContainerType,{'proc','proc_bold','proc_dti'})
          dlist = dir(sprintf('%s/*raw*qcinfo.txt',ContainerPath));
          if ~isempty(dlist)
            if parms.verbose
              fprintf('%s: NOTE: *not* removing %s Container from %s because is older than unpack, but raw QC exists\n',...
                mfilename,ContainerType,VisitID);
              fprintf('%s: NOT deleting %s ...\n',mfilename,ContainerPath);
            end;
            fprintf(parms.flog,'%s: NOTE: *not* removing %s Container from %s because is older than unpack, but raw QC exists\n',...
              mfilename,ContainerType,VisitID);
            fprintf(parms.flog,'%s: NOT deleting %s ...\n',mfilename,ContainerPath);
            return;
          end;
        end;
      end;
    
      % write messages
      if parms.verbose
        fprintf('%s: NOTE: removing %s Container from %s because is older than unpack\n',...
          mfilename,ContainerType,VisitID);
        fprintf('%s: deleting %s ...\n',mfilename,ContainerPath);
      end;
      fprintf(parms.flog,'%s: NOTE: removing %s Container from %s because is older than unpack\n',...
        mfilename,ContainerType,VisitID);
      fprintf(parms.flog,'%s: deleting %s ...\n',mfilename,ContainerPath);
      % delete container
      if parms.check_only_flag
        if parms.verbose
          fprintf('%s: NOT REALLY! just checking...\n',mfilename);
        end;
        fprintf(parms.flog,'%s: NOT REALLY! just checking...\n',mfilename);
      else
        cmd = sprintf('rm -r %s',ContainerPath);
        [s,r] = unix(cmd);
        if s, error('cmd %s failed:\n%s',cmd,r); end;
      end;
    end;
  end;
return;

