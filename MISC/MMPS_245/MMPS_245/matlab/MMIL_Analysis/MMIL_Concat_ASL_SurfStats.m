function [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_ASL_SurfStats(ProjID,varargin)
%function [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_ASL_SurfStats(ProjID,[options])
%
% Usage:
%  [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_ASL_SurfStats(ProjID,'key1', value1,...);
%
% Required Input:
%   ProjID: Project ID string
%     used to to load ProjInfo and StudyInfo from user's home
%       (e.g. '/home/{user}/ProjInfo/MMIL_ProjInfo.csv'
%             '/home/{user}/ProjInfo/{ProjID}/{ProjID}_VisitInfo.csv' )
%     may be empty if StudyInfo and RootDirs are supplied directly
%
% Optional Parameters that specify study specific information
%  'StudyInfo': struct array of study information
%       (e.g. read from csv file with MMIL_Read_StudyInfo)
%     If empty, will use ProjID to get StudyInfo
%     {default = []}
%  'RootDirs': struct that must contain the following fields:
%      proc_asl, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%
% Optional Parameters:
%  'measlist': list of ASL "measures" to extract for each ROI
%     e.g. 'CBF', cerebral blood flow
%     {default = {'CBF'}}
%  'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast maps
%     {default = 1}
%  'gwnorm_flag': [0|1] for gray/white contrast, normalize difference by mean
%     {default = 0}
%  'smoothing': surface smoothing steps (on sphere)
%     slope of FWHM vs. sqrt(N) is ~1.13 for fsaverage
%     (FWHM = full-width-half-max smoothing kernel
%         N = number of smoothing steps)
%    {default = 0}
%  'analdir': name of analysis folder relative to ContainerPath
%    {default = 'analysis'}
%  'outdir': output directory
%     full path or relative to /home/{user}/MetaData/{ProjID}
%     {default = 'ASL_SurfStats'}
%  'outstem': output file stem
%     relative to outdir unless full path given
%     {default = 'ASL'}
%  'meanflag': [0|1] calculate mean instead of concatenating
%     {default = 0}
%  'options': option string to use any of the mri_concat command line options
%     {default = []}
%  'verbose': [0|1|2] display status messages
%     0: no messages except errors
%     1: no messages except WARNING
%     2: frequent status messages
%     {default = 1}
%  'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fnames_mgz: cell array of output mgz file names
%     that contain concatenated surface stats
%     size = [nhemi,nsubs] or [nhemi,1] if seed_roinames is empty
%   fname_csv: output csv file name containing relevant info for each visit
%   fnames_log: cell array of output log file names
%     that contain lists of input file names
%     size = [nhemi,nseeds] or [nhemi,1] if seed_roinames is empty
%
% Created:  04/26/15 by Don Hagler
% Last Mod: 04/26/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
fnames_mgz = []; fname_csv = []; fnames_log = [];

parms = check_input(ProjID,varargin);
if parms.nsubs==0, return; end;

[fnames_mgz,fnames_log] = concat_files(parms);

fname_csv = write_csv(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(ProjID,options)
  parms_filter = {...
    'ProjID',ProjID,[],...
    'StudyInfo',[],[],...
    'RootDirs',[],[],...
    'qcflag',true,[false true],...
... % specify which results to compile
    'measlist',{'CBF'},[],...
    'projdist_list',[1],[-5,5],...
    'gwnorm_flag',false,[false true],...
    'smoothing',0,[0,Inf],...
    'analdir','analysis',[],...
...
    'outdir','ASL_SurfStats',[],...
    'outstem','ASL',[],...
    'meanflag',false,[false true],...
    'options',[],[],...
    'verbose',1,[0:2],...
    'forceflag',false,[false true],...
...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'intype','mgz',{'mgh','mgz'},...
    'outtype','mgz',{'mgh','mgz'},...
    'continfo_flag',false,[false true],...
    'subjinfo_flag',false,[false true],...
...
    'info_tags',{'Age','Sex','Site','Group',...
                 'Manufacturer','ManufacturersModelName',...
                 'DeviceSerialNumber','MagneticFieldStrength',...
                 'MMPS_version','ProcDate','StudyDate'},[],...
...
    'concat_tags',{'meanflag','options','forceflag'},[],...
... % hidden
    'required_containers',{'proc_asl','fsurf'},[],...
    'QC_raw',true,[false true],... % only applies if manual raw QC exists
    'QC_recon',true,[false true],...
    'contrast_outfix','contrast',[],...
    'resample_flag',true,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(StudyInfo), error('empty StudyInfo'); end;
  parms.StudyInfo = StudyInfo;
  parms.RootDirs = RootDirs;
  parms.nsubs = length(parms.StudyInfo);
  
  parms.nhemi = length(parms.hemilist);
  % whether to calculate gray-white contrast depends on if 2 projdist
  parms.npdist = length(parms.projdist_list);
  if parms.npdist==2
    parms.contrast_flag = 1;
    parms.projdist_list = sort(parms.projdist_list); % white matter first
    parms.contflag_list = [0,1];
  else
    parms.contrast_flag = 0;
    parms.contflag_list = 0;
  end;
  parms.nsmooth = length(parms.smoothing);
  parms.nmeas = length(parms.measlist);

  parms = check_files(parms);
  if parms.nsubs==0
    fprintf('%s: ERROR: no valid subjects\n',mfilename);
    return;
  end;
  
  % get container-related info
  if parms.continfo_flag
    parms.StudyInfo = MMIL_Get_ContInfo(parms.StudyInfo,parms.RootDirs);
  end;

  % get subject-related info
  if parms.subjinfo_flag
    parms.StudyInfo = ...
      MMIL_Get_SubjInfo(parms.StudyInfo,parms.RootDirs,parms.ProjID);
  end;

  if mmil_isrelative(parms.outstem)
    if mmil_isrelative(parms.outdir)
      parms.outdir = [getenv('HOME') '/MetaData/' ProjID '/' parms.outdir];
    end;
    parms.outstem = [parms.outdir '/' parms.outstem];
  else
    parms.outdir = fileparts(parms.outstem);
  end;
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runflag,fnames_mgz,fname_csv] = check_output_files(parms)
  fnames_mgz = []; fname_csv = [];
  runflag = 0;
  fname_csv = set_fname_csv(parms);
  if ~exist(fname_csv,'file'), runflag = 1; return; end;
  fnames_mgz = cell(parms.nmeas,parms.nhemi,parms.nsmooth,parms.npdist+parms.contrast_flag);
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    for c=1:length(parms.contflag_list)
      contflag = parms.contflag_list(c);
      if contflag
        projdist_list = parms.projdist_list(2);
        npdist = 1;
      else
        projdist_list = parms.projdist_list;
        npdist = parms.npdist;
      end;
      for p=1:npdist
        if contflag
          pind = parms.npdist + 1;
        else
          pind = p;
        end;
        for h=1:parms.nhemi
          for k=1:parms.nsmooth
            fname_out = set_fname_out(parms,meas,contflag,p,k,h);
            if ~exist(fname_out,'file'), runflag = 1; return; end;
            fnames_mgz{m,h,k,pind} = fname_out;
          end;
        end;
      end;
    end;
  end;
  fnames_mgz = squeeze(fnames_mgz);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_mgz,fnames_log] = concat_files(parms)
  fnames_mgz = cell(parms.nmeas,parms.nhemi,parms.nsmooth,parms.npdist+parms.contrast_flag);
  fnames_log = cell(parms.nmeas,parms.nhemi,parms.nsmooth,parms.npdist+parms.contrast_flag);
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    for c=1:length(parms.contflag_list)
      contflag = parms.contflag_list(c);
      if contflag
        projdist_list = parms.projdist_list(2);
        npdist = 1;
      else
        projdist_list = parms.projdist_list;
        npdist = parms.npdist;
      end;
      for p=1:npdist
        if contflag
          pind = parms.npdist + 1;
        else
          pind = p;
        end;
        for h=1:parms.nhemi
          for k=1:parms.nsmooth
            fname_out = set_fname_out(parms,meas,contflag,p,k,h);
            fname_log = set_fname_log(parms,meas,contflag,p,k,h);
            if ~exist(fname_out,'file') ||...
               ~exist(fname_log,'file') || parms.forceflag
              if parms.verbose==2
                fprintf('%s: concatenating data files into %s...\n',...
                  mfilename,fname_out);
              end;
              fnamelist = cell(parms.nsubs,1);
              for s=1:parms.nsubs
                fnamelist{s} = set_fname_data(parms,s,meas,contflag,p,k,h);
              end;
              args = mmil_parms2args(parms,parms.concat_tags);
              mmil_concat(fnamelist,fname_out,args{:});
              write_log(fname_log,fnamelist);
            end;
            fnames_mgz{m,h,k,pind} = fname_out;
            fnames_log{m,h,k,pind} = fname_log;
          end;
        end;
      end;
    end;
  end;
  fnames_mgz = squeeze(fnames_mgz);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_log(fname_log,fnamelist)
  fid = fopen(fname_log,'wt');
  if fid==-1
    error('failed to open %s for writing',fname_log);
  end;
  fprintf(fid,'%s\n',fnamelist{:});
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_csv = write_csv(parms)
  fname_csv = set_fname_csv(parms);
  if ~exist(fname_csv,'file') || parms.forceflag
    if parms.verbose==2
      fprintf('%s: writing csv file...\n',mfilename);
    end;
    % initialize data cell array
    data = cell(parms.nsubs,0);
    col_labels = cell(1,0);
    % add SubjIDs as first column    
    data = cat(2,data,{parms.StudyInfo.SubjID}');
    col_labels = cat(2,col_labels,'SubjID');
    % add VisitIDs as second column
    data = cat(2,data,{parms.StudyInfo.VisitID}');
    col_labels = cat(2,col_labels,'VisitID');
    % add info from StudyInfo
    for t=1:length(parms.info_tags)
      tag = parms.info_tags{t};
      if isfield(parms.StudyInfo,tag)
        info = {parms.StudyInfo.(tag)}';
        data = cat(2,data,info);
        col_labels = cat(2,col_labels,tag);
      end;
    end;
    % add column labels
    data = cat(1,col_labels,data);
    % write file
    mmil_write_csv(fname_csv,data);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_files(parms)
  if parms.verbose==2
    fprintf('%s: checking files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    for m=1:parms.nmeas
      if ~valid_flags(s), break; end;
      meas = parms.measlist{m};
      for c=1:length(parms.contflag_list)
        if ~valid_flags(s), break; end;
        contflag = parms.contflag_list(c);
        if contflag
          projdist_list = parms.projdist_list(2);
          npdist = 1;
        else
          projdist_list = parms.projdist_list;
          npdist = parms.npdist;
        end;
        for p=1:npdist
          if ~valid_flags(s), break; end;
          for h=1:parms.nhemi
            if ~valid_flags(s), break; end;
            for k=1:parms.nsmooth
              fname = set_fname_data(parms,s,meas,contflag,p,k,h);
              if ~exist(fname,'file')
                if parms.verbose
                  fprintf('%s: WARNING: file %s not found\n',mfilename,fname);
                end;
                valid_flags(s) = 0;
                break;
              end;
            end;
          end;
        end;
      end;
    end;
  end;
  ind_valid = find(valid_flags);
  if length(ind_valid)<parms.nsubs
    parms.StudyInfo = parms.StudyInfo(ind_valid);
    parms.nsubs = length(parms.StudyInfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = set_fstem(parms,s,meas)
  % set input file stem
  fstem = [];
  ContainerPath = sprintf('%s/%s',...
    parms.RootDirs.proc_asl,parms.StudyInfo(s).proc_asl);
  if ~isempty(fstem), fstem = [fstem '_']; end;
  fstem = sprintf('%s/%s/%s',ContainerPath,parms.analdir,fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_fname_data(parms,s,meas,contflag,p,k,h)
  fstem = set_fstem(parms,s,meas);
  if isempty(fstem), return; end;
  if contflag
    fname = sprintf('%s%s_cortsurf_%s_pdist%0.1f',...
      fstem,meas,parms.contrast_outfix,parms.projdist_list(2));
  else
    fname = sprintf('%s%s_cortsurf_pdist%0.1f',...
      fstem,meas,parms.projdist_list(p));
  end;
  fname = [fname '-sphere'];
  if parms.smoothing(k)~=0
    fname = sprintf('%s-sm%d',fname,parms.smoothing(k));
  end;
  fname = [fname '-' parms.hemilist{h} '.' parms.intype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem_out = set_fstem_out(parms,meas,contflag,p,k,h)
  outstem = [parms.outstem '_' meas];
  if contflag
    outstem = sprintf('%s_contrast_pdist%0.1f',outstem,parms.projdist_list(2));
  else
    outstem = sprintf('%s_pdist%0.1f',outstem,parms.projdist_list(p));
  end;
  outstem = sprintf('%s_sm%d',outstem,parms.smoothing(k));
  fstem_out = [outstem '-' parms.hemilist{h}];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_fname_out(parms,meas,contflag,p,k,h)
  fname_out = [set_fstem_out(parms,meas,contflag,p,k,h) '.' parms.outtype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_log = set_fname_log(parms,meas,contflag,p,k,h)
  fname_log = [set_fstem_out(parms,meas,contflag,p,k,h) '.log'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_csv = set_fname_csv(parms)
  fname_csv = sprintf('%s_info.csv',parms.outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

