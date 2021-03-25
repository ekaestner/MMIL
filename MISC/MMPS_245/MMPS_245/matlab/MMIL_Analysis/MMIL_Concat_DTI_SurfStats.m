function [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_DTI_SurfStats(ProjID,varargin)
%function [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_DTI_SurfStats(ProjID,[options])
%
% Usage:
%  [fnames_mgz,fname_csv,fnames_log] = MMIL_Concat_DTI_SurfStats(ProjID,'key1', value1,...);
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
%      proc_dti, fsurf
%    If both RootDirs and StudyInfo are supplied,
%      MMIL_ProjInfo.csv is not required
%     {default = []}
%  'qcflag': use QC flag in StudyInfo to determine whether to exclude subjects
%     {default = 1}
%
% Optional Parameters that determine input diffusion data:
%  'snums_flag': [0|1|2|3] which set of "snums" to use
%     0: use all available scans
%     1: use scan numbers in StudyInfo.DTIScanNums
%     2: use scan numbers in StudyInfo.DTIScanNums2
%     3: use scan numbers in StudyInfo.DTIScanNums3
%    {default = 1}
%  'snum_index': index specifying which scan number of DTIScanNums
%     (or DTIScanNums2) to use in spread sheet (must have run DT
%     calculations separately for each scan)
%     If empty, use DT measures calculated from all DTIScanNums
%     {default = []}
%  'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix: 'corr_resDTI'
%     {default = []}
%  'auto_infix_flag': [0|1] set infix automatically based on typical
%     processing and settings in ProjInfo
%     ignored if infix is not empty
%     {default = 1}
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%  'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%     if 1, multiple b-values are required
%       also, b=0 images are still used for between image scaling
%     {default = 0}
%  'min_bval': minimum b value a scan must have to be included in tensor fit
%     {default = 1}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%  'min_ndirs': minimum number of diffusion directions to be included
%     {default = 6}
%  'full_fstem_flag': [0|1] full stem included in DT measures analysis output
%     otherwise, use only meas name (e.g. 'FA' instead of 'DTI_scans_1_2_...FA')
%     if 0, options above (snums_flag, snum_index, infix, etc.) are ignored
%     {default = 0}
%  'projdist_list': vector of mm distances along normal vector to paint from
%     negative = white matter, positive = gray matter
%     if exactly two, will be used to calculate gray-white contrast
%     {default = [-1,1]}
%  'gwnorm_flag': [0|1] when calculating gray-white contrast, normalize by mean
%     {default = 0}
%  'gwcsurf_flag': [0|1] get results of gwcsurf analysis
% 	  if true, ignore projdist_list and gwnorm_flag
%     {default = 0}
%  'gwcsurf_outfix': suffix added to output files for gwcsurf analysis
%     {default = 'gwcsurf'}
%  'smoothing': surface smoothing steps (on sphere)
%     slope of FWHM vs. sqrt(N) is ~1.13 for fsaverage
%     (FWHM = full-width-half-max smoothing kernel
%         N = number of smoothing steps)
%    {default = 0}
%  'DT_analdir': name of DTI analysis folder relative to ContainerPath
%     {default = 'DTanalysis'}
%  'DT_outfix': string attached to DT calculation and analysis file names
%     {default = []}
%
% Optional Parameters:
%  'outdir': output directory
%    full path or relative to /home/{user}/MetaData/{ProjID}
%    {default = 'DTI_SurfStats'}
%  'outstem': output file stem
%    relative to outdir unless full path given
%    {default = 'DTI'}
%  'motion_flag': [0|1] compile mean_motion for each subject
%    {default = 0}
%  'max_motion': exclude subjects with mean relative motion greater than this
%    set to Inf to include all subjects
%    {default = Inf}
%  'meanflag': [0|1] calculate mean instead of concatenating
%    {default = 0}
%  'options': option string to use any of the mri_concat command line options
%    {default = []}
%  'verbose': [0|1|2] display status messages
%    0: no messages except errors
%    1: no messages except WARNING
%    2: frequent status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
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
% Created:  07/15/13 by Don Hagler
% Last Mod: 01/30/16 by Don Hagler
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
    'measlist',{'FA','MD','LD','TD','T2w','T1w'},[],...
    'scalefacts',[],[],...
    'inputlist',[],[],...
... % specify DTI data used for tensor calculations and analysis
    'snums_flag',1,[0:3],...
    'snum_index',[],[],...
    'infix',[],[],...
    'auto_infix_flag',true,[false true],...
    'revflag',0,[0:2],...
    'nob0_flag',false,[false true],...
    'min_bval',1,[],...
    'flex_flag',false,[false true],...
    'min_nb0',1,[],...
    'min_ndirs',6,[],...
    'full_fstem_flag',false,[false true],...
    'projdist_list',[-1,1],[-5,5],...
    'gwnorm_flag',true,[false true],...
	  'gwcsurf_flag',false,[false true],...
	  'gwcsurf_outfix','gwcsurf',[],...
    'smoothing',0,[0,Inf],...
    'DT_analdir','DTanalysis',[],...
    'DT_outfix',[],[],...
...
    'outdir','DTI_SurfStats',[],...
    'outstem','DTI',[],...
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
    'ProjInfo_tags',{'min_bval','flex_flag','revflag','snums_flag',...
                     'resample_flag','regT1flag'},[],...
    'concat_tags',{'meanflag','options','forceflag'},[],...
    'scaninfo_tags',{'snums','revflag','min_nb0','min_ndirs',...
                 'min_bval','flex_flag'},[],...
    'fstem_tags',{'snums','infix','revflag','min_bval','flex_flag',...
                  'min_ndirs','min_nb0','nob0_flag','outdir','outfix'},[],...
... % hidden
		'gwcsurf_infix_list',{'wm','gm','gwc'},[],...
    'DT_outdir','DTcalc',[],...
    'RSI_outdir','RSIcalc',[],...
    'RSI_analdir','RSIanalysis',[],...
    'RSI_outfix',[],[],...
    'DT_measlist',{'FA','MD','LD','TD','b0','b0N','T2w'},[],...
    'RSI_measlist',{'T','F0','N0','F2','N2','F4','N4','FD','ND',...
                    'FT','NT','Ir','Ih','If','NIr','NIh','NIf','V0','AU'},[],...
    'required_containers',{'proc_dti','fsurf'},[],...
    'QC_raw',true,[false true],... % only applies if manual raw QC exists
    'QC_DTI',true,[false true],... % only applies if manual DTIQC.csv file exists
    'QC_recon',true,[false true],...
    'contrast_outfix','contrast',[],...
    'resample_flag',true,[false true],...
    'regT1flag',1,[0:2],...
    'motion_flag',false,[false true],...
    'max_motion',Inf,[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  args = MMIL_Args(parms,'MMIL_Check_ProjID');
  [ProjInfo,StudyInfo,RootDirs] = MMIL_Check_ProjID(ProjID,args{:});
  if isempty(StudyInfo), error('empty StudyInfo'); end;

  % for arg names present in both varargin and ProjInfo,
  %   the options values will appear in merged_args
  %   i.e. command-line takes precedence, then ProjInfo, then defaults
  if ~isempty(ProjInfo)
    % get all 'DTI_' parms, strip 'DTI_'
    ProjInfo_args = MMIL_Args(ProjInfo,'DTI');
    % convert from args cell array to parms struct
    ProjInfo_parms = mmil_args2parms(ProjInfo_args,[],0);
    % keep only those in ProjInfo_tags
    ProjInfo_args = mmil_parms2args(ProjInfo_parms,parms.ProjInfo_tags);
    % merge arguments, giving varargin precedence over ProjInfo_args
    merged_args = mmil_merge_args(options,ProjInfo_args);
    % check that parameters fit allowed range, use defaults if not supplied
    parms = mmil_args2parms(merged_args,parms_filter);
  end;

  parms.StudyInfo = StudyInfo;
  parms.RootDirs = RootDirs;
  
  % set infix
  if parms.auto_infix_flag && isempty(parms.infix)
    parms.infix = 'corr';
    if parms.resample_flag
      if parms.regT1flag==2
        parms.infix = [parms.infix '_regT1'];
      else
        parms.infix = [parms.infix '_resDTI'];
      end;
    end;
  end;

  % set DTI-specific fields in StudyInfo
  [parms.StudyInfo,parms.SIflags] = DTI_MMIL_Check_StudyInfo(...
    parms.StudyInfo,parms.RootDirs,...
    'snums_flag',parms.snums_flag,... % DTIScanNums field will get changed
    'qcflag',parms.qcflag,...
    'DTI_flag',1);
  if isempty(parms.StudyInfo), error('empty StudyInfo after check for DTI'); end;
  parms.nsubs = length(parms.StudyInfo);
  
  parms.nhemi = length(parms.hemilist);
	if ~parms.gwcsurf_flag
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
	else
		parms.ngwc = length(parms.gwcsurf_infix_list);
	end;
  parms.nsmooth = length(parms.smoothing);
  parms.nmeas = length(parms.measlist);

  if parms.motion_flag
    parms = check_motion_files(parms);
    if parms.nsubs==0
      fprintf('%s: ERROR: no valid subjects\n',mfilename);
      return;
    end;
  end;
  
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
  runflag = 0;
  fname_csv = set_fname_csv(parms);
  if ~exist(fname_csv,'file'), runflag = 1; return; end;
	if parms.gwcsurf_flag
		[runflag,fnames_mgz] = check_gwcsurf_output_files(parms);
	else
		[runflag,fnames_mgz] = check_pdist_output_files(parms);
	end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runflag,fnames_mgz] = check_pdist_output_files(parms)
	runflag = 0;
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
            fname_out = set_pdist_fname_out(parms,meas,contflag,p,k,h);
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

function [runflag,fnames_mgz] = check_gwcsurf_output_files(parms)
  runflag = 0;
  fnames_mgz = cell(parms.nmeas,parms.nhemi,parms.nsmooth,parms.ngwc);
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    for c=1:parms.ngwc
      for h=1:parms.nhemi
        for k=1:parms.nsmooth
          fname_out = set_gwcsurf_fname_out(parms,meas,c,k,h);
          if ~exist(fname_out,'file'), runflag = 1; return; end;
          fnames_mgz{m,h,k,c} = fname_out;
        end;
      end;
    end;
  end;
  fnames_mgz = squeeze(fnames_mgz);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_mgz,fnames_log] = concat_files(parms)
	if parms.gwcsurf_flag
		[fnames_mgz,fnames_log] = concat_gwcsurf_files(parms);
	else
		[fnames_mgz,fnames_log] = concat_pdist_files(parms);
	end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_mgz,fnames_log] = concat_pdist_files(parms)
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
            fname_out = set_pdist_fname_out(parms,meas,contflag,p,k,h);
            fname_log = set_pdist_fname_log(parms,meas,contflag,p,k,h);
            if ~exist(fname_out,'file') ||...
               ~exist(fname_log,'file') || parms.forceflag
              if parms.verbose==2
                fprintf('%s: concatenating data files into %s...\n',...
                  mfilename,fname_out);
              end;
              fnamelist = cell(parms.nsubs,1);
              for s=1:parms.nsubs
                fnamelist{s} = set_pdist_fname_data(parms,s,meas,contflag,p,k,h);
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

function [fnames_mgz,fnames_log] = concat_gwcsurf_files(parms)
  fnames_mgz = cell(parms.nmeas,parms.nhemi,parms.nsmooth,parms.ngwc);
  fnames_log = cell(parms.nmeas,parms.nhemi,parms.nsmooth,parms.ngwc);
  for m=1:parms.nmeas
    meas = parms.measlist{m};
    for c=1:parms.ngwc
      for h=1:parms.nhemi
        for k=1:parms.nsmooth
          fname_out = set_gwcsurf_fname_out(parms,meas,c,k,h);
          fname_log = set_gwcsurf_fname_log(parms,meas,c,k,h);
          if ~exist(fname_out,'file') ||...
             ~exist(fname_log,'file') || parms.forceflag
            if parms.verbose==2
              fprintf('%s: concatenating data files into %s...\n',...
                mfilename,fname_out);
            end;
            fnamelist = cell(parms.nsubs,1);
            for s=1:parms.nsubs
              fnamelist{s} = set_gwcsurf_fname_data(parms,s,meas,c,k,h);
            end;
            args = mmil_parms2args(parms,parms.concat_tags);
            mmil_concat(fnamelist,fname_out,args{:});
            write_log(fname_log,fnamelist);
          end;
          fnames_mgz{m,h,k,c} = fname_out;
          fnames_log{m,h,k,c} = fname_log;
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
    if parms.motion_flag
      % add motion
      motion = compile_motion(parms);
      data = cat(2,data,num2cell(motion));
      col_labels = cat(2,col_labels,'motion');
    end;
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
	if parms.gwcsurf_flag
		parms = check_gwcsurf_files(parms);
	else
		parms = check_pdist_files(parms);
	end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_pdist_files(parms)
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
              fname = set_pdist_fname_data(parms,s,meas,contflag,p,k,h);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_gwcsurf_files(parms)
  if parms.verbose==2
    fprintf('%s: checking files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    for m=1:parms.nmeas
      if ~valid_flags(s), break; end;
      meas = parms.measlist{m};
      for c=1:parms.ngwc
        if ~valid_flags(s), break; end;
        for h=1:parms.nhemi
          if ~valid_flags(s), break; end;
          for k=1:parms.nsmooth
            fname = set_gwcsurf_fname_data(parms,s,meas,c,k,h);
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
  ind_valid = find(valid_flags);
  if length(ind_valid)<parms.nsubs
    parms.StudyInfo = parms.StudyInfo(ind_valid);
    parms.nsubs = length(parms.StudyInfo);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_motion_files(parms,results)
  if parms.verbose==2
    fprintf('%s: checking motion files...\n',mfilename);
  end;
  valid_flags = ones(parms.nsubs,1);
  for s=1:parms.nsubs
    fnames = set_motion_fnames(parms,s);
		if isempty(fnames)
      valid_flags(s) = 0;
    elseif parms.max_motion<Inf
      motion = load_subj_motion(parms,s);
      if motion > parms.max_motion
        fprintf('%s: WARNING: excluding visit %s with excessive motion...\n',...
          mfilename,parms.StudyInfo(s).VisitID);
        valid_flags(s) = 0;
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
    parms.RootDirs.proc_dti,parms.StudyInfo(s).proc_dti);
  parms.snums = parms.StudyInfo(s).DTIScanNums;
  args = mmil_parms2args(parms,parms.info_tags);
  [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
  if errcode || isempty(ScanInfo)
    if parms.verbose
    	fprintf('%s: WARNING: no DTI scans for %s\n',mfilename,ContainerPath);
		end;
    return;
  end;
  switch meas
    case parms.DT_measlist
      if parms.full_fstem_flag
        parms.outdir = parms.DT_outdir;
        parms.outfix = parms.DT_outfix;
        args = mmil_parms2args(parms,parms.fstem_tags);
        fstem = DTI_MMIL_Set_DT_fstem(ContainerPath,args{:});
        if isempty(fstem), return; end;
        [tmp_path,fstem,tmp_ext] = fileparts(fstem);
      else
        fstem = parms.DT_outfix;
      end;
      indir = parms.DT_analdir;
    case parms.RSI_measlist
      if parms.full_fstem_flag
        parms.outdir = parms.RSI_outdir;
        parms.outfix = parms.RSI_outfix;
        args = mmil_parms2args(parms,parms.fstem_tags);
        fstem = DTI_MMIL_Set_RSI_fstem(ContainerPath,args{:});
        if isempty(fstem), return; end;
        [tmp_path,fstem,tmp_ext] = fileparts(fstem);
      else
        fstem = parms.RSI_outfix;
      end;
      indir = parms.RSI_analdir;
    otherwise
      fstem = [];
      indir = parms.DT_analdir;
  end;
  if ~isempty(fstem), fstem = [fstem '_']; end;
  fstem = sprintf('%s/%s/%s',ContainerPath,indir,fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_pdist_fname_data(parms,s,meas,contflag,p,k,h)
  fstem = set_fstem(parms,s,meas);
  if isempty(fstem), return; end;
  if contflag
    fname = sprintf('%s%s_%s_pdist%0.1f',...
      fstem,meas,parms.contrast_outfix,parms.projdist_list(2));
  else
    fname = sprintf('%s%s_pdist%0.1f',...
      fstem,meas,parms.projdist_list(p));
  end;
  fname = [fname '-sphere'];
  if parms.smoothing(k)~=0
    fname = sprintf('%s-sm%d',fname,parms.smoothing(k));
  end;
  fname = [fname '-' parms.hemilist{h} '.' parms.intype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = set_gwcsurf_fname_data(parms,s,meas,c,k,h)
  fstem = set_fstem(parms,s,meas);
  if isempty(fstem), return; end;
  fname = sprintf('%s%s_%s_%s',...
    fstem,meas,parms.gwcsurf_outfix,parms.gwcsurf_infix_list{c});
  fname = [fname '-sphere'];
  if parms.smoothing(k)~=0
    fname = sprintf('%s-sm%d',fname,parms.smoothing(k));
  end;
  fname = [fname '-' parms.hemilist{h} '.' parms.intype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem_out = set_pdist_fstem_out(parms,meas,contflag,p,k,h)
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

function fname_out = set_pdist_fname_out(parms,meas,contflag,p,k,h)
  fname_out = [set_pdist_fstem_out(parms,meas,contflag,p,k,h) '.' parms.outtype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_log = set_pdist_fname_log(parms,meas,contflag,p,k,h)
  fname_log = [set_pdist_fstem_out(parms,meas,contflag,p,k,h) '.log'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem_out = set_gwcsurf_fstem_out(parms,meas,c,k,h)
  outstem = [parms.outstem '_' meas];
  outstem = sprintf('%s_%s_%s',...
		outstem,parms.gwcsurf_outfix,parms.gwcsurf_infix_list{c});
  outstem = sprintf('%s_sm%d',outstem,parms.smoothing(k));
  fstem_out = [outstem '-' parms.hemilist{h}];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = set_gwcsurf_fname_out(parms,meas,c,k,h)
  fname_out = [set_gwcsurf_fstem_out(parms,meas,c,k,h) '.' parms.outtype];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_log = set_gwcsurf_fname_log(parms,meas,c,k,h)
  fname_log = [set_gwcsurf_fstem_out(parms,meas,c,k,h) '.log'];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_csv = set_fname_csv(parms)
  fname_csv = sprintf('%s_info.csv',parms.outstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnames = set_motion_fnames(parms,s)
	fname = [];
  VisitID = parms.StudyInfo(s).VisitID;
  ContainerPath = sprintf('%s/%s',...
    parms.RootDirs.proc_dti,parms.StudyInfo(s).proc_dti);
  parms.snums = parms.StudyInfo(s).DTIScanNums;
  args = mmil_parms2args(parms,parms.info_tags);
  [ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
  if errcode || isempty(ScanInfo)
    if parms.verbose
	    fprintf('%s: WARNING: no DTI scans for %s\n',mfilename,ContainerPath);
		end;
    return;
  end;
  nscans = length(SessInfo.snums_DT);
	fnames = {};
	for i=1:nscans
    j = SessInfo.snums_DT(i);
    fname_qmat = sprintf('%s/DTI%d',ContainerPath,j);
    if SessInfo.revflag
      fname_qmat = [fname_qmat '_rev'];
    end;
    if ~isempty(parms.infix)
      fname_qmat = [fname_qmat '_' parms.infix];
    end;
    fname_qmat = [fname_qmat '_qmat.mat'];
	  if ~exist(fname_qmat,'file')
  	  if parms.verbose
    	  fprintf('%s: WARNING: file %s not found\n',mfilename,fname_qmat);
      end;
			continue;
		end;
		fnames{end+1} = fname_qmat;
	end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion = load_subj_motion(parms,s)
  motion = nan;
  mean_motion = 0; mean_trans = 0; mean_rot = 0;
	fnames = set_motion_fnames(parms,s);
	nscans = length(fnames);
  for i=1:nscans 
		fname_qmat = fnames{i};
    tmp = load(fname_qmat);
    mean_motion = mean_motion + tmp.mean_motion;
    mean_trans = mean_trans + tmp.mean_trans;
    mean_rot = mean_rot + tmp.mean_rot;
  end;
  % average across scans
  mean_motion = mean_motion / nscans;
  mean_trans = mean_trans / nscans;
  mean_rot = mean_rot / nscans;
	% return value
  motion = mean_motion;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion = compile_motion(parms)
  motion = nan(parms.nsubs,1);
  for s=1:parms.nsubs
    motion(s) = load_subj_motion(parms,s);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

