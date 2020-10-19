function datadir = abcd_import(fname_json,varargin)
%function datadir = abcd_import(fname_json,[options])
%
% Required Input:
%   fname_json: full path of input json file
%
% Optional Parameters: ('key', value pairs)
%   'unpackdir': output directory that will be created
%     and filled with contents of tar or tgz corresponding to json
%     {default = [pwd '/unpack'];
%   'outdir': output directory that will contain unpacked data
%     {default = [pwd '/output'];
%   'fname_sites': spreadsheet containing site names and IDs
%     {'ABCD_Sites.csv'}
%   'mb_recon_flag': [-1|0|1|2] whether to perform multi-band recon
%     -1: unpack only, this will work with fast-track data
%     0: no multi-band recon
%     1: unpack and perform multi-band recon if required
%     2: do not unpack unless multi-band recon is required
%     {default = 1}
%   'mb_tmpdir': temporary directory that will contain untar'd P-file
%     {default = '/tmp'}
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
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%   'force_unpack_err_flag': [0|1] overwrite existing output if previous unpack errors
%     {default = 0}
%   'force_recon_err_flag': [0|1] overwrite existing output if previous recon errors
%     {default = 0}
%   'force_dicoms_err_flag': [0|1] overwrite existing output if previous dicom errors
%     for multi-band recon only
%     {default = 0}
%   'force_rename_err_flag': [0|1] overwrite existing output if previous rename errors
%     {default = 0}
%
% Output:
%   datadir: data directory containing imported dicoms
%
% Created:  07/16/16 by Don Hagler
% Prev Mod: 08/14/17 by Don Hagler
% Last Mod: 10/03/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datadir = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(fname_json,varargin);

if parms.jinfo.pfile_flag
  fprintf('%s: WARNING: skipping import because %s is a P file json file\n',mfilename,parms.fname_json);
  return;
end;

% create outdir and unpackdir
mmil_mkdir(parms.outdir);
mmil_mkdir(parms.unpackdir);

% copy json and md5sum files
copy_json_file(parms.fname_json,parms);
copy_md5sum_file(parms.fname_json,parms);
if ~isempty(parms.fname_pfile_json)
  copy_json_file(parms.fname_pfile_json,parms);
  copy_md5sum_file(parms.fname_pfile_json,parms);
end;

% determine whether mb recon is required
if parms.mb_recon_flag==2 && ~parms.jinfo.recon_required_flag
  fprintf('%s: WARNING: skipping import because recon not required and mb_recon_flag = 2\n',mfilename);
  return;
elseif parms.mb_recon_flag==0 && parms.jinfo.recon_required_flag
  fprintf('%s: WARNING: skipping import because recon required and mb_recon_flag = 0\n',mfilename);
  return;
end;

% determine whether mb recon already completed
if parms.mb_dicoms_flag==2 && parms.jinfo.recon_required_flag
  fname_out = sprintf('%s/%s.reconned',parms.unpackdir,parms.fstem_json);
  if ~exist(fname_out,'file')
    fprintf('%s: WARNING: skipping import because recon not complete and mb_dicoms_flag = 2\n',mfilename);
    return;
  end;
end;

% unpack tgz file
errcode = unpack_tar_file(parms);
if errcode, return; end;

% multi-band recon (if required)
if parms.jinfo.recon_required_flag
  % perform multi-band recon
  errcode = recon_data(parms);
  if errcode, return; end;
  % write dicom files
  if parms.mb_dicoms_flag
    errcode = write_recon_dicoms(parms);
    if errcode, return; end;
  else
    return;
  end;
end;

% link dicom data to destination directory
[datadir,errcode] = rename_data(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_json,options)
  parms = mmil_args2parms(options,{...
    'fname_json',fname_json,[],...
  ...
    'indir',[],[],...
    'unpackdir',[pwd '/unpack'],[],...
    'outdir',[pwd '/output'],[],...
    'fname_sites','ABCD_Sites.csv',[],...
    'mb_recon_flag',1,[-1:2],...
    'mb_tmpdir','/tmp',[],...
    'mb_dicoms_flag',2,[0:2],...
    'mb_cleanup_flag',2,[0:2],...
    'separate_series_flag',false,[false true],...
    'forceflag',false,[false true],...
    'force_unpack_err_flag',false,[false true],...
    'force_recon_err_flag',false,[false true],...
    'force_dicoms_err_flag',false,[false true],...
    'force_rename_err_flag',false,[false true],...
  ...
    'md5sum_ext','.md5sum',[],...
    'fstem_dcm','metadata',[],...
    'mb_recon_outstem','mb_recon',[],...
  });

  % check json file, find corresponding tar file
  [parms.jinfo,parms.fstem_json,parms.fname_tar,parms.fname_pfile_json] = ...
    abcd_check_json(parms.fname_json,parms.indir);
  if isempty(parms.jinfo)
    error('check of json file %s failed\n',parms.fname_json);
  end;
  % replace spaces with underscores
  parms.fstem_json = regexprep(parms.fstem_json,' ','_');
  
  % set name of output directory for unpacking
  parms.unpackdir = sprintf('%s/%s',parms.unpackdir,parms.fstem_json);

  % set name of output file name for multi-band recon
  if parms.jinfo.recon_required_flag
    parms.mb_recon_outdir = sprintf('%s/mb_recon',parms.unpackdir);
  end;

  % open site name to id look up table
  if ~exist(parms.fname_sites,'file')
    error('file %s not found',parms.fname_sites);
  end;
  site_info = mmil_csv2struct(parms.fname_sites);
  idk = find(~cellfun(@isempty,{site_info.DSN}));
  site_info = site_info(idk);
  parms.site_names = {site_info.Abbreviation};
  parms.site_nums = {site_info.Site_No};
  parms.site_dsn = {site_info.DSN};
  for i=1:length(parms.site_dsn)
    if isnumeric(parms.site_dsn{i})
      parms.site_dsn{i} = num2str(parms.site_dsn{i});
    end;
  end;

  % set SequenceName
  if isfield(parms.jinfo,'SequenceName')
    parms.SequenceName = parms.jinfo.SequenceName;
  elseif isfield(parms.jinfo,'PulseSequenceName')
    parms.SequenceName = parms.jinfo.PulseSequenceName;
  else
    parms.SequenceName = 'unknown';
  end;
  
  % set mb_tmpdir to be a subdir of input parameter
  if strcmp(parms.mb_tmpdir,'/scratch') && ~exist(parms.mb_tmpdir,'dir')
    parms.mb_tmpdir = '/tmp';
  end;
  parms.mb_tmpdir = sprintf('%s/%s',parms.mb_tmpdir,parms.fstem_json);

  % if recon required, check if mb_tmpdir exists
  if parms.jinfo.recon_required_flag && ~parms.forceflag
    fname_out = sprintf('%s/%s.reconned',parms.unpackdir,parms.fstem_json);
    fname_err = sprintf('%s/%s.recon_error',parms.unpackdir,parms.fstem_json);
    if ~exist(fname_out,'file') &&...
      (~exist(fname_err,'file') || parms.force_recon_err_flag)
      if ~exist(parms.mb_tmpdir,'dir')
        fname_out = sprintf('%s/%s.unpacked',parms.unpackdir,parms.fstem_json);
        if exist(fname_out,'file')
          % force unpack, since output does not exist
          fprintf('%s: mb_tmpdir %s is missing, so setting forceflag = 1\n',...
            mfilename,parms.mb_tmpdir);
          parms.forceflag = 1;
        end;
      end;
    end;    
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_json_file(fname_json,parms)
  [tmp,fstem,fext] = fileparts(fname_json);
  % replace spaces with underscores
  fstem = regexprep(fstem,' ','_');
  fname_json_out = sprintf('%s/%s.json',parms.unpackdir,fstem);
  if strcmp(fname_json,fname_json_out), return; end;
  if ~exist(fname_json_out,'file') || parms.forceflag
    cmd = sprintf('cp -p %s %s',regexprep(fname_json,' ','\\ '),fname_json_out);
    fprintf('%s: copying json file %s...\n',mfilename,fname_json);
    [s,r] = unix(cmd);
    if s
      error('failed to copy json file:\n%s\n%s\n',cmd,r);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_md5sum_file(fname_json,parms)
  % find corresponding md5sum file
  [fpath,fstem,fext] = fileparts(fname_json);
  fname_md5sum_in = sprintf('%s/%s%s',fpath,fstem,parms.md5sum_ext);
  if ~exist(fname_md5sum_in)
    fprintf('%s: WARNING: md5sum file %s not found\n',...
      mfilename,fname_md5sum_in);
    return;
  end;
  % replace spaces with underscores
  fstem = regexprep(fstem,' ','_');
  fname_md5sum_out = sprintf('%s/%s%s',parms.unpackdir,fstem,parms.md5sum_ext);
  if ~exist(fname_md5sum_out,'file') || parms.forceflag
    cmd = sprintf('cp -p %s %s',...
      regexprep(fname_md5sum_in,' ','\\ '),fname_md5sum_out);
    fprintf('%s: copying md5sum file %s to %s...\n',...
      mfilename,fname_md5sum_in,fname_md5sum_out);
    [s,r] = unix(cmd);
    if s
      error('failed to copy md5sum file:\n%s\n%s\n',cmd,r);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = unpack_tar_file(parms)
  errcode = 0;
  % check if has been unpacked before
  fname_out = sprintf('%s/%s.unpacked',parms.unpackdir,parms.fstem_json);
  fname_err = sprintf('%s/%s.unpack_error',parms.unpackdir,parms.fstem_json);
  if parms.force_unpack_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing unpack despite previous errors for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
  end;
  if (~exist(fname_out,'file') && ~exist(fname_err,'file')) || parms.forceflag
    if parms.jinfo.recon_required_flag
      unpackdir = parms.mb_tmpdir;
      mmil_mkdir(unpackdir);
    else
      unpackdir = parms.unpackdir;
    end;
    [tmp,fstem,fext] = fileparts(parms.fname_tar);
    switch fext
      case '.tar'
        cmd_base = 'tar -xvf';
      case {'.tgz' '.tar.gz','.gz'}
        cmd_base = 'tar -zxvf';
    end;
    cmd = sprintf('%s %s --directory %s',...
      cmd_base,regexprep(parms.fname_tar,' ','\\ '),unpackdir);
    fprintf('%s: unpacking data from %s...\n',mfilename,parms.fname_tar);
    [s,r] = unix(cmd);
    % check for failure of unix command
    if s
      fprintf('%s: ERROR: failed to unpack:\n%s\n%s\n',mfilename,cmd,r);
      fname_out = fname_err;
      errcode = 1;
    end;    
    % create check file
    fid = fopen(fname_out,'wt');
    if fid<0
      error('failed to open file %s for writing',fname_out);
    end;
    fprintf(fid,'%s\n%s\n',cmd,r);
    fclose(fid);
  elseif exist(fname_out,'file')
    fprintf('%s: skipping unpack because already complete\n',mfilename);
  elseif exist(fname_err,'file')
    fprintf('%s: skipping unpack because of previous errors\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = recon_data(parms)
  errcode = 0; errmsg = [];
  fname_err = sprintf('%s/%s.unpack_error',parms.unpackdir,parms.fstem_json);
  if parms.force_unpack_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing recon after previous unpack error for %s...\n',mfilename,parms.fname_tar);
    parms.forceflag = 1;
  end;
  fname_out = sprintf('%s/%s.reconned',parms.unpackdir,parms.fstem_json);
  fname_err = sprintf('%s/%s.recon_error',parms.unpackdir,parms.fstem_json);
  if parms.force_recon_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing recon despite previous errors for %s...\n',mfilename,parms.fname_tar);
    parms.forceflag = 1;
  end;
  if (~exist(fname_out,'file') && ~exist(fname_err,'file')) || parms.forceflag
    % locate directory containing P file
    [parms.pdir,errcode,errmsg] = get_pdir(parms);
    if errcode, fprintf('%s: ERROR: %s\n',mfilename,errmsg); end;
    % identify P file (k-space raw data)
    if ~errcode  
      [fname_pfile,errcode,errmsg] = get_pfile(parms);
      if errcode, fprintf('%s: ERROR: %s\n',mfilename,errmsg); end;
    end;
    % copy single dicom file
    if ~errcode
      [fname_dicom,errcode,errmsg] = get_dicom_file(parms);
      if errcode
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      else
        copy_dicom_file(fname_dicom,parms);
      end;
    end;
    % copy physio data (fMRI only)
    if ~errcode
      copy_physio_data(parms);
    end;
    % check ref.dat and vrfg.dat files
    if ~errcode
      [errcode,errmsg] = check_extra_files(parms);
      if errcode
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      end;
    end;
    % run recon
    if ~errcode
      fprintf('%s: reconning data from %s...\n',mfilename,fname_pfile);
      tic;
      try
        mb_recon(fname_pfile,parms.mb_recon_outdir,...
                 parms.mb_recon_outstem,parms.forceflag);
      catch me
        fprintf('%s: ERROR: multi-band recon failed:\n%s\n',...
          mfilename,me.message);
        errcode = 1;
        errmsg = me.message;
      end;
      toc;
    end;
    % set check file name
    if errcode
      fname_out = fname_err;
    end;
    % create check file
    fid = fopen(fname_out,'wt');
    if fid<0
      error('failed to open file %s for writing',fname_out);
    end;
    if ~errcode
      fprintf(fid,'%s\n%s\n%s\n',fname_pfile,fname_dicom,parms.mb_recon_outdir);
    elseif ~isempty(errmsg)
      fprintf(fid,'ERROR: %s\n',errmsg);
    else
      fprintf(fid,'ERROR: recon failed\n');
    end;
    fclose(fid);
  elseif exist(fname_out,'file')
    fprintf('%s: skipping recon because already complete\n',mfilename);
  elseif exist(fname_err,'file')
    fprintf('%s: skipping recon because of previous errors\n',mfilename);
    errcode = 1;
  end;
  % remove pfile or entire mb_tmpdir
  switch parms.mb_cleanup_flag
    case 1
      if ~errcode, cleanup_pfile(parms); end;
    case 2
      cleanup_tmpdir(parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_physio_data(parms)
  % skip if series is not GE fMRI
  if ~strcmp(parms.SequenceName,'research_ABCD_muxepi'), return; end;
  % pattern match to physio data in p file dir
  fnames_physio = [];
  snum = parms.jinfo.SeriesNumber;
  fpat = sprintf('%s/*ex*_se%s_*',parms.pdir,snum);
  plist = dir(fpat);
  for i=1:length(plist)
    fnames_phsyio{i} = sprintf('%s/%s',parms.pdir,plist(i).name);
  end;
  if isempty(fnames_physio)
    fprintf('%s: WARNING: no physio files found in %s\n',mfilename,parms.pdir);
  end;
  if ~isempty(fnames_physio)
    outdir = sprintf('%s/physio',parms.unpackdir);
    mmil_mkdir(outdir);
    cmd = sprintf('cp -p %s %s',sprintf('%s ',fnames_physio{:}),outdir);
    [s,r] = unix(cmd);
    if s
      error('cmd %s failed:\n%s',cmd,r);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_dicom_file(fname_dicom,parms)
  cmd = sprintf('cp -p %s %s/%s.dcm',...
    fname_dicom,parms.unpackdir,parms.fstem_dcm);
  [s,r] = unix(cmd);
  if s
    error('cmd %s failed:\n%s',cmd,r);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = write_recon_dicoms(parms)
  errcode = 0; errmsg = [];
  % check for errors
  fname_err = sprintf('%s/%s.unpack_error',parms.unpackdir,parms.fstem_json);
  if parms.force_unpack_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing dicom creation after previous unpack error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
  end;
  fname_err = sprintf('%s/%s.recon_error',parms.unpackdir,parms.fstem_json);
  if parms.force_recon_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing dicom creation after previous recon error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
  end;
  fname_out = sprintf('%s/%s.dicoms',parms.unpackdir,parms.fstem_json);
  fname_err = sprintf('%s/%s.dicoms_error',parms.unpackdir,parms.fstem_json);
  if parms.force_dicoms_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing dicom write after previous error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
  end;
  if (~exist(fname_out,'file') && ~exist(fname_err,'file')) || parms.forceflag
    fprintf('%s: writing dicoms after multi-band recon...\n',mfilename);
    % use copied dicom in unpackdir
    fname_dicom = sprintf('%s/%s.dcm',parms.unpackdir,parms.fstem_dcm);
    % read dicom header
    series_info = dicominfo(fname_dicom);
    % determine outdir
    outdir = sprintf('%s/%s/%s',parms.mb_tmpdir,...
                     series_info.StudyInstanceUID,series_info.SeriesInstanceUID);
    mmil_mkdir(outdir);
    % create json file
    fname_json = sprintf('%s.json',outdir);
    savejson([],series_info,fname_json);
    % write dicom files
    mb_write_dicoms(parms.mb_recon_outdir,fname_dicom,outdir);
    % tar dicoms
    fname_tar = sprintf('%s/dicoms.tar',parms.mb_tmpdir);
    cmd = sprintf('tar --directory %s -cvf %s %s',...
      parms.mb_tmpdir,fname_tar,series_info.StudyInstanceUID);
    fprintf('%s: packaging dicoms into %s...\n',mfilename,fname_tar);
    [s,r] = unix(cmd);
    % check for failure of unix command
    if s
      errmsg = sprintf('failed to tar dicoms:\n%s\n%s',cmd,r);
      errcode = 1;
    end;
    if ~errcode
      % untar dicoms to parms.unpackdir
      cmd = sprintf('tar -xvf %s --directory %s',...
        fname_tar,parms.unpackdir);
      fprintf('%s: unpacking data from %s to %s...\n',...
        mfilename,fname_tar,parms.unpackdir);
      [s,r] = unix(cmd);
      % check for failure of unix command
      if s
        errmsg = sprintf('failed to unpack dicoms:\n%s\n%s',cmd,r);
        errcode = 1;
      end;    
    end;
    % set check file name
    if errcode
      fname_out = fname_err;
    end;
    % create check file
    fid = fopen(fname_out,'wt');
    if fid<0
      error('failed to open file %s for writing',fname_out);
    end;
    if ~errcode
      fprintf(fid,'%s\n%s\n%s/%s/%s\n',...
        parms.mb_recon_outdir,fname_dicom,parms.unpackdir,...
        series_info.StudyInstanceUID,series_info.SeriesInstanceUID);
    elseif ~isempty(errmsg)
      fprintf(fid,'ERROR: %s\n',errmsg);
    else
      fprintf(fid,'ERROR: dicom write failed\n');
    end;
    fclose(fid);
  elseif exist(fname_out,'file')
    fprintf('%s: skipping dicom write because already complete\n',mfilename);
  elseif exist(fname_err,'file')
    fprintf('%s: skipping dicom write because of previous errors\n',mfilename);
    errcode = 1;
  end;
  if parms.mb_cleanup_flag==2
    cleanup_tmpdir(parms);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pdir,errcode,errmsg] = get_pdir(parms)
  pdir = [];
  errcode = 0; errmsg = [];
  snum = parms.jinfo.SeriesNumber;
  fpat = sprintf('%s/*ex*_se%s_*',parms.mb_tmpdir,snum);
  dlist = dir(fpat);
  if ~isempty(dlist)
    dflags = [dlist.isdir];
    dlist = dlist(find(dflags));
  end;
  if isempty(dlist)
    errcode = 1;
    errmsg = sprintf('no P file directory in %s',parms.mb_tmpdir);
  elseif length(dlist)>1
    errcode = 1;
    errmsg = sprintf('multiple P file directories in %s',parms.mb_tmpdir);
  else
    pdir = sprintf('%s/%s',parms.mb_tmpdir,dlist.name);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_pfile,errcode,errmsg] = get_pfile(parms)
  fname_pfile = [];
  errcode = 0; errmsg = [];
  plist = dir(sprintf('%s/P*',parms.pdir));
  if ~isempty(plist)
    m = regexp({plist.name},'^P\d{5}\.\d$');
    idp = find(~cellfun(@isempty,m));
    plist = plist(idp);
  end;
  if ~isempty(plist)
    if length(plist)>1
      % use biggest P file (skip aborted scans)
      psz = [plist.bytes];
      [~,idp] = max(psz);
      plist = plist(idp);
    end;
    fname_pfile = sprintf('%s/%s',parms.pdir,plist.name);
  end;
  if isempty(fname_pfile)
    errcode = 1;
    errmsg = sprintf('no P files in %s\n',parms.pdir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_dicom,errcode,errmsg] = get_dicom_file(parms)
  fname_dicom = [];
  errcode = 0; errmsg = [];
  % look for dicom file with specific file name pattern
  dlist = dir(sprintf('%s/i*MRDC*',parms.pdir));
  if isempty(dlist)
    dlist = dir(sprintf('%s/*.dcm',parms.pdir));
  end;
  % look for any dicom file
  if isempty(dlist)
    dlist = dir(sprintf('%s/*',parms.pdir));
    ind = [];
    for i=1:length(dlist)
      fname_tmp = sprintf('%s/%s',parms.pdir,dlist(i).name);
      if mmil_isdicomfile(fname_tmp)
        ind = [ind,i];
      end;
    end;
    dlist = dlist(ind);
  end;
  if ~isempty(dlist)
    fname_dicom = sprintf('%s/%s',parms.pdir,dlist(1).name);
  else
    % get dicom from scanner dicom tgz
    [fpath,fstem] = fileparts(parms.fname_json);
    fpath = fileparts(parms.fname_tar);
    fname_tar = sprintf('%s/%s.tgz',fpath,fstem);
    if ~exist(fname_tar,'file')
      errcode = 1;
      errmsg = sprintf('no dicoms in %s and %s not found\n',...
        parms.pdir,fname_tar);
    else
      cmd = sprintf('tar -zxvf %s --directory %s',...
        regexprep(fname_tar,' ','\\ '),parms.mb_tmpdir);
      fprintf('%s: unpacking dicoms from %s...\n',mfilename,fname_tar);
      [s,r] = unix(cmd);
      % check for failure of unix command
      if s
        errcode = 1;
        errmsg = sprintf('no dicoms in %s and failed to unpack %s:\n%s\n%s\n',...
          parms.pdir,fname_tar,cmd,r);
      else
        % get dicom from unpacked directory
        sdir = sprintf('%s/%s/%s',parms.mb_tmpdir,...
          parms.jinfo.StudyInstanceUID,parms.jinfo.SeriesInstanceUID);
        if ~exist(sdir,'dir')
          errcode = 1;
          errmsg = sprintf('no dicoms in %s and %s not found after unpacking %s\n',...
            parms.pdir,sdir,fname_tar);
        else
          % find dicom
          dlist = dir(sprintf('%s/*',sdir));
          for i=1:length(dlist)
            fname_tmp = sprintf('%s/%s',sdir,dlist(i).name);
            if mmil_isdicomfile(fname_tmp)
              fname_dicom = fname_tmp;
              break;
            end;
          end;
          if isempty(fname_dicom)
            errcode = 1;
            errmsg = sprintf('no dicoms in %s or in %s not found after unpacking %s\n',...
              parms.pdir,sdir,fname_tar);
          end;
        end;
      end;
    end;
  end;
  if ~errcode && isempty(fname_dicom)
    errcode = 1;
    errmsg = sprintf('failed to find dicom in %s\n',...
      parms.pdir);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = check_extra_files(parms)
  errcode = 0; errmsg = [];
  extra_files = {'ref.dat','vrgf.dat'};
  alt_extra_files = {'ref_sav.dat','vrgf_sav.dat'};
  % check for rev_sav.dat and vrgf_sav.dat files
  for i=1:length(extra_files)
    fname_alt = alt_extra_files{i};
    fname_extra = extra_files{i};
    if exist(sprintf('%s/%s',parms.pdir,fname_alt),'file')
      cmd = sprintf('mv %s/%s %s/%s',...
        parms.pdir,fname_alt,parms.pdir,fname_extra);
      fprintf('%s: moving %s to %s...\n',...
        mfilename,fname_alt,fname_extra);
      [s,r] = unix(cmd);
      if s
        error('cmd %s failed:\n%s',cmd,r);
      end;
    end;
  end;  
  % check for ref.dat and vrgf.dat files
  for i=1:length(extra_files)
    fname_extra = extra_files{i};
    dlist = dir(sprintf('%s/*%s',parms.pdir,fname_extra));
    if length(dlist)>1
      fprintf('%s: WARNING: multiple %s files in %s\n',...
        mfilename,fname_extra,parms.pdir);
      SUID = parms.jinfo.SeriesInstanceUID;
      xdir = sprintf('%s/extra',parms.pdir);
      mmil_mkdir(xdir);
      for i=1:length(dlist)
        fname = dlist(i).name;
        if ~isempty(regexp(fname,SUID))
          cmd = sprintf('cp -p %s/%s %s/%s',...
            parms.pdir,fname,parms.pdir,fname_extra);
          fprintf('%s: copying %s to %s...\n',...
            mfilename,fname,fname_extra);
          [s,r] = unix(cmd);
          if s
            error('cmd %s failed:\n%s',cmd,r);
          end;
        end;
        cmd = sprintf('mv %s/%s %s',parms.pdir,fname,xdir);
        [s,r] = unix(cmd);
        if s
          error('cmd %s failed:\n%s',cmd,r);
        end;
      end;
    end;
    dlist = dir(sprintf('%s/*%s',parms.pdir,fname_extra));
    if isempty(dlist)
      errcode = 1;
      errmsg = sprintf('no %s file in %s\n',...
        fname_extra,parms.pdir);
      break;
    end;    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if ~exist(parms.mb_tmpdir,'dir')
%    fprintf('%s: WARNING: mb_tmpdir %s not found\n',mfilename,parms.mb_tmpdir);
    return;
  end;
  cmd = sprintf('rm -r %s',parms.mb_tmpdir);
  fprintf('%s: deleting %s...\n',mfilename,parms.mb_tmpdir);
  [s,r] = unix(cmd);
  if s
    error('cmd %s failed:\n%s',cmd,r);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_pfile(parms)
  % delete P file if it exists and there were no errors
  fname_pfile = get_pfile(parms);
  if ~isempty(fname_pfile)
    fprintf('%s: recon successful, deleting %s...\n',mfilename,fname_pfile);
    delete(fname_pfile);
  else
    fprintf('%s: WARNING: recon successful, but no P file found to delete\n',mfilename);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [datadir,errcode] = rename_data(parms)
  datadir = []; errcode = 0; errmsg = [];
  % read dicom header
  fname_dicom = sprintf('%s/%s.dcm',parms.unpackdir,parms.fstem_dcm);
  if exist(fname_dicom,'file')
    series_info = dicominfo(fname_dicom);
  else
    series_info = parms.jinfo;
  end;
  % check series dir
  sdir = sprintf('%s/%s/%s',parms.unpackdir,...
    series_info.StudyInstanceUID,series_info.SeriesInstanceUID);
  if ~exist(sdir)
    sdir = sprintf('%s/%s',parms.unpackdir,...
      series_info.SeriesInstanceUID);
  end;
  if ~exist(sdir)
    % check for match to alternate naming convention
    fpat = sprintf('%s/SUID_%s*_se%s_*',...
      parms.unpackdir,parms.jinfo.StudyInstanceUID,parms.jinfo.SeriesNumber);
    jlist = dir(fpat);
    if ~isempty(jlist)
      jlist = jlist([jlist.isdir]);
    end;
    if ~isempty(jlist)
      sdir = sprintf('%s/%s',parms.unpackdir,jlist.name);
    end;
  end;
  fname_json = sprintf('%s.json',sdir);
  % check for errors
  fname_err = sprintf('%s/%s.unpack_error',parms.unpackdir,parms.fstem_json);
  if parms.force_unpack_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing rename after previous unpack error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
    delete(fname_err);
  end;
  fname_err = sprintf('%s/%s.recon_error',parms.unpackdir,parms.fstem_json);
  if parms.force_recon_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing rename after previous recon error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
    delete(fname_err);
  end;
  fname_err = sprintf('%s/%s.dicoms_error',parms.unpackdir,parms.fstem_json);
  if parms.force_dicoms_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing rename after previous dicom error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
    delete(fname_err);
  end;
  fname_err = sprintf('%s/%s.rename_error',parms.unpackdir,parms.fstem_json);
  if parms.force_rename_err_flag && exist(fname_err,'file')
    fprintf('%s: forcing rename after previous rename error for %s...\n',...
      mfilename,parms.fname_tar);
    parms.forceflag = 1;
    delete(fname_err);
  end;

  fname_out = sprintf('%s/%s.renamed',parms.unpackdir,parms.fstem_json);
  fname_err = sprintf('%s/%s.rename_error',parms.unpackdir,parms.fstem_json);
  if (~exist(fname_out,'file') && ~exist(fname_err,'file')) || parms.forceflag
    if ~exist(sdir,'dir')
      errcode = 1;
      errmsg = sprintf('failed to rename: series dir %s not found',sdir);
      fname_out = fname_err;
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    else
      % get header info from one dicom
      series_info = [];
      flist = dir(sdir);
      fnames = setdiff({flist.name},{'.','..'});
      for d=1:length(fnames)
        fname = fnames{d};
        dname = sprintf('%s/%s',sdir,fname);
        if mmil_isdicomfile(dname)
          series_info = dicominfo(dname);
          break;
        end;
      end;
      if isempty(series_info)
        errcode = 1;
        errmsg = sprintf('no dicom files found in %s',sdir);
        fname_out = fname_err;
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      else
        % get SubjID from fstem_json if possible or jinfo if not
        SubjID = abcd_get_SubjID(parms.fstem_json,series_info);
        % get other info
        StudyDate = series_info.StudyDate;
        StudyTime = series_info.StudyTime;
        % if StudyTime includes ".", remove everything after and including the "."
        k = regexp(StudyTime,'^(?<Time>[^]+)\.','names');
        if ~isempty(k), StudyTime = k.Time; end;
        snum = series_info.SeriesNumber;
        % set prefix for manufacturer
        switch series_info.Manufacturer
          case 'GE MEDICAL SYSTEMS'
            manu_prefix = 'G';
          case 'Philips Medical Systems'
            manu_prefix = 'P';
          case 'SIEMENS'
            manu_prefix = 'S';
          otherwise
            manu_prefix = 'X';
            fprintf('%s: WARNING: unrecognized manufacturer: %s\n',mfilename,series_info.Manufacturer);
        end;
        % get Device Serial Number, etc.
        DSN = series_info.DeviceSerialNumber;
        InstitutionName = series_info.InstitutionName;
        StationName = series_info.StationName;
        % look up site number using DSN
        ids = find(strcmp(parms.site_dsn,DSN));
        if isempty(ids)
          errcode = 1;
          errmsg = sprintf('site with DSN %s not found',DSN);
          errmsg = sprintf('%s\n          InstitutionName = %s',errmsg,InstitutionName);
          errmsg = sprintf('%s\n          StationName = %s',errmsg,StationName);
          fname_out = fname_err;
          fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        else
          site_name = parms.site_names{ids};
          site_num = parms.site_nums{ids};
          % set VisitID
          VisitID = sprintf('%s%03d_%s_%s_%s',...
                            manu_prefix,site_num,SubjID,StudyDate,StudyTime);
          if parms.separate_series_flag
            VisitID = sprintf('%s_se%d',VisitID,snum);
          end;
          % link series dir to be a subdirectory of newdir
          newdir = sprintf('%s/%s',parms.outdir,VisitID);
          mmil_mkdir(newdir);
          % set datadir
          [tmp,fstem,fext] = fileparts(sdir);
          fstem = [fstem fext];
          datadir = sprintf('%s/%s',newdir,fstem);
          fname_json_link = sprintf('%s.json',datadir);
          % link json file
          if ~exist(fname_json,'file')
            fprintf('%s: WARNING: json file %s not found\n',mfilename,fname_json);
          else
            % remove link if it exists
            if exist(fname_json_link,'file')
              cmd = sprintf('rm %s',fname_json_link);
              [s,r] = unix(cmd);
              if s
                error('failed to remove old json file link:\n%s\n%s\n',...
                  mfilename,cmd,r);
              end;
            end;
            if ~exist(fname_json_link,'file')
              cmd = sprintf('ln -s %s %s',fname_json,fname_json_link);
              [s,r] = unix(cmd);
              if s
                fprintf('%s: WARNING: failed to link json file:\n%s\n%s\n',...
                  mfilename,cmd,r);
              end;
            end;
          end;
          % remove link if it exists
          if exist(datadir,'file')
            cmd = sprintf('rm %s',datadir);
            [s,r] = unix(cmd);
            if s
              errcode = 1;
              errmsg = sprintf('failed to remove old datadir link:\n%s\n%s',cmd,r);
              fname_out = fname_err;
              fprintf('%s: ERROR: %s\n',mfilename,errmsg);
            end;
          end;
          if ~errcode
            % create link
            if ~exist(datadir,'file')
              fprintf('%s: linking data from %s into %s...\n',mfilename,sdir,datadir);
              cmd = sprintf('ln -s %s %s',sdir,datadir);
              [s,r] = unix(cmd);
              % check for failure of unix command
              if s
                errcode = 1;
                errmsg = sprintf('failed to link:\n%s\n%s',cmd,r);
                fname_out = fname_err;
                fprintf('%s: ERROR: %s\n',mfilename,errmsg);
              end;
            end;
          end;
        end;
      end;
    end;
    % create check file
    fid = fopen(fname_out,'wt');
    if fid<0
      error('failed to open file %s for writing',fname_out);
    end;
    if ~errcode
      fprintf(fid,'%s',datadir);
    elseif ~isempty(errmsg)
      fprintf(fid,'ERROR: %s\n',errmsg);
    else
      fprintf(fid,'ERROR: rename failed\n');
    end;
    fclose(fid);
  elseif exist(fname_out,'file')
    fprintf('%s: skipping rename because already complete\n',mfilename);
    if exist(fname_out,'file')
      fid = fopen(fname_out,'rt');
      datadir = fscanf(fid,'%s');
    end;
  elseif exist(fname_err,'file')
    fprintf('%s: skipping rename because of previous errors\n',mfilename);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

