function cmd = fs_recon(subj,varargin)
%function cmd = fs_recon(subj,varargin)
%
% Purpose: create csh script for running FreeSurfer's recon-all
%
% Usage: fs_recon(subj,'key1',value1,...)
%
% Required Parameters:
%  subj: subject name (name of FreeSurfer recon output dir)
%
% Options specifying input / output:
%  'fname_T1': input T1-weighted file name
%     Required if subj dir does not already exist
%    {default = []}
%  'fname_T2': input T2-weighted file name
%     Used to refine pial surface if supplied
%    {default = []}
%  'fname_wm': input white matter file name to replace auto white matter mask
%    {default = []}
%  'fname_bm': input brain mask file name to replace FreeSurfer brain mask
%    {default = []}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    {default = $SUBJECTS_DIR}
%
% Options specifying whether to run particular steps:
%  'rescale_orig_flag': [0|1] whether to allow orig.mgz to be rescaled
%    when converting to uchar
%    if 0, orig.mgz is a uchar copy of rawavg.mgz
%       as long as highest value in input image is 255,
%       orig.mgz will not be rescaled
%    {default = 1}
%  'talairach_flag': [0|1] whether to do Talairach registration
%    {default = 1}
%  'nu_flag': [0|1[ whether to run nonuniformity intensity correction
%    if 0, nu.mgz is simply a copy of orig.mgz
%    {default = 1}
%  'nu3T_flag': [0|1] whether to use nu settings
%      optimized for 3T (version >= 500 only)
%    irrelevant if nu_flag = 0
%    {default = 0}
%  'noaseg_flag': [0|1] do not run subcortical segmentation
%    and do not use aseg results in cortical reconstruction
%    {default = 0}
%  'nocort_flag': [0|1] do not run cortical reconstruction
%    {default = 0}
%  'labelV1_flag': [0|1] whether to label V1 (version >= 400 only)
%    {default = 0}
%  'surfsegedit_flag': [0|1] whether to edit aseg.mgz with
%    cortical surfaces (version >= 500 only)
%    {default = 0}
%  'steps': cell array of the steps to run:
%    if empty, will be created to match standard set of reconstruction
%      steps for the specified FreeSurfer version
%    {default = []}
%
% Options for GCA atlas (used for aseg):
%  'gca_dir': full path containing substitute GCA atlas files (for aseg)
%    {default = []}
%  'gca': name of substitute GCA atlas file (for aseg)
%    {default = []}
%  'gca_skull': name of substitute GCA atlas file with skull (for aseg)
%    {default = []}
%
% Options for rerunning steps (i.e. after manual intervention):
%  'remake_all': [0|1] whether to run all steps (even if done already)
%    {default = 0}
%  'remake_bm': [0|1] whether to run all steps following brain mask
%    (will lose any edits to white matter segmentation)
%    {default = 0}
%  'remake_cp': [0|1] whether to run all steps starting with intensity norm 2
%    (for control points; will lose any edits to white matter segmentation)
%    {default = 0}
%  'remake_cp_quick': [0|1] whether to run intensity norm 2
%    and white matter segmentation
%    (for control points; will lose any edits to white matter segmentation)
%    {default = 0}
%  'remake_surf': [0|1] whether to remake surfaces for white matter edits
%    {default = 0}
%  'remake_final': [0|1] whether to remake final surfaces for brain mask edits
%    {default = 0}
% NOTE: These remake options direct the removal of sets of touch files,
%     causing certain steps to be re-run 
%   These options can also be set by the presence of files of the same
%     names in the touch directory of the subj directory (e.g. 'remake_all')
%   If there are multiple remake options set, the order of precedence
%     is as listed above (only one will be used)
%
% Other options:
%  'make_surfs_noauto_flag': [0|1] whether to set -noauto flag
%     when running mris_make_surfaces
%    {default = 0}
%  'make_surfs_parms_file': full path name of text file containing
%     parameter names and values (e.g. min_gray_at_white_border, etc.)
%     (one parameter per line, values separated by spaces)
%    {default = []}
%  'conform_min_flag': [0|1] resample volumes to lowest voxel size
%     for use with higher than 1mm resolution scans
%    {default = 0}
%  'version': FreeSurfer version number (e.g. 302, 305, 450, 510, 530)
%    {default = 530}
%  'touchonly_flag': [0|1] whether to create touchfiles and return empty cmd
%    {default = 0}
%  'forceflag': [0|1] whether to output recon commands even if recon is complete
%    To force overwrite of FreeSurfer output, use forceflag=1 and remake_all=1
%    {default = 0}
%
% Output:
%   cmd: csh commands to run for FreeSurfer recon
%     If recon is complete, and forceflag = 0, and remake flags are all 0,
%     cmd will return as []
%
% Created:  10/03/10 by Don Hagler
% Prev Mod: 11/11/15 by Don Hagler
% Last Mod: 05/11/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: one remake option?
%%       accept string, e.g. 'remake_all', 'remake_surf'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(subj,varargin);

% create list of steps to run
if isempty(parms.steps)
  parms.steps = set_steps(parms);
end;

% check remake options / files
parms = check_remake(parms);

% check for final touchfile to see if recon is complete
if ~parms.forceflag && isempty(parms.remake_type)
  fname_touch = [parms.touchdir '/' parms.fname_finish_all];
  if exist(fname_touch,'file'), return; end;
end;

% initialize cmd
cmd = init_cmd(parms);

% add each step to cmd
cmd = build_cmd(parms,cmd);

% finish command by writing last touch file
cmd = finish_cmd(parms,cmd);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(subj,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'subj',subj,[],...
  ... % input/output options
    'fname_T1',[],[],...
    'fname_T2',[],[],...
    'fname_wm',[],[],...
    'fname_bm',[],[],...
    'subjdir',[],[],...
  ... % step options
    'rescale_orig_flag',true,[false true],...
    'talairach_flag',true,[false true],...
    'nu_flag',true,[false true],...
    'nu3T_flag',false,[false true],...
    'noaseg_flag',false,[false true],...
    'nocort_flag',false,[false true],...
    'labelV1_flag',false,[false true],...
    'surfsegedit_flag',false,[false true],...
    'steps',[],[],...
  ... % gca options
    'gca_dir',[],[],...
    'gca',[],[],...
    'gca_skull',[],[],...
  ... % remake options
    'remake_all',false,[false true],...
    'remake_bm',false,[false true],...
    'remake_cp',false,[false true],...
    'remake_cp_quick',false,[false true],...
    'remake_surf',false,[false true],...
    'remake_final',false,[false true],...
  ... % other options
    'make_surfs_noauto_flag',false,[false true],...
    'make_surfs_parms_file',[],[],...
    'conform_min_flag',false,[false true],...
    'version',530,[1,1000],...
    'touchonly_flag',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden options
    'fname_start_all','fs.start.all.touch',[],...
    'fname_finish_all','fs.finish.all.touch',[],...
    'fname_log','recon-all.log',[],...
    'ribbon_steps',{'cortribbon','aparc2aseg','wmparc'},[],...
    'remake_types',{'all','bm','cp','cp_quick','surf','final'},[],...
  ... % lists of steps corresponding to remake flags (order not important)
    'all_steps',{'motioncor' 'nuintensitycor' 'talairach' 'tal-check' ...
      'normalization' 'skullstrip'...
      'gcareg' 'canorm' 'careg' 'careginv' ...
      'rmneck' 'skull-lta' 'calabel' 'segstats' 'normalization2' ...
      'maskbfs' 'segmentation' 'fill' 'tessellate' 'smooth1' ...
      'inflate1' 'qsphere' 'fix' 'finalsurfs' 'white' 'smooth2' 'inflate2' ...
      'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
      'cortribbon' 'contrasurfreg' 'avgcurv' 'cortparc' 'pial' 'T2pial' 'pctsurfcon'...
      'parcstats' 'cortparc2' 'parcstats2'...
      'cortparc3' 'parcstats3'...
      'aparc2aseg' 'wmparc' ...
      'label_v1' 'surfsegedit' 'label-exvivo-ec' 'ba-labels'},[],...
    'cp_steps',{'normalization' 'skullstrip'...
      'gcareg' 'canorm' 'careg' 'careginv' ...
      'rmneck' 'skull-lta' 'calabel' 'segstats' 'normalization2' ...
      'maskbfs' 'segmentation' 'fill' 'tessellate' 'smooth1' ...
      'inflate1' 'qsphere' 'fix' 'finalsurfs' 'white' 'smooth2' 'inflate2' ...
      'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
      'cortribbon' 'contrasurfreg' 'avgcurv' 'cortparc' 'pial' 'T2pial' 'pctsurfcon'...
      'parcstats' 'cortparc2' 'parcstats2'...
      'cortparc3' 'parcstats3'...
      'aparc2aseg' 'wmparc' ...
      'label_v1' 'surfsegedit' 'label-exvivo-ec' 'ba-labels'},[],...
    'cp_quick_steps',{'normalization2' 'segmentation'},[],...
    'bm_steps',{'gcareg' 'canorm' 'careg' 'careginv' ...
      'rmneck' 'skull-lta' 'calabel' 'segstats' 'normalization2' ...
      'maskbfs' 'segmentation' 'fill' 'tessellate' 'smooth1' ...
      'inflate1' 'qsphere' 'fix' 'finalsurfs' 'white' 'smooth2' 'inflate2' ...
      'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
      'cortribbon' 'contrasurfreg' 'avgcurv' 'cortparc' 'pial' 'T2pial' 'pctsurfcon'...
      'parcstats' 'cortparc2' 'parcstats2'...
      'cortparc3' 'parcstats3'...
      'aparc2aseg' 'wmparc' ...
      'label_v1' 'surfsegedit' 'label-exvivo-ec' 'ba-labels'},[],...
    'surf_steps',{'maskbfs' 'segmentation' 'fill' 'tessellate' 'smooth1' ...
      'inflate1' 'qsphere' 'fix' 'finalsurfs' 'white' 'smooth2' 'inflate2' ...
      'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
      'cortribbon' 'contrasurfreg' 'avgcurv' 'cortparc' 'pial' 'T2pial' 'pctsurfcon'...
      'parcstats' 'cortparc2' 'parcstats2'...
      'cortparc3' 'parcstats3'...
      'aparc2aseg' 'wmparc' ...
      'label_v1' 'surfsegedit' 'label-exvivo-ec' 'ba-labels'},[],...
    'final_steps',{'maskbfs' 'finalsurfs' 'white' 'smooth2' 'inflate2' ...
      'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
      'cortribbon' 'contrasurfreg' 'avgcurv' 'cortparc' 'pial' 'T2pial' 'pctsurfcon'...
      'parcstats' 'cortparc2' 'parcstats2'...
      'cortparc3' 'parcstats3'...
      'aparc2aseg' 'wmparc' ...
      'label_v1' 'surfsegedit' 'label-exvivo-ec' 'ba-labels'},[],...
  });

  % set additional parameters  
  if isempty(parms.subjdir)
    parms.subjdir = getenv('SUBJECTS_DIR');
    if isempty(parms.subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  end;
  parms.fspath = [parms.subjdir '/' parms.subj];
  parms.touchdir = [parms.fspath '/touch'];
  if isempty(parms.fname_T1)
    fname_init = [parms.fspath '/mri/orig/001.mgz'];
    if ~exist(fname_init,'dir')
      error('fname_T1 not supplied and %s does not exist',fname_init);
    end;
  end;
  if parms.version<500
    parms.nu3T_flag = false;
    parms.surfsegedit_flag = false;
    parms.conform_min_flag = false;
  end;
  if parms.version >= 500
    parms.surf_steps = cat(2,parms.surf_steps,'segstats');
    parms.final_steps = cat(2,parms.final_steps,'segstats');
  end;
  if parms.version>=400
    parms.options = '-no-isrunning';
  else
    parms.options = [];
  end;
  if parms.noaseg_flag
    parms.options = [parms.options ' -noaseg'];
  end;
  if ~isempty(parms.gca_dir)
    parms.options = sprintf('%s -gca-dir %s',parms.options,parms.gca_dir);
  end;
  if ~isempty(parms.gca)
    parms.options = sprintf('%s -gca %s',parms.options,parms.gca);
  end;
  if ~isempty(parms.gca_skull)
    parms.options = sprintf('%s -gca-skull %s',parms.options,parms.gca_skull);
  end;
  if parms.version < 500
    parms.options = sprintf('%s -force',parms.options);
  end;
  parms.make_surfs_parms = [];
  if parms.make_surfs_noauto_flag
    parms.make_surfs_parms{1} = 'noauto';
  end;
  if ~isempty(parms.make_surfs_parms_file)
    if ~exist(parms.make_surfs_parms_file,'file')
      error('make_surfs_parms_file %s not found',parms.make_surfs_parms_file);
    end;
    tmp_parms = mmil_readtext(parms.make_surfs_parms_file);
    parms.make_surfs_parms = cat(1,parms.make_surfs_parms,tmp_parms);  
  end;
  if ~isempty(parms.make_surfs_parms)
    parms.fname_xopts = '$fspath/scripts/expert_mris_make_surfaces.opts';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function steps = set_steps(parms)
  % create list of steps to run (order of steps will be followed)
  steps = {'motioncor'};
  if parms.version >= 500 && parms.talairach_flag
    steps = cat(2,steps,{'talairach','tal-check'});
  end;
  steps = cat(2,steps,'nuintensitycor');
  if parms.version < 500 && parms.talairach_flag
    steps = cat(2,steps,'talairach');
    if parms.version >= 400
      steps = cat(2,steps,'tal-check');
    end;
  end;
  steps = cat(2,steps,{...
    'normalization','skullstrip','gcareg','canorm','careg','careginv',...
    'rmneck','skull-lta','calabel'});
  if parms.version < 500 || parms.nocort_flag
    steps = cat(2,steps,'segstats');
  end;
  if ~parms.nocort_flag
    steps = cat(2,steps,{...
      'normalization2','maskbfs','segmentation','fill','tessellate',...
      'smooth1','inflate1','qsphere','fix'});
    if parms.version >= 500
      steps = cat(2,steps,'white');
    else
      steps = cat(2,steps,'finalsurfs');
    end;
    steps = cat(2,steps,{'smooth2','inflate2'});
    if parms.version >= 400
      steps = cat(2,steps,'curvstats');
    end;
    steps = cat(2,steps,{'sphere','surfreg'});
    if parms.version >= 400
      steps = cat(2,steps,'jacobian_white');
    else
      steps = cat(2,steps,'contrasurfreg');
    end;
    steps = cat(2,steps,{'avgcurv','cortparc'});
    if parms.version >= 500
      steps = cat(2,steps,'pial');
    end;  
    if parms.version >= 530
      if ~isempty(parms.fname_T2)
        steps = cat(2,steps,'T2pial');
      end;
      steps = cat(2,steps,'pctsurfcon');
    end;
    steps = cat(2,steps,{'parcstats','cortparc2','parcstats2'});
    if parms.version >= 530
      steps = cat(2,steps,{'cortparc3','parcstats3'});
    end;
    if parms.version >= 400 && parms.labelV1_flag
      steps = cat(2,steps,'label_v1');
    end;
    if parms.version >= 500 && parms.surfsegedit_flag
      steps = cat(2,steps,'surfsegedit');
    end;
    steps = cat(2,steps,'cortribbon');
    if parms.version >= 500
      steps = cat(2,steps,'segstats');
    end;
    steps = cat(2,steps,'aparc2aseg');
    if parms.version >= 400
      steps = cat(2,steps,'wmparc');
    end;
    if parms.version >= 500
      steps = cat(2,steps,{'label-exvivo-ec','ba-labels'});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_remake(parms)
  % check for remake files in touch dir
  remake_files = [];
  for r=1:length(parms.remake_types)
    remake_type = parms.remake_types{r};
    remake_file = sprintf('%s/remake_%s',...
      parms.touchdir,remake_type);
    if exist(remake_file,'file')
      remake_files{end+1} = remake_file;
      parms.(['remake_' remake_type]) = 1;
    end;
  end;

  % determine remake type
  remake_type = [];
  remake_steps = [];
  if parms.remake_all
    remake_type = 'all';
    remake_steps = parms.all_steps;
  elseif parms.remake_bm
    remake_type = 'bm';
    remake_steps = parms.bm_steps;
  elseif parms.remake_cp
    remake_type = 'cp';
    remake_steps = parms.cp_steps;
  elseif parms.remake_cp_quick
    remake_type = 'cp_quick';
    remake_steps = parms.cp_quick_steps;
  elseif parms.remake_surf
    remake_type = 'surf';
    remake_steps = parms.surf_steps;
  elseif parms.remake_final
    remake_type = 'final';
    remake_steps = parms.final_steps;
  end;

  % remove touch files
  if ~isempty(remake_type)  
    % create flist of touch files to remove
    flist = cell(2*length(remake_steps)+2,1);
    f = 1;
    flist{f} = parms.fname_start_all; f=f+1;
    flist{f} = parms.fname_finish_all; f=f+1;
    for r=1:length(remake_steps)
      flist{f} = ['fs.start.'  remake_steps{r} '.touch']; f=f+1;
      flist{f} = ['fs.finish.' remake_steps{r} '.touch']; f=f+1;
    end;

    touchfiles = dir([parms.touchdir '/fs.*.touch']);
    fnames = {touchfiles.name};
    ind = find(ismember(fnames,flist));
    fnames = fnames(ind);
    if length(fnames)>0
      fname_list = [];
      for f=1:length(fnames)
        fname_list = [fname_list ' ' parms.touchdir '/' fnames{f}];
      end;
      [status,result] = unix(['rm ' fname_list]);
      if status
        error('failed to remove touch files from %s:\n%s',...
          parms.touchdir,result);
      end;
    end;

    % create "remade" file
    remade_file = sprintf('%s/remade_%s',...
      parms.touchdir,remake_type);
    [status,result]=unix(sprintf('touch %s\n',remade_file));
    if status
      error('failed to create remade file %s:\n%s',...
        remade_file,result);
    end;

    % remove all "remake" files
    if ~isempty(remake_files)
      fname_list = [];
      for f=1:length(remake_files)
        fname_list = [fname_list ' ' remake_files{f}];
      end;
      [status,result] = unix(['rm ' fname_list]);
      if status
        error('failed to remove remake files from %s:\n%s',...
          parms.touchdir,result);
      end;
    end;

    % limit list of steps to run
    ind = find(ismember(parms.steps,remake_steps));
    parms.steps = parms.steps(ind);
  end;

  parms.remake_type = remake_type;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmd = init_cmd(parms)
  cmd = [];

  % determine whether running on cluster (so we can skip cortical ribbon)
  if parms.version>=500
    cmd = sprintf('%sset rocks_flag = 0\n',cmd);
  else
    cmd = sprintf('%sif (-f /etc/rocks-release) then\n',cmd);
    cmd = sprintf('%s  set rocks_flag = 1\n',cmd);
    cmd = sprintf('%selse\n',cmd);
    cmd = sprintf('%s  set rocks_flag = 0\n',cmd);
    cmd = sprintf('%sendif\n',cmd);
  end;

  % set subject-specific variables
  cmd = sprintf('%s\n',cmd);
  cmd = sprintf('%sset subjdir = %s\n',cmd,parms.subjdir);
  cmd = sprintf('%sset subj = %s\n',cmd,parms.subj);
  cmd = sprintf('%sset fspath = %s\n',cmd,parms.fspath);
  cmd = sprintf('%sset touchdir = %s\n',cmd,parms.touchdir);

  % if fspath does not already exists, create it and touchdir
  cmd = sprintf('%s\n',cmd);
  cmd = sprintf('%sif (! -e $touchdir) then\n',cmd);
  cmd = sprintf('%s  mkdir -p $touchdir\n',cmd);
  cmd = sprintf('%sendif\n',cmd);

  % if orig/001.mgz does not already exists, copy fname_T1
  cmd = sprintf('%s\n',cmd);
  cmd = sprintf('%sif (! -e $fspath/mri/orig/001.mgz) then\n',cmd);
  cmd = sprintf('%s  mkdir -p $fspath/mri/orig\n',cmd);
  cmd = sprintf('%s  mri_convert %s $fspath/mri/orig/001.mgz\n',...
    cmd,parms.fname_T1);
  cmd = sprintf('%sendif\n',cmd);

  % if orig/T1raw.mgz does not already exists, copy fname_T2
  if ~isempty(parms.fname_T2)
    cmd = sprintf('%s\n',cmd);
    cmd = sprintf('%sif (! -e $fspath/mri/orig/T2raw.mgz) then\n',cmd);
    cmd = sprintf('%s  mkdir -p $fspath/mri/orig\n',cmd);
    cmd = sprintf('%s  mri_convert %s $fspath/mri/orig/T2raw.mgz\n',...
      cmd,parms.fname_T2);
    cmd = sprintf('%sendif\n',cmd);
  end;

  % create expert options file if make_surfs_noauto_flag
  if ~isempty(parms.make_surfs_parms)
    cmd = sprintf('%s\n',cmd);
    cmd = sprintf('%sif (! -e %s) then\n',cmd,parms.fname_xopts);
    cmd = sprintf('%s  mkdir -p $fspath/scripts\n',cmd);
    cmd = sprintf('%s  echo mris_make_surfaces',cmd);
    for i=1:length(parms.make_surfs_parms)
      cmd = sprintf('%s -%s',cmd,parms.make_surfs_parms{i});
    end;
    cmd = sprintf('%s > %s\n',cmd,parms.fname_xopts);
    cmd = sprintf('%sendif\n',cmd);
  end;
  cmd = sprintf('%s\n',cmd);
  cmd = sprintf('%secho "started all on `date`" > $touchdir/%s\n',...
    cmd,parms.fname_start_all);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmd = build_cmd(parms,cmd)
  % loop over steps, add to cmd
  for s=1:length(parms.steps)
    step = parms.steps{s};
    fname_start =  ['fs.start.'  step '.touch'];
    fname_finish = ['fs.finish.' step '.touch'];
    cmd = sprintf('%s\n',cmd);
    cmd = sprintf('%sif (! -e $touchdir/%s) then\n',cmd,fname_finish);
    if ismember(step,parms.ribbon_steps)
      % for earlier versions, ribbon_steps cannot run on cluster
      cmd = sprintf('%s  if ($rocks_flag == 0) then\n',cmd);
    end;
    cmd = sprintf('%s    echo "started %s on `date`" > $touchdir/%s\n',...
      cmd,step,fname_start);
    if ~parms.touchonly_flag
      if strcmp(step,'nuintensitycor') && parms.nu3T_flag
        step = 'nuintensitycor-3T';
      end;
      if strcmp(step,'nuintensitycor') && ~parms.nu_flag
        tmp_cmd ='cp -p $fspath/mri/orig.mgz $fspath/mri/nu.mgz';
        cmd = sprintf('%s    echo %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
        cmd = sprintf('%s    %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
      else
        cmd = sprintf('%s    recon-all -sd $subjdir -subject $subj %s -%s',...
          cmd,parms.options,step);
        if ismember(step,{'white','pial','finalsurfs'}) &&...
           parms.make_surfs_noauto_flag
          % pass expert options file with -experts flag
          cmd = sprintf('%s -expert %s -xopts-overwrite',cmd,parms.fname_xopts);
        end;
      end;
      if parms.conform_min_flag
        cmd = sprintf('%s -cm',cmd);
      end;
      cmd = sprintf('%s\n',cmd);
      cmd = sprintf('%s    if (! -e $fspath/scripts/%s) exit\n',...
        cmd,parms.fname_log);
      cmd = sprintf('%s    set errstr = `tail -5 $fspath/scripts/%s | grep ERROR`\n',...
        cmd,parms.fname_log);
      cmd = sprintf('%s    if ($#errstr > 0) exit\n',cmd);
      % special cases
      switch step
        case 'motioncor'
          if ~parms.rescale_orig_flag
            tmp_cmd = 'mri_convert -odt uchar -ns 1 $fspath/mri/rawavg.mgz $fspath/mri/orig.mgz';
            cmd = sprintf('%s    echo %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
            cmd = sprintf('%s    %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
            tmp_cmd = 'mri_add_xform_to_header -c $fspath/mri/transforms/talairach.xfm $fspath/mri/orig.mgz $fspath/mri/orig.mgz';
            cmd = sprintf('%s    echo %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
            cmd = sprintf('%s    %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
          end;
        case 'skullstrip'
          if ~isempty(parms.fname_bm)
            tmp_cmd = sprintf('mri_convert %s $fspath/mri/brainmask.mgz',parms.fname_bm);
            cmd = sprintf('%s    echo %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
            cmd = sprintf('%s    %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
          end;
        case 'segmentation'
          if ~isempty(parms.fname_wm)
            tmp_cmd = sprintf('mri_convert %s $fspath/mri/wm.mgz',parms.fname_wm);
            cmd = sprintf('%s    echo %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
            cmd = sprintf('%s    %s >> $fspath/scripts/recon-all.log\n',cmd,tmp_cmd);
          end;
      end;
    end;
    cmd = sprintf('%s    echo "finished %s on `date`" > $touchdir/%s\n',...
      cmd,step,fname_finish);
    if ismember(step,parms.ribbon_steps)
      cmd = sprintf('%s  endif\n',cmd);
    end;
    cmd = sprintf('%sendif\n',cmd);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmd = finish_cmd(parms,cmd)
  cmd = sprintf('%s\n',cmd);
  cmd = sprintf('%secho "finished all on `date`" > $touchdir/%s\n',...
    cmd,parms.fname_finish_all);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



