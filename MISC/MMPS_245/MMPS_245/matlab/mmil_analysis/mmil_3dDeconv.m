function [fname_out,fname_info,fnames_iresp] = mmil_3dDeconv(fname,stim_fnames,varargin)
%function [fname_out,fname_info,fnames_iresp] = mmil_3dDeconv(fname,stim_fnames,varargin)
%
% Purpose: wrapper around AFNI's 3dDeconvolve
%   input mgh file, get back mgh file as output
%
% Usage:
%  mmil_3dDeconv(fname,stim_fnames,'key1', value1,...);
%
% Required Parameters:
%   fname: full or relative path name of input mgh file
%     This can be a cell array of input file names to be concatenated
%   stim_fnames: string or cell array of strings with full or relative path
%     names of "1D" stimulus time course files
%     If "fname" is a cell array, "stim_fnames" should be a nested cell array
%
% Optional parameters:
%  'fname_motion': name of 1D file containing motion estimates from 3dvolreg
%    to be used as regressors
%    If "fname" is a cell array, "fname_motion" should be matching cell array
%    {default = []}
%  'stim_labels': string or cell array of strings with stimulus condition names
%    If empty, will use labels such as "cond1", "cond2", etc.
%    {default = []}
%   'stim_times_flag': [0|1|2] use stimulus time files
%     0 = .1D, 1 = .txt and 2 = _block.txt 
%     {default = 1}
%  'stim_times_model': name of stimulus model
%    model names in 3dDeconvolve for -stim_times
%     allowed:  'SPMG','TENT','CSPLIN','GAM','BLOCK'
%     {default = 'SPMG'}
%  'stim_times_nbasis': max number of basis functions for HRF
%     only applies for TENT or CSPLIN
%    {default = 10}
%  'contrasts_flag': [0|1] calculate glt contrasts between each condition
%    {default = 0}
%  'iresp_flag': [0|1] output impulse response functions
%    {default = 0}
%  'censor_flag': [0|1] whether to censor 
%    {default = true}
%  'contigous_censor_nTRs': censor TRs in which there are 
%    not contigous_censor_nTRs contiguous, non-censored TRs.     
%    {default = 0}
%  'deriv_flag': [0|1] whether to include motion derivatives in the analysis 
%    {default = true}
%  'reml_flag': [0|1] whether to use AFNI's 3dREMLfit in the analysis 
%    {default = false}
%  'outdir': output directory
%    If empty, will write output files to directory containing input fname
%    {default = []}
%  'outstem': output file stem
%    If empty, will write output files with same stem as fname
%    {default = []}
%  'skipTRs': number of initial repetitions to remove from data and motion files
%    NOTE: it is assumed that stim files do NOT include these extra repetitions
%    {default = 0}
%  'minlag': number of TRs for minimum "lag" between stimulus and response
%    {default = 0}
%  'maxlag': number of TRs for maximum "lag" between stimulus and response
%    {default = 4}
%  'norm_flag': [0|1] whether to normalize input timeseries
%    by mean for each voxel (new mean = 100) before doing calculations
%    {default = 1}
%  'detrend': [0|1|2|3] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    3: cubic detrend
%    {default = 2}
%  'out_ext': output file extension ('.mgh' or '.mgz')
%    if empty, will use input file extension
%    {default = []}
%  'forceflag': [0|1] overwrite existing output files
%    {default = 0}
%  'goforit_flag': [0|1] whether to use -GOFORIT option in 
%    AFNI's 3dDeconvolve 
%    {default = true}
%  'goforit_val': Threshold for -GOFORIT option in 
%    AFNI's 3dDeconvolve 
%    {default = 10}
%
% Output:
%   fname_out: output file name
%   fname_info: name of mat file containing stats_info struct array
%      with information about each frame in fname_out
%      including name, type, and dofs
%   fnames_iresp: cell array of iresp file names (only if iresp_flag=1)
%
%   Output file names will be created with name based on input fname
%
% Created:  09/01/08 Don Hagler
% Prev Mod: 09/06/17 Don Hagler
% Last Mod: 10/27/17 Dani Cornejo 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];

if ~mmil_check_nargs(nargin, 2), return; end;
parms = mmil_args2parms(varargin, { ...
  'fnames',fname,[],...
  'stim_fnames',stim_fnames,[],...
  'stim_labels',[],[],...
  'stim_times_flag',1,0:2,...
  'stim_times_model','SPMG',{'SPMG','TENT','CSPLIN','GAM','BLOCK'},...
  'stim_times_nbasis',10,[1,100],...
  'contrasts_flag',false,[false true],...
  'iresp_flag',false,[false true],...
  'censor_flag',true,[false true],... 
  'contigous_censor_nTRs',0,[0 100],... 
  'deriv_flag',true,[false true],... 
  'reml_flag',false,[false true],...
  'censor_thresh',0.9,[0,10],... 
  'outdir',[],[],...
  'outstem',[],[],...
  'skipTRs',0,[0 Inf],...
  'TR',[],[0 Inf],...
  'norm_flag',true,[false true],...
  'detrend',3,0:3,...
  'thresh',10,[0 Inf],...
  'out_ext','.mgh',{'.mgh','.mgz'},...
  'minlag',0,[0,10],...
  'maxlag',4,[0,30],...
  'fname_motion',[],[],...
  'motion_labels',{'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'},[],...
  'glt_fnames',[],[],...
  'glt_labels',[],[],...
  'forceflag',false,[false true],...
...
  'goforit_flag',true,[false true],... 
  'goforit_val',10,1:20,... 
...
  'conds_contrast',[],[],...
...
  'glt_tags',{'nruns','stim_labels',...
              'stim_times_flag','stim_times_model','stim_times_nbasis',...
              'contrasts_flag','minlag','maxlag',...
              'detrend','motion_flag','forceflag',...
              'conds_contrast'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check/set parameters

fprintf('%s: checking input..\n',mfilename);

% check stim_times_model if reml_flag
if parms.reml_flag && ~ismember(parms.stim_times_model,{'GAM','SPGM'})
  error('unsupported stim times model with REML');
end;

% check input file(s)
if isempty(parms.fnames), error('input fname is empty'); end;
if ~iscell(parms.fnames), parms.fnames = {parms.fnames}; end;
parms.nfiles = length(parms.fnames);
fpath = [];
for f=1:parms.nfiles
  fname = parms.fnames{f};
  if ~exist(fname,'file'), error('file %s not found',fname); end;
  if isempty(fpath)
    [fpath,fstem,fext] = fileparts(fname);
    if ~ismember(fext,{'.mgh','.mgz'})
      error('input file must be mgh or mgz file type (has %s extension)',...
        fext);
    end;
  end;
end;
if isempty(parms.out_ext), parms.out_ext = fext; end;
if isempty(parms.outdir), parms.outdir = fpath; end;
if isempty(parms.outstem), parms.outstem = fstem; end;

% check motion file(s)
if ~isempty(parms.fname_motion)
  if ~iscell(parms.fname_motion), parms.fname_motion = {parms.fname_motion}; end;
  if parms.nfiles ~= length(parms.fname_motion)
    error('number of input files (%d) does not match number of motion files (%d)',...
      parms.nfiles,length(parms.fname_motion));
  end;
  parms.nmotion = length(parms.motion_labels);
  for f=1:parms.nfiles
    fname_motion = parms.fname_motion{f};
    if ~exist(fname_motion,'file')
      error('file %s not found',fname_motion);
    end;
  end;
  
else
  parms.nmotion = 0;
end;

% check stim files
if isempty(parms.stim_fnames), error('no stimulus files specified'); end;
if ~iscell(parms.stim_fnames) ||...
  (iscell(parms.stim_fnames) && ~iscell(parms.stim_fnames{1}))
  parms.stim_fnames = {parms.stim_fnames};
end;
if parms.nfiles ~= length(parms.stim_fnames)
  error('number of input files (%d) does not match number of sets of stimulus files (%d)',...
    parms.nfiles,length(parms.stim_fnames));
end;
parms.nconds = [];

for f=1:parms.nfiles
  stim_fnames = parms.stim_fnames{f};
  if ~iscell(stim_fnames), stim_fnames = {stim_fnames}; end;
  if isempty(parms.nconds)
    parms.nconds = length(stim_fnames);
  elseif parms.nconds ~= length(stim_fnames)
    error('number of conditions (%d) does not match number of conditions (%d) for run %d',...
      parms.nconds,length(stim_fnames),f);
  end;
  for c=1:parms.nconds
    if ~exist(stim_fnames{c},'file')
      error('file %s not found',stim_fnames{c});
    end;
  end;
  if isempty(parms.stim_labels)
    for c=1:parms.nconds
      parms.stim_labels{c} = sprintf('cond%d',c);
    end;
  end;
end;

parms.nstims = parms.nconds + parms.nmotion; % motion regressors
if parms.deriv_flag 
    parms.nstims = parms.nconds + 2*parms.nmotion;
end 

if isempty(parms.glt_fnames)
  outdir = sprintf('%s/gltmat',parms.outdir);
  mmil_mkdir(outdir);
  outstem = sprintf('%s/contrast',outdir);
  parms.nruns = parms.nfiles;
  if parms.deriv_flag 
    parms.motion_flag = 2; %(parms.nmotion>0);
  end 
  args = mmil_parms2args(parms,parms.glt_tags);
  [parms.glt_fnames,parms.glt_labels] = mmil_create_glt_files(outstem,args{:});
else
  if isempty(parms.glt_labels)
    for g=1:length(parms.glt_fnames)
      parms.glt_labels{g} = sprintf('glt%d',g);
    end;
  end;
end;
parms.nglt = length(parms.glt_fnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct output file names, check if they exist

fprintf('%s: checking output..\n',mfilename);

fname_out = sprintf('%s/%s_3dDeconv%s',...
  parms.outdir,parms.outstem,parms.out_ext);
fname_xmat = sprintf('%s/%s_3dDeconv.xmat.1D',...
  parms.outdir,parms.outstem);
fname_info = sprintf('%s/%s_3dDeconv.mat',...
  parms.outdir,parms.outstem);

missing_iresp_flag = 0;
if parms.iresp_flag
  fnames_iresp = cell(1,parms.nconds);
  for c=1:parms.nconds
    fnames_iresp{c} = sprintf('%s/%s_3dDeconv_iresp_%s%s',...
      parms.outdir,parms.outstem,parms.stim_labels{c},parms.out_ext);
    if ~exist(fnames_iresp{c},'file'), missing_iresp_flag = 1; end;
  end;
else
  fnames_iresp = [];
end;

if exist(fname_out,'file') &&...
   exist(fname_xmat,'file') && exist(fname_info,'file') &&...
   ~missing_iresp_flag && ~parms.forceflag
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: getting data file info..\n',mfilename);

% find number of TRs for each file
numTRs = zeros(1,parms.nfiles);
for f=1:parms.nfiles
  fname = parms.fnames{f};
  [M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);
  numTRs(f) = volsz(4);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: motion and censoring..\n',mfilename);

% if skipTRs>0, strip that many entries from start of motion file(s)
if parms.skipTRs>0 && parms.nmotion
  for f=1:parms.nfiles
    fname = parms.fnames{f};
    fname_motion = parms.fname_motion{f};
    motion_data = load(fname_motion);
    num_motion_TRs = size(motion_data,1);
    if num_motion_TRs~=numTRs(f)
      error('length of motion file %s (%d) does not match number of TRs (%d) in %s',...
        fname_motion,num_motion_TRs,numTRs(f),fname);
    end;
    if parms.skipTRs<numTRs
      [fpath,fstem,fext] = fileparts(fname_motion);
      tmp_fname_motion = sprintf('%s/scan%d_%s_skip%dTRs%s',...
        parms.outdir,f,fstem,parms.skipTRs,fext);
      write_motion_data(motion_data,tmp_fname_motion,parms.skipTRs, parms.forceflag);
      parms.fname_motion{f} = tmp_fname_motion;
      tmp_fname_censor = sprintf('%s/scan%d_%s_censor_skip%dTRs%s',...
        parms.outdir,f,fstem,parms.skipTRs,fext);
      write_censor_data(fname_motion,tmp_fname_censor,parms.skipTRs,...
          parms.censor_flag,parms.censor_thresh,...
          parms.contigous_censor_nTRs,parms.forceflag);
      tmp_fname_deriv = sprintf('%s/scan%d_%s_deriv_skip%dTRs%s',...
        parms.outdir,f,fstem,parms.skipTRs,fext);
      write_deriv_data(fname_motion,tmp_fname_deriv,parms.skipTRs,...
          parms.deriv_flag,parms.forceflag);
      parms.fname_motion_deriv{f} = tmp_fname_deriv;
    end;
  end;
end;

fprintf('%s: preparing stimulus files..\n',mfilename);

% if nfiles>1, concat motion and stim files
if parms.nfiles>1
  if parms.nmotion
    tmp_fname_motion = [parms.outdir '/concat_motion.1D'];
    tmp_fname_censor = [parms.outdir '/concat_censor.1D'];
    tmp_fname_deriv = [parms.outdir '/concat_deriv.1D'];
    if ~exist(tmp_fname_motion,'file')
      motion_data = [];
      for f=1:parms.nfiles
        tmp_motion_data = load(parms.fname_motion{f});
        motion_data = cat(1,motion_data,tmp_motion_data);
      end;
      write_motion_data(motion_data,tmp_fname_motion,0,parms.forceflag);
    end;
    write_censor_data(tmp_fname_motion,tmp_fname_censor,0,...
        parms.censor_flag,parms.censor_thresh,...
        parms.contigous_censor_nTRs,parms.forceflag);
    write_deriv_data(tmp_fname_motion,tmp_fname_deriv,0,...
        parms.deriv_flag,parms.forceflag);
    parms.fname_motion = tmp_fname_motion;
    parms.fname_motion_deriv = tmp_fname_deriv;
  end;
  
  tmp_stim_fnames = cell(parms.nconds,1);
  for c=1:parms.nconds
    outdir = sprintf('%s/stim',parms.outdir);
    mmil_mkdir(outdir);
    if parms.stim_times_flag == 1  
      tmp_fname_stim = sprintf('%s/concat_stim_%s.txt',...
        outdir,parms.stim_labels{c});
    elseif parms.stim_times_flag == 2
         tmp_fname_stim = sprintf('%s/concat_stim_%s_block.txt',...
        outdir,parms.stim_labels{c}); 
         tmp_fname_stim_dur = sprintf('%s/concat_stim_%s_block_dur.txt',...
        outdir,parms.stim_labels{c});
    else
        tmp_fname_stim = sprintf('%s/concat_stim_%s.1D',...
        outdir,parms.stim_labels{c});
    end;
    if ~exist(tmp_fname_stim,'file')
      if parms.stim_times_flag == 0
        stim_vec = [];
        for f=1:parms.nfiles
          fname_stim = parms.stim_fnames{f}{c} 
          tmp_stim_vec = load(fname_stim);
          num_stim_TRs = length(tmp_stim_vec);
          if num_stim_TRs < numTRs(f) - parms.skipTRs
            error('length of stim file %s (%d) is less than number of TRs (%d) in %s',...
              fname_stim,num_stim_TRs,numTRs(f),fname);
          elseif num_stim_TRs > numTRs(f) - parms.skipTRs
            fprintf('%s: removing %d extra TRs from stim file %s...\n',...
              mfilename,num_stim_TRs - numTRs(f),fname_stim);
            tmp_stim_vec = tmp_stim_vec(1:numTRs(f));
          end;
          stim_vec = cat(1,stim_vec,tmp_stim_vec);
        end;
        write_stim_vec(stim_vec,tmp_fname_stim,parms.forceflag);
      elseif parms.stim_times_flag == 1
        stim_mat = cell(parms.nfiles,1);
        for f=1:parms.nfiles
          fname_stim = parms.stim_fnames{f}{c}; 
          try
            tmp_stim_vec = load(fname_stim);
          catch
            fid = fopen(fname_stim,'rt');
            if fid<0
              error('failed to open %s for reading',fname_stim);
            end;
            tmp_stim_vec = fscanf(fid,'%s');
          end
          stim_mat{f} = tmp_stim_vec;
        end;
        write_stim_mat(stim_mat,tmp_fname_stim,parms.forceflag);
      else
        stim_mat = cell(parms.nfiles,1);
        stim_mat_dur = cell(parms.nfiles,1); 
        for f=1:parms.nfiles
          fname_stim = parms.stim_fnames{f}{c};   
          [block_filepath, block_name, block_ext] = fileparts(fname_stim); 
          fname_stim_dur = sprintf('%s/%s_dur.txt',block_filepath,block_name);
          try
            tmp_stim_vec = load(fname_stim);
            tmp_stim_vec_dur = load(fname_stim_dur); 
          catch
            fid = fopen(fname_stim,'rt');
            if fid<0
              error('failed to open %s for reading',fname_stim);
            end;
            tmp_stim_vec = fscanf(fid,'%s');
            fid = fopen(fname_stim_dur,'rt');
            if fid<0
              error('failed to open %s for reading',fname_stim_dur);
            end;
            tmp_stim_vec_dur = fscanf(fid,'%s');
          end
          stim_mat{f} = tmp_stim_vec; 
          stim_mat_dur{f} = tmp_stim_vec_dur;  
        end;    
        if ~std(cell2mat(stim_mat_dur))
          write_stim_mat(stim_mat,tmp_fname_stim,parms.forceflag);
          write_stim_mat(stim_mat_dur(1),tmp_fname_stim_dur,parms.forceflag)
        else 
          error('blocks have different lenghts among the runs, BLOCK is not supported yet.Check file: %s',...
            fname_stim);
        end
      end %parms.stim_times_flag 
    end
    tmp_stim_fnames{c} = tmp_fname_stim;
  end
  parms.stim_fnames = tmp_stim_fnames;
else
  if parms.nmotion
    parms.fname_motion = parms.fname_motion{1};
    parms.fname_motion_deriv = parms.fname_motion_deriv{1};
  end;
  parms.stim_fnames = parms.stim_fnames{1};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert files to nii
fprintf('%s: preparing data..\n',mfilename);
for f=1:parms.nfiles
  fname = parms.fnames{f};

  % preprocess: skip TRs, normalize
  if parms.skipTRs>0 || parms.norm_flag
    if parms.skipTRs && parms.norm_flag
      fprintf('%s: removing %d dummy TRs from data, normalizing...\n',...
        mfilename,parms.skipTRs);
    elseif parms.skipTRs
      fprintf('%s: removing %d dummy TRs from data...\n',...
        mfilename,parms.skipTRs);
    else
      fprintf('%s: normalizing data...\n',mfilename);
    end;
    [fpath,fstem,fext] = fileparts(fname);
    tmp_fname = sprintf('%s/scan%d_%s.mgh',parms.outdir,f,fstem);
    tmp_fname = mmil_normdetrend(fname,'fname_out',tmp_fname,...
      'skipTRs',parms.skipTRs,...
      'detrend',0,'norm_flag',parms.norm_flag,...
      'forceflag',parms.forceflag);
    fname = tmp_fname;
    [fpath,fstem,fext] = fileparts(fname);
  else 
      % 
  end;
  
  % convert to nii
  tmp_fname = sprintf('%s/%s.nii',parms.outdir,fstem);
  fs_mri_convert(fname,tmp_fname,'forceflag',parms.forceflag,'TR',parms.TR);
  
  parms.fnames{f} = tmp_fname;
end;

% delete output files if they exist already
tmp_fstem_out = sprintf('%s/%s_3dDeconv',parms.outdir,parms.outstem);
tmp_fname_out = sprintf('%s+orig.BRIK',tmp_fstem_out);
tmp_fname_hdr = sprintf('%s+orig.HEAD',tmp_fstem_out);
if exist(tmp_fname_out,'file') || exist(tmp_fname_hdr,'file')
  delete(tmp_fname_out);
  delete(tmp_fname_hdr);
end;

% delete iresp files if they exist already
if parms.iresp_flag
  for c=1:parms.nconds
    tmp_iresp_fstem_out = sprintf('%s/%s_3dDeconv_iresp_%s',...
      parms.outdir,parms.outstem,parms.stim_labels{c});
    tmp_iresp_fname_out = sprintf('%s+orig.BRIK',tmp_iresp_fstem_out);
    tmp_iresp_fname_hdr = sprintf('%s+orig.HEAD',tmp_iresp_fstem_out);
    if exist(tmp_iresp_fname_out,'file') || exist(tmp_iresp_fname_hdr,'file')
      delete(tmp_iresp_fname_out);
      delete(tmp_iresp_fname_hdr);
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: creating shell script..\n',mfilename);

% create shell script to run 3dDeconvolve
deconvscript = sprintf('%s/%s-deconv.sh',parms.outdir,parms.outstem);
fid = fopen(deconvscript,'wt');
if fid==-1
  error('unable to open file %s for writing');
end;
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'3dDeconvolve \\\n');
fprintf(fid,'  -input \\\n');
for f=1:parms.nfiles
  fprintf(fid,'    ''%s'' \\\n',parms.fnames{f});
end;
fprintf(fid,'  \\\n');
fprintf(fid,'  -polort %d \\\n',parms.detrend);
if parms.reml_flag
  fprintf(fid,'  -x1D_stop \\\n');
end
if parms.censor_flag 
  fprintf(fid,'  -censor %s \\\n',tmp_fname_censor);
end 
fprintf(fid,'  -num_stimts %d \\\n',parms.nstims);
fprintf(fid,'  -allzero_OK \\\n');
fprintf(fid,'  -local_times \\\n');
s=1;

for c=1:parms.nconds
  if parms.stim_times_flag
    switch parms.stim_times_model
      case 'GAM'
        stim_model = parms.stim_times_model;
      case 'SPMG'
        stim_model = parms.stim_times_model;  
      case {'TENT','CSPLIN'}
        stim_model = sprintf('%s(%d,%d,%d)',...
                       parms.stim_times_model,parms.minlag,parms.maxlag,parms.stim_times_nbasis);
      case 'BLOCK'
        [block_filepath,block_name,block_ext] = fileparts(parms.stim_fnames{c}); 
        block_dur_path = sprintf('%s/%s_dur.txt',block_filepath,block_name);
        block_dur = load(block_dur_path); 
        stim_model =  sprintf('BLOCK(%0.4f)',block_dur);  
      otherwise
        error('unsupported stim times model');
    end;
    fprintf(fid,'  -stim_times %d ''%s'' ''%s'' -stim_label %d %s \\\n',...
      s,parms.stim_fnames{c},stim_model,s,parms.stim_labels{c});
  else
    fprintf(fid,'  -stim_file %d ''%s''  -stim_label %d %s \\\n',...
      s,parms.stim_fnames{c},s,parms.stim_labels{c});
    fprintf(fid,'  -stim_minlag %d %d  -stim_maxlag %d %d \\\n',...
      s,parms.minlag,s,parms.maxlag);
  end;
  if parms.iresp_flag
    fprintf(fid,'  -iresp %d %s_iresp_%s \\\n',...
      s,tmp_fstem_out,parms.stim_labels{c});
  end;
  s=s+1;
end;

for m=1:parms.nmotion
  fprintf(fid,'  -stim_file %d ''%s[%d]'' -stim_label %d %s -stim_base %d \\\n',...
    s,parms.fname_motion,m,s,parms.motion_labels{m},s);
  s=s+1;
end;

if parms.deriv_flag 
  for m=1:parms.nmotion
    mm1 = m-1;
    fprintf(fid,'  -stim_file %d ''%s[%d]'' -stim_label %d %s_d -stim_base %d \\\n',...
    s,parms.fname_motion_deriv,mm1,s,parms.motion_labels{m},s);
    s=s+1;
  end;
end 

if parms.nglt>0
  fprintf(fid,'  -num_glt %d \\\n',parms.nglt);
  for g=1:parms.nglt
    fprintf(fid,'  -glt 1 ''%s'' -glt_label %d "%s" \\\n',...
      parms.glt_fnames{g},g,parms.glt_labels{g});
  end;
end;
 
% number of stims to chop 
if parms.reml_flag 
  switch parms.stim_times_model
    case 'GAM'
      reml_stims = size(parms.stim_labels,2)*2+1;
    case 'SPMG'
      reml_stims = size(parms.stim_labels,2)*3+1;
    otherwise
      error('unsupported stim times model with REML');
  end;
end 

if parms.reml_flag 
  fprintf(fid,'  -xjpeg ''%s'' \\\n',tmp_fstem_out);
  fprintf(fid,'  -x1D ''%s'' \n',tmp_fstem_out);
  fprintf(fid,'\n');
  fprintf(fid,'3dREMLfit -matrix ''%s.xmat.1D'' \\\n',tmp_fstem_out);
  fprintf(fid,'  -input \\\n');
  fprintf(fid,'  ''');
  for f=1:parms.nfiles
    fprintf(fid,'%s ',parms.fnames{f});
  end;
  fprintf(fid,'  ''');
  fprintf(fid,'  \\\n');
  fprintf(fid,'  -fout \\\n');
  fprintf(fid,'  -Rbuck ''%s_stims'' \\\n',tmp_fstem_out);
  if parms.goforit_flag
    fprintf(fid,'  -verb \\\n');
    fprintf(fid,'  -GOFORIT ''%s''\n',num2str(parms.goforit_val));
  else 
    fprintf(fid,'  -verb \n');
  end 
  fprintf(fid,'\n');
  fprintf(fid,'3dTcat -prefix ''%s'' \\\n',tmp_fstem_out);
  fprintf(fid,'''%s_stims+orig[0,%d..$]'' \n',tmp_fstem_out,reml_stims);
else
  fprintf(fid,'  -nocout -tout -float \\\n');
  fprintf(fid,'  -xjpeg ''%s'' \\\n',tmp_fstem_out);
  fprintf(fid,'  -x1D ''%s'' \\\n',tmp_fstem_out);
  if parms.goforit_flag
    fprintf(fid,'  -bucket ''%s'' \\\n',tmp_fstem_out);
    fprintf(fid,'  -GOFORIT ''%s''\n',num2str(parms.goforit_val));
  else 
    fprintf(fid,'  -bucket ''%s''\n',tmp_fstem_out);
  end 
end

fclose(fid);

% run script
fprintf('%s: running 3dDeconvolve with %s...\n',mfilename,deconvscript);
[status,result] = unix(sprintf('source %s',deconvscript));
if status
  error('failed to run shell script %s:\n%s',deconvscript,result);
else
  disp(result);
end;

% read HEAD file to get labels and degrees of freedom
stats_info = mmil_read_HEAD_statsinfo(tmp_fname_out);
save(fname_info,'stats_info');

% convert output to mgh
mmil_BRIK2mgh(tmp_fname_out,fname_out);

% convert iresp
if parms.iresp_flag
  for c=1:parms.nconds
    tmp_iresp_fstem_out = sprintf('%s/%s_3dDeconv_iresp_%s',...
      parms.outdir,parms.outstem,parms.stim_labels{c});
    tmp_iresp_fname_out = sprintf('%s+orig.BRIK',tmp_iresp_fstem_out);
    mmil_BRIK2mgh(tmp_iresp_fname_out,fnames_iresp{c});
  end;
end;

% remove REML_cmd file
fname_reml_cmd = sprintf('%s.REML_cmd',regexprep(tmp_fstem_out,'\..+',''));
if exist(fname_reml_cmd,'file')
  delete(fname_reml_cmd);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_motion_data(motion_data,fname_motion,skipTRs,forceflag)
  if ~exist(fname_motion,'file') || forceflag
    numTRs = size(motion_data,1);
    ncols = size(motion_data,2);
    if ncols==7
      motion_data = cat(2,motion_data,zeros(numTRs,2));
    end;
    fid = fopen(fname_motion,'w');
    k = 1;
    for t=skipTRs+1:numTRs
      fprintf(fid,'  %d  %s %s\n',k,...
        sprintf('%0.4f  ',squeeze(motion_data(t,2:7))),...
        sprintf('    %0.1f',squeeze(motion_data(t,8:9))));
      k = k + 1;
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_censor_data(fname_in,fname_out,skipTRs,censor_flag,thrs,contigous_trs,forceflag)

 if ~exist(fname_out,'file') || forceflag || censor_flag 
      
   motion_data = mmil_load_motion_1D(fname_in);  
   [stats, fd] = mmil_motion_stats(motion_data,[],[],[],1);
   numTRs = size(motion_data,1);
   fd_x = fd(skipTRs+1:numTRs); 
   censor = logical(fd_x > thrs);  
   
   % code from Eric Earl, earl@ohsu.edu
   contiguous = zeros(size(censor)); 
   contiguous_groups = bwlabel(logical((censor-1)*-1)); 
   for group = 1:max(contiguous_groups)
     if sum(contiguous_groups == group) < contigous_trs
        contiguous(contiguous_groups == group) = true;
    else
        contiguous(contiguous_groups == group) = false;
     end
   end
   double_censor = ~logical(censor + contiguous);      
        
   fid = fopen(fname_out,'w');
   fprintf(fid,'%d\n',double_censor);
   fclose(fid);
 end;
 
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function write_deriv_data(fname_in,fname_out,skipTRs,deriv_flag,forceflag)
 
 if ~exist(fname_out,'file') || forceflag || deriv_flag 
     
   motion_data = mmil_load_motion_1D(fname_in);
   [stats, fd, deriv] = mmil_motion_stats(motion_data,[],[],[],0);
   deriv = [zeros(1,6) ; deriv];
   numTRs = size(motion_data,1);
   deriv_x = deriv(skipTRs+1:numTRs,:);
   
   fid = fopen(fname_out,'w');
   for i=1:length(deriv_x)
     fprintf(fid,'  %s\n',sprintf('%0.4f  ',deriv_x(i,:)));
   end 
   fclose(fid);
 end;
 
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stim_vec(stim_vec,fname_stim,forceflag)
  if ~exist(fname_stim,'file') || forceflag
    fid = fopen(fname_stim,'w');
    fprintf(fid,'%d\n',stim_vec);
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stim_mat(stim_mat,fname_stim,forceflag)
  if ~exist(fname_stim,'file') || forceflag
    fid = fopen(fname_stim,'w');
    for i=1:length(stim_mat)
      stim_vec = stim_mat{i};
      if isempty(stim_vec) || strcmp(stim_vec,'*')
        fprintf(fid,'*\n');
      else
        fprintf(fid,'%s\n',sprintf('%0.3f ',stim_vec));
      end;
    end;
    fclose(fid);
  end;
return;

