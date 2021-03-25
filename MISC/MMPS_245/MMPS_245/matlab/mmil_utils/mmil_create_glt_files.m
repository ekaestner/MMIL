function [glt_fnames,glt_labels] = mmil_create_glt_files(outstem,varargin)
%function [glt_fnames,glt_labels] = mmil_create_glt_files(outstem,[options])
%
% Purpose: create general linear test files for mmil_3dDeconv
%
% Required Input:
%   outstem: output file stem
%
% Optional Parameters:
%  'nruns': number of separate runs
%     {default = 1}
%  'nconds': number of conditions
%     if stim_labels supplied, nconds determined from number of elements
%     {default = 1}
%  'contrasts_flag': whether to create contrasts for combinations of conditions
%     {default = 0}
%  'stim_labels': string or cell array of strings with stimulus condition names
%    If empty, will use labels such as "cond1", "cond2", etc.
%    {default = []}
%  'stim_times_flag': [0|1|2] use stimulus time files
%    0 = 1D, 1 = txt and 2 = block 
%    {default = 1}
%  'stim_times_model': name of stimulus model
%    model names in 3dDeconvolve for -stim_times
%    allowed: 'SPMG','GAM','TENT','CSPLIN', 'BLOCK'
%    {default = 'SPMG'}
%  'stim_times_nbasis': max number of basis functions for HRF
%    {default = 10}
%  'out_ext': output file extension
%    {default = '.gltmat'}
%  'minlag': number of TRs for minimum "lag" between stimulus and response
%    {default = 0}
%  'maxlag': number of TRs for maximum "lag" between stimulus and response
%    {default = 4}
%  'detrend': [0|1|2|3] number of zeros include for detrending
%    {default = 1}
%  'motion_flag' [0|1|2] include spaces for motion regressors
%    0 : no spaces included for motion
%    1 : 6 spaces (motion estimates)
%    2: 12 spaces (motion estimates + derivatives)
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  11/01/12 by Don Hagler
% Prev Mod: 09/06/17 by Don Hagler 
% Last Mod: 10/24/17 by Dani Cornejo 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% based on code from mmil_3dDeconv, created 09/01/08 by Don Hagler

if ~mmil_check_nargs(nargin, 1), return; end;
parms = mmil_args2parms(varargin, { ...
  'nruns',1,[1,1000],...
  'nconds',1,[1,1000],...
  'contrasts_flag',false,[false true],...
  'stim_labels',[],[],...
  'stim_times_flag',1,0:2,...
  'stim_times_model','SPMG',{'SPMG','TENT','CSPLIN','GAM','BLOCK'},...
  'stim_times_nbasis',10,[1,100],...
  'out_ext','.gltmat',[],...
  'minlag',0,[0,10],...
  'maxlag',4,[0,30],...
  'detrend',1,[0:3],...
  'motion_flag',1,[0:2],...
  'forceflag',false,[false true],...
...
  'conds_contrast',[],[],...
});

glt_fnames = [];
glt_labels = [];

% determine how many parameters for each condition
if parms.stim_times_flag
  switch parms.stim_times_model
    case 'SPMG'
      nbasis = 2;
      ncoef = 1;
    case 'GAM'
      nbasis = 1;
      ncoef = nbasis;
    case 'BLOCK'      
      nbasis = 1;
      ncoef = nbasis;
    case {'TENT','CSPLIN'}
      nbasis = parms.stim_times_nbasis;
      ncoef = nbasis;
  end;
  parms.minlag = 1;
  parms.maxlag = nbasis;
else
  ncoef = parms.maxlag;
end;
parms.nlags = parms.maxlag - parms.minlag + 1;
parms.nregs = parms.nconds * parms.nlags;

% prepare general linear tests for each condition
%   (area under curve of hemodynamic response)
if ~isempty(parms.stim_labels)
  if ~iscell(parms.stim_labels), parms.stim_labels = {parms.stim_labels}; end;
  parms.nconds = length(parms.stim_labels);
end;

switch parms.motion_flag
  case 0
    parms.nmotion = 0;    
  case 1
    parms.nmotion = 6;    
  case 2
    parms.nmotion = 12;            
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare general linear tests for each condition
g = 1;
for c=1:parms.nconds
  if isempty(parms.stim_labels)
    stim_label = sprintf('cond%d',c);
  else
    stim_label = parms.stim_labels{c};
  end;
  fname_out = sprintf('%s_%s%s',outstem,stim_label,parms.out_ext);
  if ~exist(fname_out,'file') || parms.forceflag
    fid = fopen(fname_out,'wt');
    if fid==-1
      error('failed to open file %s for writing',fname_out);
    end;
    % if multiple runs, need separate baseline regs for each
    for r=1:parms.nruns
      for p=0:parms.detrend % number of baseline regs depends on detrend (detrend)
        fprintf(fid,'0 ');
      end;
    end;
    % a '1' for each "lag" for this particular condition
    for tmpc=1:parms.nconds
      for lag = parms.minlag:parms.maxlag
        if tmpc==c && lag<=ncoef
          val = 1;
        else
          val = 0;
        end;
        fprintf(fid,'%d ',val);
      end;
    end;
    % 6 extra zeros at end for motion regressors
    if parms.motion_flag
      for i=1:parms.nmotion
        fprintf(fid,'0 ');
      end;
    end;
    fprintf(fid,'\n');
    fclose(fid);
  end;
  glt_fnames{g} = fname_out;
  glt_labels{g} = stim_label;
  g=g+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare general linear tests for each contrast
%   (one condition vs. another)
if parms.contrasts_flag && parms.nconds>1
  for c1=1:parms.nconds
    if isempty(parms.stim_labels)
      stim_label1 = sprintf('cond%d',c1);
    else
      stim_label1 = parms.stim_labels{c1};
    end;
    for c2=1:parms.nconds
      if c1==c2, continue; end;
      if isempty(parms.stim_labels)
        stim_label2 = sprintf('cond%d',c2);
      else
        stim_label2 = parms.stim_labels{c2};
      end;
      fname_out = sprintf('%s_%sVS%s%s',outstem,...
        stim_label1,stim_label2,parms.out_ext);
      if ~exist(fname_out,'file') || parms.forceflag
        fid = fopen(fname_out,'wt');
        if fid==-1
          error('failed to open file %s for writing',fname_out);
        end;
        % if multiple runs, need separate baseline regs for each
        for r=1:parms.nruns
          for p=0:parms.detrend % number of baseline regs depends on detrend
            fprintf(fid,'0 ');
          end;
        end;
        % a '1' for each "lag" for this particular condition
        for tmpc=1:parms.nconds
          for lag = parms.minlag:parms.maxlag
            if tmpc==c1 && lag<=ncoef
              val = 1;
            elseif tmpc==c2 && lag<=ncoef
              val = -1;
            else
              val = 0;
            end;
            fprintf(fid,'%d ',val);
          end;
        end;
        % 6 extra zeros at end for motion regressors
        if parms.motion_flag
          for i=1:parms.nmotion
            fprintf(fid,'0 ');
          end;
        end;
        fprintf(fid,'\n');
        fclose(fid);
      end;
      glt_fnames{g} = fname_out;
      glt_labels{g} = sprintf('%sVS%s',stim_label1,stim_label2);
      g=g+1;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare general linear tests for contrasts between combinations of conditions (condition names)
if ~isempty(parms.conds_contrast) && ~isempty(parms.stim_labels)
  ncontrasts = length(parms.conds_contrast);
  for i=1:ncontrasts
    groupA = parms.conds_contrast(i).groupA;
    groupB = parms.conds_contrast(i).groupB;
    stim_label = parms.conds_contrast(i).name;
    fname_out = sprintf('%s_%s%s',outstem,stim_label,parms.out_ext);
    if ~exist(fname_out,'file') || parms.forceflag
      fid = fopen(fname_out,'wt');
      if fid==-1
        error('failed to open file %s for writing',fname_out);
      end;
      % if multiple runs, need separate baseline regs for each
      for r=1:parms.nruns
        for p=0:parms.detrend % number of baseline regs depends on detrend
          fprintf(fid,'0 ');
        end;
      end;
      % a '1' for each "lag" for this particular condition
      for c=1:parms.nconds
        % if part of groupA, set to 1
        % if part of groupB, set to -1
        [tmp,i_condA] = intersect(parms.stim_labels,groupA);
        [tmp,i_condB] = intersect(parms.stim_labels,groupB);
        for lag = parms.minlag:parms.maxlag
          if ~isempty(i_condA) && ismember(c,i_condA) && lag<=ncoef
            val = 1/length(i_condA);
          elseif ~isempty(i_condB) && ismember(c,i_condB) && lag<=ncoef
            val = -1/length(i_condB);
          else
            val = 0;
          end;
          if val==0
            fprintf(fid,'%d ',val);
          else
            fprintf(fid,'%0.2f ',val);
          end;
        end;
      end;
      % 6 extra zeros at end for motion regressors
      if parms.motion_flag
        for i=1:parms.nmotion
          fprintf(fid,'0 ');
        end;
      end;
      fprintf(fid,'\n');
      fclose(fid);
    end;
    glt_fnames{g} = fname_out;
    glt_labels{g} = sprintf('%s',stim_label);
    g=g+1;
  end;

end;























    
    
