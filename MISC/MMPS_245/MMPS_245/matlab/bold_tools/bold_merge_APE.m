function [fname_out,fname_motion_out] = bold_merge_APE(fname_rev,fname_for,varargin)
%function [fname_out,fname_motion_out] = bold_merge_APE(fname_rev,fname_for,[options])
%
% Required Input:
%   fname_rev: full path for reverse phase-encode image file
%   fname_for: full path for forward phase-encode image file
%
% Optional Input:
%   'pepolar': phase-encode polarity
%      2 = odd frames are rev, even frames are for
%      3 = odd frames are for, even frames are rev
%     {default = 2}
%   'mc_flag' - [0|1] also merge motion.1D files created by AFNI's 3dvolreg
%     {default = 1}
%   'infix': string included at end of file names
%     used to find corresponding motion.1D files if mc_flag = 1
%     {default = []}
%   'fname_out': output file name
%     if empty, name will be generated from input
%     {default = []}
%   'fname_motion_out': output motion.1D file name
%     if empty, name will be generated from input
%     {default = []}
%   'normflag': [0|1] output volume will contain norm of
%     for and rev values (sqrt of sum of squares)
%     {default = 1}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_out: output file name
%   fname_motion_out: output motion.1D file name
%
% Created:  01/13/11 by Josh Kuperman
% Last Mod: 03/23/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'pepolar',2,[2,3],...
  'mc_flag',true,[false true],...
  'infix',[],[],...
  'fname_out',[],[]...
  'fname_motion_out',[],[]...
  'normflag',true,[false true],...
  'verbose',true,[false true],...
  'forceflag', false,[false true],...
});

fname_out = [];
fname_motion_out = [];

[tpath_rev,tstem_rev,text_rev] = fileparts(fname_rev);
[tpath_for,tstem_for,text_for] = fileparts(fname_for);

if isempty(parms.fname_out)
  if ~isempty(regexp(fname_rev,'rev'))
    parms.fname_out = regexprep(fname_rev,'rev','ape');
  else
    parms.fname_out = [tpath_rev '/' tstem_rev '_ape.mgz'];
  end;
end;

if parms.mc_flag
  pathlist = {tpath_rev,tpath_for};
  stemlist = {tstem_rev,tstem_for};
  for i=1:length(stemlist)
    tpath = pathlist{i};
    tstem = stemlist{i};
    fname_motion = [tpath '/' tstem '_motion.1D'];
    if ~exist(fname_motion,'file')
      fprintf('%s: WARNING: motion 1D file %s not found\n',...
        mfilename,fname_motion);
      fname_motion = [];
    end;
    if i==1
      parms.fname_motion_rev = fname_motion;
    else
      parms.fname_motion_for = fname_motion;
    end;
  end;
  if isempty(parms.fname_motion_out) && ~isempty(parms.fname_motion_rev)
    if ~isempty(regexp(parms.fname_motion_rev,'rev'))
      parms.fname_motion_out = regexprep(parms.fname_motion_rev,'rev','ape');
    else
      parms.fname_motion_out = [tpath_rev '/' tstem_rev '_ape_motion.1D'];
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% merge rev and for image volumes
if ~exist(parms.fname_out) || parms.forceflag
  if parms.verbose
    fprintf('%s: merging data files %s and %s...\n',...
      mfilename,fname_rev,fname_for);
  end;
  [vol_rev,M_rev,mr_parms_rev,volsz_rev] = fs_load_mgh(fname_rev);
  [vol_for,M_for,mr_parms_for,volsz_for] = fs_load_mgh(fname_for);
  if any(M_rev(:)~=M_for(:)) || any(volsz_rev(1:3)~=volsz_for(1:3))
    error('mismatch between rev and for volumes');
  end;
  nf_rev = size(vol_rev,4);
  nf_for = size(vol_for,4);
  if parms.normflag
    nf_ape = min(nf_rev,nf_for);
    vol_ape = hypot(vol_rev(:,:,:,1:nf_ape),...
                    vol_for(:,:,:,1:nf_ape));
  else
    nf_ape = nf_rev + nf_for;
    switch parms.pepolar
      case 2
        ind_rev = [1:2:nf_ape]; % odd frames are "reverse"
        ind_for = [2:2:nf_ape]; % even frames are "forward"
      case 3
        ind_for = [1:2:nf_ape]; % odd frames are "forward"
        ind_rev = [2:2:nf_ape]; % even frames are "reverse"
    end;
    volsz_ape = [volsz_rev(1:3),nf_ape];
    vol_ape = zeros(volsz_ape);
    vol_ape(:,:,:,ind_rev) = vol_rev;
    vol_ape(:,:,:,ind_for) = vol_for;
  end;
  fs_save_mgh(vol_ape,parms.fname_out,M_rev,mr_parms_rev);
end;

% merge rev and for motion.1D files
if parms.mc_flag && (~exist(parms.fname_motion_out) || parms.forceflag)
  if ~isempty(parms.fname_motion_rev) && ~isempty(parms.fname_motion_for)
    if parms.verbose
      fprintf('%s: merging motion files %s and %s...\n',...
        mfilename,parms.fname_motion_rev,parms.fname_motion_for);
    end;
    motion_data_rev = mmil_load_motion_1D(parms.fname_motion_rev);
    motion_data_for = mmil_load_motion_1D(parms.fname_motion_for);
    nmotion = size(motion_data_rev,2);
    nf_rev = size(motion_data_rev,1);
    nf_for = size(motion_data_for,1);
    if parms.normflag
      nf_ape = min(nf_rev,nf_for);
      motion_data_ape = (motion_data_rev(1:nf_ape,:) +...
                         motion_data_for(1:nf_ape,:))/2;
    else
      nf_ape = nf_rev + nf_for;
      switch parms.pepolar
        case 2
          ind_rev = [1:2:nf_ape]; % odd frames are "reverse"
          ind_for = [2:2:nf_ape]; % even frames are "forward"
        case 3
          ind_for = [1:2:nf_ape]; % odd frames are "forward"
          ind_rev = [2:2:nf_ape]; % even frames are "reverse"
      end;
      motion_data_ape = zeros(nf_ape,nmotion);
      motion_data_ape(ind_rev,:) = motion_data_rev;
      motion_data_ape(ind_for,:) = motion_data_for;
    end;
    fid = fopen(parms.fname_motion_out,'w');
    if fid==-1
      error('failed to open %s for writing',parms.fname_motion_out);
    end;
    for t=1:nf_ape
      fprintf(fid,'%4d%s\n',t-1,sprintf('%10.4f',motion_data_ape(t,:)));
    end;
    fclose(fid);
  end;
end;

fname_out = parms.fname_out;
fname_motion_out = parms.fname_motion_out;

