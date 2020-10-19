function [ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,varargin)
%function [ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,[options])
%
% Purpose: load BOLD scan information, select reference scans
%   for between-scan registration and B0 unwarping, and select valid
%   scans for processing
%
% Required Input:
%   ContainerPath: full path of BOLDPROC Container
%
% Optional Parameters:
%   'snums': list of scan numbers to process
%     if empty, use all BOLD scans in container
%     {default = []}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     0: use only forward phase-encode polarity data
%     1: use only reverse phase-encode polarity data
%     2: use both forward and reverse data
%     {default = 2}
%
% Output:
%   ScanInfo: struct array containing info about each BOLD scan
%     from ContainerInfo
%   SessInfo: struct containing summary info about session, including:
%     revflag (1 if most scans are rev, 0 otherwise)
%     MagneticFieldStrength
%     nscans
%     snums_valid
%     snums_for
%     snums_rev
%     B0uw_refs_for
%     B0uw_refs_rev
%     B0reg_refs_for
%     B0reg_refs_rev
%     reg_ref_for
%     reg_ref_rev
%     regT1_ref
%
%   errcode: 0 if successful, 1 if error
%
% Created:  04/29/10 by Don Hagler
% Last Mod: 07/23/17 by Don Hagler
%

%% todo: check that reference scans have same EchoSpacing, AcquisitionRows,
%%         and AcquisitionColumns as scan for which they serve as reference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScanInfo = []; SessInfo = []; errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'revflag',2,[0,1,2],...
...
  'ContainerPath',ContainerPath,[],...
  'fnamestem','BOLD',[],...
  'copy_tags',{'VisitID','StudyDate','StudyTime','StudyInstanceUID',...
    'MagneticFieldStrength','Manufacturer',...
    'ManufacturersModelName','MagneticFieldStrength'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ContainerInfo, get info for BOLD scans, and identify scans of each type
[ScanInfo,SessInfo,errcode] = get_info(parms);
if errcode | isempty(ScanInfo), return; end;

% set revflag to indicate which scan type is more prevalent
SessInfo = set_revflag(ScanInfo,SessInfo,parms);

% choose pairs of forward and reverse scans for B0 unwarp
SessInfo = set_B0uw_refs(SessInfo);

% choose reference scans for pre-B0unwarp inter-scan registration
SessInfo = set_B0reg_refs(SessInfo);

% select references for post-B0unwarp inter-scan registration
SessInfo = set_reg_refs(SessInfo);

% select reference for registration to T1
SessInfo = set_regT1_ref(SessInfo);

% select B0uw and B0reg reference scans for each scan
ScanInfo = set_scan_refs(ScanInfo,SessInfo);

% set file stems for each scan according to type of scan
[ScanInfo,SessInfo] = set_scan_fstems(ScanInfo,SessInfo,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ContainerInfo, get info for BOLD scans, and identify scans of each type
function [ScanInfo,SessInfo,errcode] = get_info(parms)
  ScanInfo = []; SessInfo = []; errcode = 0;

  % load ContainerInfo
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(parms.ContainerPath);
  if errcode~=0, return; end;

  % get info for BOLD scans
  ScanInfo = mmil_getfield(ContainerInfo.ScanInfo,...
    parms.fnamestem,[]);
  if isempty(ScanInfo), return; end;

  % copy fields from ContainerInfo
  for t=1:length(parms.copy_tags)
    tag = parms.copy_tags{t};
    SessInfo.(tag) = ContainerInfo.(tag);
  end;

  % identify scans of each type
  SessInfo.nscans = length(ScanInfo);
  if isempty(parms.snums), parms.snums = [1:SessInfo.nscans]; end;
  SessInfo.snums_undef = [];
  SessInfo.snums_for = [];
  SessInfo.snums_rev = [];
  SessInfo.snums_SEpep = [];
  SessInfo.snums_valid = [];
  SessInfo.snums_SE_for = [];
  SessInfo.snums_SE_rev = [];
  for s=1:SessInfo.nscans
    if ~ismember(s,parms.snums), continue; end;

    % exclude scans that failed to convert
    if ~ScanInfo(s).valid
      fprintf('%s: WARNING: skipping invalid BOLD scan %d\n',...
        mfilename,s);
      continue;
    end;

    % exclude scans with empty pepolar
    ScanInfo(s).pepolar = mmil_getfield(ScanInfo(s),'pepolar',0);  
    % exclude scans without pepolar info
    if isempty(ScanInfo(s).pepolar)
      fprintf('%s: WARNING: skipping BOLD scan %d: missing pepolar info\n',...
        mfilename,s);
      continue;
    end;

    % exclude scans without ScanType
    ScanInfo(s).ScanType = mmil_getfield(ScanInfo(s),'ScanType',[]);
    if isempty(ScanInfo(s).ScanType)
      fprintf('%s: WARNING: skipping BOLD scan %d: missing ScanType info\n',...
        mfilename,s);
      continue;
    end;

    % categorize scans
    try
      switch ScanInfo(s).ScanType
        case 'BOLD_bu' % bottom-up "forward"
          SessInfo.snums_for = [SessInfo.snums_for,s];
        case 'BOLD_td' % top-down  "reverse"
          SessInfo.snums_rev = [SessInfo.snums_rev,s];
        case 'BOLD_ipp' % integrated pepolar
          if ismember(ScanInfo(s).pepolar,[0,2,3])
            SessInfo.snums_for = [SessInfo.snums_for,s];
          end;
          if ismember(ScanInfo(s).pepolar,[1,2,3])
            SessInfo.snums_rev = [SessInfo.snums_rev,s];
          end;
        case 'BOLD_ape' % alternating pepolar
          if ismember(ScanInfo(s).pepolar,[0,2,3])
            SessInfo.snums_for = [SessInfo.snums_for,s];
          end;
          if ismember(ScanInfo(s).pepolar,[1,2,3])
            SessInfo.snums_rev = [SessInfo.snums_rev,s];
          end;
        case 'BOLD_pep' % spin-echo pepolar
          SessInfo.snums_SEpep = [SessInfo.snums_SEpep,s];
        case 'BOLD_SE' %SE field maps for Siemens/Philips
          switch ScanInfo(s).pepolar
            case 0
              SessInfo.snums_SE_for = [SessInfo.snums_SE_for,s];
            case 1
              SessInfo.snums_SE_rev = [SessInfo.snums_SE_rev,s]; 
          end;
        otherwise
          SessInfo.snums_undef = [SessInfo.snums_undef,s];
      end;
    catch
      fprintf('%s: WARNING: skipping BOLD scan %d: invalid scan type:\n%s\n',...
        mfilename,parms.fnamestem,s,lasterr);
      continue;
    end;
  end;
  SessInfo.nscans_undef = length(SessInfo.snums_undef);
  SessInfo.nscans_for = length(SessInfo.snums_for);
  SessInfo.nscans_rev = length(SessInfo.snums_rev);
  SessInfo.nscans_SEpep = length(SessInfo.snums_SEpep);
  SessInfo.nscans_SE_for = length(SessInfo.snums_SE_for);
  SessInfo.nscans_SE_rev = length(SessInfo.snums_SE_rev);

  % determine set of scans valid for processing and analysis
  switch parms.revflag
    case 0
      SessInfo.snums_valid = sort([SessInfo.snums_for,...
                                   SessInfo.snums_SEpep,SessInfo.snums_SE_for]);
    case 1
      SessInfo.snums_valid = sort([SessInfo.snums_rev,...
                                   SessInfo.snums_SEpep,SessInfo.snums_SE_rev]);
    case 2
      SessInfo.snums_valid = unique([SessInfo.snums_for,SessInfo.snums_rev,...
                                    SessInfo.snums_SEpep,...
                                    SessInfo.snums_SE_for,SessInfo.snums_SE_rev]);
  end;
  if SessInfo.nscans_undef>0 && isempty(SessInfo.snums_valid)
    SessInfo.snums_valid = SessInfo.snums_undef;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set revflag to indicate which scan type is more prevalent
function SessInfo = set_revflag(ScanInfo,SessInfo,parms)
  % calculate total number of forward and reverse reps
  SessInfo.nreps_rev = 0;
  SessInfo.nreps_for = 0;
  for i=1:SessInfo.nscans_rev
    s = SessInfo.snums_rev(i);
    switch ScanInfo(s).ScanType
      case 'BOLD_ipp'
        switch ScanInfo(s).pepolar
          case 2
            SessInfo.nreps_rev = SessInfo.nreps_rev + 1;
          case 3
            SessInfo.nreps_rev = SessInfo.nreps_rev + ScanInfo(s).nreps - 1;
          otherwise
            SessInfo.nreps_rev = SessInfo.nreps_rev + ScanInfo(s).nreps;
        end;
      case 'BOLD_ape'
        SessInfo.nreps_rev = SessInfo.nreps_rev + ScanInfo(s).nreps/2;
      otherwise
        SessInfo.nreps_rev = SessInfo.nreps_rev + ScanInfo(s).nreps;
    end;
  end;
  for i=1:SessInfo.nscans_for
    s = SessInfo.snums_for(i);
    switch ScanInfo(s).ScanType
      case 'BOLD_ipp'
        switch ScanInfo(s).pepolar
          case 2
            SessInfo.nreps_for = SessInfo.nreps_for + ScanInfo(s).nreps - 1;
          case 3
            SessInfo.nreps_for = SessInfo.nreps_for + 1;
          otherwise
            SessInfo.nreps_for = SessInfo.nreps_for + ScanInfo(s).nreps;
        end;
      case 'BOLD_ape'
        SessInfo.nreps_for = SessInfo.nreps_for + ScanInfo(s).nreps/2;
      otherwise
        SessInfo.nreps_for = SessInfo.nreps_for + ScanInfo(s).nreps;
    end;
  end;

  % set revflag to indicate which scan type is more prevalent
  if ismember(parms.revflag,[0,1])
    % override with user-choice
    SessInfo.revflag = parms.revflag;
  elseif (SessInfo.nreps_rev >= SessInfo.nreps_for)
    SessInfo.revflag = 1;
  else
    SessInfo.revflag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose pairs of forward and reverse scans for B0 unwarp
function SessInfo = set_B0uw_refs(SessInfo)
  SessInfo.B0uw_refs_for = [];
  SessInfo.B0uw_refs_rev = [];
  snums_for = [];
  snums_rev = [];
  if SessInfo.nscans_SEpep>0
    snums_for = SessInfo.snums_SEpep;
    snums_rev = SessInfo.snums_SEpep;
  elseif SessInfo.nscans_SE_for>0 && SessInfo.nscans_SE_rev >0
    snums_for = SessInfo.snums_SE_for;
    snums_rev = SessInfo.snums_SE_rev;
  elseif SessInfo.nscans_for>0 && SessInfo.nscans_rev>0
    snums_for = SessInfo.snums_for;
    snums_rev = SessInfo.snums_rev;
  end;
  if ~isempty(snums_for) && ~isempty(snums_rev)
    if SessInfo.revflag
      snumsA = snums_for;
      snumsB = snums_rev;
    else
      snumsA = snums_rev;
      snumsB = snums_for;
    end;
    for i=1:length(snumsA)
      sA = snumsA(i);
      tmp_diff = abs(sA - snumsB);
      [tmp,ind] = min(tmp_diff);
      % if tie, choose later scan
      ind_close = find(tmp_diff==tmp);
      if length(ind_close)>1, ind = ind_close(2); end;
      sB = snumsB(ind);    
      if SessInfo.revflag
        SessInfo.B0uw_refs_for = [SessInfo.B0uw_refs_for,sA];
        SessInfo.B0uw_refs_rev = [SessInfo.B0uw_refs_rev,sB];
      else
        SessInfo.B0uw_refs_rev = [SessInfo.B0uw_refs_rev,sA];
        SessInfo.B0uw_refs_for = [SessInfo.B0uw_refs_for,sB];
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose reference scans for pre-B0unwarp inter-scan registration
function SessInfo = set_B0reg_refs(SessInfo)
  SessInfo.B0reg_refs_for = [];
  SessInfo.B0reg_refs_rev = [];
  for i=1:length(SessInfo.B0uw_refs_for)
    % find nearest forward scan
    if SessInfo.nscans_for>0
      tmp = SessInfo.snums_for - SessInfo.B0uw_refs_for(i);
      tmp(tmp<0) = Inf; % prefer scans after the B0uw ref
      [mindiff,ind]=min(tmp);
      if ~isinf(mindiff)
        SessInfo.B0reg_refs_for(i) = SessInfo.snums_for(ind);
      else % no scans after, so look before
        tmp = abs(SessInfo.snums_for - SessInfo.B0uw_refs_for(i));
        [mindiff,ind]=min(tmp);
        SessInfo.B0reg_refs_for(i) = SessInfo.snums_for(ind);
      end;
    end;
    % repeat for reverse scans
    if SessInfo.nscans_rev>0
      tmp = SessInfo.snums_rev - SessInfo.B0uw_refs_rev(i);
      tmp(tmp<0) = Inf; % prefer scans after the B0uw ref
      [mindiff,ind]=min(tmp);
      if ~isinf(mindiff)
        SessInfo.B0reg_refs_rev(i) = SessInfo.snums_rev(ind);
      else % no scans after, so look before
        tmp = abs(SessInfo.snums_rev - SessInfo.B0uw_refs_rev(i));
        [mindiff,ind]=min(tmp);
        SessInfo.B0reg_refs_rev(i) = SessInfo.snums_rev(ind);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select references for post-B0unwarp inter-scan registration
function SessInfo = set_reg_refs(SessInfo)
  if ~isempty(SessInfo.B0reg_refs_for)
    if SessInfo.revflag
      if SessInfo.snums_for>0
        ind = floor(median(1:length(SessInfo.B0reg_refs_for)));
      elseif SessInfo.snums_rev>0
        ind = floor(median(1:length(SessInfo.B0reg_refs_rev)));
      else
        ind = [];
      end;
    else
      if SessInfo.snums_rev>0
        ind = floor(median(1:length(SessInfo.B0reg_refs_rev)));
      elseif SessInfo.snums_for>0
        ind = floor(median(1:length(SessInfo.B0reg_refs_for)));
      else
        ind = [];
      end;
    end;
    if SessInfo.snums_rev>0
      SessInfo.reg_ref_rev = SessInfo.B0reg_refs_rev(ind);
    else
      SessInfo.reg_ref_rev = [];
    end;
    if SessInfo.snums_for>0
      SessInfo.reg_ref_for = SessInfo.B0reg_refs_for(ind);
    else
      SessInfo.reg_ref_for = [];
    end;
  else
    if SessInfo.snums_rev>0
      ind = floor(median(1:SessInfo.nscans_rev));
      SessInfo.reg_ref_rev = SessInfo.snums_rev(ind);
    else
      SessInfo.reg_ref_rev = [];
    end;
    if SessInfo.snums_for>0
      ind = floor(median(1:SessInfo.nscans_for));
      SessInfo.reg_ref_for = SessInfo.snums_for(ind);
    elseif SessInfo.snums_undef>0
      ind = floor(median(1:SessInfo.nscans_undef));
      SessInfo.reg_ref_for = SessInfo.snums_undef(ind);
    else
      SessInfo.reg_ref_for = [];
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select reference for registration to T1
function SessInfo = set_regT1_ref(SessInfo)
  if (SessInfo.nscans_SEpep>0 ||...
     (SessInfo.nscans_SE_rev > 0 && SessInfo.nscans_SE_rev > 0)) && ...
     (~isempty(SessInfo.B0uw_refs_for) && ~isempty(SessInfo.B0uw_refs_rev))
    if SessInfo.revflag
      ind = floor(median(1:length(SessInfo.B0uw_refs_rev)));
      SessInfo.regT1_ref = SessInfo.B0uw_refs_rev(ind);
    else
      ind = floor(median(1:length(SessInfo.B0uw_refs_for)));
      SessInfo.regT1_ref = SessInfo.B0uw_refs_for(ind);
    end;
    SessInfo.regT1_ref_rev = SessInfo.B0uw_refs_rev(ind);
    SessInfo.regT1_ref_for = SessInfo.B0uw_refs_for(ind);
  else
    if SessInfo.revflag
      SessInfo.regT1_ref = SessInfo.reg_ref_rev;
    else
      SessInfo.regT1_ref = SessInfo.reg_ref_for;
    end;
    SessInfo.regT1_ref_rev = SessInfo.reg_ref_rev;
    SessInfo.regT1_ref_for = SessInfo.reg_ref_for;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select B0uw and B0reg reference scans for each scan
function ScanInfo = set_scan_refs(ScanInfo,SessInfo)
  for s=1:SessInfo.nscans
    ScanInfo(s).B0uw_ref_for = [];
    ScanInfo(s).B0uw_ref_rev = [];
    ScanInfo(s).B0reg_ref_for = [];
    ScanInfo(s).B0reg_ref_rev = [];
  end;
  if ~isempty(SessInfo.B0uw_refs_for)
    % choose B0uw_ref for each scan 
    for s=1:SessInfo.nscans
      if ismember(s,SessInfo.snums_undef)
        continue;
      elseif ismember(s,SessInfo.snums_rev)
        B0uw_refs = SessInfo.B0uw_refs_rev;
      else
        B0uw_refs = SessInfo.B0uw_refs_for;
      end;
      tmp = B0uw_refs - s;
      if any(tmp<0) % use most recent reference
        tmp(tmp>0) = Inf; % ignore future scans
        [minval,ind]=min(abs(tmp));
      else % ref scan is here or ahead, find closest
        [minval,ind]=min(tmp);
      end;
      ScanInfo(s).B0uw_ref_for = SessInfo.B0uw_refs_for(ind);
      ScanInfo(s).B0uw_ref_rev = SessInfo.B0uw_refs_rev(ind);
      if ~isempty(SessInfo.B0reg_refs_for) ||...
         ~isempty(SessInfo.B0reg_refs_rev)
        switch ScanInfo(s).ScanType
          case 'BOLD_bu' % bottom-up "forward"
            ScanInfo(s).B0reg_ref_for = SessInfo.B0reg_refs_for(ind);
          case 'BOLD_td' % top-down  "reverse"
            ScanInfo(s).B0reg_ref_rev = SessInfo.B0reg_refs_rev(ind);
          case 'BOLD_ipp' % integrated pepolar
            if ismember(ScanInfo(s).pepolar,[0,2,3])
              ScanInfo(s).B0reg_ref_for = SessInfo.B0reg_refs_for(ind);
            end;
            if ismember(ScanInfo(s).pepolar,[1,2,3])
              ScanInfo(s).B0reg_ref_rev = SessInfo.B0reg_refs_rev(ind);
            end;
          case 'BOLD_ape' % alternating pepolar
            if ismember(ScanInfo(s).pepolar,[0,2,3])
              ScanInfo(s).B0reg_ref_for = SessInfo.B0reg_refs_for(ind);
            end;
            if ismember(ScanInfo(s).pepolar,[1,2,3])
              ScanInfo(s).B0reg_ref_rev = SessInfo.B0reg_refs_rev(ind);
            end;
          case 'BOLD_pep' % spin-echo pepolar
            ScanInfo(s).B0reg_ref_for = SessInfo.B0uw_refs_for(ind);
            ScanInfo(s).B0reg_ref_rev = SessInfo.B0uw_refs_rev(ind);
          case 'BOLD_SE'
            switch ScanInfo(s).pepolar
              case 0
                ScanInfo(s).B0reg_ref_for = SessInfo.B0uw_refs_for(ind);
              case 1
                ScanInfo(s).B0reg_ref_rev = SessInfo.B0uw_refs_rev(ind);
            end;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set file stems for each scan according to type of scan
function [ScanInfo,SessInfo] = set_scan_fstems(ScanInfo,SessInfo,parms)
  for s=1:SessInfo.nscans
    switch ScanInfo(s).ScanType
      case 'BOLD_bu' % "forward"
        ScanInfo(s).fstem = sprintf('%s%d_for',parms.fnamestem,s);
        ScanInfo(s).fstemlist = {ScanInfo(s).fstem};
      case 'BOLD_td' % "reverse"
        ScanInfo(s).fstem = sprintf('%s%d_rev',parms.fnamestem,s);
        ScanInfo(s).fstemlist = {ScanInfo(s).fstem};
      case 'BOLD_ape' % alternating phase-encode polarity
        if ismember(ScanInfo(s).pepolar,[0,1])
          revflag = ScanInfo(s).pepolar;
        else
          revflag = parms.revflag;
        end;
        switch revflag
          case 0
            ScanInfo(s).fstem = sprintf('%s%d_for',parms.fnamestem,s);
          case 1
            ScanInfo(s).fstem = sprintf('%s%d_rev',parms.fnamestem,s);
          case 2
            ScanInfo(s).fstem = sprintf('%s%d_ape',parms.fnamestem,s);
        end;
        ScanInfo(s).fstemlist = {ScanInfo(s).fstem};
      case 'BOLD_ipp' % integrated pepolar
        if ismember(ScanInfo(s).pepolar,[0,1])
          revflag = ScanInfo(s).pepolar;
        else
          revflag = SessInfo.revflag;
        end;
        if ~revflag
          ScanInfo(s).fstem = sprintf('%s%d_for',parms.fnamestem,s);
        else
          ScanInfo(s).fstem = sprintf('%s%d_rev',parms.fnamestem,s);
        end;
        if ismember(ScanInfo(s).pepolar,[0,1])
          ScanInfo(s).fstemlist = {ScanInfo(s).fstem};
        else
          ScanInfo(s).fstemlist = {...
            sprintf('%s%d_for',parms.fnamestem,s);
            sprintf('%s%d_rev',parms.fnamestem,s);
          };
        end;
      case 'BOLD_pep' % spin-echo pepolar
        if ~SessInfo.revflag
          ScanInfo(s).fstem = sprintf('%s%d_for',parms.fnamestem,s);
        else
          ScanInfo(s).fstem = sprintf('%s%d_rev',parms.fnamestem,s);
        end;
        ScanInfo(s).fstemlist = {...
          sprintf('%s%d_for',parms.fnamestem,s);
          sprintf('%s%d_rev',parms.fnamestem,s);
        };        
      case 'BOLD_SE' % spin-echo field maps for Siemens/Philips
        switch ScanInfo(s).pepolar
          case 0
            ScanInfo(s).fstem = sprintf('%s%d_for',parms.fnamestem,s);
          case 1
            ScanInfo(s).fstem = sprintf('%s%d_rev',parms.fnamestem,s);
        end;
        ScanInfo(s).fstemlist = {ScanInfo(s).fstem};
      otherwise      % "neutral"
        ScanInfo(s).fstem = sprintf('%s%d',parms.fnamestem,s);
        ScanInfo(s).fstemlist = {ScanInfo(s).fstem};
    end;
  end;
  
  % set fstem for regT1 reference scan
  SessInfo.fstem_regT1_ref = [];
  refsnum = SessInfo.regT1_ref;
  if ~isempty(refsnum)
    fstem_ref = ScanInfo(refsnum).fstem;
    if strcmp(ScanInfo(refsnum).ScanType,'BOLD_ape')
      if SessInfo.revflag
        fstem_ref = regexprep(fstem_ref,'ape','rev');
      else
        fstem_ref = regexprep(fstem_ref,'ape','for');
      end;      
    end;
    fname = sprintf('%s/%s.mgz',parms.ContainerPath,fstem_ref);
    [M,volsz] = mmil_load_mgh_info(fname);
    if volsz(4)>1, fstem_ref = [fstem_ref '_f0']; end;
    SessInfo.fstem_regT1_ref = fstem_ref;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
