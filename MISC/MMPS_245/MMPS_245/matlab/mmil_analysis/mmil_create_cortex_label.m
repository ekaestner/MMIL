function mmil_create_cortex_label(fspath,varargin)
%function mmil_create_cortex_label(fspath,[options])
%
% Purpose: create cortex label files from aparc.annot files
%
% Required Input:
%   fspath: FreeSurfer recon container path
%
% Created:  09/29/13 by Don Hagler
% Last Mod: 09/29/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'fspath',fspath,[],...
  'cortex_codes',[1001:1003,1005:1034,2001:2003,2005:2034],[],...
  'forceflag',false,[false true],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
});

parms.nhemi = length(parms.hemilist);

% check FreeSurfer recon exists
if ~exist(parms.fspath,'dir')
  error('FreeSurfer recon dir %s not found',parms.fspath);
end;
[parms.subjdir,parms.subj,text] = fileparts(parms.fspath);
parms.subj = [parms.subj text];

% convert aparc to cortex label
for h=1:parms.nhemi
  hemi = parms.hemilist{h};
  fname_aparc = sprintf('%s/label/%s.aparc.annot',parms.fspath,hemi);  
  fname_label = sprintf('%s/label/%s.cortex.label',parms.fspath,hemi);
  if ~exist(fname_label,'file') || parms.forceflag
    if ~exist(fname_aparc)
      error('aparc file %s not found',fname_aparc);
    end;
    % match roilabels to roicodes
    [roinums,roilabels] = fs_read_annotation(fname_aparc);
    [roicodes,roinames] = fs_colorlut;
    [tmp,ind_cortex] = intersect(roicodes,parms.cortex_codes);
    roicodes = roicodes(ind_cortex);
    roinames = roinames(ind_cortex);
    ind_cortex = [];
    for i=1:length(roilabels)
      roiname = sprintf('ctx-%s-%s',hemi,roilabels{i});
      ind = find(strcmp(roiname,roinames));
      if ~isempty(ind)
        ind_cortex = union(ind_cortex,i);
      end;
    end;
    % set v for cortex vertices
    v = find(ismember(roinums,ind_cortex));
    % write label file
    fs_write_label(v,fname_label,parms.subj);
  end;
end;

