function fnames = fs_tal2labels(fname_points,varargin)
%function fnames = fs_tal2labels(fname_points,[options])
%
% Purpose: sample Talairach coordinates in a csv file
%   to a subject's cortical surfaces and create an annotation file
%
% Required Input:
%   fname_points: csv file containing Talairach coordinates
%     with column headers: 'roiname', 'hemi', 'x', 'y', 'z'
%     hemi must be 'lh', or 'rh'
%     additional columns will be ignored
%     
% Optional Input:
%   'outdir': output directory
%     {default = [pwd '/labels']}
%   'subj': subject name
%     {default = 'fsaverage'}
%   'MNIflag': [0|1] whether coordinates are MNI or "true" Talairach
%     {default = 1}
%   'vol_sm': smoothing sigma in volume
%     {default = 0}
%   'surf_sm': number of smoothing iterations on surface
%     {default = 50}
%   'thresh': threshold applied after smoothing
%     {default = 0.1}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%     
% Output:
%   fnames: cell array of label files created
%
% Created:  11/03/11 by Don Hagler
% Last Mod: 09/05/12 by Don Hagler
%

%% todo: get MNIflag from fname_points
%%  process MNIflag=0 and MNIflag=1 separately

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir',[pwd '/labels'],[],...
  'subj','fsaverage',[],...
  'MNIflag',true,[false true],...
  'vol_sm',0,[0,1000],...
  'surf_sm',50,[0,1000],...
  'thresh',0.1,[0,1],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
...
  'surf_tags',{'subj','subjdir','MNIflag','surfname','vol_sm','surf_sm',...
               'thresh','binarize_flag','verbose','tmpdir'},[],...
});
fnames = {};

if isempty(parms.tmpdir)
  parms.tmpdir = [parms.outdir '/tmp'];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);
roi_info = mmil_csv2struct(fname_points);

% check if expected output files exist
if ~parms.forceflag
  runflag = 0;
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    ind = find(strcmp({roi_info.hemi},hemi));
    npoints = length(ind);
    for i=1:npoints
      j = ind(i);
      roiname = roi_info(j).roiname;
      fname_out = sprintf('%s/%s.%s.label',parms.outdir,hemi,roiname);
      if ~exist(fname_out,'file')
        runflag = 1;
      else
        fnames{end+1} = fname_out;
      end;
    end;
  end;
  if ~runflag
    return;
  else
    fnames = {};
  end;
end;  

% sample points to labels for each hemisphere
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  ind = find(strcmp({roi_info.hemi},hemi));
  npoints = length(ind);
  points = zeros(npoints,3);
  for i=1:length(ind)
    j = ind(i);
    points(i,1) = roi_info(j).x;
    points(i,2) = roi_info(j).y;
    points(i,3) = roi_info(j).z;
  end;
  args = mmil_parms2args(parms,parms.surf_tags);
  vals = fs_tal2surf(points,'hemilist',hemi,args{:});
  for i=1:length(ind)
    j = ind(i);
    roiname = roi_info(j).roiname;
    fname_out = sprintf('%s/%s.%s.label',parms.outdir,hemi,roiname);
    if ~exist(fname_out,'file') || parms.forceflag
      v = find(vals(:,i));
      fs_write_label(v,fname_out,parms.subj);
    end;
    fnames{end+1} = fname_out;
  end;
end;

