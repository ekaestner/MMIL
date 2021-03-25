function fs_write_ctab(roinames,fname_ctab,roicolors)
%function fs_write_ctab(roinames,fname_ctab,roicolors)
%
% Purpose: create custom color table file for a set of ROIs
%   using randomly selected, unique colors for each ROI
%
% Required Input:
%   roinames: cell array of ROI names (e.g. corresponding to label files)
%
% Optional Input:
%   fname_ctab: output file name
%     {default = 'custom.ctab'}
%   roicolors: matrix of rgb color vectors for each ROI (nroi x 3)
%     If empty, will randomly select unique colors for each ROI
%     {default = []}
% 
% Created:  12/01/11 by Don Hagler
% Last Mod: 09/11/15 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('fname_ctab','var') || isempty(fname_ctab)
  fname_ctab = 'custom.ctab';
end;
if ~exist('roicolors','var'), roicolors = []; end;

if ~iscell(roinames), roinames = {roinames}; end;
roinames = mmil_colvec(roinames);

nroi = length(roinames);
% check roicolors has right number of entries
if ~isempty(roicolors)
  if numel(roicolors) ~= nroi*3 || any(size(roicolors)~=[nroi,3])
    error('size of roicolors must be [%d,%d], but is [%d,%d]',...
      nroi,3,size(roicolors,1),size(roicolors,2));
  end;
end;

% find uniq roinames without changing order
[uniq_roinames,ind_uniq] = unique(roinames,'first');
if length(uniq_roinames)<length(roinames)
  ind_uniq = sort(ind_uniq);
  roinames = mmil_colvec(roinames(ind_uniq));
  if ~isempty(roicolors)
    roicolors = roicolors(ind_uniq,:);
  end;
  nroi = length(roinames);
end;

if isempty(roicolors)
  nsteps = max(2,ceil(nroi^(1/3))+1);
  x = round(linspace(0,255,nsteps));
  [R,G,B] = meshgrid(x,x,x); % create 3D matrices of colors
  color_order = randperm(numel(R)); % randomly shuffle colors
  % do not use gray (reserve for unlabeled)
  color_order = color_order(R~=128 | G~=128 | B~=128);
  % do not use black (treated as transparent, non-label)
  color_order = color_order(color_order~=1);
  roicolors = zeros(nroi,3);
  for i=1:nroi
    j = color_order(i);
    %% todo: require color components to be different from each other
    %%       to avoid gray
    roicolors(i,:) = [R(j),G(j),B(j)];
  end;
end;

if ~ismember('unknown',roinames)
  roinames = cat(1,'unknown',roinames);
  roicolors = cat(1,[0,0,0],roicolors);
  nroi = length(roinames);
end;

if ismember('unlabeled',roinames)
  i = find(strcmp('unlabeled',roinames));
  roicolors(i,:) = [128,128,128];
end;

fid = fopen(fname_ctab,'wt');
if fid<0, error('failed to open %s for writing',fname_ctab); end;
for i=1:nroi
  rgb = sprintf(' %3d',roicolors(i,:));
  fprintf(fid,'%4d  %-40s %s  0\n',i-1,roinames{i},rgb);
end;
fclose(fid);

