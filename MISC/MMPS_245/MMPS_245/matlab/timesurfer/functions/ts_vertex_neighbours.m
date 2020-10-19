function cfg = ts_vertex_neighbours (cfg,subjdir,subj,hemi)
%function cfg = ts_vertex_neighbours (cfg,subjdir,subj)
%
% Required input:
%  cfg  is a structure specifying the configuration for MC clustering
%  subj is a string specifying the subject name
%  hemi is a string specifying the hemisphere, i.e. 'lh' or 'rh'
%
% Output:
%   cfg is a structure containing:
%     neighbours.label       -- a string specifying the label for each vertex
%     neighbours.neighblabel -- a cell array of strings specifying the
%     labels for each neighbour associated with the vertex defined in
%     neighbours.label
%
% created:        04/27/10 Andrei Irimia
% last modified:  04/27/10 Andrei Irimia

if nargin < 4, help(mfilename); return; end;

if ~ismember(hemi,{'lh','rh'}), error('hemi must be lh or rh (is %s)',hemi); end;

if ~exist('subjdir','var') | isempty(subjdir)
    subjdir = getenv('SUBJECTS_DIR');
    if isempty(subjdir), error('SUBJECTS_DIR not defined as an environment variable'); end;
end;

surf = fs_load_subj(subj,hemi,'white',0,subjdir);
labels = cellfun(@(x)(sprintf('vertex%d',x)),num2cell(1:surf.nverts),'uniformoutput',false);
for k = 1:surf.nverts
    x = surf.nbrs(k,find(surf.nbrs(k,:) ~= 0));
    cfg.neighbours{k}.label = labels{k};
    cfg.neighbours{k}.neighblabel = labels(x);
end
return