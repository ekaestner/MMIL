function networks = abcd_define_rsBOLD_networks(roi_type)
%function networks = abcd_define_rsBOLD_networks(roi_type)
%
% Purpose: create struct containing cell arrays listing
%   members of cortical "networks" for use in resting-state
%   BOLD analysis
%
% Optional Input:
%   roi_type: type of cortical ROI
%     allowed: 'gordon','fparc', 'aparc', 'talrc'
%     {default = 'fparc'}
%
% Created:  12/17/12 by Don Hagler
% Prev Mod: 08/11/17 by Don Hagler
% Last Mod: 10/09/17 by Don Hagler
%

networks = [];
if ~exist('roi_type','var') || isempty(roi_type)
  roi_type = 'fparc';
end;

switch roi_type
  case 'gordon'
    indir = [getenv('HOME') '/ProjInfo/network_gordon'];
    instem = 'gordon_networks';
    % read network info from csv
    fname_in = sprintf('%s/%s.csv',indir,instem);
    if ~exist(fname_in,'file')
      error('network file %s not found',fname_in);
    end;
    rois = mmil_csv2struct(fname_in);
    % find ROIs in each unique network
    roi_names = {rois.roiname};
    each_network = {rois.network};
    network_names = unique(each_network);
    % define network struct
    networks = [];
    for i=1:length(network_names)
      network_name = network_names{i};
      ind = find(strcmp(each_network,network_name));
      networks.(network_name) = roi_names(ind);
    end;
  case 'fparc'
    networks.default = {'mPFC','pCC','pIPL','aSFC','aMTG'}; % 'PHC'
    networks.attention = {'FEF','IPS','iPCS'}; % 'MT'
    networks.control = {'aCC','aIPL','aPFC','dlPFC','pMTG'};
    networks.reference = {'AC','MC','VC'};
  case 'talrc'
    networks.default = {'mPFC','pCC','pIPL'};
    networks.attention = {'FEF','IPS','MT'};
    networks.control = {'aCC','aINS','aIPL','aPFC','dlPFC'};
    networks.reference = {'AC','MC','VC'};
  case 'aparc'
    networks.default = {...
      'medialorbitofrontal',...
      'isthmuscingulate',...%'precuneus',...
      'inferiorparietal'};
    networks.attention = {...
      'precentral',...%'superiorfrontal',...
      'superiorparietal',...
      'inferiortemporal'};%,'middletemporal'};
    networks.control = {...
      'caudalanteriorcingulate',...
      'insula',...%'lateralorbitofrontal',...
      'supramarginal',...
      'rostralmiddlefrontal',...
      'caudalmiddlefrontal'};
    networks.reference = {...
      'transversetemporal',...
      'precentral',...
      'lateraloccipital'};
  otherwise
    error('invalid roi_type %s',roi_type);
end;

