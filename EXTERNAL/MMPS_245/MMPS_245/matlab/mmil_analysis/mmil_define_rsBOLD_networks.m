function networks = mmil_define_rsBOLD_networks(roi_type)
%function networks = mmil_define_rsBOLD_networks(roi_type)
%
% Purpose: create struct containing cell arrays listing
%   members of cortical "networks" for use in resting-state
%   BOLD analysis
%
% Optional Input:
%   roi_type: type of cortical ROI
%     allowed: 'fparc', 'aparc', 'talrc'
%     {default = 'fparc'}
%
% Created:  12/17/12 by Don Hagler
% Prev Mod: 08/28/14 by Don Hagler
% Last Mod: 08/11/17 by Don Hagler
%

networks = [];
if ~exist('roi_type','var') || isempty(roi_type)
  roi_type = 'fparc';
end;

switch roi_type
  case 'fparc'
    networks.default = {'mPFC','pCC','pIPL','aSFC','aMTG'}; % 'PHC'
%    networks.defaultA = {'mPFC','pCC','pIPL'};
%    networks.defaultB = {'pCC','aSFC','aMTG'};
    networks.attention = {'FEF','IPS','iPCS'}; % 'MT'
    networks.control = {'aCC','aIPL','aPFC','dlPFC','pMTG'};
%    networks.controlA = {'aCC','aIPL','aPFC'};
%    networks.controlB = {'aCC','dlPFC','pMTG'};
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

