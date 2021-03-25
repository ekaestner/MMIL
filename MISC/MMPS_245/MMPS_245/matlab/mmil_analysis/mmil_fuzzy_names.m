function fuzzy_names = mmil_fuzzy_names(varargin)
%function fuzzy_names = mmil_fuzzy_names([options])
%
% Purpose:
%  return cell array of anatomical names for fuzzy clusters
%   derived from genetic correlation of cortical thickness
%
% Optional Parameters:
%   'fuzzy_fstem': file stem for fuzzy cluster ROIs
%     {default = 'fuzzy'}
%   'fuzzy_order': [2|4|12|18] number of fuzzy cluster ROIs
%     note: set of 18 includes combined sets of 2, 4, and 12
%     {default = 18}
%
% Created:  01/09/13 by Don Hagler
% Last Mod: 01/09/13 by Don Hagler
%

fuzzy_names = [];
parms = mmil_args2parms(varargin,{...
  'fuzzy_fstem','fuzzy',[],...
  'fuzzy_order',18,[2,4,12,18],...
...
  'order2_labels',{'frontal','posterior'},[],...
  'order4_labels',{'frontal','occipital','temporal','parietal'},[],...
  'order12_labels',{'central','occipital','posterolateraltemporal',...
    'superiorparietal','orbitofrontal','superiortemporal','inferiorparietal',...
    'dorsomedialfrontal','anteromedialtemporal','precuneus',...
    'dorsolateralprefrontal','parsopercularis'},[],...
});

fuzzy_names = cell(1,parms.fuzzy_order);

switch parms.fuzzy_order
  case {2,4,12}
    fuzzy_labels = parms.(sprintf('order%d_labels',parms.fuzzy_order));
    for i=1:parms.fuzzy_order
      fuzzy_names{i} = sprintf('%s%d_%s',...
          parms.fuzzy_fstem,parms.fuzzy_order,fuzzy_labels{i});
    end;
  otherwise
    fuzzy_labels = cat(2,parms.order2_labels,parms.order4_labels,parms.order12_labels);
    for i=1:parms.fuzzy_order
      if i<=2
        fuzzy_order = 2;
      elseif i<=6
        fuzzy_order = 4;
      else
        fuzzy_order = 12;
      end;
      fuzzy_names{i} = sprintf('%s%d_%s',...
        parms.fuzzy_fstem,fuzzy_order,fuzzy_labels{i});
    end;
end;

