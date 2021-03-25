function [hdr] = ts_get_header(data,varargin)

% Created by Jason Sherfey on 15-Oct-2008

datafield = ts_objecttype(data);  % (epochs, averages, timefreq)

parms = mmil_args2parms(varargin,{'condition',1,[],'zparam',[],[]},false);
if parms.condition>length(data.(datafield))
    fprintf('WARNING: condition specified is too large.  Returning header for condition 1.\n');
    parms.condition=1; 
end

if isempty(parms.zparam)
    fnames        = fieldnames(data.(datafield));
    parms.zparam    = intersect(fnames,{'data','power','cmplx','stdev','std_dev','plv','coherence','coh'});
end

hdr = rmfield(data,datafield);
hdr.(datafield) = rmfield(data.(datafield)(parms.condition),parms.zparam);
