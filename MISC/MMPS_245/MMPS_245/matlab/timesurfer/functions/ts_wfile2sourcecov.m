function R = ts_wfile2sourcecov(wfile, varargin)
% R = ts_wfile2sourcecov(wfile,[options])
%
% Usage:
%  R = ts_wfile2sourcecov (wfile, 'key1', value1,...);
%
% Required input:
%  wfile - full or relative path name of w file (FreeSurfer surface data format)
%
% Optional parameters:
%  'decdips'     - vector of 0's and 1's indicating which vertices should be included
%                  in source covariance matrix
%  'thresh'      - threshold value
%    {default = 0}
%  'threshabs'   - [0 | 1] toggle whether threshold is applied to absolute values
%    {default = 1}
%  'maxvar'      - value assigned to diagonal element of covariance matrix for
%                  suprathreshold vertices
%    {default = 0.9}
%  'minvar'      - value assinged to diagonal for subthreshold vertices
%    {default = 0.09}
%  'ncomponents' - number of components per dipole (e.g. x,y,z)
%    {default = 3}
%
%  created:       03/23/06   by Don Hagler
%  last modified: 08/07/06   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse options
if nargin < 1
   help(mfilename)
   return
end

R=[];

if(~exist(wfile,'file'))
  fprintf('%s error: wfile %s not found\n',mfilename,wfile);
  return;
end
try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end
  end
  if ~isempty( varargin ), g=struct(options{:}); 
  else g= []; end
catch
  disp([mfilename ': calling convention {''key'', value, ... } error']); return;
end    

% set defaults if not already set by user
try, g.decdips;      catch, g.decdips     = [];   end
try, g.thresh;       catch, g.thresh      = 0;    end
try, g.threshabs;    catch, g.threshabs   = 1;    end
try, g.maxvar;       catch, g.maxvar      = 0.9;  end
try, g.minvar;       catch, g.minvar      = 0.09; end
try, g.ncomponents;  catch, g.ncomponents = 3;    end


% catch unrecognized options
gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'decdips' 'thresh' 'threshabs' 'maxvar' 'minvar' 'ncomponents'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end
end

% catch bad param values
g.ncomponents = max(g.ncomponents,1);


% open w file
[w,v]=fs_read_wfile(wfile);

% open dec file
if isempty(g.decdips)
  g.decdips = ones(length(v));
end
n_dips = length(g.decdips);
v_dec = find(g.decdips);
n_decdips = length(v_dec);

% threshold values
[w,v]=ts_thresh_weights(w,v,g.thresh,g.threshabs);

% expand to full number of vertices
w_dec = zeros(n_dips,1);
w_dec(v) = w;

% reduce to dec dips
w_dec = w_dec(v_dec);
v_nonzero = find(w_dec);
v_zero = find(~w_dec);

% create sparse matrix
nc = g.ncomponents;
ncomps=nc*n_decdips;
R=sparse(ncomps,ncomps);
c_nonzero=[];
c_zero=[];
for i=1:nc
  c_nonzero = [c_nonzero,v_nonzero*nc-(nc-i)];
  c_zero = [c_zero,v_zero*nc-(nc-i)];
end;

% set diagonals of R for suprathreshold dip components
i_diag = sub2ind(size(R),c_nonzero,c_nonzero);
R(i_diag) = g.maxvar;
i_diag = sub2ind(size(R),c_zero,c_zero);
R(i_diag) = g.minvar;

