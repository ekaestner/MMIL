function [status] = fs_joinmasks(masklist,hemi,varargin)
%function [status] = fs_joinmasks(masklist,hemi,[options])
%
% Usage:
%  [status] = fs_joinmasks(hemi,masklist,'key1', value1,...);
%
% Required input:
%  hemi should be either 'lh' for left hemisphere or 'rh' for right hemi
%
%  masklist should be a cell array containing a list of file stems
%    corresponding to w mask files defining ROIs (regions of interest)
%    (e.g. maskstem1-lh.w, maskstem2-lh.w, etc.)
%    hemi and .w extension should be omitted
%
% Optional parameters:
%  'maskdir' - directory containing mask files
%    {default: '.'}
%  'outstem'     - file stem for output joined mask
%    {default = joinedmasks}
%  'outdir' - directory for output file
%    {default: '.'}
%
% Output:
%  status is returned with value of 1 if no errors occured, 0 otherwise
%
%  created:       03/29/06   by Don Hagler
%  last modified: 08/20/07   by Don Hagler
%
% see also: fs_read_wfile(), fs_write_wfile
%

status = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options
if nargin < 2, help(mfilename); return; end;
try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end
  end
  if ~isempty( varargin ), g=struct(options{:}); 
  else g= []; end
catch
  error('calling convention {''key'', value, ... } error');
end    

% set defaults if not already set by user
try
  if ~isstr(g.outstem)
    g.outstem = 'joinedmasks';
  end;
catch
  g.outstem = 'joinedmasks';
end;
try
  if ~isstr(g.maskdir)
    g.maskdir = '.';
  end;
catch
  g.maskdir = '.';
end;
try
  if ~isstr(g.outdir)
    g.outdir = '.';
  end;
catch
  g.outdir = '.';
end;

% catch unrecognized options
gfields = fieldnames(g);
for index=1:length(gfields)
  switch gfields{index}
  case {'outstem' 'maskdir' 'outdir'},;
  otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
  end
end
maskdir  = g.maskdir;
outstem  = g.outstem;
outdir  = g.outdir;

[success,msg] = mkdir(outdir);
if ~success
  error('failed to create outdir %s: %s',outdir,msg);
end;

% catch bad inputs
Nroi = length(masklist);
if ~Nroi
  fprintf('%s: masklist is empty... quitting\n',mfilename);
  return;
end;
for r=1:Nroi
  try
    maskstem = masklist{r};
    if ~ischar(maskstem)
      fprintf('%s: masklist:\n',mfilename);
      disp(masklist);
      error('masklist should be a cell array of strings');
    else
      infile = sprintf('%s/%s-%s.w',maskdir,maskstem,hemi);
      if ~exist(infile,'file')
        error(sprintf('maskfile %s not found',infile);
      end;
    end;
  catch
    fprintf('%s: masklist:\n',mfilename);
    disp(masklist);
    error('masklist should be a cell array of strings');
  end;  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_mask=[];
for r=1:Nroi
  maskstem = masklist{r};
  fprintf('%s: loading %s...\n',mfilename,maskstem);

  % read mask wfile
  infile = sprintf('%s/%s-%s.w',maskdir,maskstem,hemi);
  [mask,v_tmp] = fs_read_wfile(infile);

  v_mask = union(v_mask,v_tmp);
end;
w = ones(size(v_mask));

% write wfile
outfile = sprintf('%s/%s-%s.w',outdir,outstem,hemi);
fprintf('%s: writing output %s...\n',mfilename,outfile);
fs_write_wfile(outfile,w,v_mask);

return;
