function rc_get_offsets(varargin)
%function rc_get_offsets[options])
%
% 'prefix_stem': stem of RCSE output prefix
%   {default = 'RCSE'}
% 'areas': vector of area numbers
%   If empty, assume idependent offset fits
%     for different visual areas were not done
%   {default = []}
% 'hemis':  vector of hemifield indices (1=right hemifield, 2=left hemifield)
%   If empty, assume idependent offset fits
%     for different hemifields were not done
%   {default = []}
% 'uplows':  vector of upper/lower field indices (1=upper field, 2=lower field)
%   If empty, assume idependent offset fits
%     for different fields were not done
%   {default = []}
% 'eccs': vector of eccentricity indices
%   If empty, assume idependent offset fits
%     for different eccentricities were not done
%   {default = []}
% 'contvec': vector of contrast indices used simultaneously for fits
%   If empty, assume all contrast levels were used
%   {default = []}
% 'polarity_penalty': penalty for positive or negative polarity source waweforms
%   {default = 0}
% 'fname_out': name of output file
%   If empty, write to stdout
%   {default = []}
% 'forceflag': [0|1[ whether to overwrite output file
%   {default = 1}
%
% Created:  02/05/09 by Don Hagler
% Last Mod: 02/19/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin, { ...
  'prefix_stem','RCSE',[],...
  'areas',[],[],...
  'hemis',[],[],...
  'uplows',[],[],...
  'eccs',[],[],...
  'contvec',[],[],...
  'polarity_penalty',0,[-Inf,Inf],...
  'fname_out',[],[],...
  'forceflag',true,[false true],...
...
  'cond_offsets_flag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.areas), parms.areas = 0; end;
if isempty(parms.hemis), parms.hemis = 0; end;
if isempty(parms.uplows), parms.uplows = 0; end;
if isempty(parms.eccs), parms.eccs = 0; end;

if ~isempty(parms.fname_out) && exist(parms.fname_out,'file') && ~parms.forceflag
  return;
end;

% open csv file for output
if isempty(parms.fname_out)
  fid = 1; % write to stdout instead
else
  fid = fopen(parms.fname_out,'wt');
  if fid==-1
    error('failed to open %s for writing',parms.fname_out);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'"area", "hemi", "uplow", "ecc", "r_offset", "th_offset", "min_error"\n');

for a=parms.areas
  for h=parms.hemis
    for u=parms.uplows
      for e=parms.eccs
        prefix = parms.prefix_stem;
        if h~=0
          prefix = sprintf('%s_hemi%d',prefix,h);
        end;
        if u~=0
          prefix = sprintf('%s_uplow%d',prefix,u);
        end;
        if e~=0
          prefix = sprintf('%s_ecc%d',prefix,e);
        end;
        if length(parms.contvec)==1
          prefix = sprintf('%s_cont%d',...
            prefix,parms.contvec);
        elseif ~isempty(parms.contvec)
          prefix = sprintf('%s_cont%s',...
            prefix,sprintf('_%d',parms.contvec));
        end;
        if a~=0
          prefix = sprintf('%s_areas_%d',prefix,a);
        end;
        if parms.cond_offsets_flag
          prefix = sprintf('%s_condoffsets',prefix);
        end;
        if parms.polarity_penalty~=0
          prefix = sprintf('%s_polpenalty%0.0f',...
            prefix,parms.polarity_penalty);
        end;
        matfile = sprintf('matfiles/%s_ret_forward.mat',prefix);
        if ~exist(matfile,'file')
          fprintf('%s: WARNING: file %s not found\n',mfilename,matfile);
          continue;
        end;
        clear results;
        load(matfile);
        fprintf(fid,'%d, %d, %d, %d, %0.2f, %0.2f, %0.4f\n',...
          a,h,u,e,...
          retforward.retmap.r_offset,...
          retforward.retmap.th_offset,...
          retforward.min_error);
      end;
    end;
  end;
end;

if fid~=1, fclose(fid); end;

